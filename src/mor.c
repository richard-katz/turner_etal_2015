#include "mor.h"

/*-----------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
/*-----------------------------------------------------------------------*/
{
  AppCtx         user;
  PetscErrorCode ierr;
  PetscInt       cs;
  PetscFunctionBegin;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize PETSc and options.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  ierr = PetscInitialize(&argc,&argv,(char *)0,PETSC_NULL); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&cs);
  PetscOptionsSetValue(NULL,"-ksp_gmres_restart","360");
  PetscOptionsSetValue(NULL,"-snes_monitor","");
  PetscOptionsSetValue(NULL,"-snes_converged_reason","");
  if (cs==1) {
    PetscOptionsSetValue(NULL,"-pc_type","lu");
  } else {
    if (OptionsHasName("-mumps")) {
      PetscOptionsSetValue(NULL,"-pc_factor_mat_solver_type","mumps"); 
      PetscOptionsSetValue(NULL,"-pc_type","lu");    
    } else {
      PetscOptionsSetValue(NULL,"-pc_type","asm");
      PetscOptionsSetValue(NULL,"-sub_pc_type","lu");
    }
  }
  user.comm = PETSC_COMM_WORLD;

 /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set up data structures.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  ierr = InitializeViscosityModule(user.comm,&user.VM); CHKERRQ(ierr);
  ierr = SetParams(&user);CHKERRQ(ierr);
  ierr = CreateDataStructures(&user); CHKERRQ(ierr);
  ierr = CreateInitialGuess(&user); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve with continuation
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */ 
  ierr = ContinuationSolve(&user);CHKERRQ(ierr);
  ierr = DoOutput(&user);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = FinalizeViscosityModule(&user.VM); CHKERRQ(ierr);
  ierr = DestroyDataStructures(&user);CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "ContinuationSolve"
PetscErrorCode ContinuationSolve(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *p = user->param;
  PetscReal      atol, rtol, stol;
  PetscReal      contin_old = 0, contin_new = 0;
  SNESConvergedReason conv;
  PetscFunctionBegin;

  if (p->const_visc) { // Constant viscosity solve

    p->contin = 0;
    ierr = PetscPrintf(user->comm,"v Solving with constant viscosity.\n");CHKERRQ(ierr);
    ierr = SNESSolve(user->snes,PETSC_NULL,user->X);CHKERRQ(ierr);
    PetscFunctionReturn(0);

  } else {             // Variable viscosity solve
    
      p->contin = 0;
      atol = 1e-1;  rtol = 5e-6; stol = 1e-16;
      ierr = SNESSetTolerances(user->snes,atol,rtol,stol,3,1000); 
      ierr = PetscPrintf(user->comm,"v Iterating with lagged viscosity, continuation parameter = %.3f\n",p->contin);CHKERRQ(ierr);
      ierr = SolveWithLaggedViscosity(5e-1,75,user); CHKERRQ(ierr);
      ierr = VecCopy(user->X,user->Xps);
      contin_old = p->contin;
      p->contin = 0.45;  
      PetscPrintf(user->comm,"-----------------------------------------------------------------\n");  
      while (p->contin < 1) {
      ierr = PetscPrintf(user->comm,"v Iterating with lagged viscosity, continuation parameter = %.3f\n",p->contin);CHKERRQ(ierr);
      ierr = SolveWithLaggedViscosity(5e-1,75,user); CHKERRQ(ierr);
         if (p->contin + 0.45 < 1) {
	   contin_new = p->contin + 0.45;
	   ierr = ContinuationCorrection(contin_old, contin_new, user); 
	   contin_old = p->contin;
	   p->contin += 0.45;
	 } else {
	   contin_new = 1;
	   p->contin += 0.45;
	 }

	 PetscPrintf(user->comm,"-----------------------------------------------------------------\n");
      }

    p->contin = 1;
    atol = 4e-4;  rtol = 5e-6; stol = 1e-16;
    ierr = SNESSetTolerances(user->snes,atol,rtol,stol,3,1000); 
    ierr = PetscPrintf(user->comm,"v Iterating with lagged viscosity, continuation parameter = %.3f\n",p->contin);CHKERRQ(ierr);
    ierr = SolveWithLaggedViscosity(p->FinalTolerance,350,user); CHKERRQ(ierr);
    
    if (p->do_unlagged_final_solve) {     
      PetscPrintf(user->comm,"-----------------------------------------------------------------\n");      
      p->lag_visc = PETSC_FALSE;
      ierr = PetscPrintf(user->comm,"v Solving with non-lagged viscosity\n");CHKERRQ(ierr);
      atol = 1e-6;  rtol = 1e-15; stol = 1e-16;
      ierr = SNESSetTolerances(user->snes,atol,rtol,stol,10,1000);                        
      ierr = VecCopy(user->X,user->Xo);
      ierr = SNESSolve(user->snes,PETSC_NULL,user->X);CHKERRQ(ierr);
      ierr = SNESGetConvergedReason(user->snes,&conv);CHKERRQ(ierr);
      if (conv<0) {
	ierr = PetscPrintf(user->comm,"^ Solve failed, using solution from last lagged iterate\n");CHKERRQ(ierr);
	ierr = VecCopy(user->Xo,user->X);
      }
    }
    PetscFunctionReturn(0);
  }
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "SolveWithLaggedViscosity"
PetscErrorCode SolveWithLaggedViscosity(PetscReal tol, PetscInt maxit,
					AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscInt       n=0;
  PetscReal      norm=1e10;
  PetscBool      lag_visc_current = user->param->lag_visc;
  user->param->lag_visc = PETSC_TRUE;
  while (norm>tol && n<maxit) {                                           
    ierr = VecCopy(user->X,user->Xo);                            
    ierr = SNESSolve(user->snes,PETSC_NULL,user->X);CHKERRQ(ierr);
    ierr = CheckSNESResidualNoLag(user->snes,user->X,&norm,user);CHKERRQ(ierr);
    n++;
  }
  ierr = PetscPrintf(user->comm,"^ Non-lagged SNES solution norm = %.3e (tolerance %.3e)\n",norm,tol);CHKERRQ(ierr);
  user->param->lag_visc = lag_visc_current;
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DoOutput"
PetscErrorCode DoOutput(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter      *p = user->param;
  PetscViewer    viewer;
  PetscErrorCode ierr;
  PetscFunctionBegin;  
  if (p->output_to_file) { 
    ierr = PetscPrintf(user->comm," Generating output file: \"%s\"\n",p->output_filename);
    ierr = PetscViewerBinaryOpen(user->comm,p->output_filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_BINARY_MATLAB);CHKERRQ(ierr);
    ierr = PetscBagView(user->bag,viewer);CHKERRQ(ierr);
    ierr = PetscBagView(user->VM.bag,viewer);CHKERRQ(ierr);
    ierr = VecView(user->X,viewer);CHKERRQ(ierr);
    ierr = UpdateAuxilliaryFields(user);CHKERRQ(ierr);
    ierr = VecView(user->aux,viewer);CHKERRQ(ierr);
    ierr = VecView(user->R,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "UpdateAuxilliaryFields"
PetscErrorCode UpdateAuxilliaryFields(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter      *p = user->param;
  PetscErrorCode ierr;
  PetscInt       i, is, ie, j, js, je;
  PetscReal      T, terms[4], dominant[4];
  //  PetscReal      WATER;
  Auxfield       **u;
  Field          **x;
  PetscReal      h, z;  //Needed to calculate lithostatic pressure
  Vec            v;

  PetscFunctionBegin;
  /* get the residual array WITHOUT ghost points */
  ierr = DMDAVecGetArray(user->da_aux,user->aux,&u);CHKERRQ(ierr);

  /* get the residual array WITH ghost points */
  ierr = DMDAGetGhostedArray(user->da,user->X,&v,&x);CHKERRQ(ierr);

  /* get dimensions of LOCAL piece of grid (in parallel, this is not full grid) */
  ierr = DMDAGetGridInfo(user->da,&is,&js,0,&ie,&je,0,0,0,0,0);CHKERRQ(ierr);
  

  h = pow((p->nj-4),-1);
  for (j=js;j<je;j++) {
    z = PetscMax(h*(j-1.5),0);
    for (i=is;i<ie;i++) {
      /* don't compute viscosity in buffer points */
      if (i<2 || j<2 || i>p->ni-3 || j>p->nj-3) { u[j][i].eta=1; u[j][i].e_2=1; u[j][i].beta=1; 
	u[j][i].growth=1; u[j][i].decay=1; u[j][i].advn=1;                // use with aux tests
	u[j][i].diffusion=1; u[j][i].plastic=1;    // use with aux tests
	u[j][i].dislocation=1; u[j][i].GBS=1;
	u[j][i].Press = 1;
	continue; }
      if (p->const_visc) { T = p->Tp + 273; }
      else               { T = x[j][i].T*p->Tp + 273; } 
      u[j][i].Press = TotalPressure(x,i,j,GRID_CENTER,user);
      u[j][i].e_2 = StrainRateSecInv(x,i,j,GRID_CENTER,user)*(p->U0/1e2/SEC_PER_YR)/(p->height*1e3);
      u[j][i].Water = WaterContent(j,p);
      //      WATER = WaterContent(j,p);
      //      ierr = ComputeViscosityModule(T,u[j][i].Press,u[j][i].e_2,p->grain_size,WATER,&u[j][i].eta,&u[j][i].beta,dominant,&user->VM);CHKERRQ(ierr);
      ierr = ComputeViscosityModule(T,u[j][i].Press,u[j][i].e_2,p->grain_size*exp(x[j][i].A),u[j][i].Water,&u[j][i].eta,&u[j][i].beta,dominant,&user->VM);CHKERRQ(ierr);
      u[j][i].eta = pow(p->eta0/u[j][i].eta + p->eta0/p->etamax,-1);
      ierr = GrainSizeEvolution(x,x,i,j,terms,user);
      u[j][i].decay  = terms[0];
      u[j][i].growth = terms[1];
      u[j][i].advn   = terms[2];
      u[j][i].diffusion = dominant[0];
      u[j][i].plastic = dominant[1];
      u[j][i].dislocation = dominant[2];
      u[j][i].GBS = dominant[3];
    }
  }

  ierr = DMDARestoreGhostedArray(user->da,user->X,&v,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da_aux,user->aux,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CreateInitialGuess"
PetscErrorCode CreateInitialGuess(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter      *p = user->param;
  PetscErrorCode ierr;
  PetscInt       i, is, ie, j, js, je;
  PetscReal      h, x, z;
  Field          **u;
  PetscFunctionBegin; 
  /* get the residual array WITHOUT ghost points */
  ierr = DMDAVecGetArray(user->da,user->X,&u);CHKERRQ(ierr);

  /* get dimensions of LOCAL piece of grid (in parallel, this is not full grid) */
  ierr = DMDAGetGridInfo(user->da,&is,&js,0,&ie,&je,0,0,0,0,0);CHKERRQ(ierr);

  h = pow((p->nj-4),-1);
  for (j=js;j<je;j++) {
    z = PetscMax(h*(j-1.5),0);
    for (i=is;i<ie;i++) {
      x = h*(i-1.5) - p->xridge/p->height;
      u[j][i].U = CornerFlowU(x+0.5*h,z);
      u[j][i].W = CornerFlowW(x,z+0.5*h);
      u[j][i].P = CornerFlowP(x,z);
      u[j][i].T = HalfSpaceCooling(x,z,p);
      u[j][i].A = log(1);                    // All grains start at a_o (p->grain_size)
      if (i==1 && j<=1) {
	u[j][i].U = 0;
	u[j][i].W = 0;
	u[j][i].P = 0;
	u[j][i].T = 0;
	u[j][i].A = 0;
      }
    }
  }

  ierr = DMDAVecRestoreArray(user->da,user->X,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}



/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "ContinuationCorrection"
PetscErrorCode ContinuationCorrection(PetscReal contin_old, PetscReal contin_new, AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter      *p = user->param;
  PetscErrorCode ierr;
  PetscInt       i, is, ie, j, js, je;
  PetscReal      alpha_step, delta_alpha;
  Field          **x, **x0, **xps;
  PetscFunctionBegin; 

  alpha_step = contin_new-p->contin;
  delta_alpha = p->contin - contin_old;

  ierr = VecCopy(user->X,user->Xo); 

  /* get the residual array WITHOUT ghost points */
  ierr = DMDAVecGetArray(user->da,user->X,&x);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,user->Xo,&x0);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(user->da,user->Xps,&xps);CHKERRQ(ierr);

  /* get dimensions of LOCAL piece of grid (in parallel, this is not full grid) */
  ierr = DMDAGetGridInfo(user->da,&is,&js,0,&ie,&je,0,0,0,0,0);CHKERRQ(ierr);


  for (j=js;j<je;j++) {
    for (i=is;i<ie;i++) {
      x[j][i].U = x0[j][i].U + (alpha_step/delta_alpha)*(x0[j][i].U - xps[j][i].U);
      x[j][i].W = x0[j][i].W + (alpha_step/delta_alpha)*(x0[j][i].W - xps[j][i].W);
    }
  }

  ierr = DMDAVecRestoreArray(user->da,user->X,&x);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,user->Xo,&x0);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(user->da,user->Xps,&xps);CHKERRQ(ierr);

  ierr = VecCopy(user->Xo,user->Xps);
  ierr = VecCopy(user->X,user->Xo);


  PetscFunctionReturn(0);
}
