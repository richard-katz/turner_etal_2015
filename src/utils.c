#include "mor.h"

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CheckSNESResidualNoLag"
PetscErrorCode CheckSNESResidualNoLag(SNES snes, Vec X, PetscReal *norm, 
				      AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  Parameter      *p = user->param;
  PetscBool      lag_visc_current = p->lag_visc;
  PetscFunctionBegin;  
  p->lag_visc = PETSC_FALSE;
  ierr = SNESComputeFunction(user->snes, X, user->R); CHKERRQ(ierr);
  ierr = VecNorm(user->R, NORM_2, norm); CHKERRQ(ierr);
  p->lag_visc = lag_visc_current;
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CornerFlowU"
PetscReal CornerFlowU(PetscReal x, PetscReal z)
/*---------------------------------------------------------------------*/
{
  PetscReal r, st, ct, th;
  r = sqrt(x*x+z*z);th = atan(x/z); st = x/r;  ct = z/r;  
  return 2.0*(th - ct*st)/PETSC_PI;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CornerFlowW"
PetscReal CornerFlowW(PetscReal x, PetscReal z)
/*---------------------------------------------------------------------*/
{
  PetscReal r, ct;
  r = sqrt(x*x+z*z); ct = z/r;
  return -2.0*ct*ct/PETSC_PI;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CornerFlowP"
PetscReal CornerFlowP(PetscReal x, PetscReal z)
/*---------------------------------------------------------------------*/
{
  PetscReal r, ct;
  r = sqrt(x*x+z*z); ct = z/r;
  return (-4.0/PETSC_PI/r + r)*ct;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "HalfSpaceCooling"
PetscReal HalfSpaceCooling(PetscReal x, PetscReal z, Parameter *p)
/*---------------------------------------------------------------------*/
{
  PetscReal  V;
  V = p->U0/1e2/SEC_PER_YR;      // velocity in units of m/s
  x *= p->height*1e3;            // Distance in m
  z *= p->height*1e3;            // Depth in m  
  return erf(z/(2*sqrt(p->K*abs(x)/V)));
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DestroyDataStructures"
PetscErrorCode DestroyDataStructures(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscFunctionBegin; 
  ierr = VecDestroy(&user->X);CHKERRQ(ierr);
  ierr = VecDestroy(&user->Xo);CHKERRQ(ierr);
  ierr = VecDestroy(&user->Xps);CHKERRQ(ierr);
  ierr = VecDestroy(&user->R);CHKERRQ(ierr);
  ierr = VecDestroy(&user->aux);CHKERRQ(ierr);
  ierr = SNESDestroy(&user->snes);CHKERRQ(ierr);
  ierr = DMDestroy(&user->da_aux);CHKERRQ(ierr);
  ierr = DMDestroy(&user->da);CHKERRQ(ierr);
  ierr = PetscBagDestroy(&user->bag);CHKERRQ(ierr); 
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "CreateDataStructures"
PetscErrorCode CreateDataStructures(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode ierr;
  PetscInt       dofs;
  Parameter      *p = user->param;
  PetscFunctionBegin;  

  /* set up solution and residual vectors */
  dofs = (PetscInt)(sizeof(Field)/sizeof(PetscReal));
  ierr = DMDACreate2d(user->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
		      DMDA_STENCIL_BOX,p->ni,p->nj,PETSC_DECIDE,PETSC_DECIDE,
		      dofs,2,PETSC_NULL,PETSC_NULL,&user->da);CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da,0,"U");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da,1,"W");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da,2,"P");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da,3,"T");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da,4,"A");CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(user->da,&user->X);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->X,"grid");CHKERRQ(ierr);
  ierr = VecDuplicate(user->X,&user->R);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->R,"residual");CHKERRQ(ierr);
  ierr = VecDuplicate(user->X,&user->Xo);CHKERRQ(ierr);
  ierr = VecDuplicate(user->X,&user->Xps);CHKERRQ(ierr);

  /* set up SNES */
  ierr = SNESCreate(user->comm,&user->snes);CHKERRQ(ierr);
  ierr = SNESSetDM(user->snes,user->da);CHKERRQ(ierr);
  ierr = SNESSetFunction(user->snes,user->R,FormResidual,user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(user->snes);CHKERRQ(ierr);

  /* set up auxilliary vector for output */
  dofs = (PetscInt)(sizeof(Auxfield)/sizeof(PetscReal));
  ierr = DMDACreate2d(user->comm,DMDA_BOUNDARY_NONE,DMDA_BOUNDARY_NONE,
		      DMDA_STENCIL_BOX,p->ni,p->nj,PETSC_DECIDE,PETSC_DECIDE,
		      dofs,0,PETSC_NULL,PETSC_NULL,&user->da_aux);CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da_aux,0,"eta");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da_aux,1,"e_2");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da_aux,2,"beta");CHKERRQ(ierr);
  ierr = DMDASetFieldName(user->da_aux,3,"decay");CHKERRQ(ierr);           // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,4,"growth");CHKERRQ(ierr);          // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,5,"advn");CHKERRQ(ierr);            // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,6,"dislocation");CHKERRQ(ierr);     // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,7,"diffusion");CHKERRQ(ierr);       // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,8,"plastic");CHKERRQ(ierr);         // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,9,"GBS");CHKERRQ(ierr);             // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,10,"Press");CHKERRQ(ierr);          // use with aux tests
  ierr = DMDASetFieldName(user->da_aux,11,"Water");CHKERRQ(ierr);          // use with aux tests
  ierr = DMCreateGlobalVector(user->da_aux,&user->aux);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)user->aux,"grid");CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetGhostedArray"
/* Gets array from a vector associated with a DMDA, with ghost points */
PetscErrorCode DMDAGetGhostedArray(DM da, Vec globvec, Vec *locvec, void *arr)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMGetLocalVector(da,locvec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,globvec,INSERT_VALUES,*locvec);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd  (da,globvec,INSERT_VALUES,*locvec);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,*locvec,arr); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "DMDARestoreGhostedArray"
/* Restores array from a vector associated with a DMDA, with ghost points */
PetscErrorCode DMDARestoreGhostedArray(DM da, Vec globvec, Vec *locvec, void *arr)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = DMDAVecRestoreArray(da,*locvec,arr); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,locvec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "DMDAGetGridInfo"
/* Get useful information from the DMDA */
PetscErrorCode DMDAGetGridInfo(DM da, int *is, int *js, int *ks, int *ie, 
			       int *je, int *ke, int *ni, int *nj, int *nk, 
			       int *dim)
/* ------------------------------------------------------------------- */
{
  PetscErrorCode ierr;
  PetscInt       im, jm, km;
  PetscFunctionBegin;
  ierr = DMDAGetCorners(da,is,js,ks,&im,&jm,&km); CHKERRQ(ierr);
  ierr = DMDAGetInfo(da,dim,ni,nj,nk,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  if (ie) *ie = *is + im; 
  if (je) *je = *js + jm; 
  if (ke) *ke = *ks + km;
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "OptionsHasName"
/* check whether an option has been set by the user on the command line */
PetscBool OptionsHasName(const char name[])
/* ------------------------------------------------------------------- */
{
  PetscBool retval; 
  PetscFunctionBegin;
  PetscOptionsHasName(PETSC_NULL,name,&retval);
  PetscFunctionReturn(retval);
}


