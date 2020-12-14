#include "mor.h"
#include "advectionSchemes.h"
PetscReal StrainRateTensor(Field**, int, int, int, int, AppCtx*);
PetscReal DeviatoricStressTensor(Field**, Field**, int, int, int, int, AppCtx*);
PetscReal StrainRateSecInv(Field**, int, int, int, AppCtx*);
PetscReal ViscosityTensor(Field**, Field**, int, int, int, int, int, AppCtx*);
PetscReal ForceBalance_X(Field**, Field**, PetscInt, PetscInt, AppCtx*);
PetscReal ForceBalance_Z(Field**, Field**, PetscInt, PetscInt, AppCtx*);
PetscReal Continuity(Field**, PetscInt, PetscInt, AppCtx*);
PetscReal EnergyEquation(Field**, PetscInt, PetscInt, AppCtx*);
PetscReal PlateSpeed(PetscInt, Parameter*);
PetscReal WaterContent(PetscInt, Parameter*);
PetscReal TotalPressure(Field**, PetscInt, PetscInt, PetscInt, AppCtx*);
PetscReal Temperature(Field**, PetscInt, PetscInt, PetscInt, AppCtx*);

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormResidual"
PetscErrorCode FormResidual(SNES snes, Vec X, Vec R, void *ptr)
/* ------------------------------------------------------------------- */
{
  AppCtx         *user = (AppCtx*)ptr;
  Parameter      *p = user->param;
  PetscErrorCode ierr;
  PetscReal      x, z, h = pow(p->nj-4,-1);
  PetscInt       i, is, ie, j, js, je;
  Field          **r, **u, **u0;
  Vec            Xloc, Xoloc;
  PetscFunctionBegin;
  /* get the residual array WITHOUT ghost points */
  ierr = DMDAVecGetArray(user->da,R,&r);CHKERRQ(ierr);

  /* get the current guess of solution WITH ghost points */
  ierr = DMDAGetGhostedArray(user->da,X,&Xloc,&u);CHKERRQ(ierr);

  if (p->lag_visc) {
    /* get the last solution WITH ghost points */
    ierr = DMDAGetGhostedArray(user->da,user->Xo,&Xoloc,&u0);
  } else { u0 = u; }

  /* get dimensions of LOCAL piece of grid (in parallel, this is not full grid) */
  ierr = DMDAGetGridInfo(user->da,&is,&js,0,&ie,&je,0,0,0,0,0);CHKERRQ(ierr);

  /* z = 0 top boundary */
  if (js==0) {
    js += 2;
    for (j=0;j<2;j++) {
      for (i=is;i<ie;i++) {
	r[j][i].U = u[j][i].U + u[j+1][i].U - 2*PlateSpeed(i,p);
	r[j][i].W = u[j][i].W;
	r[j][i].P = u[j][i].P - u[j+1][i].P;
	r[j][i].T = u[j][i].T + u[j+1][i].T;
	r[j][i].A = u[j][i].A - u[j+1][i].A;
      }
    }
  }

  /* bottom boundary */
  if (je==p->nj) {
    je -= 2;
    for (j=p->nj-2;j<p->nj;j++) {
      for (i=is;i<ie;i++) {
	x = h*(i-1) - p->xridge/p->height; z = h*(j-1.5);
	if (p->const_visc) {
	  r[j][i].U = u[j][i].U - CornerFlowU(x,z);
	  r[j][i].W = u[j][i].W - CornerFlowW(x,z+0.5*h);
	} else {
	  /* non-constant viscosity */
	  r[j][i].U = u[j][i].U - u[j-1][i].U;
	  r[j][i].W = u[j][i].W - u[j-1][i].W;
	}
	r[j][i].P = u[j][i].P - u[j-1][i].P;
	r[j][i].T = u[j][i].T + u[j-1][i].T - 2;
	r[j][i].A = u[j][i].A + u[j-1][i].A - 2*log(p->BB_grain_size/p->grain_size); // set the inflow size of grain size with BB_grain_size
      }
    }
  }

  /* left boundary */
  if (is==0) {
    is += 2;
    for (j=js;j<je;j++) {
      for (i=0;i<2;i++) {
	if (i == 1) {
	  r[j][i].U = u[j][i].U;
	  r[j][i].W = u[j][i].W - u[j][i+1].W;
	  r[j][i].T = u[j][i].T - u[j][i+1].T;
	  r[j][i].P = u[j][i].P - u[j][i+1].P;
	  r[j][i].A = u[j][i].A - u[j][i+1].A;
	} else {
	  r[j][i].U = u[j][i].U + u[j][i+2].U;
	  r[j][i].W = u[j][i].W - u[j][i+3].W;
	  r[j][i].T = u[j][i].T - u[j][i+3].T;
	  r[j][i].P = u[j][i].P - u[j][i+3].P;
	  r[j][i].A = u[j][i].A - u[j][i+3].A;
	}
      }
    }
  }

  /* right boundary */
  if (ie==p->ni) {
    ie -= 2;
    for (j=js;j<je;j++) {
      for (i=p->ni-2;i<p->ni;i++) {
	x = h*(i-1) - p->xridge/p->height; z = h*(j-1); 
	if (p->const_visc) {
	  r[j][i].U = u[j][i].U - CornerFlowU(x,z-0.5*h);
	  r[j][i].W = u[j][i].W - CornerFlowW(x-0.5*h,z);
	} else {
	  /* non-constant viscosity */
	  r[j][i].U = u[j][i].U - u[j][i-1].U;
	  r[j][i].W = u[j][i].W - u[j][i-1].W;
	}
	r[j][i].P = u[j][i].P;
	r[j][i].T = u[j][i].T - u[j][i-1].T;
	r[j][i].A = u[j][i].A - u[j][i-1].A;
      }
    }
  }

  /* interior of domain */
  for (j=js;j<je;j++) {
    for (i=is;i<ie;i++) {
      r[j][i].U = ForceBalance_X(u,u0,i,j,user);
      r[j][i].W = ForceBalance_Z(u,u0,i,j,user);
      r[j][i].P = Continuity(u,i,j,user);
      r[j][i].T = EnergyEquation(u,i,j,user);
      r[j][i].A = GrainSizeEvolution(u,u0,i,j,NULL,user);
    }
  }

  /* clean up */
  ierr = DMDARestoreGhostedArray(user->da,X,&Xloc,&u);CHKERRQ(ierr);
  if (p->lag_visc) { ierr = DMDARestoreGhostedArray(user->da,user->Xo,&Xoloc,&u0);CHKERRQ(ierr); }
  ierr = DMDAVecRestoreArray(user->da,R,&r);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "EnergyEquation"
PetscReal EnergyEquation(Field **x, PetscInt i, PetscInt j, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  PetscReal h = pow(user->param->nj-4,-1);
  PetscReal advn, difn, v[4], q[9];
  difn = (x[j][i-1].T - 2*x[j][i].T + x[j][i+1].T + 
	  x[j-1][i].T - 2*x[j][i].T + x[j+1][i].T)/h/h;
  v[GRID_N] = x[j][i].W; v[GRID_S] = x[j-1][i].W;
  v[GRID_E] = x[j][i].U; v[GRID_W] = x[j][i-1].U;
  q[GRID_N] = x[j+1][i].T; q[GRID_NN] = x[j+2][i].T;
  q[GRID_S] = x[j-1][i].T; q[GRID_SS] = x[j-2][i].T;
  q[GRID_E] = x[j][i+1].T; q[GRID_EE] = x[j][i+2].T;
  q[GRID_W] = x[j][i-1].T; q[GRID_WW] = x[j][i-2].T;
  q[GRID_C] = x[j][i].T;
  advn = FrommAdvection(v,q,h,h);
  return advn - difn/user->param->Pe;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ForceBalance_X"
PetscReal ForceBalance_X(Field **x, Field **x0, PetscInt i, PetscInt j, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  PetscReal h = pow(user->param->nj-4,-1);
  PetscReal dPdx, sxxp, sxxm, szxp, szxm;
  dPdx = (x[j][i+1].P - x[j][i].P)/h;
  sxxp = DeviatoricStressTensor(x,x0,i+1,j  ,XX,GRID_CENTER,user);
  sxxm = DeviatoricStressTensor(x,x0,i  ,j  ,XX,GRID_CENTER,user);
  szxp = DeviatoricStressTensor(x,x0,i  ,j  ,ZX,GRID_IJPLUS,user);
  szxm = DeviatoricStressTensor(x,x0,i  ,j-1,ZX,GRID_IJPLUS,user);
  return dPdx - (sxxp-sxxm + szxp-szxm)/h;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ForceBalance_Z"
PetscReal ForceBalance_Z(Field **x, Field **x0, PetscInt i, PetscInt j, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  PetscReal h = pow(user->param->nj-4,-1);
  PetscReal dPdz, szzp, szzm, szxp, szxm;
  dPdz = (x[j+1][i].P - x[j][i].P)/h;
  szzp = DeviatoricStressTensor(x,x0,i  ,j+1,ZZ,GRID_CENTER,user);
  szzm = DeviatoricStressTensor(x,x0,i  ,j  ,ZZ,GRID_CENTER,user);
  szxp = DeviatoricStressTensor(x,x0,i  ,j  ,ZX,GRID_IJPLUS,user);
  szxm = DeviatoricStressTensor(x,x0,i-1,j  ,ZX,GRID_IJPLUS,user);
  return dPdz - (szzp-szzm + szxp-szxm)/h;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Continuity"
PetscReal Continuity(Field **x, PetscInt i, PetscInt j, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  PetscReal h = pow(user->param->nj-4,-1);
  return (x[j][i].U - x[j][i-1].U + x[j][i].W - x[j-1][i].W)/h;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DeviatoricStressTensor"
PetscReal DeviatoricStressTensor(Field **x, Field **x0, int i, int j, int comp, 
				 int pos, AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscReal sigma=0, dir[3] = {XX, ZZ, ZX};
  PetscInt  n;
  for (n=0;n<3;n++) {
    sigma += 2*ViscosityTensor(x,x0,i,j,comp,dir[n],pos,user)*StrainRateTensor(x,i,j,dir[n],pos,user);
  }
  return sigma;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "ViscosityTensor"
PetscReal ViscosityTensor(Field **x, Field **x0, int i, int j, int C1, int C2, 
			  int pos, AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter       *p = user->param;
  PetscErrorCode  ierr;
  PetscReal       eta, T, P, e_2 = 1.0, beta = 0.0, a = 0.0, Water = 0.0;  

  if (p->const_visc) { 
    /* constant viscosity */
    P   = 0;
    T   = p->Tp + 273; 
    e_2 = 1e-15;
    a   = p->grain_size;
  } else { 
    /* variable viscosity */
    P   = TotalPressure(x,i,j,pos,user);                                                // Pascals
    T   = Temperature(x,i,j,pos,user);                                              // Kelvin
    e_2 = StrainRateSecInv(x0,i,j,pos,user)*(p->U0/1e2/SEC_PER_YR)/(p->height*1e3); // sec^{-1}
    /* dimensional grain size, m */
    if (pos==GRID_CENTER)      { a = p->grain_size*exp(x0[j][i].A); } 
    else if (pos==GRID_IPLUS)  { a = p->grain_size*(exp(x0[j][i].A) + exp(x0[j][i+1].A))*0.5; } 
    else if (pos==GRID_JPLUS)  { a = p->grain_size*(exp(x0[j][i].A) + exp(x0[j+1][i].A))*0.5; }
    else if (pos==GRID_IJPLUS) { a = p->grain_size*(exp(x0[j][i].A) + exp(x0[j][i+1].A) + exp(x0[j+1][i].A) + exp(x0[j+1][i+1].A))*0.25; } 
  } 

  if (p->WetComp){
    Water = WaterContent(j,p);
  } else {
    Water = 0;
  }

  ierr = ComputeViscosityModule(T, P, e_2, a, Water, &eta, &beta, NULL, &user->VM);CHKERRQ(ierr);
  eta = pow(p->eta0/eta + p->eta0/p->etamax,-p->contin);                         // ??? factor of two comes from the harmonic mean equation

  /* entries in viscosity tensor */
  if      (C1==XX && C2==XX) { return eta; }
  else if (C1==ZZ && C2==ZZ) { return eta; }
  else if (C1==XX && C2==ZZ) { return 0; }
  else if (C1==ZZ && C2==XX) { return 0; }
  else if (C1==ZX && C2==ZX) { return eta; }
  else                       { return 0; }
  return 0;
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "StrainRateTensor"
PetscReal StrainRateTensor(Field **x, int i, int j, int com, 
			   int pos, AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscReal e=1e10, h = pow(user->param->nj-4,-1);
  if (pos==GRID_CENTER) {
    if (com==XX) { e = (x[j][i].U - x[j][i-1].U)/h; } // chk
    if (com==ZZ) { e = (x[j][i].W - x[j-1][i].W)/h; } // chk
    if (com==ZX) { e = (0.25*(x[j][i+1].W+x[j-1][i+1].W-x[j][i-1].W-x[j-1][i-1].W)/h +
			0.25*(x[j+1][i].U+x[j+1][i-1].U-x[j-1][i].U-x[j-1][i-1].U)/h)/2; } // chk
  } else if (pos==GRID_IPLUS) {
    if (com==XX) { e = 0.5*(x[j][i+1].U-x[j][i-1].U)/h; } // chk
    if (com==ZZ) { e = 0.5*(x[j][i].W+x[j][i+1].W-x[j-1][i].W-x[j-1][i+1].W)/h; } // chk
    if (com==ZX) { e = (0.5*(x[j][i+1].W+x[j-1][i+1].W-x[j][i].W-x[j-1][i].W)/h +
			0.5*(x[j+1][i].U-x[j-1][i].U)/h)/2; } // chk
  } else if (pos==GRID_JPLUS) {
    if (com==XX) { e = 0.5*(x[j][i].U+x[j+1][i].U-x[j][i-1].U-x[j+1][i-1].U)/h; } // chk
    if (com==ZZ) { e = 0.5*(x[j+1][i].W-x[j-1][i].W)/h; }
    if (com==ZX) { e = (0.5*(x[j][i+1].W-x[j][i-1].W)/h +
			0.5*(x[j+1][i].U+x[j+1][i-1].U-x[j][i].U-x[j][i-1].U)/h)/2; } // chk
  } else if (pos==GRID_IJPLUS) {
    if (com==XX) { e = 0.25*(x[j][i+1].U+x[j+1][i+1].U-x[j][i-1].U-x[j+1][i-1].U)/h; } // chk
    if (com==ZZ) { e = 0.25*(x[j+1][i].W+x[j+1][i+1].W-x[j-1][i].W-x[j-1][i+1].W)/h; } // chk
    if (com==ZX) { e = ((x[j][i+1].W - x[j][i].W)/h + 
			(x[j+1][i].U - x[j][i].U)/h)/2; } // chk
  } 
  return e;
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "PlateSpeed"
PetscReal PlateSpeed(PetscInt i, Parameter *p)
/* ------------------------------------------------------------------- */
{
  PetscReal x, h = pow(p->nj-4,-1);
  x = h*(i-1) - p->xridge/p->height;
  return tanh(2*x/(p->wridge/p->height));
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "StrainRateSecInv"
/* 
   calculate the second invariant of the strain rate 
*/
PetscReal StrainRateSecInv(Field **x, int i, int j, int pos, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  PetscReal exx = 1.0, ezz = 1.0, ezx = 1.0;
  
  exx = StrainRateTensor(x,i,j,XX,pos,user);
  ezz = StrainRateTensor(x,i,j,ZZ,pos,user);
  ezx = StrainRateTensor(x,i,j,ZX,pos,user);
  return sqrt(0.5*(exx*exx + ezz*ezz + ezx*ezx + ezx*ezx));
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "TotalPressure"
/* 
   Calculate the total pressure (lithostatic plus dynamic) in Pascals.
 */
PetscReal TotalPressure(Field **x, int i, int j, int pos, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  Parameter       *p = user->param;
  PetscReal        h, p_l, p_d;

  /* Lithostatic Pressure */
  h = 1.0/(p->nj-1);

  if (pos==GRID_CENTER)      { p_l = p->rho*p->g*PetscMax(h*(j-1.5),0)*p->height*1e3; } 
  else if (pos==GRID_IPLUS)  { p_l = p->rho*p->g*PetscMax(h*(j-1.5),0)*p->height*1e3; } 
  else if (pos==GRID_JPLUS)  { p_l = p->rho*p->g*(PetscMax(h*(j-1.5),0) + PetscMax(h*(j-0.5),0))*0.5*p->height*1e3; }
  else if (pos==GRID_IJPLUS) { p_l = p->rho*p->g*(PetscMax(h*(j-1.5),0) + PetscMax(h*(j-0.5),0))*0.5*p->height*1e3; } 

  /* Dynamic Pressure */

  if (pos==GRID_CENTER)      { p_d = x[j][i].P*p->eta0*p->U0/(1e5*SEC_PER_YR*p->height); } 
  else if (pos==GRID_IPLUS)  { p_d = (x[j][i].P + x[j][i+1].P)*0.5*p->eta0*p->U0/(1e5*SEC_PER_YR*p->height); } 
  else if (pos==GRID_JPLUS)  { p_d = (x[j][i].P + x[j+1][i].P)*0.5*p->eta0*p->U0/(1e5*SEC_PER_YR*p->height); }
  else if (pos==GRID_IJPLUS) { p_d = (x[j][i].P + x[j][i+1].P + x[j+1][i].P + x[j+1][i+1].P)*0.25*p->eta0*p->U0/(1e5*SEC_PER_YR*p->height); } 

  //  p_d = x[j][i].P*p->eta0*p->U0/(1e5*SEC_PER_YR*p->height);    // Pascals

  /* Total pressure (force non-negative) */
  if(p->TotalPressure == 1) {  return PetscMax(p_l+p_d, 0 ); }
  else {  return p_l; }
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Temperature"
/* 
   Calculate interpolated Temperature in Kelvin
*/
PetscReal Temperature(Field **x, int i, int j,int pos, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  Parameter       *p = user->param;
  PetscReal       T = 0;

  if (pos==GRID_CENTER)      { T = x[j][i].T*p->Tp + 273; } 
  else if (pos==GRID_IPLUS)  { T = (x[j][i].T + x[j][i+1].T)*0.5*p->Tp + 273; } 
  else if (pos==GRID_JPLUS)  { T = (x[j][i].T + x[j+1][i].T)*0.5*p->Tp + 273; }
  else if (pos==GRID_IJPLUS) { T = (x[j][i].T + x[j][i+1].T + x[j+1][i].T + x[j+1][i+1].T)*0.25*p->Tp + 273; } 
  return T;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "GrainSizeEvolution"
/* 
   Calculate residual of grain size evolution equation
*/
PetscReal GrainSizeEvolution(Field **x, Field **x0, int i, int j, 
			     PetscReal *terms, AppCtx *user)
/* ------------------------------------------------------------------- */
{
  Parameter       *p = user->param;
  PetscErrorCode  ierr;
  PetscReal       eta, beta, GrowthRate, DecayRate, Advection, Water;
  PetscReal       P=0, T=0, e_2=0, a;
  PetscReal       h = pow(user->param->nj-4,-1);
  PetscReal       v[4], q[9];
  if (p->const_visc) { return x[j][i].A; }          // if user has selected const visc, don't evolve grain size.
  if (p->constant_grain_size) { return x[j][i].A; } // if user has selected constant grain size, don't evolve grain size.

  P   = TotalPressure(x,i,j,GRID_CENTER,user);                // Pascals
  T   = Temperature(x,i,j,GRID_CENTER,user);      // Kelvin
  e_2 = StrainRateSecInv(x0,i,j,GRID_CENTER,user); // dim'less
  a   = p->grain_size*exp(x[j][i].A);              // meters

  if (p->WetComp){
    Water = WaterContent(j,p);
  } else {
    Water = 0;
  }

  ierr = ComputeViscosityModule(T, P, e_2*(p->U0/1e2/SEC_PER_YR)/(p->height*1e3),
				a, Water, &eta, &beta, NULL, &user->VM);CHKERRQ(ierr);
 
  if (beta < 0) { return x[j][i].A; }                    // if user has selected diffusion only, don't evolve grain size.

  eta = pow(p->eta0/eta + p->eta0/p->etamax,-p->contin); // dimensionless eta

  DecayRate  = p->dest_coef*e_2*e_2*eta*beta*exp(x[j][i].A);

  GrowthRate = p->growth_coef
    *    exp(-(p->Q_g + (P*p->V_g))/(p->R*T))   // thermodynamic exponential
    *    exp(-p->grain_expon*x[j][i].A);        // grain size exponential 

  /* Advection term */
  v[GRID_N] = x0[j][i].W; v[GRID_S] = x0[j-1][i].W;
  v[GRID_E] = x0[j][i].U; v[GRID_W] = x0[j][i-1].U;
  q[GRID_N] = x[j+1][i].A; q[GRID_NN] = x[j+2][i].A;
  q[GRID_S] = x[j-1][i].A; q[GRID_SS] = x[j-2][i].A;
  q[GRID_E] = x[j][i+1].A; q[GRID_EE] = x[j][i+2].A;
  q[GRID_W] = x[j][i-1].A; q[GRID_WW] = x[j][i-2].A;
  q[GRID_C] = x[j][i].A;
  Advection = FrommAdvection(v,q,h,h);

  if (terms) {
    terms[0] = GrowthRate;
    terms[1] = DecayRate;
    terms[2] = Advection;
  }

  return GrowthRate - DecayRate - Advection;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "WaterContent"
PetscReal WaterContent(PetscInt j, Parameter *p)
/* ------------------------------------------------------------------- */
{
  PetscReal z, h = pow((p->nj-4),-1);

  z = PetscMax(h*(j-1.5),0)*p->height;     //current depth, km
  
  /*   Use a model for water simlar in functional form to batch melting */
    if (z > p->Z_w){
      return p->COH;
    } else if (z < p->Z_d){
      return 0;
    } else {
      return p->COH*(z - p->Z_d)/(p->Z_w - p->Z_d);
    }

}
