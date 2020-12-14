#include "mor.h"
PetscReal DiffusionCreepViscosity(PetscReal T, PetscReal PR, PetscReal grain_size, PetscReal COH, ViscosityParameters*);
PetscReal DislocationCreepViscosity(PetscReal T, PetscReal PR, PetscReal e_2, PetscReal COH, ViscosityParameters*);
PetscReal GrainBoundarySlidingViscosity(PetscReal T, PetscReal PR, PetscReal grain_size, PetscReal e_2, PetscReal COH, ViscosityParameters*);
PetscReal Plasticity(PetscReal PR, PetscReal e_2, ViscosityParameters*);
#define REG_REAL(A,B,C,D,E)   ierr = PetscBagRegisterReal(A,B,C,D,E);CHKERRQ(ierr)
#define REG_INTG(A,B,C,D,E)   ierr = PetscBagRegisterInt(A,B,C,D,E);CHKERRQ(ierr)

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "InitializeViscosityModule"
PetscErrorCode InitializeViscosityModule(MPI_Comm comm, ViscosityModule *VM)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode       ierr;
  PetscBag             *bag = &VM->bag;
  ViscosityParameters  *p = VM->param;
  PetscFunctionBegin;
  VM->comm = comm;
  /* Create bag */
  ierr = PetscBagCreate(comm,sizeof(ViscosityParameters),bag);CHKERRQ(ierr);
  ierr = PetscBagSetName(*bag,"par_vm","--- params for diffusion and dislocation viscosity");CHKERRQ(ierr);
  ierr = PetscBagSetOptionsPrefix(*bag,"vm_");CHKERRQ(ierr);
  ierr = PetscBagGetData(*bag,(void**)&VM->param);CHKERRQ(ierr); p = VM->param;

  /* PARAMS */
  /* Values from Q. Wang, Lithos 120 (2010) 30-41 - A review of water contents and ductile deformation mechanisms 
     of olivine: Implications for the lithosphere-asthenosphere boundary of continents. */
  REG_REAL(*bag,&p->R,8.314462      ,"R","J/mol-K, Universal gas constant");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Qd,3.75e5       ,"Qd","J/mol, Activation Energy, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Ad,1.5e9        ,"Ad","(mu m)^3 Sec^{-1} (MPa)^{-1}, Prefactor, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Vd,4e-6         ,"Vd","m^3/mol, Activation volume, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Ql,5.3e5        ,"Ql","J/mol, Activation Energy, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Al,1.1e5        ,"Al","sec^{-1} (MPa)^{-n}, Prefactor, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Vl,1.6e-5       ,"Vl","m^3/mol, Activation volume, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->ndisl,3.5       ,"ndisl","Non-Newtonian exponent, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->mdisl,0         ,"mdisl","grain exponent, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->mdiff,3         ,"mdiff","grain exponent, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->ndiff,1         ,"ndiff","Non-Newtonian exponent, dffusion creep");CHKERRQ(ierr);

  /* wet pararmeters */
  /* Below are essential for wet composition */
  REG_REAL(*bag,&p->Ad_w,2.5e7      ,"Ad_w" ,"(mu m)^3 Sec^{-1} (MPa)^{-1}, Prefactor, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->rdiff,1         ,"rdiff","Water exponent, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Al_w,1.6e3      ,"Al_w" ,"sec^{-1} (MPa)^{-n}, Prefactor, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->rdisl,1.2       ,"rdisl","water exponent, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->rgbs,0          ,"rgbs" ,"water exponent, grain boundary sliding");CHKERRQ(ierr);

  /* Below are optional parameters than may be changed for wet composition */
  /* By default these parameters are the same as for dry */
  REG_REAL(*bag,&p->Qd_wet,3.75e5       ,"Qd_wet","J/mol, Activation Energy, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Vd_wet,4e-6         ,"Vd_wet","m^3/mol, Activation volume, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Ql_wet,5.3e5        ,"Ql_wet","J/mol, Activation Energy, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Vl_wet,1.6e-5       ,"Vl_wet","m^3/mol, Activation volume, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->ndisl_wet,3.5       ,"ndisl_wet","Non-Newtonian exponent, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->mdisl_wet,0         ,"mdisl_wet","grain exponent, dislocation creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->mdiff_wet,3         ,"mdiff_wet","grain exponent, diffusion creep");CHKERRQ(ierr);
  REG_REAL(*bag,&p->ndiff_wet,1         ,"ndiff_wet","Non-Newtonian exponent, dffusion creep");CHKERRQ(ierr);

  /*GBS parameters */

  REG_REAL(*bag,&p->Qgbs,4.4e5      ,"Qgbs","J/mol, Activation Energy, Grain boundary sliding");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Agbs,6.3096e4   ,"Agbs","(mu m)^3 Sec^{-1} (MPa)^{-1}, Prefactor, Grain boundary sliding");CHKERRQ(ierr);
  REG_REAL(*bag,&p->Vgbs,1.6e-5     ,"Vgbs","m^3/mol, Activation volume, Grain boundary sliding");CHKERRQ(ierr);
  REG_REAL(*bag,&p->ngbs,2.9        ,"ngbs","Non-Newtonian exponent, Grain boundary sliding");CHKERRQ(ierr);
  REG_REAL(*bag,&p->mgbs,0.7        ,"mgbs","grain exponent, Grain boundary sliding");CHKERRQ(ierr);

  /*Plastic stress limiter parameters */

  REG_REAL(*bag,&p->cohesion,1e8    ,"cohesion","Cohesion for Drucker-Prager plasticity, Pa");CHKERRQ(ierr);
  REG_REAL(*bag,&p->fric_angle,30   ,"fric_angle","Friction angle for Drucker-Prager plasticity, degrees");CHKERRQ(ierr);
  
  /*Rheology combination */

  REG_INTG(*bag,&p->ViscType, 0     ,"ViscType","Control the type of viscosity. Number from 0-8."); CHKERRQ(ierr);

  /* report parameter values */
  PetscPrintf(comm,"-----------------------------------------------------------------\n");
  ierr = PetscBagView(*bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscPrintf(comm,"-----------------------------------------------------------------\n");
  PetscFunctionReturn(0);
}


/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "FinalizeViscosityModule"
PetscErrorCode FinalizeViscosityModule(ViscosityModule *VM)
/*---------------------------------------------------------------------*/
{
  PetscErrorCode  ierr;
  PetscFunctionBegin;
  ierr = PetscBagDestroy(&VM->bag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "ComputeViscosityModule"
/*
   Input:  
           Temperature  [Kelvin]
	   PRessure     [Pa]

   Return: Viscosity    [Pa-sec]
           Beta         [dimensionless]
	   dominant determines the relative contribution to eta from the physical process available
*/
PetscErrorCode ComputeViscosityModule(PetscReal T, PetscReal PR, PetscReal e_2,
				      PetscReal grain_size, PetscReal COH, PetscReal *eta, 
				      PetscReal *beta, PetscReal *dominant, ViscosityModule *VM)
/*---------------------------------------------------------------------*/
{
  ViscosityParameters  *p = VM->param;
  PetscReal            difn = 0.0, disl = 0.0 , plst = 0.0, gbs = 0.0;
  PetscFunctionBegin;

  if (p->ViscType == 0){
    /* DIFFUSION + DISLOCATION + PLASTICITY + GRAIN-BOUNDARY SLIDING */
    difn = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
    disl = DislocationCreepViscosity(T,PR,e_2,COH,p);
    plst = Plasticity(PR,e_2,p);
    gbs  = GrainBoundarySlidingViscosity(T,PR,grain_size,e_2,COH,p);
      *eta  = pow(1/difn + 1/disl + 1/gbs + 1/plst, -1);
      *beta = (*eta)/pow(1/disl + 1/gbs,-1);
     
      
    if (dominant) {
      dominant[0] = difn;
      dominant[1] = plst;  
      dominant[2] = disl;
      dominant[3] = gbs;  
    }

  } else if (p->ViscType == 1){
    /* DIFFUSION + DISLOCATION + GRAIN-BOUNDARY SLIDING */
      difn = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
      disl = DislocationCreepViscosity(T,PR,e_2,COH,p);
      gbs  = GrainBoundarySlidingViscosity(T,PR,grain_size,e_2,COH,p);
      *eta  = pow(1/difn + 1/disl + 1/gbs, -1);
      *beta = (*eta)/pow(1/disl + 1/gbs,-1);
     
      
    if (dominant) {
      dominant[0] = difn;
      dominant[1] =0;  
      dominant[2] = disl;
      dominant[3] = gbs;  
    }

  } else if (p->ViscType == 2){
    /* DIFFUSION + DISLOCATION + PLASTICITY */
      difn = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
      disl = DislocationCreepViscosity(T,PR,e_2,COH,p);
      plst = Plasticity(PR,e_2,p);
      *eta  = pow(1/difn + 1/disl + 1/plst, -1);
      *beta = (*eta)/disl;

    if (dominant) {
      dominant[0] = difn;
      dominant[1] = plst;
      dominant[2] = disl;
      dominant[3] = 0;  
    }

  } else if (p->ViscType == 3) {
    /* DIFFUSION + DISLOCATION */
      difn = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
      disl = DislocationCreepViscosity(T,PR,e_2,COH,p);
      *eta  = pow(1/difn + 1/disl, -1);
      *beta = (*eta)/disl;

    if (dominant) {
      dominant[0] = difn;
      dominant[1] = 0;
      dominant[2] = disl;
      dominant[3] = 0;  
    }

  } else if (p->ViscType == 4) {
    /* DIFFUSION + GBS */
      difn = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
      gbs  = GrainBoundarySlidingViscosity(T,PR,grain_size,e_2,COH,p);
      *eta  = pow(1/difn + 1/gbs, -1);
      *beta = (*eta)/gbs;

    if (dominant) {
      dominant[0] = difn;
      dominant[1] = 0;
      dominant[2] = 0;
      dominant[3] = gbs;  
    }

  } else if (p->ViscType == 5){
    /* DIFFUSION */
      *eta = DiffusionCreepViscosity(T,PR,grain_size,COH,p);
      *beta = -1;

    if (dominant) {
      dominant[0] = *eta;
      dominant[1] = 0;
      dominant[2] = 0;
      dominant[3] = 0;  
    }


  } else if (p->ViscType == 6){
    /* DISLOCATION */
      *eta = DislocationCreepViscosity(T,PR,e_2,COH,p);
      *beta = 1;

    if (dominant) {
      dominant[0] = 0;
      dominant[1] = 0;
      dominant[2] = *eta;
      dominant[3] = 0;  
    } 

  } else if (p->ViscType == 7){
    /* Grain-boundary Sliding */
      *eta  = GrainBoundarySlidingViscosity(T,PR,grain_size,e_2,COH,p);
      *beta = 1;

    if (dominant) {
      dominant[0] = 0;
      dominant[1] = 0;
      dominant[2] = 0;
      dominant[3] = *eta;  
    }

  } else {

  *eta = 1e20;
  *beta = 0;

  }
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------
  Private functions 
  ---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DiffusionCreepViscosity"
/* 
   Input:  
           Temperature  [Kelvin]
	   PRessure     [Pa]
	   grain_size   [m]

   Return: Viscosity    [Pa-sec]
*/
PetscReal DiffusionCreepViscosity(PetscReal T, PetscReal PR, PetscReal a,  
				  PetscReal COH, ViscosityParameters *p)
/*---------------------------------------------------------------------*/
{
  PetscReal A0, A_prime, A_prime_wet, eta_wet, eta_dry, A0_w;

  if (COH <=30){

    /*return dry diffusion creep */

  A_prime = (pow(3,(p->ndiff+1)/2)/pow(2,1-p->ndiff))*pow(1e-6, p->mdiff + p->ndiff)*p->Ad;
  A0      = 0.5*pow(1/A_prime,1/p->ndiff);

  return A0*pow(a,p->mdiff)*exp((p->Qd + PR*p->Vd)/(p->R*T));

  } else {

    /*return the minimum viscosity from dry and wet creep */

  A_prime = (pow(3,(p->ndiff+1)/2)/pow(2,1-p->ndiff))*pow(1e-6, p->mdiff + p->ndiff)*p->Ad;
  A0      = 0.5*pow(1/A_prime,1/p->ndiff);

  eta_dry = A0*pow(a,p->mdiff)*exp((p->Qd + PR*p->Vd)/(p->R*T));

  A_prime_wet = (pow(3,(p->ndiff_wet+1)/2)/pow(2,1-p->ndiff_wet))*pow(1e-6, p->mdiff_wet + p->ndiff_wet)*p->Ad_w;
  A0_w      = 0.5*pow(1/A_prime_wet,1/p->ndiff_wet);

  eta_wet = A0_w*pow(a,p->mdiff_wet)*pow(COH,-1.0*p->rdiff/p->ndiff_wet)*exp((p->Qd_wet + PR*p->Vd_wet)/(p->R*T));

  return PetscMin(eta_dry, eta_wet );
  }
}



/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "GrainBoundarySlidingViscosity"
/* 
   Input:  
           Temperature  [Kelvin]
	   PRessure     [Pa]
	   grain_size   [m]

   Return: Viscosity    [Pa-sec]
*/
 PetscReal GrainBoundarySlidingViscosity(PetscReal T, PetscReal PR, PetscReal a, PetscReal e_2,  
					 PetscReal COH, ViscosityParameters *p)
/*---------------------------------------------------------------------*/
{
  PetscReal A0, E, A_prime;

  A_prime = (pow(3,(p->ngbs+1)/2)/pow(2,1-p->ngbs))*pow(1e-6, p->mgbs + p->ngbs)*p->Agbs;
  A0 = 0.5*pow(1/A_prime,1/p->ngbs);
  E = pow(e_2,(1-p->ngbs)/p->ngbs);

  return A0*pow(a,p->mgbs/p->ngbs)*E*exp((p->Qgbs + PR*p->Vgbs)/(p->R*T*p->ngbs));
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "DislocationCreepViscosity"
/* 
   Input:  
           Temperature  [Kelvin]
	   PRessure     [Pa]

   Return: Viscosity    [Pa-sec]
*/
PetscReal DislocationCreepViscosity(PetscReal T, PetscReal PR, PetscReal e_2, 
				    PetscReal COH, ViscosityParameters *p)
/*---------------------------------------------------------------------*/
{
  PetscReal E = 1.0, A0, A_prime, A_prime_wet, eta_wet, eta_dry, A0_wet;
  PetscReal E_wet;

  if (COH <=30)  {

    /*return dry diffusion creep */

  A_prime = (pow(3,(p->ndisl+1)/2)/pow(2,1-p->ndisl))*pow(1e-6, p->mdisl + p->ndisl)*p->Al;
  A0      = 0.5*pow(1/A_prime,1/p->ndisl);
  E       = pow(e_2,(1-p->ndisl)/p->ndisl);

  return E*A0*exp((p->Ql + PR*p->Vl)/(p->R*T*p->ndisl));

  } else {

    /*return the minimum viscosity from dry and wet creep */

  A_prime = (pow(3,(p->ndisl+1)/2)/pow(2,1-p->ndisl))*pow(1e-6, p->mdisl + p->ndisl)*p->Al;
  A0      = 0.5*pow(1/A_prime,1/p->ndisl);
  E       = pow(e_2,(1-p->ndisl)/p->ndisl);

  eta_dry = E*A0*exp((p->Ql + PR*p->Vl)/(p->R*T*p->ndisl));

  A_prime_wet = (pow(3,(p->ndisl_wet+1)/2)/pow(2,1-p->ndisl_wet))*pow(1e-6, p->mdisl_wet + p->ndisl_wet)*p->Al_w;
  A0_wet      = 0.5*pow(1/A_prime_wet,1/p->ndisl_wet);
  E_wet       = pow(e_2,(1-p->ndisl_wet)/p->ndisl_wet);
  eta_wet     = E_wet*A0_wet*pow(COH,-1.0*p->rdisl/p->ndisl_wet)*exp((p->Ql_wet + PR*p->Vl_wet)/(p->R*T*p->ndisl_wet));

  return PetscMin(eta_dry, eta_wet );
  }

}


/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "Plasticity"
/* 
   Input:  
           PRessure     [Pa]

   Return: Viscosity    [Pa-sec]
*/
PetscReal Plasticity(PetscReal PR, PetscReal e_2, ViscosityParameters *p)
/*---------------------------------------------------------------------*/
{
  return (p->cohesion*cos(p->fric_angle*PETSC_PI/180) 
  	  + PR*sin(p->fric_angle*PETSC_PI/180))/e_2; 
}
