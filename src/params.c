#include "mor.h"
#define REG_REAL(A,B,C,D,E)   ierr=PetscBagRegisterReal(A,B,C,D,E);CHKERRQ(ierr)
#define REG_INTG(A,B,C,D,E)   ierr=PetscBagRegisterInt(A,B,C,D,E);CHKERRQ(ierr)
#define REG_IARR(A,B,C,D,E)   ierr=PetscBagRegisterIntArray(A,B,C,D,E);CHKERRQ(ierr)
#define REG_TRUE(A,B,C,D,E)   ierr=PetscBagRegisterBool(A,B,C,D,E);CHKERRQ(ierr)
#define REG_STRG(A,B,C,D,E,F) ierr=PetscBagRegisterString(A,B,C,D,E,F);CHKERRQ(ierr)
#define REG_ENUM(A,B,C,D,E,F) ierr=PetscBagRegisterEnum(A,B,C,D,E,F);CHKERRQ(ierr)
PetscErrorCode ReportParams(AppCtx*);

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "SetParams"
PetscErrorCode SetParams(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  PetscBag       bag;
  Parameter      *p;
  PetscErrorCode ierr;
  PetscReal      P_ref;           // Reference pressure
  PetscReal      e_2_ref = 1e-15; // Reference strain-rate, sec^-1
  PetscReal      beta_ref;        // Reference beta for grain size deformation, unitless.
  PetscFunctionBegin;

  /* create bag */
  ierr = PetscBagCreate(PETSC_COMM_WORLD,sizeof(Parameter),&user->bag); CHKERRQ(ierr);
  ierr = PetscBagGetData(user->bag,(void**)&user->param);CHKERRQ(ierr);
  bag = user->bag; p = user->param; 
  ierr = PetscBagSetName(bag,"par","--- params for mor simulation");CHKERRQ(ierr);

  /* grid parameters */
  REG_INTG(bag,&p->ndims,2   ,"ndims","[GRID] Grid dimensionality");
  REG_INTG(bag,&p->ni,124     ,"ni","[GRID] Grid points in x-direction, incl buffer points");
  REG_INTG(bag,&p->nj,84      ,"nj","[GRID] Grid points in z-direction, incl buffer points");

  /* domain parameters */
  REG_REAL(bag,&p->height,80      ,"height","[DOMAIN] Height of the box, km");
  p->width = p->height*(PetscReal)(p->ni-4)/(PetscReal)(p->nj-4);
  REG_REAL(bag,&p->width ,p->width,"width","[DOMAIN] <DO NOT SET> Width of the box, km");
  REG_REAL(bag,&p->xridge,0        ,"xridge","[DOMAIN] x-position of ridge, km");
  REG_REAL(bag,&p->wridge,4        ,"wridge","[DOMAIN] half-width of ridge spreading zone, km");

  /* Grain evolution parameters - subset of physical parameters */
  REG_REAL(bag,&p->lambda, 1             ,"lambda","[PHYSICAL] Fraction of dislocation work to grain-size reduction, dimensionless");
  REG_REAL(bag,&p->grain_size,9.5e-3     ,"grain_size","[PHYSICAL] Default constant grain-size, m");
  REG_REAL(bag,&p->geometrical_coef,3    ,"geometrical_coef","[PHYSICAL] Grain size geometric factor, dimensionless");
  REG_REAL(bag,&p->gamma,1               ,"gamma","[PHYSICAL] surface energy at grain-grain contacts, Kg/s^2");
  REG_REAL(bag,&p->grain_expon,3         ,"grain_expon","[PHYSICAL] grain size growth exponent, dimensionless");
  REG_REAL(bag,&p->K_g,1e-5              ,"K_g","[PHYSICAL] Grain growth prefactor, m^p/s");
  REG_REAL(bag,&p->R,8.314462            ,"R","[PHYSICAL] Universal gas constant, J/mol-K");
  REG_REAL(bag,&p->Q_g,3.5e5             ,"Q_g","[PHYSICAL] grain growth activation Energy, J/mol");
  REG_REAL(bag,&p->V_g,8e-6              ,"V_g","[PHYSICAL] Grain growth activation volume, m^3/mol");
  REG_REAL(bag,&p->BB_grain_size,2e-2   ,"BB_grain_size","[PHYSICAL] Inflow grain-size at the bottom boundary, m");

  /* physical parameters */
  REG_TRUE(bag,&p->const_visc,PETSC_FALSE,"constant_viscosity","[PHYSICAL] isoviscous model, true/false");
  REG_REAL(bag,&p->rho,3300         ,"rho","[PHYSICAL] mantle density, kg/m^3");
  REG_REAL(bag,&p->g,9.8            ,"g","[PHYSICAL] accelleration of gravity, m^2/sec");
  REG_REAL(bag,&p->U0,4             ,"U0","[PHYSICAL] half-spreading rate, cm/yr");
  REG_REAL(bag,&p->Tp,1350          ,"Tp","[PHYSICAL] potential temperature, deg C");
  REG_REAL(bag,&p->K,1e-6           ,"kappa","[PHYSICAL] thermal diffusivity, m^2/sec");
  P_ref = p->rho*p->g*10e3;  // 10 km depth, Pascals
  ierr = ComputeViscosityModule(p->Tp+273,P_ref,e_2_ref,p->grain_size,0,&p->eta0,&beta_ref,NULL,&user->VM);CHKERRQ(ierr); //Dry
  REG_REAL(bag,&p->eta0,p->eta0,"eta0","[PHYSICAL] reference viscosity, Pa-sec");
  REG_REAL(bag,&p->etamax,1e24 ,"etamax","[PHYSICAL] maximum viscosity, Pa-sec");
 
  /*water concerntration parameters */
  REG_REAL(bag,&p->COH,0 ,"COH","[PHYSICAL] max water content, H/10^6 Si");
  REG_REAL(bag,&p->Z_d,57 ,"Z_d","[PHYSICAL] Depth at which mantle is dry, km");
  REG_REAL(bag,&p->Z_w,160 ,"Z_w","[PHYSICAL] Depth of wet solidus, km");

  /* nondimensional parameters */
  p->Pe = p->U0/100/SEC_PER_YR /* m/sec */
    *     p->height*1000       /* m */
    /     p->K;                /* m^2/sec */
  REG_REAL(bag,&p->Pe,p->Pe,"Peclet","[NONDIMENSIONAL] <DO NOT SET> Peclet number");

  p->dest_coef = p->lambda              /* dimensionless */
    *            p->eta0                /* Pa sec */
    *            (p->U0/100/SEC_PER_YR) /* m/sec */
    *            p->grain_size          /* m */
    /            p->geometrical_coef    /* dimensionless */
    /            p->gamma               /* kg/sec^2      */
    /            (p->height*1000) ;     /* m */
  REG_REAL(bag,&p->dest_coef,p->dest_coef,"dest_coef","[NONDIMENSIONAL] <DO NOT SET> Grain destruction coef");

  p->growth_coef = (p->K_g*(p->grain_expon/3)*pow(p->grain_size,p->grain_expon - 3))         /* m^p/s, adjusted for grain_expon */
    *              (p->height*1000)                                                          /* m */
    /              (p->U0/100/SEC_PER_YR)                                                    /* m/sec */
    /              pow(p->grain_size, p->grain_expon)                                        /* m^p */
    /              p->grain_expon;                                                           /*dimensionless */
  REG_REAL(bag,&p->growth_coef,p->growth_coef,"growth_coef","[NONDIMENSIONAL] <DO NOT SET> Grain growth coef");

  /* solver params */
  REG_TRUE(bag,&p->constant_grain_size,PETSC_FALSE,"constant_grain_size","[SOLVER} Use a constant grain-size, true/false");
  REG_TRUE(bag,&p->do_unlagged_final_solve,PETSC_FALSE,"unlagged_final_solve","[SOLVER] Do a final solve with unlagged viscosity, true/false");
  REG_TRUE(bag,&p->lag_visc,PETSC_FALSE,"lag_viscosity","[SOLVER] Compute residual (and do solve) with viscosity from previous solve <DO NOT SET>");
  REG_TRUE(bag,&p->TotalPressure,PETSC_TRUE,"TotalPressure","[SOLVER] Compute with total presseure i.e. lithostatic and dynamic");
  REG_TRUE(bag,&p->WetComp,PETSC_FALSE,"WetComp","[SOLVER] Wet or dry composition, false/true");
  REG_TRUE(bag,&p->LinearWater,PETSC_TRUE,"LinearWater","[SOLVER] Do you want a piecewise linear water profile, true/false");
  REG_REAL(bag,&p->FinalTolerance,7.5e-4,"FinalTolerance","[SOLVER] Final Tolerance");

  /* output params */
  REG_STRG(bag,&p->output_filename,FNAME_LENGTH,"null","output_file","[I/O] Name base for output files, set with: -output_file <filename>");
  REG_TRUE(bag,&p->output_to_file,PETSC_FALSE,"do_output","[I/O] <DO NOT SET> Flag will be true if you specify an output file name");
  p->output_to_file  = (PetscBool)(strcmp(p->output_filename,"null")!=0);
  REG_TRUE(bag,&p->test,PETSC_FALSE,"test","[I/O] Do parameter i/o only (T/F)");

  ierr = ReportParams(user); CHKERRABORT(PETSC_COMM_WORLD,ierr);
  PetscFunctionReturn(0);
}

/*---------------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "ReportParams"
PetscErrorCode ReportParams(AppCtx *user)
/*---------------------------------------------------------------------*/
{
  Parameter      *p = user->param;
  PetscInt       cs;
  char           date[30], *options_string;
  MPI_Comm       comm = PETSC_COMM_WORLD;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = MPI_Comm_size(user->comm,&cs);
  ierr = PetscGetDate(date,30);CHKERRQ(ierr);
  ierr = PetscOptionsGetAll(&options_string);CHKERRQ(ierr);

  PetscPrintf(comm,"-------------------mor: %s-----------------\n",&(date[0]));
  PetscPrintf(comm,"Number of MPI processes: %d\n",cs);
  PetscPrintf(comm,"Command line: %s\n",options_string);
  PetscPrintf(comm,"-----------------------------------------------------------------\n");
  ierr = PetscBagView(user->bag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscPrintf(comm,"-------------------------END PARAM REPORT------------------------\n");
  ierr = PetscFree(options_string); CHKERRQ(ierr);
  if ( p->test ) { PetscEnd(); }
  PetscFunctionReturn(0);
}
