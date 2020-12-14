#if !defined(__VISCOSITYMODULE_H)
#define __VISCOSITYMODULE_H

#if !defined(__PETSCBAG_H)
#include <petscbag.h>
#endif

#if !defined(__PETSCVIEWER_H)
#include <petscviewer.h>
#endif

typedef struct _VisMod_par {
  /* PARAMETERS */
  PetscReal  R;          /* units, Gas constant */
  PetscReal  Qd;         /* J/mol, Activiation energy, diffusion creep */
  PetscReal  Ad;         /* Pa-sec, Prefactor, diffusion creep */
  PetscReal  Vd;         /* m^3/mol, Activation volume, diffusion creep */
  PetscReal  Ql;         /* J/mol, Activiation energy, dislocation creep */
  PetscReal  Al;         /* Pa-sec, Prefactor, dislocation creep */
  PetscReal  Vl;         /* m^3/mol, Activation volume, dislocation creep */
  PetscReal  mdisl;      /* grain power-law, dslocation creep */
  PetscReal  ndisl;      /* non-Newtonain power, dislocation creep */
  PetscReal  mdiff;      /* grain power-law, diffusion creep */
  PetscReal  ndiff;      /* non-Newtonain power, diffusion creep */
  PetscReal  Qgbs;       /* J/mol, Activiation energy, Grain-boundary sliding */
  PetscReal  Agbs;       /* Pa-sec, Prefactor,  Grain-boundary sliding */
  PetscReal  Vgbs;       /* m^3/mol, Activation volume,  Grain-boundary sliding */
  PetscReal  ngbs;       /* non-Newtonain power,  Grain-boundary sliding */
  PetscReal  mgbs;       /* grain power-law,  Grain-boundary sliding */
  /* Wet parameters */
  PetscReal  Qd_wet;         /* J/mol, Activiation energy, diffusion creep */
  PetscReal  Ad_wet;         /* Pa-sec, Prefactor, diffusion creep */
  PetscReal  Vd_wet;         /* m^3/mol, Activation volume, diffusion creep */
  PetscReal  Ql_wet;         /* J/mol, Activiation energy, dislocation creep */
  PetscReal  Al_wet;         /* Pa-sec, Prefactor, dislocation creep */
  PetscReal  Vl_wet;         /* m^3/mol, Activation volume, dislocation creep */
  PetscReal  mdisl_wet;      /* grain power-law, dslocation creep */
  PetscReal  ndisl_wet;      /* non-Newtonain power, dislocation creep */
  PetscReal  mdiff_wet;      /* grain power-law, diffusion creep */
  PetscReal  ndiff_wet;      /* non-Newtonain power, diffusion creep */
  PetscReal  Ad_w;       /* Wet prefactor, diffusion creep */
  PetscReal  rdiff;      /* water exponent, diffusion creep */
  PetscReal  Al_w;       /* Wet prefactor, dislocation creep */
  PetscReal  rdisl;      /* water exponent, dislocation creep */
  PetscReal  rgbs;       /* water exponent, gbs creep */
  /* **** */
  PetscReal  cohesion;   /* Pa, plasticity parameter */
  PetscReal  fric_angle; /* degrees, plasticity parameter */
  PetscInt   ViscType;  /* Type of visc: [0] diff + disl + plasticity [1] disl + diff [2] diff [3] disl */
  //  PetscInt   GBSType;   /* Type of Grain-boundary sliding: [0] harmonic with dislocation [1] serial sum */
} ViscosityParameters;

typedef struct _VisMod {
  MPI_Comm              comm;
  PetscBag              bag;
  ViscosityParameters   *param;
} ViscosityModule;

PetscErrorCode InitializeViscosityModule(MPI_Comm, ViscosityModule*);
PetscErrorCode FinalizeViscosityModule(ViscosityModule*);
PetscErrorCode ComputeViscosityModule(PetscReal, PetscReal, PetscReal,PetscReal, PetscReal, PetscReal*, PetscReal*, PetscReal*, ViscosityModule*);

#endif /*__VISTEMP_H*/
