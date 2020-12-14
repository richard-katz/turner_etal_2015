#include <petsc.h>
#include "ViscosityModule.h"

/* -------- PREPROCESSOR MACROS -------- */

#define FNAME_LENGTH  120
#define GRID_CENTER   900
#define GRID_IPLUS    901
#define GRID_JPLUS    902
#define GRID_IJPLUS   903
#define XX            10
#define ZZ            11
#define ZX            12
#define SEC_PER_YR    31557600.0

/* -------- FIELD DATA STRUCTURES -------- */

/* PRIMARY HYPERBOLIC VARIABLES */
typedef struct {
  PetscReal U,W;  /* mantle flow velocity */
  PetscReal P;    /* mantle dynamic pressure */
  PetscReal T;    /* mantle potential temperature */
  PetscReal A;    /* log(a/a_0). Dimensionless grain size */
} Field;

/* AUXILIARY (OUTPUT) VARIABLES */
typedef struct {
  PetscReal eta;             /* mantle viscosity */
  PetscReal e_2;             /* second strain-rate invariant */
  PetscReal beta;            /* percentage of work resulting in dislocation creep */
  PetscReal decay;           /* grain decay rate */             // use with aux tests
  PetscReal growth;          /* grain growth rate */            // use with aux tests
  PetscReal advn;            /* grain rate of chnge due to advection */            // use with aux tests
  PetscReal dislocation;     /* grain decay rate */             // use with aux tests
  PetscReal diffusion;       /* grain growth rate */            // use with aux tests
  PetscReal plastic;         /* grain rate of chnge due to advection */            // use with aux tests
  PetscReal GBS;             /* grain growth rate */            // use with aux tests
  PetscReal Press;           /* Total pressure */
  PetscReal Water;           /* Distribution of water */
} Auxfield;

/* -------- PARAMETER DATA STRUCTURES -------- */

/* PHYSICAL AND NUMERICAL PARAMETERS */
typedef struct {
  PetscInt       ndims, ni, nj;                     // grid params
  PetscReal      width, height, xridge, wridge;     // domain params
  PetscReal      U0, Tp, K, Pe, eta0, etamax, g;    // physical params
  PetscReal      FinalTolerance;                    // Solver parameter
  PetscReal      rho;
  PetscReal      contin;                            // continuation parameter
  PetscBool      test, output_to_file, const_visc, TotalPressure;
  PetscBool      lag_visc, do_unlagged_final_solve, constant_grain_size;
  PetscBool      WetComp, LinearWater;                           // switch between wet and dry composition
  char           output_filename[FNAME_LENGTH];
  /* Grain size parameters */
  PetscReal      lambda, grain_size, geometrical_coef;  // grain parameters
  PetscReal      grain_expon, gamma, BB_grain_size;     // grain parameters
  PetscReal      K_g, R, Q_g, V_g;                      // grain thermodynamic parameters
  PetscReal      dest_coef, growth_coef;                // dimensionless parameters for grain size change
  PetscReal      COH, Z_d, Z_w;                                   // Water profile parameters
} Parameter;

/* HIGHEST-LEVEL STRUCTURE FOR PARAMETERS & PROBLEM DATA */
typedef struct {
  PetscBag          bag;
  Parameter         *param;
  ViscosityModule   VM;
  MPI_Comm          comm;
  SNES              snes;
  DM                da, da_aux;
  Vec               X, Xo, Xps, R, aux;
} AppCtx;


PetscErrorCode SetParams(AppCtx*);
PetscErrorCode ContinuationSolve(AppCtx*);
PetscErrorCode CreateDataStructures(AppCtx*);
PetscErrorCode DestroyDataStructures(AppCtx*);
PetscErrorCode CreateInitialGuess(AppCtx*);
PetscErrorCode DoOutput(AppCtx*);
PetscErrorCode UpdateAuxilliaryFields(AppCtx*);
PetscErrorCode FormResidual(SNES, Vec, Vec, void*);
PetscErrorCode CheckSNESResidualNoLag(SNES, Vec, PetscReal*, AppCtx*);
PetscErrorCode SolveWithLaggedViscosity(PetscReal, PetscInt, AppCtx*);
PetscErrorCode ContinuationCorrection(PetscReal, PetscReal, AppCtx*);

PetscReal CornerFlowU(PetscReal, PetscReal);
PetscReal CornerFlowW(PetscReal, PetscReal);
PetscReal CornerFlowP(PetscReal, PetscReal);
PetscReal HalfSpaceCooling(PetscReal, PetscReal, Parameter*);
PetscReal StrainRateSecInv(Field**, int, int, int, AppCtx*);
PetscReal GrainSizeEvolution(Field**, Field**, PetscInt, PetscInt, PetscReal*, AppCtx*); 
PetscReal TotalPressure(Field**, PetscInt, PetscInt,PetscInt, AppCtx*);
PetscReal WaterContent(PetscInt, Parameter*);

PetscErrorCode DMDAGetGhostedArray(DM, Vec, Vec*, void*);
PetscErrorCode DMDARestoreGhostedArray(DM, Vec, Vec*, void*);
PetscErrorCode DMDAGetGridInfo(DM, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);

PetscBool OptionsHasName(const char*);
