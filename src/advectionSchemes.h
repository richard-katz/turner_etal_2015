#ifndef __ADVECTIONSCHEMES_H
#define __ADVECTIONSCHEMES_H

#ifndef __PETSCSYS_H
#include "petscsys.h"
#endif

/* PUBLIC */
#define GRID_C        4
#define GRID_N        0
#define GRID_S        1
#define GRID_E        2
#define GRID_W        3
#define GRID_NN       5
#define GRID_SS       6
#define GRID_EE       7
#define GRID_WW       8

PetscReal TVDAdvection    (PetscReal*, PetscReal*, PetscReal, PetscReal);
PetscReal FrommAdvection  (PetscReal*, PetscReal*, PetscReal, PetscReal);
PetscReal FVAdvection     (PetscReal*, PetscReal*, PetscReal, PetscReal);

#endif /*__ADVECTIONSCHEMES_H*/
