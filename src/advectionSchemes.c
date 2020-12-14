#include "advectionSchemes.h"

#define SUPERBEE(R) PetscMax(0,PetscMax(PetscMin(1,2*(R)),PetscMin(2,(R))))
#define MINMOD(R)   PetscMax(0,PetscMin(1,(R)))
#define MC(R)       PetscMax(0,PetscMin(0.5*PetscMin((1+(R)),2),2*(R)))
#define TVDR(A,B,C,D,E,F) ( ((A)-(B)) / ((C)-(D)) * (E) / (F) )
#define Np  0
#define Nm  1
#define Sp  2
#define Sm  3
#define Ep  4
#define Em  5
#define Wp  6
#define Wm  7
#define NNm 8
#define SSp 9
#define EEm 10
#define WWp 11

/*------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "TVDAdvection"
PetscReal TVDAdvection(PetscReal *v, PetscReal *q, PetscReal dx, 
		       PetscReal dz)
/*------------------------------------------------------------*/
{
  PetscReal AV[12],ax,bx,cx,az,bz,cz;
  PetscReal REp,REm,RWp,RWm,RNp,RNm,RSp,RSm;
  PetscReal LEp,LEm,LWp,LWm,LNp,LNm,LSp,LSm;

  AV[Ep] = 0.5*(v[GRID_E] + fabs(v[GRID_E]));  AV[Em] = 0.5*(v[GRID_E] - fabs(v[GRID_E]));
  AV[Wp] = 0.5*(v[GRID_W] + fabs(v[GRID_W]));  AV[Wm] = 0.5*(v[GRID_W] - fabs(v[GRID_W]));
  AV[WWp]= 0.5*(v[GRID_WW]+ fabs(v[GRID_WW])); AV[EEm]= 0.5*(v[GRID_EE]- fabs(v[GRID_EE]));

  REp = TVDR(q[GRID_C] ,q[GRID_W] ,q[GRID_E],q[GRID_C],AV[Wp] ,AV[Ep]);
  REm = TVDR(q[GRID_EE],q[GRID_E] ,q[GRID_E],q[GRID_C],AV[EEm],AV[Em]); 
  RWp = TVDR(q[GRID_W] ,q[GRID_WW],q[GRID_C],q[GRID_W],AV[WWp],AV[Wp]);
  RWm = TVDR(q[GRID_E] ,q[GRID_C] ,q[GRID_C],q[GRID_W],AV[Em] ,AV[Wm]);

  LEp = SUPERBEE(REp);
  LEm = SUPERBEE(REm);
  LWp = SUPERBEE(RWp);
  LWm = SUPERBEE(RWm);

  ax = - AV[Wp]*(1 - LWp*0.5) - AV[Wm]*LWm*0.5;
  bx =   AV[Ep]*(1 - LEp*0.5) + AV[Em]*LEm*0.5
       - AV[Wm]*(1 - LWm*0.5) - AV[Wp]*LWp*0.5;
  cx =   AV[Em]*(1 - LEm*0.5) + AV[Ep]*LEp*0.5;

  AV[Np] = 0.5*(v[GRID_N] + fabs(v[GRID_N]));  AV[Nm] = 0.5*(v[GRID_N] - fabs(v[GRID_N]));
  AV[Sp] = 0.5*(v[GRID_S] + fabs(v[GRID_S]));  AV[Sm] = 0.5*(v[GRID_S] - fabs(v[GRID_S]));
  AV[SSp]= 0.5*(v[GRID_SS]+ fabs(v[GRID_SS])); AV[NNm]= 0.5*(v[GRID_NN]- fabs(v[GRID_NN]));

  RNp = TVDR(q[GRID_C] ,q[GRID_S] ,q[GRID_N],q[GRID_C],AV[Sp] ,AV[Np]);
  RNm = TVDR(q[GRID_NN],q[GRID_N] ,q[GRID_N],q[GRID_C],AV[NNm],AV[Nm]);
  RSp = TVDR(q[GRID_S] ,q[GRID_SS],q[GRID_C],q[GRID_S],AV[SSp],AV[Sp]);
  RSm = TVDR(q[GRID_N] ,q[GRID_C] ,q[GRID_C],q[GRID_S],AV[Nm] ,AV[Sm]);

  LNp = SUPERBEE(RNp);
  LNm = SUPERBEE(RNm);
  LSp = SUPERBEE(RSp);
  LSm = SUPERBEE(RSm);

  az = - AV[Sp]*(1 - LSp*0.5) - AV[Sm]*LSm*0.5;
  bz =   AV[Np]*(1 - LNp*0.5) + AV[Nm]*LNm*0.5
       - AV[Sm]*(1 - LSm*0.5) - AV[Sp]*LSp*0.5;
  cz =   AV[Nm]*(1 - LNm*0.5) + AV[Np]*LNp*0.5;

  return (ax*q[GRID_W] + bx*q[GRID_C] + cx*q[GRID_E])/dx + 
         (az*q[GRID_S] + bz*q[GRID_C] + cz*q[GRID_N])/dz; 
}

/*------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "FrommAdvection"
PetscReal FrommAdvection(PetscReal *v, PetscReal *q, PetscReal dx,
			 PetscReal dz)
/*------------------------------------------------------------*/
{
  PetscReal vN,vS,vE,vW;
  PetscReal qC,qN,qNN,qS,qSS,qE,qEE,qW,qWW;
  vN = v[GRID_N];  vE = v[GRID_E];  vS = v[GRID_S];  vW = v[GRID_W]; 
  qC = q[GRID_C]; 
  qN = q[GRID_N]; qNN = q[GRID_NN]; qS = q[GRID_S]; qSS = q[GRID_SS]; 
  qE = q[GRID_E]; qEE = q[GRID_EE]; qW = q[GRID_W]; qWW = q[GRID_WW];
  return ((vE *(-qEE + 5*(qE+qC)-qW)/8  - fabs(vE)*(-qEE + 3*(qE-qC)+qW)/8) -
	  (vW *(-qE  + 5*(qC+qW)-qWW)/8 - fabs(vW)*(-qE  + 3*(qC-qW)+qWW)/8))/dx
       + ((vN *(-qNN + 5*(qN+qC)-qS)/8  - fabs(vN)*(-qNN + 3*(qN-qC)+qS)/8) -
          (vS *(-qN  + 5*(qC+qS)-qSS)/8 - fabs(vS)*(-qN  + 3*(qC-qS)+qSS)/8))/dz;
}

/*------------------------------------------------------------*/
#undef __FUNCT__
#define __FUNCT__ "FVAdvection"
PetscReal FVAdvection(PetscReal *v, PetscReal *q, PetscReal dx, 
		      PetscReal dz)
/*------------------------------------------------------------*/
{
  PetscReal qE, qW, qN, qS;
  qE = 0.5*(q[GRID_E]+q[GRID_C]); qW = 0.5*(q[GRID_W]+q[GRID_C]);
  qN = 0.5*(q[GRID_N]+q[GRID_C]); qS = 0.5*(q[GRID_S]+q[GRID_C]);
  return (v[GRID_E]*qE - v[GRID_W]*qW)/dx + 
         (v[GRID_N]*qN - v[GRID_S]*qS)/dz;
}
