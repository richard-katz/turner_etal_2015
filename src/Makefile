CFLAGS   =
FFLAGS   = 
CPPFLAGS =
FPPFLAGS =
PROG     = mor
OBJECTS  = ${PROG}.o params.o residual.o utils.o \
	advectionSchemes.o ViscosityModule.o

all:	${PROG}
	${RM} ${OBJECTS}

debug:  ${PROG}

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


${PROG}: ${OBJECTS} chkopts
	${RM} ${PROG} ${CFLAGS}
	-${CLINKER} -o ${PROG} ${OBJECTS} ${PETSC_SNES_LIB} ${PETSC_CONTRIB}
