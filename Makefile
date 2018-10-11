#  User-editable options
#
# 1. choose parallel (PARALLEL) or sequential (SEQ) build

#COMP=PARALLEL
 COMP=SEQ

# 2. choose compiler (INTEL or gfortran)

 COMPILER=INTEL
#COMPILER=gfortran

# 3. choose SIMINT base path

 SIMINT_BASE=/home/someone/simint    

# 4. choose MPI library (full path) - optional

#MPI_LIBRARY=/home/my_libs/mpi  


######   WARNING: some expertise required to edit below here!


 GIVENS=USE_FORTRAN
#GIVENS=USE_FORTRAN_KNL
#GIVENS=USE_C
#GIVENS=USE_LAPACK

#PRINT_MATRIX=yes
 PRINT_MATRIX=no
#PRINT_TIMING=yes
 PRINT_TIMING=no
#PRINT_COUNTERS=yes
 PRINT_COUNTERS=no

#INTEGRALS=  #place holder
 INTEGRALS=USE_SIMINT
ifeq ($(INTEGRALS), USE_SIMINT)
SIMINT_LIBRARY=$(SIMINT_BASE)/lib/
SIMINT_INCLUDE=$(SIMINT_BASE)/include/
endif

ifeq ($(COMP), SEQ)
  FC=ifort
  F77=ifort
else
  FC=mpiifort
  #FC=mpif90
  77=mpiifort
  #F77=mpif77
  LINK_FLAGS+=MPI_LIBRARY
endif

#CC=clang++
CC=gcc

ifeq ($(findstring theta,$(HOSTNAME)),theta)
  COMP=PARALLEL
  FC=ftn
  F77=ftn
  CC=cc
endif

ifeq ($(findstring mira,$(HOSTNAME)),mira)
 COMP=PARALLEL
 FC=mpixlf90_r    # bgq
 F77=mpixlf77_r        # bgq
 CC=mpixlc
endif


TARGET=valence.exe
TARGETLIB=libvalence.a

CFLAGS=-O3 -g
FFLAGS=-O3 -g

OBJS=givens.o rsg.o moduletools.o modulevalence_simint.o moduledensity.o moduleintegrals.o valence.o timing_flops.o xm.o valence_api.o valence_initialize.o  valence_finalize.o

ifeq ($(COMPILER), gfortran)
  FC=g95
  F77=g95
  CC=gcc
endif

ifeq ($(COMPILER), INTEL)
 FFLAGS+=-align array64byte -mkl
endif

ifeq ($(GIVENS), USE_C)
  FFLAGS=-O3 -g -DUSE_C_VERSION
  OBJS += givens_in_c.o
endif

ifeq ($(GIVENS), USE_LAPACK)
  FFLAGS+=-DUSE_LAPACK
endif

# alignment here is necessary if simint is vectorized
# CCA order changes the ordering of integrals to Common Component Architecture order 
ifeq ($(INTEGRALS), USE_SIMINT)
  FFLAGS+=-DSIMINT_INT -DCCA_ORDER -align array64byte -I$(SIMINT_INCLUDE)
  LINK_FLAGS+=$(SIMINT_LIBRARY)libsimint.a
endif

ifeq ($(GIVENS), USE_FORTRAN)
  FFLAGS+=-DUSE_FORTRAN_VERSION
endif

ifeq ($(GIVENS), USE_FORTRAN_KNL)
  FFLAGS+=-DGIVENS_KNL
endif

ifeq ($(PRINT_MATRIX), yes)
  FFLAGS+=-DPRINT_MATRIX
endif

ifeq ($(PRINT_TIMING), yes)
  FFLAGS+=-DPRINT_TIMING
endif

ifeq ($(PRINT_COUNTERS), yes)
  FFLAGS+=-DPRINT_COUNTERS
endif

all: libvalence.a valence.exe

libvalence.a: $(OBJS)
	ar rcs $(TARGETLIB) $(OBJS)
valence.exe $(OBJS) valence_driver.o
	$(FC) -o $(TARGET) $(FFLAGS) valence_driver.o $(TARGETLIB) $(LINK_FLAGS)
valence_driver.o : valence_driver.F90 moduletools.o  tools.mod
	$(FC) $(FFLAGS) -c $< -o $@
givens.o : givens.F90 moduletools.o  tools.mod
	$(FC) $(FFLAGS) -c $< -o $@
rsg.o : rsg.F
	$(F77) $(FFLAGS) -c $< -o $@
moduletools.o tools.mod : moduletools.F90
	$(FC) $(FFLAGS) -c $< -o moduletools.o
modulevalence_simint.o valence_simint.mod : modulevalence_simint.F90 moduleintegrals.o integrals.mod
	$(FC) $(FFLAGS) -c $< -o modulevalence_simint.o
moduledensity.o densitywork.mod : moduledensity.F90 moduletools.o tools.mod
	$(FC) $(FFLAGS) -c $< -o moduledensity.o
moduleintegrals.o integrals.mod : moduleintegrals.F90 moduletools.o tools.mod
	$(FC) $(FFLAGS) -c $< -o moduleintegrals.o
timing_flops.o timing_flops.mod : timing_flops.F90
	$(FC) $(FFLAGS) -c $< -o timing_flops.o
valence.o : valence.F90 moduletools.o moduledensity.o moduleintegrals.o tools.mod integrals.mod densitywork.mod timing_flops.o timing_flops.mod modulevalence_simint.o valence_simint.mod xm.o xm.mod valence_initialize.o valence_finalize.o valence_finit.mod
	$(FC) $(FFLAGS) -c $< -o $@ 
givens_in_c.o : givens_in_c.cpp
	$(CC) $(CFLAGS) -c $< -o $@
valence_api.o : valence_api.F90
	$(FC) $(FFLAGS) -c $< -o $@
valence_initialize.o : valence_initialize.F90 
	$(FC) $(FFLAGS) -c $< -o $@
valence_finalize.o : valence_finalize.F90 
	$(FC) $(FFLAGS) -c $< -o $@

ifeq ($(COMP), PARALLEL) 
  FFLAGS+=-DMPI
endif

xm.o : xm.F90 moduletools.o moduledensity.o moduleintegrals.o integrals.mod densitywork.mod timing_flops.o timing_flops.mod tools.mod
	$(FC) $(FFLAGS) -c $< -o $@

doc:
	doxygen Doxyfile
clean:
	rm -f *.o *.mod $(TARGETLIB)
test:
	../testing/test_script_small $(TARGET)

test-large:
	../testing/test_script $(TARGET)

.PHONY: all clean doc test
