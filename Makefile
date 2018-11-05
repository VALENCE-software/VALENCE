#  User-editable options
#
# 1. choose parallel or sequential build

 SEQUENTIAL=true
#SEQUENTIAL=false

ifeq ($(SEQUENTIAL), false) 
  FFLAGS+=-DVALENCE_MPI
endif

# 2. choose compiler (INTEL or GNU)

#COMPILER_VENDOR=INTEL
 COMPILER_VENDOR=GNU

ifeq ($(COMPILER_VENDOR), GNU)
  ifeq ($(SEQUENTIAL), false) 
    FC=mpifort
    F77=mpifort
    CC=mpicc
  else
    FC=gfortran
    F77=gfortran
    CC=gcc
  endif
endif

ifeq ($(COMPILER_VENDOR), INTEL)
  ifeq ($(SEQUENTIAL), false) 
    FC=mpiifort
    F77=mpiifort
    CC=mpiicc
  else
    FC=ifort
    F77=ifort
    CC=icc
  endif
  FFLAGS+=-align array64byte -mkl
endif

# 3. choose SIMINT base path

SIMINT_BASE=./simint
INTEGRALS=USE_SIMINT
ifeq ($(INTEGRALS), USE_SIMINT)
  SIMINT_LIBRARY=$(SIMINT_BASE)/lib/
  SIMINT_INCLUDE=$(SIMINT_BASE)/include/simint

# alignment here is necessary if simint is vectorized
# CCA order changes the ordering of integrals to Common Component Architecture order 

  FFLAGS+=-DSIMINT_INT -DCCA_ORDER -I$(SIMINT_INCLUDE)
  LINK_FLAGS+=$(SIMINT_LIBRARY)libsimint.a
endif


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

ifeq ($(findstring theta,$(HOSTNAME)),theta)
  FC=ftn
  F77=ftn
  CC=cc
  FFLAGS+=-align array64byte -mkl
endif

ifeq ($(findstring mira,$(HOSTNAME)),mira)
  FC=mpixlf90_r
  F77=mpixlf77_r
  CC=mpixlc
endif


CFLAGS=-O3 -g
FFLAGS+=-O3 -g

OBJS=givens.o rsg.o tools_module.o valence_simint_module.o density_module.o\
	 integrals_module.o valence.o timing_flops_module.o xm_module.o \
	 valence_api_nitrogen.o valence_initialize_module.o  \
	 valence_finalize_module.o

ifeq ($(GIVENS), USE_C)
  FFLAGS=-O3 -g -DUSE_C_VERSION
  OBJS += givens_in_c.o
endif

ifeq ($(GIVENS), USE_LAPACK)
  FFLAGS+=-DUSE_LAPACK
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

# change these to proper directories where each file should be
SRCDIR   = src
BINDIR   = bin
LIBDIR   = lib

TARGET=valence
TARGETMK=modelkit
TARGETLIB=libvalence.a

all: $(TARGETLIB) $(TARGET) $(TARGETMK) clean
$(TARGETLIB): $(OBJS)
	mkdir -p $(LIBDIR)
	ar rcs $(LIBDIR)/$(TARGETLIB) $(OBJS)
$(TARGET): $(OBJS) valence_driver.o
	mkdir -p $(BINDIR)
	$(FC) -o $(BINDIR)/$(TARGET) $(FFLAGS) valence_driver.o $(LIBDIR)/$(TARGETLIB) $(LINK_FLAGS)
$(TARGETMK) : 
	$(FC) $(FFLAGS) $(SRCDIR)/modelkit.F90 -o $(BINDIR)/$(TARGETMK) 
valence_driver.o : $(SRCDIR)/valence_driver.F90 tools_module.o  tools.mod
	$(FC) $(FFLAGS) -c $< -o $@
givens.o : $(SRCDIR)/givens.F90 tools_module.o  tools.mod
	$(FC) $(FFLAGS) -c $< -o $@
rsg.o : $(SRCDIR)/rsg.F
	$(F77) $(FFLAGS) -c $< -o $@
tools_module.o tools.mod : $(SRCDIR)/tools_module.F90
	$(FC) $(FFLAGS) -c $< -o tools_module.o
valence_simint_module.o valence_simint.mod : $(SRCDIR)/valence_simint_module.F90 \
	integrals_module.o integrals.mod
	$(FC) $(FFLAGS) -c $< -o valence_simint_module.o
density_module.o densitywork.mod : $(SRCDIR)/density_module.F90 tools_module.o tools.mod
	$(FC) $(FFLAGS) -c $< -o density_module.o
integrals_module.o integrals.mod : $(SRCDIR)/integrals_module.F90 tools_module.o tools.mod
	$(FC) $(FFLAGS) -c $< -o integrals_module.o
timing_flops_module.o timing_flops.mod : $(SRCDIR)/timing_flops_module.F90
	$(FC) $(FFLAGS) -c $< -o timing_flops_module.o
valence.o : $(SRCDIR)/valence.F90 tools_module.o density_module.o integrals_module.o \
	tools.mod integrals.mod densitywork.mod timing_flops_module.o \
	timing_flops.mod valence_simint_module.o valence_simint.mod \
	 xm_module.o xm.mod valence_initialize_module.o \
	 valence_finalize_module.o valence_finit.mod
	$(FC) $(FFLAGS) -c $< -o $@ 
givens_in_c.o : $(SRCDIR)/givens_in_c.cpp
	$(CC) $(CFLAGS) -c $< -o $@
valence_api_nitrogen.o : $(SRCDIR)/valence_api_nitrogen.F90
	$(FC) $(FFLAGS) -c $< -o $@
valence_initialize_module.o : $(SRCDIR)/valence_initialize_module.F90 xm_module.o xm.mod
	$(FC) $(FFLAGS) -c $< -o $@
valence_finalize_module.o : $(SRCDIR)/valence_finalize_module.F90 
	$(FC) $(FFLAGS) -c $< -o $@
xm_module.o xm.mod : $(SRCDIR)/xm_module.F90 tools_module.o density_module.o \
	integrals_module.o integrals.mod densitywork.mod timing_flops_module.o \
	 timing_flops.mod tools.mod
	$(FC) $(FFLAGS) -c $< -o $@
doc:
	doxygen doc/Doxyfile
clean:
	rm -f *.o *.mod 
veryclean:
	rm -f *.o *.mod $(LIBDIR)/$(TARGETLIB) $(BINDIR)/$(TARGET) $(BINDIR)/$(TARGETMK)
test:
	testing/test_script $(PWD)/$(BINDIR)/$(TARGET) $(SEQUENTIAL) small
test-large:
	testing/test_script $(PWD)/$(BINDIR)/$(TARGET) $(SEQUENTIAL) large
.PHONY: all clean doc test
