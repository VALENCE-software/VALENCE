FC=gfortran #ifort
FFLAGS=-O3 -g 

all: modelkit

modelkit: modelkit.f90
	$(FC) $(FFLAGS) -o modelkit $<

clean:
	rm -f *.o *.mod 
