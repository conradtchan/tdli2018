FC = gfortran -c
FL = gfortran

FLTFL   = -fdefault-real-8

FC_OPTS = -O3 -fopenmp #-g -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace -Wall
FL_OPTS = $(FC_OPTS)

FC_STR  = $(FC) $(FLTFL) $(FC_OPTS)
FL_STR  = $(FL) $(FLTFL) $(FL_OPTS)

hydro: hydro.o boundary.o eos.o primitives.o riemann.o recon.o timestep.o sweep.o rk_step.o setup.o output.o grid.o
	$(FL_STR) -o hydro hydro.o boundary.o eos.o primitives.o riemann.o recon.o timestep.o sweep.o rk_step.o setup.o output.o grid.o

hydro.o: hydro.f90 primitives.o rk_step.o setup.o output.o grid.o
	$(FC_STR) -o hydro.o hydro.f90

test: test.o boundary.o eos.o primitives.o riemann.o recon.o timestep.o sweep.o rk_step.o
	$(FL_STR) -o test test.o boundary.o eos.o primitives.o riemann.o recon.o timestep.o sweep.o rk_step.o

test.o: test.f90 boundary.o eos.o primitives.o riemann.o recon.o timestep.o sweep.o rk_step.o
	$(FC_STR) -o test.o test.f90

boundary.o: boundary.f90
	$(FC_STR) -o boundary.o boundary.f90

eos.o: eos.f90
	$(FC_STR) -o eos.o eos.f90

primitives.o: primitives.f90 eos.o
	$(FC_STR) -o primitives.o primitives.f90

riemann.o: riemann.f90
	$(FC_STR) -o riemann.o riemann.f90

recon.o: recon.f90
	$(FC_STR) -o recon.o recon.f90

timestep.o: timestep.f90
	$(FC_STR) -o timestep.o timestep.f90

sweep.o: sweep.f90 boundary.o eos.o recon.o riemann.o
	$(FC_STR) -o sweep.o sweep.f90

rk_step.o: rk_step.f90 primitives.o timestep.o sweep.o
	$(FC_STR) -o rk_step.o rk_step.f90

setup.o: setup.f90 eos.o
	$(FC_STR) -o setup.o setup.f90

output.o: output.f90
	$(FC_STR) -o output.o output.f90

grid.o: grid.f90
	$(FC_STR) -o grid.o grid.f90

clean:
	rm -f *.o F*.f *.mod output/out*

arx:
	tar -cf hydro.tar *.F Makefile
