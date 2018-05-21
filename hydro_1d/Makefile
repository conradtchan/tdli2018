FC = gfortran -c
FL = gfortran

FLTFL   = -fdefault-real-8

FC_OPTS = -g -O2 -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace
FL_OPTS = $(FC_OPTS)

FC_STR  = $(FC) $(FLTFL) $(FC_OPTS)
FL_STR  = $(FL) $(FLTFL) $(FL_OPTS)

test: test.o boundary.o eos.o primitives.o
	$(FL_STR) -o test test.o boundary.o eos.o primitives.o

test.o: test.f90 boundary.o eos.o primitives.o
	$(FC_STR) -o test.o test.f90

boundary.o: boundary.f90
	$(FC_STR) -o boundary.o boundary.f90

eos.o: eos.f90
	$(FC_STR) -o eos.o eos.f90

primitives.o: primitives.f90
	$(FC_STR) -o primitives.o primitives.f90

clean:
	rm -f *.o F*.f *.mod

arx:
	tar -cf hydro.tar *.F Makefile
