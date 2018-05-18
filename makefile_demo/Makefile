
# helloworld is the final product - we want to link together all of the
# object files that we've made. The first line specifies the depdendencies:
# if program.o, hello.o, and world.o do not exist yet, then it will
# automatically first build those recipes before starting this one.
# The second line is the command that is run. In this case, we use gfortran
# to link together all of the .o files, and make 'helloworld'
helloworld: program.o hello.o world.o
	gfortran -o helloworld program.o hello.o world.o

# These are the objects for each .f90 file
program.o: program.f90 hello.o world.o
	gfortran -c -o program.o program.f90

hello.o: hello.f90
	gfortran -c -o hello.o hello.f90

world.o: world.f90
	gfortran -c -o world.o world.f90

# The 'make clean' command is a nice shortcut to clean up the mess you've made
clean:
	rm -f *.o *.mod helloworld
