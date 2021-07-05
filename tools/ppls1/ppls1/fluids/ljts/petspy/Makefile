# compiler
FC  = gfortran
#FPC = f2py
FPC = python3 -m numpy.f2py


# compile flags
#FCFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0
FCFLAGS  = -shared -fdefault-real-8 -fdefault-double-8 -ffree-line-length-0 -fPIC
FPCFLAGS = -c --fcompiler=gnu95 --f90flags="-fdefault-real-8 -fdefault-double-8 -ffree-line-length-0"

PYMODULE=petspy
SOURCE  =pets.F90
LIB     =libpets.so

petspy.so: pywrappets.F90 
	$(FC) $(FCFLAGS) -I. $(SOURCE) -o $(LIB)
	$(FPC) $(FPCFLAGS) -L. -I. -lpets -m $(PYMODULE) $^



clean: 
		rm libpets.so petspy.so pets.mod 
