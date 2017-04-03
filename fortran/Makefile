LIBS = /usr/local/lib/liblapack.a  /usr/local/lib/libblas.a /usr/local/lib/libma57.a /usr/local/lib/libmetis.a
OBJECTS = crs.o read_file.o integration.o face.o shapefunction.o material.o externalforce.o internalforce.o \
	      tangentstiffness.o symmetric_solver.o output.o mass.o \
		  solidsolver.o
MODULES = read_file.mod integration.mod face.mod shapefunction.mod material.mod externalforce.mod internalforce.mod \
          tangentstiffness.mod symmetric_solver.mod output.mod mass.mod
FC = gfortran
CPP = g++
FLAG = -lblas -lstdc++

main: $(OBJECTS)
	$(FC) -o main $(OBJECTS) $(LIBS) $(FLAG)

pre: pre_processing.cpp
	$(CPP) pre_processing.cpp -o pre -std=c++11

crs.o: crs.cpp
	$(CPP) -c crs.cpp -std=c++11

symmetric_solver.o: symmetric_solver.f90 read_file.mod
	$(FC) -c symmetric_solver.f90

internalforce.o: internalforce.f90 read_file.mod shapefunction.mod integration.mod material.mod
	$(FC) -c internalforce.f90

mass.o: mass.f90 read_file.mod shapefunction.mod integration.mod
	$(FC) -c mass.f90

output.o: output.f90 read_file.mod shapefunction.mod integration.mod material.mod
	$(FC) -c output.f90

read_file.o: read_file.f90 integration.mod
	$(FC) -c read_file.f90

tangentstiffness.o: tangentstiffness.f90 read_file.mod shapefunction.mod integration.mod material.mod
	$(FC) -c tangentstiffness.f90

%.mod: %.f90 %.o
	@true

%.o: %.f90
	$(FC) -c $<

clean:
	rm -f *.o *.mod main pre
