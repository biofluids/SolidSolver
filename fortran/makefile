LIBS = /usr/local/lib/liblapack.a  /usr/local/lib/libblas.a
OBJECTS = read_file.o integration.o face.o shapefunction.o material.o externalforce.o internalforce.o \
	      tangentstiffness.o mgmres.o output.o mass.o solidsolver.o
MODULES = read_file.mod integration.mod face.mod shapefunction.mod material.mod externalforce.mod internalforce.mod \
          tangentstiffness.mod mgmres.mod output.mod mass.mod
FC = gfortran

main.exe: $(OBJECTS)
	$(FC) -fcheck=all -o main $(OBJECTS) $(LIBS)
		
#read_file.o: read_file.f90
#	$(FC) -c $(switch) read_file.f90
	
solidsolver.o: $(MODULES) solidsolver.f90
	$(FC) -c solidsolver.f90

%.mod: %.f90
	$(FC) -c $<
	
%.o: %.f90
	$(FC) -c $<
	
clean:
	rm -f *.o *.mod main
