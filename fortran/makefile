LIBS = /usr/local/lib/liblapack.a  /usr/local/lib/libblas.a /usr/local/lib/libma57.a /usr/local/lib/libmetis.a
OBJECTS = read_file.o integration.o face.o shapefunction.o material.o externalforce.o internalforce.o \
	      tangentstiffness.o mgmres.o directsolver.o output.o mass.o volumecheck.o \
		  full_tangentstiffness.o full_internalforce.o solidsolver.o
MODULES = read_file.mod integration.mod face.mod shapefunction.mod material.mod externalforce.mod internalforce.mod \
          tangentstiffness.mod mgmres.mod directsolver.mod output.mod mass.mod volumecheck.mod \
		  full_tangentstiffness.mod full_internalforce.mod
FC = gfortran

main.exe: $(OBJECTS)
	$(FC) -o main.exe $(OBJECTS) $(LIBS)
			
solidsolver.o: $(MODULES) solidsolver.f90
	$(FC) -c solidsolver.f90

%.mod: %.f90
	$(FC) -c $<
	
%.o: %.f90
	$(FC) -c $<
	
clean:
	rm -f *.o *.mod main.exe
