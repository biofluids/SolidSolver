OBJECTS = read_file.o integration.o solidsolver.o
FC = gfortran

# any source file which uses a module needs to be stated explicitly
# those source files which don't use modules could be included in the macro
main.exe: $(OBJECTS)
	$(FC) -o main $(OBJECTS)
		
#read_file.o: read_file.f90
#	$(FC) -c $(switch) read_file.f90
	
solidsolver.o: read_file.mod integration.mod solidsolver.f90
	$(FC) -c solidsolver.f90
	
read_file.mod: read_file.o read_file.f90
	$(FC) -c read_file.f90
	
integration.mod: integration.o integration.f90
	$(FC) -c integration.f90
	
%.o: %.f90
	$(FC) -c $<
	
clean:
	rm -f *.o *.mod main
