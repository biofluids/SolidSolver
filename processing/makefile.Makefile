OBJECTS = read_file.o solidsolver.o

MODULES = read_file.mod

.PHONY: clean
	
output.txt: main.exe
	./main.exe > output.txt
	
main.exe: $(MODULES) $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe
	
%.o: %.f90
	gfortran -c $<
	
%.mod: %.f90
	gfortran -c $<
	
clean:
	rm -f $(OBJECTS) $(MODULES) main.exe