# Makefile for kinetic Monte Carlo surface simulations

# Variables 
FC = /usr/bin/gfortran
FIL ?= case_1   # set FIL=name if different case

objects = $(FIL).o kinetic_tools_$(FIL).o main_$(FIL).o spline.o  
obj_mod = $(FIL).mod kinetic_tools_$(FIL).mod

LDFLAGS = -fopenmp -g -Wall -pedantic -fcheck=all -fbacktrace

# Executable file
walks: $(objects)
	$(FC) $(LDFLAGS) -o walks $(objects)

# Compile to .o files
%.o: %.f90
	$(FC) $(LDFLAGS) -c $<
%.o: %.f
	$(FC) $(LDFLAGS) -c $<

# Clean up .o files
clean:
	rm $(objects), $(obj_mod)

valgrind:
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-origins=yes \
	--num-callers=20 --track-fds=yes ./walks

valgrind_parallel:
	valgrind --tool=helgrind ./walks

run:
	./walks &
