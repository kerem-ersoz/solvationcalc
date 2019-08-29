#OBJS specifies which files to compile as part of the project
OBJS = calc.c readopt.c dg_ld.c smol.c lgvnx.c #gen_gridx.c

#CC specifies which compiler we're using
CC = gcc

#COMPILER_FLAGS specifies the additional compilation options we're using
# -w suppresses all warnings
COMPILER_FLAGS = -g

#LINKER_FLAGS specifies the libraries we're linking against
LINKER_FLAGS = -lm

#OBJ_NAME specifies the name of our executable
OBJ_NAME = calc

#This is the target that compiles our executable
all : $(OBJS)
	$(CC) $(OBJS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

#############################################################################


FC = gfortran
FCFLAGS = -g

mods:	
	$(FC) $(FCFLAGS) -fdump-core -c push_array.f90 chemsol.f90
cs21:
	$(FC) $(FCFLAGS) chemsol.o push_array.o main.f90 -o cs21

.PHONY: clean test diff

clean: 
	rm -f *.o *.mod *.MOD cs21
31:
	rm cs21