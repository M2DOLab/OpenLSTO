CC = g++
CFLAGS = -std=c++14 -Wall -g -fopenmp -w -O3
# CHANGE TO directory where FEA and LSM include folders are
# parent directory for FEA files
# parent directory for FEA files

FEA_DIR = ./../../M2DO_FEA/
LSM_DIR = ./../../M2DO_LSM/

# CHANGE TO the directory where you want your .o files stored
M2DO_FEA_OBJ_DIR = ./../../M2DO_FEA/bin/

# Include directories
FEA_INCL = $(FEA_DIR)/include
LSM_INCL = $(LSM_DIR)/include
INCLUDES = -I$(FEA_INCL) -I$(LSM_INCL)
LIBS = -l lapack -l blas #-l metis -L /usr/local/gfortran/lib -l gfortran #-l ma57 -l metis -l nlopt #-lglpk
# Write the directory where your M2DO_FEA source files are

SRC_DIR = $(FEA_DIR)/src
# List of .o files
# Comment out uncessary object files (if needed)
M2DO_OBJ = $(M2DO_FEA_OBJ_DIR)/quadrature.o \
$(M2DO_FEA_OBJ_DIR)/linear_shape.o \
$(M2DO_FEA_OBJ_DIR)/node.o \
$(M2DO_FEA_OBJ_DIR)/solid_element.o \
$(M2DO_FEA_OBJ_DIR)/solid_material.o \
$(M2DO_FEA_OBJ_DIR)/mesh.o \
$(M2DO_FEA_OBJ_DIR)/boundary_conditions.o \
$(M2DO_FEA_OBJ_DIR)/stationary_study.o \
$(M2DO_FEA_OBJ_DIR)/sensitivity.o \
# Compile all the M2DO FEA object files
all: $(M2DO_OBJ)
# Make commands for each individual object file
$(M2DO_FEA_OBJ_DIR)/quadrature.o: $(SRC_DIR)/quadrature.cpp $(FEA_INCL)/quadrature.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/quadrature.o $(SRC_DIR)/quadrature.cpp
$(M2DO_FEA_OBJ_DIR)/linear_shape.o: $(SRC_DIR)/linear_shape_function.cpp $(FEA_INCL)/linear_shape_function.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/linear_shape.o $(SRC_DIR)/linear_shape_function.cpp
$(M2DO_FEA_OBJ_DIR)/node.o: $(SRC_DIR)/node.cpp $(FEA_INCL)/node.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/node.o $(SRC_DIR)/node.cpp
$(M2DO_FEA_OBJ_DIR)/solid_element.o: $(SRC_DIR)/solid_element.cpp $(FEA_INCL)/solid_element.h $(FEA_INCL)/quadrature.h $(FEA_INCL)/linear_shape_function.h $(FEA_INCL)/mesh.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/solid_element.o $(SRC_DIR)/solid_element.cpp
$(M2DO_FEA_OBJ_DIR)/solid_material.o: $(SRC_DIR)/solid_material.cpp $(FEA_INCL)/solid_material.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/solid_material.o $(SRC_DIR)/solid_material.cpp
$(M2DO_FEA_OBJ_DIR)/mesh.o: $(SRC_DIR)/mesh.cpp $(FEA_INCL)/mesh.h $(FEA_INCL)/quadrature.h $(FEA_INCL)/linear_shape_function.h $(FEA_INCL)/node.h $(FEA_INCL)/solid_element.h $(FEA_INCL)/solid_material.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/mesh.o $(SRC_DIR)/mesh.cpp
$(M2DO_FEA_OBJ_DIR)/boundary_conditions.o: $(SRC_DIR)/boundary_conditions.cpp $(FEA_INCL)/boundary_conditions.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/boundary_conditions.o $(SRC_DIR)/boundary_conditions.cpp
$(M2DO_FEA_OBJ_DIR)/sensitivity.o: $(SRC_DIR)/sensitivity.cpp $(FEA_INCL)/sensitivity.h $(FEA_INCL)/stationary_study.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/sensitivity.o $(SRC_DIR)/sensitivity.cpp
$(M2DO_FEA_OBJ_DIR)/stationary_study.o: $(SRC_DIR)/stationary_study.cpp $(FEA_INCL)/stationary_study.h $(FEA_INCL)/mesh.h $(FEA_INCL)/boundary_conditions.h #$(FEA_INCL)/ma57_solver.h
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $(M2DO_FEA_OBJ_DIR)/stationary_study.o $(SRC_DIR)/stationary_study.cpp
###############################################################################
# Add your own makes here
###############################################################################
# List of targets for executables you write
TARGETS = 
main: main.cpp $(M2DO_OBJ)
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCLUDES) -o ./bin/a.out main.cpp $(M2DO_OBJ) $(LIBS)
.PHONY: clean clean_obj clean_output
clean:
	rm -f bin/a.out *~
clean_obj:
	rm -f $(M2DO_OBJ)
clean_results:
	find ./results -name '*.vtk' -delete
	find ./results -name '*.txt' -delete
