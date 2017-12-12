#
# Makefile for writeCellsPolyMesh program
#         

F90FLAGS = -O2 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      qsort_c_module.f90 \
      cells_module.f90 \
      writeCellsPolyMesh.f90



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: writeCellsPolyMesh

writeCellsPolyMesh: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/writeCellsPolyMesh

.PHONY: clean
clean:
	@rm  *.o *.mod ../../bin/writeCellsPolyMesh

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
