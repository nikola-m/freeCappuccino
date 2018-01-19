#
# Makefile for cellConnectivity program
#         

F90FLAGS = -O2 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      utils.f90 \
      cellConnectivity.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: cellConnectivity

cellConnectivity: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o ../../bin/cellConnectivity

.PHONY: clean
clean:
	@rm  *.o *.mod ../../bin/cellConnectivity

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
