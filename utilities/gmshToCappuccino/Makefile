#
# Makefile for gmshToCappuccino program
#         

F90FLAGS = -O2 -Wall

# Compiler:
F90 = gfortran

F90FILES=\
      util.f90 \
      qsort_c_module.f90 \
      gmshToCappuccino.f90 



F90OBJS = ${F90FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

gmshToCappuccino: ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} -o gmshToCappuccino 
clean:
	@rm  *.o *.mod gmshToCappuccino

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c  ${@:.o=.f90}

%.o: %.mod
