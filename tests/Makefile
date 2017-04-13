#
# Makefile for Cappuccino
#         

F90FLAGS = -O3 -Wall
LFLAGS = 

OBJ_DIR = ../


# Compiler:
F90 = gfortran

# F90FILES=\
#       types.f90 \
#       utils.f90 \
#       mesh_geometry_and_topology.f90 \
#       tensor_fields.f90 \
#       output.f90 

# PoissonFILES=\
#       poissonCappuccino.f90

# GradFILES=\
#       testGaussGradient.f90

testFieldFILES=\
        test_tensor_field_operations.f90


# F90OBJS = ${F90FILES:.f90=.o}
# PoissonOBJS = ${PoissonFILES:.f90=.o}
# GradOBJS = ${GradFILES:.f90=.o}
testFieldOBJS = ${testFieldFILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################
# all: heatCappuccino testGaussGrad testFieldOperations
all: testFieldOperations

# heatCappuccino: ${F90OBJS} ${PoissonOBJS}
# 	@echo  "Linking" $@ "... "
# 	${F90} ${F90OBJS} ${PoissonOBJS} ${LFLAGS} -o heatCappuccino 

# testGaussGrad: ${F90OBJS} ${GradOBJS}
# 	@echo  "Linking" $@ "... "
# 	${F90} ${F90OBJS} ${GradOBJS} ${LFLAGS} -o testGaussGrad 


testFieldOperations: ${testFieldOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${testFieldOBJS} ${LFLAGS} -o testFieldOperations



clean:
	@rm  *.o *.mod

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -I${OBJ_DIR} -c ${@:.o=.f90}

%.o: %.mod
