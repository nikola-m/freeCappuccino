#
# Makefile for caffa3d program
#

# Compiler flags:
LFLAGS = -llapack 
F90FLAGS = -Wall -O2 

# Compiler:
F90 = gfortran

MOD_FILES=\
    modules_allocatable.f90 \
    utils.f90 \
    matrix.f90 \
    mesh_geometry_and_topology.f90 \
    sparse_matrix.f90 \
    gradients.f90 \
    output.f90 \
    allocate.f90

LINEAR_SOLVER_FILES=\
    iccg.f90 \
    bicgstab.f90
    # sipsol.f90 \
    # cgstab_sip.f90 \
    # pcg-jacobi.f90 \
    # pcg-sip.f90

TURBULENCE_FILES=\
    temperature.f90 \
    k_epsilon_std.f90

SRCS=\
    set_parameters.f90 \
    asm_stress_terms.f90 \
    asm_heatflux_terms.f90 \
    adjustMassFlow.f90 \
    bcin.f90 \
    bpres.f90 \
    fieldManipulation.f90 \
    calcheatflux.f90 \
    calcp-multiple_correction_SIMPLE.f90 \
    calcstress.f90 \
    calcuvw.f90 \
    calc_statistics.f90 \
    correctBoundaryConditionsVelocity.f90 \
    correct_turbulence.f90 \
    correct_turbulence_inlet.f90 \
    facefluxmass.f90 \
    facefluxmass_piso.f90 \
    facefluxMassCorr.f90 \
    facefluxuvw.f90 \
    fvm_laplacian.f90 \
    find_strain_rate.f90 \
    get_rAU_x_UEqnH.f90 \
    init.f90 \
    openfiles.f90 \
    PISO_multiple_correction.f90 \
    PIMPLE_multiple_correction.f90 \
    readfiles.f90 \
    random_seed.f90 \
    writefiles.f90 \
    write_restart_files.f90 \
    writehistory.f90 \
    main.f90 

RK4FILES=\
    assemble_pressure_eq_rk4Projection.f90 \
    fluxmass_plain.f90 \
    fluxuvw-explicit.f90 \
    main_rk4Projection.f90 

#
# How to create object files:
# 
MODS = ${MOD_FILES:.f90=.o}
TURBULENCE = ${TURBULENCE_FILES:.f90=.o}
LINEAR_SOLVERS = ${LINEAR_SOLVER_FILES:.f90=.o}
F90OBJS = ${SRCS:.f90=.o}
OBJS = ${RK4FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: caffa3d #rk4Projection

caffa3d: ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${F90OBJS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o caffa3d 

rk4ProjectionChannel: ${F90OBJS} ${RK4OBJS} ${LINEAR_SOLVERS}
	@echo  "Linking" $@ "... "
	${F90} ${F90OBJS} ${RK4OBJS} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o rk4ProjectionChannel 

.PHONY: clean
clean:
	@rm  *.o *.mod

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 

.f90.o:
	${F90} ${F90FLAGS} -c ${INCS}  ${@:.o=.f90}
