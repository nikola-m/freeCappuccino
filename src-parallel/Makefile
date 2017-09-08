#
# Makefile for caffa3d program
#

# Compiler flags:
LFLAGS = -llapack
F90FLAGS = -Wall -O2

L95FLAGS = -fopenmp -llis -llapack
F95FLAGS = -Wall -O2 -cpp

LIS_DIR = /usr/local/include 

# Compiler:
F90 = mpif90

MOD_FILES=\
    modules_allocatable.f90 \
    utils.f90 \
    matrix.f90 \
    mesh_geometry_and_topology.f90 \
    tensor_fields.f90 \
    sparse_matrix.f90 \
    LIS_linear_solver_library.f90 \
    gradients.f90 \
    output.f90 \
    interpolation.f90 \
    scalar_fluxes.f90

MPI_FILES=\
    global_max_mpi.f90 \
    global_sum_mpi.f90 \
    abort_mission_mpi.f90 \
    synchronize_mpi.f90 \
    mpi_exchange.f90

LINEAR_SOLVER_FILES=\
    iccg.f90 \
    bicgstab.f90 \
    pcg-jacobi.f90 \
    gauss-seidel.f90

TURBULENCE_FILES=\
    temperature.f90 \
    concentration.f90 \
    k_epsilon_std.f90 \
    k_omega_sst.f90 \
    k_epsilon_rng.f90 \
    spalart_allmaras.f90 \
    k_eqn_eddy.f90

CAFFA_FILES=\
    allocate.f90 \
    asm_stress_terms.f90 \
    asm_heatflux_terms.f90 \
    adjustMassFlow.f90 \
    bcin.f90 \
    bpres.f90 \
    fieldManipulation.f90 \
    faceflux_velocity.f90 \
    facefluxmass.f90 \
    facefluxmass_piso.f90 \
    facefluxMassCorr.f90 \
    calcheatflux.f90 \
    calcp-multiple_correction_SIMPLE.f90 \
    calcstress.f90 \
    calcuvw.f90 \
    calc_statistics.f90 \
    correctBoundaryConditionsVelocity.f90 \
    correct_turbulence.f90 \
    correct_turbulence_inlet.f90 \
    fvm_laplacian.f90 \
    find_strain_rate.f90 \
    get_rAU_x_UEqnH.f90 \
    init.f90 \
    openfiles.f90 \
    PISO_multiple_correction.f90 \
    PIMPLE_multiple_correction.f90 \
    readfiles.f90 \
    read_input.f90 \
    random_seed.f90 \
    writefiles.f90 \
    write_restart_files.f90 \
    writehistory.f90 \
    main.f90 

RK4_FILES=\
    assemble_pressure_eq_rk4Projection.f90 \
    fluxmass_plain.f90 \
    fluxuvw-explicit.f90 \
    main_rk4Projection.f90 

POISSON_FILES=\
    fvm_laplacian.f90 \
    poisson.f90

#
# How to create object files:
# 
MODS = ${MOD_FILES:.f90=.o}
MPI = ${MPI_FILES:.f90=.o}
TURBULENCE = ${TURBULENCE_FILES:.f90=.o}
LINEAR_SOLVERS = ${LINEAR_SOLVER_FILES:.f90=.o}
CAFFAOBJS = ${CAFFA_FILES:.f90=.o}
RK4OBJS = ${RK4_FILES:.f90=.o}
POISSONOBJS = ${POISSON_FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: caffa3d-mpi poisson-mpi #rk4ProjectionCaffa

caffa3d-mpi: ${MODS} ${MPI} ${TURBULENCE} ${LINEAR_SOLVERS} ${CAFFAOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${CAFFAOBJS} ${MODS} ${MPI} ${TURBULENCE} ${LINEAR_SOLVERS} ${LFLAGS} ${L95FLAGS} ${INCS} -o ../bin/caffa3d-mpi 

rk4ProjectionCaffa: ${CAFFAOBJS} ${RK4OBJS} ${LINEAR_SOLVERS} ${MPI} 
	@echo  "Linking" $@ "... "
	${F90} ${CAFFAOBJS} ${MPI} ${RK4OBJS} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o rk4ProjectionChannel 

poisson-mpi: ${MODS} ${MPI} ${LINEAR_SOLVERS} ${POISSONOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${POISSONOBJS} ${MODS} ${MPI} ${LINEAR_SOLVERS} ${L95FLAGS} ${INCS} -o poisson-mpi

.PHONY: clean
clean:
	@rm  *.o *.mod caffa3d-mpi poisson-mpi

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90 .f95

.f90.o:
	${F90} ${F90FLAGS} -c ${INCS}  ${@:.o=.f90}

.f95.o:
	${F90} ${F95FLAGS} -I${LIS_DIR} -c ${INCS}  ${@:.o=.f95}
