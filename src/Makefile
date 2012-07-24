# Compile/link options which you will need to set:
# ------------------------------------------------
# FC = Command which invokes your Fortran compiler (e.g. f90, ifort, gfortran, ...)
# FFLAGS = List of flags for your Fortran compiler (optimization/debugging/profiling)
# DEBUGS = List of strict debugging flags for your Fortran compiler (for testing)
# LAPACK = Compiler flag for picking up system copy of LAPACK & BLAS
#
# NCLINK = Linker flags for picking up the netCDF library
#
# MAKEFILE = Makefile on which one may want the make process to depend,
#            (e.g. to ensure rebuild following preprocessor option changes)


# ANU NCI system (VAYU):
# ----------------------
#FC = ifort
## SEQUENTIAL:
#FFLAGS = -O3 -mp1 -auto -p ### FOR PROFILING
## OPENMP
###LAPACK = -lscs_mp
#FFLAGS = -ip -O3 -xHost -warn
#FFLAGS_OPENMP = -openmp -DOPENMP
#FFLAGS_DEBUG = -debug -g -check all -traceback
#FFLAGS = -O3 -mp1 -p -openmp   ### FOR PROFILING
#NCLINK = -lnetcdf -lnetcdff
#LAPACK = -L$(MKL)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_em64t -lguide -lpthread
#FFTW = -lfftw3

## Linux gfortran
FC := gfortran
FFLAGS = -march=native -O3 -Wall -I/usr/local/include -I/usr/include
FFLAGS_DEBUG = -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero
FFLAGS_OPENMP = -fopenmp -DOPENMP
FFLAGS_PROFILE = -pg -fno-omit-frame-pointer -g
NCLINK = -lnetcdff
LAPACK = -llapack
FFTW = -L/usr/local/lib -lfftw3

MAKEFILE = Makefile

OBJS =	openqg.o radsubs.o valsubs.o subsampling.o forcing.o \
	tavsubs.o areasubs.o covsubs.o vorsubs.o linalg.o \
	omlsubs.o core.o state.o print.o ml_monitor.o \
	amlsubs.o coupler.o windstress.o driver.o basin.o \
	ncutils.o bicubic.o box.o glam.o grid.o mesh.o numerics.o topog.o qg_monitor.o \
	util.o intsubs.o constants.o units.o mixed.o qg.o modes.o homog.o inhomog.o \
	ekman.o clock.o vorticity.o mixed_layer.o constraint.o pressure.o

EXEC := openqg


default: parallel

debug-single: FFLAGS += ${FFLAGS_DEBUG}
debug-single: EXEC :=${EXEC}-debug-single
debug-single: executable
profile-single: FFLAGS += ${FFLAGS_PROFILE}
profile-single: EXEC := ${EXEC}-profile-single
profile-single: executable
single: EXEC := ${EXEC}-single
single: executable

debug-parallel: FFLAGS += ${FFLAGS_DEBUG}
debug-parallel: EXEC :=${EXEC}-debug-parallel
debug-parallel: parallel
profile-parallel: FC := kinst-ompp-papi ${FC}
profile-parallel: EXEC := ${EXEC}-profile-parallel
profile-parallel: parallel
parallel: FFLAGS += ${FFLAGS_OPENMP}
parallel: executable
	cp ${EXEC} ${EXEC}-parallel

executable: ${OBJS}
	${FC} ${FFLAGS} -o ${EXEC} ${OBJS} ${LAPACK} ${NCLINK} ${FFTW}

## Util

util.o: util/util.F90
	${FC} ${FFLAGS} -c util/util.F90
util.mod: util.o

ncutils.o: util/ncutils.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c util/ncutils.F90

intsubs.o: util/intsubs.F90
	${FC} ${FFLAGS} -c util/intsubs.F90

constants.o: util/constants.F90
	${FC} ${FFLAGS} -c util/constants.F90

units.o: util/units.F90
	${FC} ${FFLAGS} -c util/units.F90

linalg.o: util/linalg.F90
	${FC} ${FFLAGS} -c util/linalg.F90

## GLAM

mesh.o: glam/mesh.F90 ncutils.o constants.o units.o
	${FC} ${FFLAGS} -c glam/mesh.F90

box.o: glam/box.F90 ncutils.o mesh.o units.o
	${FC} ${FFLAGS} -c glam/box.F90

grid.o: glam/grid.F90 ncutils.o box.o units.o
	${FC} ${FFLAGS} -c glam/grid.F90

topog.o: glam/topog.F90 ncutils.o util.o intsubs.o box.o grid.o units.o
	${FC} ${FFLAGS} -c glam/topog.F90

bicubic.o: glam/bicubic.F90 box.o grid.o
	${FC} ${FFLAGS} -c glam/bicubic.F90

glam.o: glam/glam.F90 $ ncutils.o mesh.o
	${FC} ${FFLAGS} -c glam/glam.F90

numerics.o: glam/numerics.F90 intsubs.o box.o constants.o
	${FC} ${FFLAGS} -c glam/numerics.F90

#######################

## QG 

modes.o: qg/modes.F90 box.o topog.o
	${FC} ${FFLAGS} -c qg/modes.F90

inhomog.o: qg/inhomog.F90 box.o constants.o modes.o
	${FC} ${FFLAGS} -c qg/inhomog.F90

constraint.o: qg/constraint.F90 box.o modes.o intsubs.o
	${FC} ${FFLAGS} -c qg/constraint.F90

homog.o: qg/homog.F90 box.o modes.o intsubs.o numerics.o constraint.o inhomog.o linalg.o
	${FC} ${FFLAGS} -c qg/homog.F90

qg.o: qg/qg.F90 util.o ncutils.o box.o numerics.o intsubs.o topog.o modes.o homog.o inhomog.o
	${FC} ${FFLAGS} -c qg/qg.F90

pressure.o: qg/pressure.F90 intsubs.o box.o topog.o modes.o homog.o inhomog.o qg.o 
	${FC} ${FFLAGS} -c qg/pressure.F90

vorsubs.o: qg/vorsubs.F90 ${MAKEFILE} box.o topog.o modes.o numerics.o qg.o
	${FC} ${FFLAGS} -c qg/vorsubs.F90

vorticity.o: qg/vorticity.F90 intsubs.o qg.o numerics.o qg.o
	${FC} ${FFLAGS} -c qg/vorticity.F90

core.o: qg/core.F90 ${MAKEFILE} qg.o vorticity.o pressure.o vorsubs.o
	${FC} ${FFLAGS} -c qg/core.F90

## Mixed layer

mixed.o: mixed/mixed.F90 ncutils.o box.o numerics.o
	${FC} ${FFLAGS} -c mixed/mixed.F90

mixed_layer.o: mixed/mixed_layer.F90 box.o numerics.o mixed.o
	${FC} ${FFLAGS} -c mixed/mixed_layer.F90

#---

omlsubs.o: mixed/omlsubs.F90 ${MAKEFILE} box.o mixed.o numerics.o mixed_layer.o
	${FC} ${FFLAGS} -c mixed/omlsubs.F90

radsubs.o: mixed/radsubs.F90 ${MAKEFILE} box.o grid.o ncutils.o intsubs.o mixed.o constants.o linalg.o
	${FC} ${FFLAGS} -c mixed/radsubs.F90

amlsubs.o: mixed/amlsubs.F90 ${MAKEFILE} ncutils.o box.o topog.o mixed.o radsubs.o numerics.o mixed_layer.o
	${FC} ${FFLAGS} -c mixed/amlsubs.F90

forcing.o: mixed/forcing.F90 ${MAKEFILE} grid.o box.o topog.o mixed.o radsubs.o amlsubs.o
	${FC} ${FFLAGS} -c mixed/forcing.F90

## Coupler

ekman.o: coupler/ekman.F90 intsubs.o box.o numerics.o ncutils.o
	${FC} ${FFLAGS} -c coupler/ekman.F90

coupler.o: coupler/coupler.F90 ${MAKEFILE} ekman.o box.o
	${FC} ${FFLAGS} -c coupler/coupler.F90

## Windstress

windstress.o: windstress/windstress.F90 intsubs.o box.o grid.o bicubic.o topog.o numerics.o ncutils.o qg.o ekman.o
	${FC} ${FFLAGS} -c windstress/windstress.F90

#####################

## Types

clock.o: types/clock.F90 ncutils.o constants.o
	${FC} ${FFLAGS} -c types/clock.F90

## Diagnostics

ml_monitor.o: diagnostics/ml_monitor.F90 ${MAKEFILE} box.o qg.o grid.o ekman.o omlsubs.o amlsubs.o intsubs.o 
	${FC} ${FFLAGS} -c diagnostics/ml_monitor.F90

valsubs.o: diagnostics/valsubs.F90 ${MAKEFILE} mixed.o qg.o ekman.o omlsubs.o amlsubs.o
	${FC} ${FFLAGS} -c diagnostics/valsubs.F90

areasubs.o: diagnostics/areasubs.F90 ${MAKEFILE} ncutils.o mixed.o box.o units.o
	${FC} ${FFLAGS} -c diagnostics/areasubs.F90

print.o: diagnostics/print.F90 ${MAKEFILE} qg.o omlsubs.o amlsubs.o
	${FC} ${FFLAGS} -c diagnostics/print.F90

subsampling.o: diagnostics/subsampling.F90 ${MAKEFILE}  box.o mixed.o qg.o ekman.o ncutils.o clock.o units.o
	${FC} ${FFLAGS} -c diagnostics/subsampling.F90

covsubs.o: diagnostics/covsubs.F90 ${MAKEFILE} ncutils.o mixed.o qg.o clock.o box.o
	${FC} ${FFLAGS} -c diagnostics/covsubs.F90

tavsubs.o: diagnostics/tavsubs.F90 ${MAKEFILE} ncutils.o box.o mixed.o qg.o ekman.o clock.o units.o
	${FC} ${FFLAGS} -c diagnostics/tavsubs.F90

qg_monitor.o: diagnostics/qg_monitor.F90 box.o qg.o ncutils.o ekman.o numerics.o intsubs.o units.o
	${FC} ${FFLAGS} -c diagnostics/qg_monitor.F90

####

basin.o: basin.F90 ${MAKEFILE} omlsubs.o ekman.o box.o qg.o amlsubs.o covsubs.o tavsubs.o qg_monitor.o coupler.o print.o ml_monitor.o subsampling.o
	${FC} ${FFLAGS} -c basin.F90

##

state.o: state.F90 ${MAKEFILE} qg.o mixed.o box.o radsubs.o util.o ncutils.o vorsubs.o
	${FC} ${FFLAGS} -c state.F90

driver.o: driver.F90 ${MAKEFILE} coupler.o windstress.o omlsubs.o amlsubs.o core.o forcing.o basin.o constants.o
	${FC} ${FFLAGS} -c driver.F90

openqg.o: openqg.F90 ${MAKEFILE}  \
	ncutils.o glam.o omlsubs.o forcing.o mixed.o constants.o units.o \
	qg.o modes.o homog.o inhomog.o subsampling.o clock.o vorsubs.o tavsubs.o qg_monitor.o \
	valsubs.o covsubs.o areasubs.o vorticity.o radsubs.o constraint.o pressure.o covsubs.o \
	amlsubs.o core.o state.o coupler.o windstress.o driver.o basin.o print.o ml_monitor.o
	${FC} ${FFLAGS} -c openqg.F90


##
nuke: clean
	rm -f openqg openqg-single openqg-profile-single openqg-debug-single openqg-debug-parallel openqg-debug-profile


clean:
	rm -f *.o *.l *.lst *.w2f.* *.fln *.list *.chek *.mod *genmod* *opari.inc *.mod.F90 opari.rc opari.tab.c
	find . -name *.mod.F90 | xargs -r rm

tidy:
	rm -f *.o *.l *.lst *.w2f.* *.fln *.list *.chek

