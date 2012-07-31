FC := gfortran
FFLAGS = -march=native -O3 -Wall -I/usr/local/include -I/usr/include
FFLAGS_DEBUG = -fbacktrace -g -fbounds-check -ffpe-trap=invalid,zero
FFLAGS_OPENMP = -fopenmp -DOPENMP
FFLAGS_PROFILE = -pg -fno-omit-frame-pointer -g
NCLINK = -lnetcdff
LAPACK = -llapack
FFTW = -L/usr/local/lib -lfftw3

MAKEFILE = Makefile

OBJS = ncutils.o intsubs.o constants.o units.o linalg.o mesh.o box.o numerics.o inhomog.o test_inhomog.o

default: test_inhomog
EXEC = test_inhomog

test_inhomog: ${OBJS}
	${FC} ${FFLAGS} -o ${EXEC} ${OBJS} ${LAPACK} ${NCLINK} ${FFTW}


ncutils.o: ../util/ncutils.F90 ${MAKEFILE}
	${FC} ${FFLAGS} -c ../util/ncutils.F90

intsubs.o: ../util/intsubs.F90
	${FC} ${FFLAGS} -c ../util/intsubs.F90

constants.o: ../util/constants.F90
	${FC} ${FFLAGS} -c ../util/constants.F90

units.o: ../util/units.F90
	${FC} ${FFLAGS} -c ../util/units.F90

linalg.o: ../util/linalg.F90
	${FC} ${FFLAGS} -c ../util/linalg.F90

mesh.o: ../glam/mesh.F90 ncutils.o constants.o units.o
	${FC} ${FFLAGS} -c ../glam/mesh.F90

box.o: ../glam/box.F90 ncutils.o mesh.o units.o
	${FC} ${FFLAGS} -c ../glam/box.F90

numerics.o: ../glam/numerics.F90 intsubs.o box.o constants.o
	${FC} ${FFLAGS} -c ../glam/numerics.F90

inhomog.o: ../qg/inhomog.F90 box.o constants.o
	${FC} ${FFLAGS} -c ../qg/inhomog.F90

test_inhomog.o: test_inhomog.F90 inhomog.o
	${FC} ${FFLAGS} -c test_inhomog.F90

clean:
	rm -f *.o *.l *.lst *.w2f.* *.fln *.list *.chek *.mod *genmod* *opari.inc *.mod.F90 opari.rc opari.tab.c
	find . -name *.mod.F90 | xargs -r rm
