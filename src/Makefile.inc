#############
# Compiler
#############

## The Intel(R) fortran compiler (ifort)
#ifeq ($(DEBUG),1)
#  FC = ifort
#  FCFLAGS = -g -CU -C -traceback -fpe0 -debug -openmp
#  LDFLAGS = -openmp
#else
#  FC = ifort
#  FCFLAGS = -O2 -openmp
#  LDFLAGS = -O2 -openmp
#endif

# The GNU fortran compiler (gfortran)
ifeq ($(DEBUG),1)
  FC=gfortran
  FCFLAGS=-O3 -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core -fopenmp
  LDFLAGS=-fopenmp
else
  FC=gfortran
  FCFLAGS=-O3 -fopenmp
  LDFLAGS=-fopenmp
endif
