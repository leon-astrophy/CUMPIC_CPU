#***************************************************************
#
# make.inc file incluing compilation flags
#
#***************************************************************
# section for Fortran compiler #

# Link to MPI fortran, prefer h5pfc #
F90 = h5pfc

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# section for parallel computing 

# use MPI? #
#PARALLEL = -DMPI

# use OPENMP? #
MULTICORE = -fopenmp

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# section for math packages #

# need extra math packages? #
#MATHPACK = -llapack -lblas

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# these are for standard gcc compiler flags

# Define F90 flags #
F90FLAGS= -cpp -Ofast -O4 -march=native -mtune=native -flto -funroll-loops -ftree-vectorize -pipe -ffree-line-length-999 $(PARALLEL) $(MULTICORE) $(MATHPACK) 
LDFLAGS= -cpp -Ofast -O4 -march=native -mtune=native -flto -funroll-loops -ftree-vectorize -pipe -ffree-line-length-999 $(PARALLEL) $(MULTICORE) $(MATHPACK)
          
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
