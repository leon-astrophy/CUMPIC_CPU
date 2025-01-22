#***************************************************************
#
# make.inc file incluing compilation flags
#
#***************************************************************
# section for Fortran compiler #

# MPI fortran #
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

# compiler flags
#HPC = -fast -O4 -Mconcur -Munroll -Mdalign -Mframe -Mlre -Mstack_arrays -Mvect=short -Mvect=prefetch -Mipa=fast,inline -Mvect=levels:4 -Mvect=fuse -Mvect=gather -Mvect=idiom -Mvect=altcode -Mvect=simd 
ifort = -traceback

# Define F90 flags #
F90FLAGS= -cpp $(ifort) $(HPC) $(PARALLEL) $(MULTICORE) $(MATHPACK) $(COMPILE) 
LDFLAGS= -cpp $(ifort) $(HPC) $(PARALLEL) $(MULTICORE) $(MATHPACK) $(COMPILE)   

#F90FLAGS= -cpp -Ofast -O4 -march=native -mtune=native -flto -funroll-loops -ftree-vectorize -pipe -ffree-line-length-999 $(PARALLEL) $(MULTICORE) $(MATHPACK) 
#LDFLAGS= -cpp -Ofast -O4 -march=native -mtune=native -flto -funroll-loops -ftree-vectorize -pipe -ffree-line-length-999 $(PARALLEL) $(MULTICORE) $(MATHPACK)
          
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
