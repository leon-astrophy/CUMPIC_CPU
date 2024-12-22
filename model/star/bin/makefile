#include header#
include $(CUMPIC_DIR)/make.inc

#specify HDF5#
EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS +=  $(HDF5LIBS)

#specify directories#
SOLVER  := $(CUMPIC_DIR)/solver
MODEL   := $(CUMPIC_DIR)/model/star/src
CUSTOM  := $(CUMPIC_DIR)/model/star/custom-solver

vpath %.f90 $(CUSTOM):$(SOLVER):$(MODEL)

#source list#
SOURCE = $(shell cat $(SOLVER)/source.list) $(shell cat $(MODEL)/source.list)
SRC = $(notdir $(SOURCE))
OBJ = $(patsubst %.f90,%.o,$(SRC))

#core object files#
DEF_OBJ = definition.o
CUS_MOD = def.o

OTHER_OBJ = $(filter-out $(SIMPLE_MOD) $(DEF_OBJ), $(OBJ))

#you can use make -j to build in parallel#
ALL: $(OBJ)
	$(F90) $(LDFLAGS) -o ../CUMC3D $(OBJ) $(EXTRAOBJECTS)

$(CUS_MOD): $(DEF_OBJ)
$(OTHER_OBJ) : $(DEF_OBJ) $(CUS_MOD)

%.o: %.f90
	$(F90) $(F90FLAGS) $(EXTRAINCS) -c $< -o $@
