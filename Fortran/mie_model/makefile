# Define some variables

MOD_PATH  = mod
OBJ_PATH  = obj
RTSMIE_PATH = ./Sourcecode

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

#SHELL = /bin/bash

# Define FORTRAN90 compiler to use 
#  gfortran or Ifort

FC = gfortran
#FC = ifort

# Define FORTRAN90 compiler flags

FFLAGS = -c

# Additional flags for gfortran
ifeq ($(FC), gfortran)
	FFLAGS_TRACE = -Wall -C -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS := $(FFLAGS) $(FFLAGS_TRACE)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
endif

# Additional flags for ifort
ifeq ($(FC), ifort)
	FFLAGS_TRACE = -warn all -check all -traceback
	FFLAGS := $(FFLAGS) $(FFLAGS_TRACE)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
endif

#  Additional Optimization and Debug flags (Uncomment as required)
#	FFLAGS_OPTIM :=  -O3
#	FFLAGS_DEBUG := -g $(FFLAGS_TRACE)
#	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG) 
#	FFLAGS := $(FFLAGS) $(FFLAGS_OPTIM) 

# Define list of source files
# (Note: ordering is important because of dependencies)

SOURCES_TM1 =
SOURCES_TM1 +=   \
   $(RTSMIE_PATH)/RTSMie_parameters.f90          \
   $(RTSMIE_PATH)/RTSMie_distribution.f90        \
   $(RTSMIE_PATH)/RTSMie_sourcecode_2.f90        \
   $(RTSMIE_PATH)/RTSMie_master_TRIMODAL.f90      \
   MIE_MODEL.f90


# Define pattern rules for creating object files:

.SUFFIXES:
$(OBJ_PATH)/%.o : $(RTSMIE_PATH)/%.f90
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

# Define object files

F90SOURCES_TM1 := $(notdir $(filter %.f90, $(SOURCES_TM1)))
F90OBJECTS_TM1 := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(F90SOURCES_TM1)))

ALL:  MIE_MODEL.exe

MIE_MODEL.exe: $(F90OBJECTS_TM1)
	$(FC) $^ -o $@

.PHONY: clean
clean:
	rm -f *.o $(OBJ_FILES) *.mod $(MOD_FILES) *.log *.exe *~ Ice_SSprops_V2/mie3* *.res
