# MAKEFILE for VQS Mark3 Experimental Code

#  ================================================================
#  9/21/23. Upgrade. Introduce new test for LC software (LC93AerBB)
#           -- Original tests unchanged, but now ACTOOL test is #4
#           -- In this makefile, introduce LC and LCS sources
#  ================================================================

# Define some variables

.DEFAULT_GOAL := set2

UTIL_PATH           = util
VSUP_PATH           = vsup
VQS_DEF_PATH        = vqs_def
VQS_MAIN_PATH       = vqs_main
VQS_MASTER_PATH     = $(VQS_MAIN_PATH)/masters
VQS_MASTER_AUX_PATH = $(VQS_MAIN_PATH)/masters_aux
VQS_FOMS_PATH       = $(VQS_MAIN_PATH)/foms
VQS_MSSEG_PATH1     = $(VQS_MAIN_PATH)/msseg/regular
VQS_MSSEG_PATH2     = $(VQS_MAIN_PATH)/msseg/linearized

SUPPORT_PATH  = ../driver_support

MOD_PATH = mod
OBJ_PATH = obj

MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

# Define default shell make will use

SHELL = cmd.exe

# Define FORTRAN90 compiler to use (can be defined here, but usually defined on makefile command line)

#  Intel
#FC = ifort

#  gfortran
FC = gfortran

#  g95
#FC = g95

#  NAG
#FC = f95

#
# Define FORTRAN90 compiler flags
#

FFLAGS = -c

# Additional flags for Intel
ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -module $(MOD_PATH)
	FFLAGS_DEBUG = -g -warn all -check all -traceback
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

# Additional flags for gfortran
ifeq ($(FC), gfortran)
        FFLAGS := $(FFLAGS) -std=legacy -mcmodel=large -O3
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -J$(MOD_PATH)
	FFLAGS_DEBUG  = -g -C -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_Catch  =       -Wall -fbounds-check -frange-check -ffpe-trap=invalid,zero,overflow
	FFLAGS_OPT    = -O3
	FFLAGS_OPENMP = -fopenmp -frecursive
endif

# Additional flags for g95
ifeq ($(FC), g95)
#      older g95
	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check
	FFLAGS_DEBUG = -g -Wall -Wuninitialized -fbounds-check
	FFLAGS_OPT = -O3
#      g95 v0.92
#	FFLAGS := $(FFLAGS) -I$(MOD_PATH) -fmod=$(MOD_PATH)
#	FFLAGS_DEBUG = -g -Wall -fbounds-check -ftrace=full
endif

# Additional flags for NAG
#ifeq ($(FC), f95)
ifeq ($(FC), $(filter $(FC), f95 nagfor))
	FFLAGS := $(FFLAGS) -mdir $(MOD_PATH) -I$(MOD_PATH)
	FFLAGS_DEBUG = -w=obs -w=unused -C=all -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -C=array -C=undefined -gline
#	FFLAGS_DEBUG = -w=obs -w=unused -C=all -C=undefined -gline
#	FFLAGS_DEBUG = -g -C=all -C=undefined -gline -mtrace=all -nan
#	FFLAGS_DEBUG = -w=obs -w=unused -g -C=all -C=undefined -gline -mtrace=all -nan
	FFLAGS_OPT = -O3
        FFLAGS_OPENMP = -openmp
endif

#  for compilation detailed
#	FFLAGS := $(FFLAGS) $(FFLAGS_Catch)

# For debug build, use "make DEBUG=t"
ifeq ($(DEBUG), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_DEBUG)
endif

# For optimized build, use "make OPT=t"
ifeq ($(OPT), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPT)
endif

# For parallel build using OpenMP, use "make OPENMP=t"
ifeq ($(OPENMP), t)
	FFLAGS := $(FFLAGS) $(FFLAGS_OPENMP)
endif

.SUFFIXES:

#
#  Define list of source files
#  (Note: ordering is important because of dependencies)

#  VQS Sources
#  -----------

#  Regular MSSeg

SOURCES_MSSEG =
SOURCES_MSSEG += \
   $(VQS_DEF_PATH)/MSSeg_pars.f90	           \
   $(VQS_MSSEG_PATH1)/lapack_tools.f90		   \
   $(VQS_MSSEG_PATH1)/MSSeg_aux.f90		   \
   $(VQS_MSSEG_PATH1)/MSSeg_HPTypes.f90		   \
   $(VQS_MSSEG_PATH1)/MSSeg_Homog_miscsetups.f90   \
   $(VQS_MSSEG_PATH1)/MSSeg_Homog_solutions.f90	   \
   $(VQS_MSSEG_PATH1)/MSSeg_Homog_bvproblem.f90	   \
   $(VQS_MSSEG_PATH1)/MSSeg_Homog_Fourier.f90	   \
   $(VQS_MSSEG_PATH1)/MSSeg_Partic_miscsetups.f90  \
   $(VQS_MSSEG_PATH1)/MSSeg_Partic_solutions.f90   \
   $(VQS_MSSEG_PATH1)/MSSeg_Partic_bvproblem.f90   \
   $(VQS_MSSEG_PATH1)/MSSeg_Partic_PostProcess.f90 \
   $(VQS_MSSEG_PATH1)/MSSeg_Partic_Fourier.f90

#  Linearized MSSeg

SOURCES_MSSEG_LA = 
SOURCES_MSSEG_LA +=   \
   $(SOURCES_MSSEG) \
   $(VQS_MSSEG_PATH2)/MSSeg_L_HPTypes.f90	\
   $(VQS_MSSEG_PATH2)/MSSeg_LA_miscsetups.f90	\
   $(VQS_MSSEG_PATH2)/MSSeg_LA_solutions.f90	\
   $(VQS_MSSEG_PATH2)/MSSeg_LA_Fourier.f90

SOURCES_MSSEG_LC = 
SOURCES_MSSEG_LC +=   \
   $(SOURCES_MSSEG_LA) \
   $(VQS_MSSEG_PATH2)/MSSeg_LC_miscsetups.f90  \
   $(VQS_MSSEG_PATH2)/MSSeg_LC_solutions.f90   \
   $(VQS_MSSEG_PATH2)/MSSeg_LC_bvproblem.f90   \
   $(VQS_MSSEG_PATH2)/MSSeg_LC_PostProcess.f90 \
   $(VQS_MSSEG_PATH2)/MSSeg_LC_Fourier.f90  

SOURCES_MSSEG_LP = 
SOURCES_MSSEG_LP +=   \
   $(SOURCES_MSSEG_LA) \
   $(VQS_MSSEG_PATH2)/MSSeg_LP_miscsetups.f90  \
   $(VQS_MSSEG_PATH2)/MSSeg_LP_solutions.f90   \
   $(VQS_MSSEG_PATH2)/MSSeg_LP_bvproblem.f90   \
   $(VQS_MSSEG_PATH2)/MSSeg_LP_PostProcess.f90 \
   $(VQS_MSSEG_PATH2)/MSSeg_LP_Fourier.f90  

#  FOMS General

SOURCES_FOMS =
SOURCES_FOMS += \
   $(VQS_DEF_PATH)/FOMS_pars.f90 \
   $(VQS_FOMS_PATH)/FOMS_Geometry_Routines.f90 \
   $(VQS_FOMS_PATH)/FOMS_Geometry_Masters.f90  \
   $(VQS_FOMS_PATH)/FOMS_VectorSpherfuncs.f90

#  FOMS Primary

SOURCES_VQS_VectorRT =
SOURCES_VQS_VectorRT += \
   $(SOURCES_MSSEG) \
   $(SOURCES_FOMS) \
   $(VQS_FOMS_PATH)/FOMS_FO_VectorRT.f90         \
   $(VQS_FOMS_PATH)/FOMS_MS_VectorRT_Mk2_Alt.f90 \
   $(VQS_FOMS_PATH)/FOMS_MS_VectorRT_Mk3_Exp.f90

SOURCES_VQS_LP_VectorRT =
SOURCES_VQS_LP_VectorRT += \
   $(SOURCES_VQS_VectorRT) \
   $(SOURCES_MSSEG_LP) \
   $(VQS_FOMS_PATH)/FOMS_FO_LP_VectorRT.f90  	\
   $(VQS_FOMS_PATH)/FOMS_MS_LP_VectorRT_Mk2_Alt.f90 \
   $(VQS_FOMS_PATH)/FOMS_MS_LP_VectorRT_Mk3_Exp.f90

SOURCES_VQS_LC_VectorRT =
SOURCES_VQS_LC_VectorRT += \
   $(SOURCES_VQS_VectorRT) \
   $(SOURCES_MSSEG_LC) \
   $(VQS_FOMS_PATH)/FOMS_FO_LC_VectorRT.f90  	\
   $(VQS_FOMS_PATH)/FOMS_MS_LC_VectorRT_Mk2_Alt.f90 \
   $(VQS_FOMS_PATH)/FOMS_MS_LC_VectorRT_Mk3_Exp.f90

#  Source Sets (Set #2)

SOURCES =
SOURCES += \
   $(SOURCES_VQS_VectorRT) \
   $(VQS_DEF_PATH)/vqs_inputs_def.f90	   \
   $(VQS_DEF_PATH)/vqs_outputs_def.f90	   \
   $(VQS_DEF_PATH)/vqs_io_defs.f90	   \
   $(VQS_MASTER_AUX_PATH)/vqs_aux.f90          \
   $(VQS_MASTER_AUX_PATH)/vqs_inputs.f90       \
   $(VQS_MASTER_AUX_PATH)/vqs_writemodules.f90 \
   $(VQS_MASTER_PATH)/vqs_masters_V1p4.f90

LG_SOURCES =
LG_SOURCES += \
   $(SOURCES) \
   $(VQS_DEF_PATH)/vqs_lin_inputs_def.f90   \
   $(VQS_DEF_PATH)/vqs_lin_outputs_def.f90  \
   $(VQS_DEF_PATH)/vqs_lin_io_defs.f90	     \
   $(VQS_MASTER_AUX_PATH)/vqs_l_inputs.f90       \
   $(VQS_MASTER_AUX_PATH)/vqs_l_writemodules.f90

LCLP_SOURCES =
LCLP_SOURCES += \
   $(SOURCES_VQS_LC_VectorRT) \
   $(SOURCES_VQS_LP_VectorRT) \
   $(LG_SOURCES) \
   $(VQS_MASTER_PATH)/vqs_lc_masters_V1p4.f90 \
   $(VQS_MASTER_PATH)/vqs_lp_masters_V1p4.f90

#  MAIN DRIVER SOURCES
#  -------------------

DRIVER_SOURCES = 
DRIVER_SOURCES += \
   $(SOURCES) \
   $(LCLP_SOURCES) \
   $(SUPPORT_PATH)/DRIVER_pars.f90 \
   $(SUPPORT_PATH)/MainDriverUtility.f90 \
   $(SUPPORT_PATH)/ReadInputData.f90 \
   $(SUPPORT_PATH)/RayleighScattering.f90 \
   DRIVER.f90


# Define pattern rules for creating object files:

# Use the same automatic directory setup as the Standard and Mie builds.
$(MOD_PATH):
	@if not exist "$(MOD_PATH)" mkdir "$(MOD_PATH)"

$(OBJ_PATH):
	@if not exist "$(OBJ_PATH)" mkdir "$(OBJ_PATH)"

$(OBJ_PATH)/%.o : $(SUPPORT_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_DEF_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_FOMS_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_MSSEG_PATH1)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_MSSEG_PATH2)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_MASTER_AUX_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : $(VQS_MASTER_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@
$(OBJ_PATH)/%.o : %.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@


# Define object files
# -------------------

DRIVER_SOURCE_NAMES := $(notdir $(filter %.f90, $(DRIVER_SOURCES)))
DRIVER_OBJECTS := $(patsubst %.f90, %.o, $(addprefix $(OBJ_PATH)/, $(DRIVER_SOURCE_NAMES)))

set2: DRIVER.exe

DRIVER.exe: $(DRIVER_OBJECTS)
	$(FC) $^ -static -o $@

.PHONY: clean
clean:
	powershell -NoProfile -Command "Remove-Item -Force -ErrorAction SilentlyContinue -Path '*.o','*.mod','Debug_FOMS*.*','*.log','*.exe','*~','*.dat','$(OBJ_FILES)','$(MOD_FILES)'"
