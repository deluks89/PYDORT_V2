# Build for copied Standard VLIDORT 2.8.3 and VLIDORT-AC sources.

FC = gfortran
.DEFAULT_GOAL := all

MOD_PATH = mod
OBJ_PATH = obj

DEF_PATH = vlidort_def
FO_PATH = vlidort_focode
REG_PATH = vlidort_main/regular
LIN_PATH = vlidort_main/linearized
SUPPORT_PATH = ../driver_support

FFLAGS = -c -O2
LDFLAGS =

ifeq ($(FC), gfortran)
	# Conservative Rayleigh LPS derivatives lose precision in binary64.
	# Promote EVERY compilation unit, including the bundled LAPACK and driver,
	# so all module interfaces use the same binary128 arithmetic.
	FFLAGS := $(FFLAGS) -Wall -fno-automatic -freal-8-real-16 -I$(MOD_PATH) -J$(MOD_PATH)
	LDFLAGS := $(LDFLAGS) -Wl,--stack,536870912
endif

ifeq ($(FC), ifort)
	FFLAGS := $(FFLAGS) -warn all -I$(MOD_PATH) -module $(MOD_PATH)
endif

.NOTPARALLEL:

REGULAR_SOURCES = \
	$(DEF_PATH)/vlidort_pars.f90 \
	$(DEF_PATH)/vlidort_inputs_def.f90 \
	$(DEF_PATH)/vlidort_outputs_def.f90 \
	$(DEF_PATH)/vlidort_sup_sleave_def.f90 \
	$(DEF_PATH)/vlidort_sup_brdf_def.f90 \
	$(DEF_PATH)/vlidort_sup_ss_def.f90 \
	$(DEF_PATH)/vlidort_sup_def.f90 \
	$(DEF_PATH)/vlidort_work_def.f90 \
	$(DEF_PATH)/vlidort_io_defs.f90 \
	$(FO_PATH)/FO_geometry_Generic.f90 \
	$(FO_PATH)/FO_WPgeometry_Routines.f90 \
	$(FO_PATH)/FO_DTWPgeometry_master.f90 \
	$(FO_PATH)/FO_SSWPgeometry_master.f90 \
	$(FO_PATH)/FO_VectorSS_Spherfuncs.f90 \
	$(FO_PATH)/FO_VectorSS_RTCalcs_I.f90 \
	$(FO_PATH)/FO_Thermal_RTCalcs_I.f90 \
	$(FO_PATH)/VFO_Master.f90 \
	$(REG_PATH)/lapack_tools.f90 \
	$(REG_PATH)/vlidort_aux.f90 \
	$(REG_PATH)/vlidort_inputs.f90 \
	$(REG_PATH)/vlidort_Taylor.f90 \
	$(REG_PATH)/vlidort_geometry.f90 \
	$(REG_PATH)/vlidort_getplanck.f90 \
	$(REG_PATH)/vlidort_solutions.f90 \
	$(REG_PATH)/vlidort_miscsetups.f90 \
	$(REG_PATH)/vlidort_bvproblem.f90 \
	$(REG_PATH)/vlidort_PostProcessing.f90 \
	$(REG_PATH)/vlidort_intensity.f90 \
	$(REG_PATH)/vlidort_converge.f90 \
	$(REG_PATH)/vlidort_thermalsup.f90 \
	$(REG_PATH)/vlidort_writemodules.f90 \
	$(REG_PATH)/vlidort_mediaprops.f90 \
	$(REG_PATH)/vlidort_pack.f90 \
	$(REG_PATH)/vlidort_unpack.f90 \
	$(REG_PATH)/vlidort_vfo_interface.f90 \
	$(REG_PATH)/vlidort_masters_V2p8p3.f90 \
	$(REG_PATH)/vlidort_scmaster_V2p8p3.f90 \
	$(REG_PATH)/vlidort_ac_scmaster_V2p8p3.f90

LINEARIZED_SOURCES = \
	$(DEF_PATH)/vlidort_lin_inputs_def.f90 \
	$(DEF_PATH)/vlidort_lin_outputs_def.f90 \
	$(DEF_PATH)/vlidort_lin_sup_sleave_def.f90 \
	$(DEF_PATH)/vlidort_lin_sup_brdf_def.f90 \
	$(DEF_PATH)/vlidort_lin_sup_ss_def.f90 \
	$(DEF_PATH)/vlidort_lin_sup_def.f90 \
	$(DEF_PATH)/vlidort_lin_work_def.f90 \
	$(DEF_PATH)/vlidort_lin_io_defs.f90 \
	$(FO_PATH)/FO_VectorSS_RTCalcs_ILCS.f90 \
	$(FO_PATH)/FO_Thermal_RTCalcs_ILCS.f90 \
	$(FO_PATH)/FO_VectorSS_RTCalcs_ILPS.f90 \
	$(FO_PATH)/FO_Thermal_RTCalcs_ILPS.f90 \
	$(FO_PATH)/VFO_LinMasters.f90 \
	$(LIN_PATH)/vlidort_lpc_solutions.f90 \
	$(LIN_PATH)/vlidort_lpc_bvproblem.f90 \
	$(LIN_PATH)/vlidort_la_miscsetups.f90 \
	$(LIN_PATH)/vlidort_l_inputs.f90 \
	$(LIN_PATH)/vlidort_l_pack.f90 \
	$(LIN_PATH)/vlidort_l_unpack.f90 \
	$(LIN_PATH)/vlidort_l_writemodules.f90 \
	$(LIN_PATH)/vlidort_l_thermalsup.f90 \
	$(LIN_PATH)/vlidort_lbbf_jacobians_vector.f90 \
	$(LIN_PATH)/vlidort_lc_solutions.f90 \
	$(LIN_PATH)/vlidort_lc_mediaprops.f90 \
	$(LIN_PATH)/vlidort_lc_miscsetups.f90 \
	$(LIN_PATH)/vlidort_lc_bvproblem.f90 \
	$(LIN_PATH)/vlidort_lc_pack.f90 \
	$(LIN_PATH)/vlidort_lc_unpack.f90 \
	$(LIN_PATH)/vlidort_lc_PostProcessing.f90 \
	$(LIN_PATH)/vlidort_lc_wfatmos.f90 \
	$(LIN_PATH)/vlidort_lcs_converge.f90 \
	$(LIN_PATH)/vlidort_lp_solutions.f90 \
	$(LIN_PATH)/vlidort_lp_mediaprops.f90 \
	$(LIN_PATH)/vlidort_lp_miscsetups.f90 \
	$(LIN_PATH)/vlidort_lp_bvproblem.f90 \
	$(LIN_PATH)/vlidort_lp_pack.f90 \
	$(LIN_PATH)/vlidort_lp_unpack.f90 \
	$(LIN_PATH)/vlidort_lp_PostProcessing.f90 \
	$(LIN_PATH)/vlidort_lp_wfatmos.f90 \
	$(LIN_PATH)/vlidort_lps_converge.f90 \
	$(LIN_PATH)/vlidort_ls_wfsleave.f90 \
	$(LIN_PATH)/vlidort_ls_wfsurface.f90 \
	$(LIN_PATH)/vlidort_vfo_lcs_interface.f90 \
	$(LIN_PATH)/vlidort_vfo_lps_interface.f90 \
	$(LIN_PATH)/vlidort_lcs_masters_V2p8p3.f90 \
	$(LIN_PATH)/vlidort_lps_masters_V2p8p3.f90 \
	$(LIN_PATH)/vlidort_lcs_scmaster_V2p8p3.f90 \
	$(LIN_PATH)/vlidort_lps_scmaster_V2p8p3.f90 \
	$(LIN_PATH)/vlidort_ac_lcs_scmaster_V2p8p3.f90 \
	$(LIN_PATH)/vlidort_ac_lps_scmaster_V2p8p3.f90

REGULAR_OBJECTS = $(patsubst %.f90,$(OBJ_PATH)/%.o,$(notdir $(REGULAR_SOURCES)))
LINEARIZED_OBJECTS = $(patsubst %.f90,$(OBJ_PATH)/%.o,$(notdir $(LINEARIZED_SOURCES)))
DRIVER_SOURCES = \
	$(SUPPORT_PATH)/DRIVER_pars.f90 \
	$(SUPPORT_PATH)/MainDriverUtility.f90 \
	$(SUPPORT_PATH)/ReadInputData.f90 \
	$(SUPPORT_PATH)/RayleighScattering.f90 \
	DRIVER.f90
DRIVER_OBJECTS = $(patsubst %.f90,$(OBJ_PATH)/%.o,$(notdir $(DRIVER_SOURCES)))
MOD_FILES = $(MOD_PATH)/*.mod
OBJ_FILES = $(OBJ_PATH)/*.o

.PHONY: all regular linearized driver clean

all: regular linearized driver

regular: $(REGULAR_OBJECTS)

linearized: regular $(LINEARIZED_OBJECTS)

driver: linearized $(DRIVER_OBJECTS)
	$(FC) $(OBJ_FILES) $(LDFLAGS) -o DRIVER.exe

# Dimension/type or build-option changes require a consistent full rebuild.
# For command-line FFLAGS overrides, use -B as well.
$(REGULAR_OBJECTS) $(LINEARIZED_OBJECTS) $(DRIVER_OBJECTS): makefile $(wildcard $(DEF_PATH)/*.f90) $(SUPPORT_PATH)/DRIVER_pars.f90

$(MOD_PATH):
	@if not exist "$(MOD_PATH)" mkdir "$(MOD_PATH)"

$(OBJ_PATH):
	@if not exist "$(OBJ_PATH)" mkdir "$(OBJ_PATH)"

$(OBJ_PATH)/%.o: $(DEF_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o: $(FO_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o: $(REG_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o: $(LIN_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o: $(SUPPORT_PATH)/%.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

$(OBJ_PATH)/%.o: %.f90 | $(MOD_PATH) $(OBJ_PATH)
	$(FC) $(FFLAGS) $< -o $@

clean:
	powershell -NoProfile -Command "Remove-Item -Force -ErrorAction SilentlyContinue -Path '*.o','*.mod','*.log','*.exe','*~','$(OBJ_FILES)','$(MOD_FILES)'"
