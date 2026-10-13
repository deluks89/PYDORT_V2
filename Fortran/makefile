# Common makefile for the buildable Fortran pieces in VLIDORT_RT.

FC = gfortran

.PHONY: all rt qs mie standard ac clean

all: rt mie

rt: qs standard

qs:
	$(MAKE) -C vlidort_qs FC=$(FC) set2

mie:
	$(MAKE) -C mie_model FC=$(FC) ALL

standard:
	$(MAKE) -C vlidort_standard FC=$(FC) all

ac: standard

clean:
	powershell -NoProfile -Command "Remove-Item -Force -ErrorAction SilentlyContinue -Path 'vlidort_qs/*.o','vlidort_qs/obj/*.o','vlidort_qs/*.mod','vlidort_qs/mod/*.mod','vlidort_qs/Debug_FOMS*.*','vlidort_qs/*.log','vlidort_qs/*.exe','vlidort_qs/*~','vlidort_qs/*.dat'"
	-$(MAKE) -C vlidort_standard clean
	powershell -NoProfile -Command "Remove-Item -Force -ErrorAction SilentlyContinue -Path 'mie_model/*.o','mie_model/obj/*.o','mie_model/*.mod','mie_model/mod/*.mod','mie_model/*.log','mie_model/*.exe','mie_model/*~','mie_model/*.res'"
