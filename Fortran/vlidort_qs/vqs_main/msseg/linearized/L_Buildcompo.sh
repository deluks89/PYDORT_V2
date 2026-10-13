#rm *.o *.mod
#exit

# HOLDING PATTERNS

#gfortran -c -Wall MSSeg_LA_Fourier.f90
#gfortran -c -Wall MSSeg_LC_Fourier.f90
#gfortran -c -Wall MSSeg_LP_Fourier.f90
#exit

#  General and Types

gfortran -c -Wall ../../../vqs_def/FOMS_pars.f90
gfortran -c -Wall ../../../vqs_def/MSSeg_pars.f90
gfortran -c -Wall ../regular/MSSeg_aux.f90
gfortran -c -Wall ../regular/lapack_tools.f90
gfortran -c -Wall ../regular/MSSeg_HPTypes.f90
gfortran -c -Wall MSSeg_L_HPTypes.f90
#exit

#  Homogeneous, regular

gfortran -c -Wall ../regular/MSSeg_Homog_miscsetups.f90
gfortran -c -Wall ../regular/MSSeg_Homog_solutions.f90
gfortran -c -Wall ../regular/MSSeg_Homog_bvproblem.f90
gfortran -c -Wall ../regular/MSSeg_Homog_Fourier.f90
#exit

# Particular, regular

gfortran -c -Wall ../regular/MSSeg_Partic_miscsetups.f90
gfortran -c -Wall ../regular/MSSeg_Partic_solutions.f90
gfortran -c -Wall ../regular/MSSeg_Partic_bvproblem.f90
gfortran -c -Wall ../regular/MSSeg_Partic_PostProcess.f90
gfortran -c -Wall ../regular/MSSeg_Partic_Fourier.f90
#exit

#  linearized

gfortran -c -Wall MSSeg_LA_miscsetups.f90
gfortran -c -Wall MSSeg_LA_solutions.f90
gfortran -c -Wall MSSeg_LP_miscsetups.f90
gfortran -c -Wall MSSeg_LP_bvproblem.f90
gfortran -c -Wall MSSeg_LP_solutions.f90
gfortran -c -Wall MSSeg_LP_PostProcess.f90
gfortran -c -Wall MSSeg_LC_miscsetups.f90
gfortran -c -Wall MSSeg_LC_bvproblem.f90
gfortran -c -Wall MSSeg_LC_solutions.f90
gfortran -c -Wall MSSeg_LC_PostProcess.f90
exit

gfortran -c -Wall MSSeg_LA_Fourier.f90
gfortran -c -Wall MSSeg_LC_Fourier.f90
gfortran -c -Wall MSSeg_LP_Fourier.f90

exit


