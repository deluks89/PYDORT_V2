#rm *.o *.mod
#exit
#  general
 gfortran -c -Wall ../../../vqs_def/MSSeg_pars.f90
 gfortran -c -Wall MSSeg_aux.f90
 gfortran -c -Wall lapack_tools.f90
exit
#  Type structures
 gfortran -c -Wall MSSeg_HPTypes.f90
exit

#  Homogeneous

gfortran -c -Wall MSSeg_Homog_miscsetups.f90
gfortran -c -Wall MSSeg_Homog_solutions.f90
gfortran -c -Wall MSSeg_Homog_bvproblem.f90
gfortran -c -Wall MSSeg_Homog_Fourier.f90
exit

# Particular

gfortran -c -Wall MSSeg_Partic_miscsetups.f90
gfortran -c -Wall MSSeg_Partic_solutions.f90
gfortran -c -Wall MSSeg_Partic_bvproblem.f90
gfortran -c -Wall MSSeg_Partic_PostProcess.f90
gfortran -c -Wall MSSeg_Partic_Fourier.f90

exit
