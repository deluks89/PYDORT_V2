rm *.o *.mod
exit

gfortran -c -Wall ../../vlidort_def/vlidort_pars.f90
#exit

gfortran -c -Wall vfzmat_Numerical_Mk2.f90
gfortran -c -Wall vfzmat_Rotation.f90
gfortran -c -Wall vfzmat_DevelopCoeffs.f90
gfortran -c -Wall vfzmat_UseCoeffs.f90
gfortran -c -Wall vfzmat_PhasMat.f90
gfortran -c -Wall vfzmat_Rayleigh_Master.f90
gfortran -c -Wall vfzmat_General_Master.f90
#exit

gfortran -c -Wall vfzmat_L_DevelopCoeffs.f90
gfortran -c -Wall vfzmat_L_UseCoeffs.f90
gfortran -c -Wall vfzmat_L_PhasMat.f90
gfortran -c -Wall vfzmat_L_General_Master.f90
exit
