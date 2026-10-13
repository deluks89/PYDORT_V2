
! ###############################################################
! #                                                             #
! #                     VLIDORT_QS_1p4                          #
! #                      ("VQS_1p4")                            #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                   (Quasi - Spherical)                       #
! #                    -       -                                #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VQS_1p4                                   #
! #  Release Date :   30 June 2023                              #
! #                                                             #
! #  Previous VQS Versions under Standard GPL 3.0:              #
! #  --------------------------------------------               #
! #                                                             #
! #      1.0   F90, limited release May    2012                 #
! #      1.1   F90, limited release June   2015                 #
! #      1.2   F90, limited release March  2016                 #
! #      1.3   F90, limited release March  2020                 #
! #                                                             #
! #  Features Summary of Recent VQS Versions:                   #
! #  ---------------------------------------                    #
! #                                                             #
! #      FO + spher-enhanced MS for Stokes vectors     (1.0)    #
! #      Profile Jacobians added                       (1.1)    #
! #      MS code upgraded to VLIDORT 2.7 standards     (1.2)    #
! #      FO/MS codes upgraded to VLIDORT 2.8 standards (1.3)    #
! #                                                             #
! #  Features Summary of This VQS Version                       #
! #  ------------------------------------                       #
! #                                                             #
! #  1.4, released 30 June 2023.                                #
! #    ==> Formal VQS package release                           #
! #    ==> Formal VQS regular and linearized prof. master subs  #
! #    ==> New Fortran Type Structures for I/O                  #
! #    ==> MS computation performance enhancements              #
! #    ==> Source code made thread-safe for OpenMP              #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 1.4 of the VQS_1p4 software library.            #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VQS Copyright (c) 2012-2023.                               #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VQS_1p4 ( Version 1.4 )                    #
! #                                                                 #
! # VQS_1p4 is free software: you can redistribute it               #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VQS_1p4 is distributed in the hope that it will be              #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VQS_1p4        #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

      MODULE VQS_LinInputs_def_m

!  This Module contains the following VQS input Structures, with Intents :

!        VQS_Fixed_LinControl    nested in VQS_Fixed_LinInputs
!        VQS_Fixed_LinOptical    nested in VQS_Fixed_LinInputs
!         VQS_Fixed_LinInputs    Intent(In)

!     VQS_Modified_LinControl    nested in VQS_Modified_LinInputs
!      VQS_Modified_LinInputs    Intent(InOut)

      USE MSSeg_PARS_m, Only : fpk, MAX_ATMOSWFS, MAX_GEOMETRIES, MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_LinControl

!  Control for atmospheric linearizations, layer by layer

      LOGICAL, dimension(MAXLAYERS)  :: LAYER_VARY_FLAG
      !INTEGER, dimension(MAXLAYERS)  :: LAYER_VARY_NUMBER

!  Total number of profile Jacobians
	  
	   INTEGER  :: N_PROFILE_WFS

!  9/21/23. Total number of column Jacobians

      INTEGER  :: N_COLUMN_WFS

      END TYPE VQS_Fixed_LinControl

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_LinOptical

!  Optical property linearizations
!  Layer linearization (bulk property variation) input
!  Layer linearization (phase function variation) input

      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: L_DELTAU_VERT_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, MAXLAYERS) :: L_OMEGA_TOTAL_INPUT
      REAL(fpk), dimension(MAX_ATMOSWFS, 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ) :: L_GREEKMAT_TOTAL_INPUT

!  Linearized Fmatrix inputs for FO calculations

      !REAL(fpk), dimension( MAX_ATMOSWFS, MAXLAYERS, 6 ) :: L_FMATRIX_UP 
      !REAL(fpk), dimension( MAX_ATMOSWFS, MAXLAYERS, 6 ) :: L_FMATRIX_DN 

      END TYPE VQS_Fixed_LinOptical

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_LinInputs

      TYPE(VQS_Fixed_LinControl)    :: Cont
      TYPE(VQS_Fixed_LinOptical)    :: Optical

      END TYPE VQS_Fixed_LinInputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VQS_Fixed_LinControl,    &
                VQS_Fixed_LinOptical,    &
                VQS_Fixed_LinInputs

   END MODULE VQS_LinInputs_def_m

