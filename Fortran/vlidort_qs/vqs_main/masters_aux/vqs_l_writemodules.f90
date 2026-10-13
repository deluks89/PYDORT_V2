
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

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #              VLIDORT_WRITE_LIN_INPUT                        #
! #                                                             #
! ###############################################################

      MODULE vqs_l_writemodules_m

      PRIVATE
      PUBLIC :: VQS_WRITE_LIN_INPUT

      CONTAINS

      SUBROUTINE VQS_WRITE_LIN_INPUT ( VQS_FixIn, VQS_ModIn, VQS_LinFixIn )

      USE MSSeg_PARS_m, Only : MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ, &
                               MAX_ATMOSWFS, &
                               DWFL, DWFL1, DWFI, DWFI1, DWFR, DWFR2, DWFR4
      USE VQS_Inputs_def_m
      USE VQS_LinInputs_def_m

      IMPLICIT NONE

!  VQS I/O
!  =======

!  VQS input structures

      TYPE(VQS_Fixed_Inputs), INTENT (IN)    :: VQS_FixIn
      TYPE(VQS_Modified_Inputs), INTENT (IN) :: VQS_ModIn

!  VQS linearized input structures

      TYPE(VQS_Fixed_LinInputs), INTENT (IN) :: VQS_LinFixIn

!  Local variables

      INTEGER :: NLAYERS
      INTEGER :: NSTOKES
      INTEGER :: NGREEK_MOMENTS_INPUT

      INTEGER :: OUTUNIT
      INTEGER :: GEO,K,LAY,MOM,NGREEK,SS,WF
      INTEGER :: GMASK(8)
      INTEGER :: N_WFS
      CHARACTER (LEN=4) :: AI
      !CHARACTER (LEN=9) :: DWFC1S

      DATA GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

!  Open output file

      OUTUNIT = 111
      !OPEN (OUTUNIT,file = 'VQS_WRITE_LIN_INPUT.dbg',status = 'replace')
      WRITE(AI,'(i4.4)') VQS_FixIn%Optical%ATMOS_INDEX
      OPEN (OUTUNIT,file = 'Debug\VQS_WRITE_LIN_INPUT_'//AI//'.dbg',status = 'replace')

!  Define local variables

      NLAYERS = VQS_FixIn%Cont%NLAYERS
      NSTOKES = VQS_FixIn%Cont%NSTOKES
      NGREEK_MOMENTS_INPUT = VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT

!  9/21/23. Add clause for N_COLUMN_WFS

      if ( VQS_LinFixIn%Cont%N_PROFILE_WFS > 0 ) then
         N_WFS = VQS_LinFixIn%Cont%N_PROFILE_WFS
      else  if ( VQS_LinFixIn%Cont%N_COLUMN_WFS > 0 ) then
         N_WFS = VQS_LinFixIn%Cont%N_COLUMN_WFS
      endif

!  Only output non-zero entries of GreekMat

      NGREEK = 1
      IF ( NSTOKES.eq.3 ) NGREEK = 5
      IF ( NSTOKES.eq.4 ) NGREEK = 8

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-------------------------'
      WRITE(OUTUNIT,'(A)') 'Linearized Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-------------------------'

!  9/21/23. Add clause for N_COLUMN_WFS

      WRITE(OUTUNIT,*)
      if ( VQS_LinFixIn%Cont%N_PROFILE_WFS > 0 ) then
        WRITE(OUTUNIT,DWFI)  'N_PROFILE_WFS      = ',VQS_LinFixIn%Cont%N_PROFILE_WFS
      else  if ( VQS_LinFixIn%Cont%N_COLUMN_WFS > 0 ) then
        WRITE(OUTUNIT,DWFI)  'NCOLUMN_WFS        = ',VQS_LinFixIn%Cont%N_COLUMN_WFS
      endif

      !WRITE(OUTUNIT,*)
      !DO LAY=1,NLAYERS
      !  WRITE(OUTUNIT,DWFL1)  'LAY = ',LAY,&
      !    ' LAYER_VARY_FLAG(LAY)   = ',VQS_LinFixIn%Cont%LAYER_VARY_FLAG(LAY)
      !END DO

      !WRITE(OUTUNIT,*)
      !DO LAY=1,NLAYERS
      !  WRITE(OUTUNIT,DWFI1)  'LAY = ',LAY,&
      !    ' LAYER_VARY_NUMBER(LAY) = ',VQS_LinFixIn%Cont%LAYER_VARY_NUMBER(LAY)
      !END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' WF = ',WF,&
            ' L_DELTAU_VERT_INPUT(WF,LAY) = ',VQS_LinFixIn%Optical%L_DELTAU_VERT_INPUT(WF,LAY)
        END DO
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        DO WF=1,N_WFS
          WRITE(OUTUNIT,DWFR2)  'LAY = ',LAY,' WF = ',WF,&
            ' L_OMEGA_TOTAL_INPUT(WF,LAY) = ',VQS_LinFixIn%Optical%L_OMEGA_TOTAL_INPUT(WF,LAY)
        END DO
      END DO

      DO K=1,NGREEK
        SS = GMASK(K) ; WRITE(OUTUNIT,*)
        DO LAY=1,NLAYERS  
          DO MOM=0,NGREEK_MOMENTS_INPUT
            DO WF=1,N_WFS
              WRITE(OUTUNIT,DWFR4) &
                'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,' WF = ',WF,&
                ' L_GREEKMAT_TOTAL_INPUT(WF,MOM,LAY,SS) = ',VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT(WF,MOM,LAY,SS)
            END DO
          END DO
        END DO
      END DO

      !DO K=1,6
      !  WRITE(OUTUNIT,*)
      !  DO LAY=1,NLAYERS
      !    DO WF=1,N_WFS
      !     WRITE(OUTUNIT,DWFR4) 'K = ',K,' LAY = ',LAY,' WF = ',WF, &
      !        ' L_FMATRIX_UP(WF,LAY,K) = ',VQS_LinFixIn%Optical%L_FMATRIX_UP(WF,LAY,K)
      !    END DO
      !  END DO
      !END DO

      !DO K=1,6
      !  WRITE(OUTUNIT,*)
      !  DO LAY=1,NLAYERS
      !    DO WF=1,N_WFS
      !      WRITE(OUTUNIT,DWFR4) 'K = ',K,' LAY = ',LAY,' WF = ',WF, &
      !        ' L_FMATRIX_DN(WF,LAY,K) = ',VQS_LinFixIn%Optical%L_FMATRIX_DN(WF,LAY,K)
      !    END DO
      !  END DO
      !END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VQS_WRITE_LIN_INPUT

      END MODULE vqs_l_writemodules_m
