
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
! #              VQS_WRITE_STD_INPUT                            #
! #                                                             #
! ###############################################################

      MODULE vqs_writemodules_m

      PRIVATE
      PUBLIC :: VQS_WRITE_STD_INPUT

      CONTAINS

      SUBROUTINE VQS_WRITE_STD_INPUT ( VQS_FixIn, VQS_ModIn )

      USE MSSeg_PARS_m, Only : MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ, &
                               DWFC, DWFL, DWFI, DWFI1, DWFR, DWFR1, DWFR2, DWFR3, DWFR4, DWFR1_3
      USE VQS_Inputs_def_m

      IMPLICIT NONE

!  VQS I/O
!  =======

!  VQS input structures

      TYPE(VQS_Fixed_Inputs), INTENT (IN)       :: VQS_FixIn
      TYPE(VQS_Modified_Inputs), INTENT (INOUT) :: VQS_ModIn

!  Local variables

      INTEGER :: NLAYERS
      INTEGER :: NSTOKES
      INTEGER :: NGREEK_MOMENTS_INPUT

      INTEGER :: OUTUNIT
      INTEGER :: GEO,K,LAY,MOM,NGREEK,SS
      INTEGER :: GMASK(8)
      CHARACTER (LEN=4) :: AI

      DATA GMASK / 1, 2, 5, 6, 11, 12, 15, 16 /

!  Open output file

      OUTUNIT = 101
      !OPEN (OUTUNIT,file = 'VQS_WRITE_STD_INPUT.dbg',status = 'replace')
      WRITE(AI,'(i4.4)') VQS_FixIn%Optical%ATMOS_INDEX
      OPEN (OUTUNIT,file = 'Debug\VQS_WRITE_STD_INPUT_'//AI//'.dbg',status = 'replace')

!  Define local variables

      NLAYERS = VQS_FixIn%Cont%NLAYERS
      NSTOKES = VQS_FixIn%Cont%NSTOKES
      NGREEK_MOMENTS_INPUT = VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT

!  Only output non-zero entries of GreekMat

      NGREEK = 1
      IF ( NSTOKES.eq.3 ) NGREEK = 5
      IF ( NSTOKES.eq.4 ) NGREEK = 8

!  Write all input to file

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '-----------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Fixed'
      WRITE(OUTUNIT,'(A)') '-----------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_MSRT                = ',VQS_FixIn%Bool%DO_MSRT
      WRITE(OUTUNIT,DWFL)  'DO_CALCULATE_FLUXES    = ',VQS_FixIn%Bool%DO_CALCULATE_FLUXES
      WRITE(OUTUNIT,DWFL)  'DO_SUNLIGHT            = ',VQS_FixIn%Bool%DO_SUNLIGHT
      WRITE(OUTUNIT,DWFL)  'DO_PLANE_PARALLEL      = ',VQS_FixIn%Bool%DO_PLANE_PARALLEL
      WRITE(OUTUNIT,DWFL)  'DO_LAMBERTIAN_SURFACE  = ',VQS_FixIn%Bool%DO_LAMBERTIAN_SURFACE

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_FO_GEOMETRY   = ',VQS_FixIn%Bool%DO_DEBUG_FO_GEOMETRY
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_FO_INPUT      = ',VQS_FixIn%Bool%DO_DEBUG_FO_INPUT
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_MS_INPUT      = ',VQS_FixIn%Bool%DO_DEBUG_MS_INPUT
      WRITE(OUTUNIT,DWFL)  'DO_DEBUG_FLUXES        = ',VQS_FixIn%Bool%DO_DEBUG_FLUXES

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'METHOD_INDEX     = ',VQS_FixIn%Cont%METHOD_INDEX
      WRITE(OUTUNIT,DWFI)  'NLAYERS          = ',VQS_FixIn%Cont%NLAYERS
      WRITE(OUTUNIT,DWFI)  'NSTOKES          = ',VQS_FixIn%Cont%NSTOKES
      WRITE(OUTUNIT,DWFI)  'NSTREAMS         = ',VQS_FixIn%Cont%NSTREAMS
      WRITE(OUTUNIT,DWFI)  'NFINELAYERS      = ',VQS_FixIn%Cont%NFINE_INPUT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DOCRIT           = ',VQS_FixIn%Cont%DOCRIT
      WRITE(OUTUNIT,DWFR)  'ACRIT            = ',VQS_FixIn%Cont%ACRIT

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'FLUX_FACTOR      = ',VQS_FixIn%Sunrays%FLUX_FACTOR

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'HEIGHT_CAMERA    = ',VQS_FixIn%UserVal%HEIGHT_CAMERA
      WRITE(OUTUNIT,DWFR)  'THETA_CAMERA     = ',VQS_FixIn%UserVal%THETA_CAMERA
      WRITE(OUTUNIT,DWFR)  'RELAZ_CAMERA     = ',VQS_FixIn%UserVal%RELAZ_CAMERA
      WRITE(OUTUNIT,DWFR)  'ALPHA_CAMERA     = ',VQS_FixIn%UserVal%ALPHA_CAMERA

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFC)  'GEO_FILENAME     = ',trim(VQS_FixIn%UserVal%GEO_FILENAME)

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'PLANET_RADIUS    = ',VQS_FixIn%Geophys%PLANET_RADIUS

      WRITE(OUTUNIT,*)
      DO LAY=0,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' HEIGHT_GRID(LAY)      = ',VQS_FixIn%Geophys%HEIGHT_GRID(LAY)
      END DO

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,&
          ' DELTAU_VERT_INPUT(LAY) = ',VQS_FixIn%Optical%DELTAU_VERT_INPUT(LAY)
      END DO

      DO K=1,NGREEK
        SS = GMASK(K) ; WRITE(OUTUNIT,*)
        DO LAY=1,NLAYERS
          DO MOM=0,NGREEK_MOMENTS_INPUT
            WRITE(OUTUNIT,DWFR3)  'SS = ',SS,' LAY = ',LAY,' MOM = ',MOM,&
              ' GREEKMAT_TOTAL_INPUT(MOM,LAY,SS) = ',VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT(MOM,LAY,SS)
          END DO
        END DO
      END DO

      !DO K=1,6
      !  WRITE(OUTUNIT,*)
      !  DO LAY=1,NLAYERS
      !    WRITE(OUTUNIT,DWFR3)  'K = ',K,' LAY = ',LAY,&
      !     ' FMATRIX_UP(LAY,K) = ',VQS_FixIn%Optical%FMATRIX_UP(LAY,K)
      !  END DO
      !END DO

      !DO K=1,6
      !  WRITE(OUTUNIT,*)
      !  DO LAY=1,NLAYERS
      !    WRITE(OUTUNIT,DWFR3)  'K = ',K,' LAY = ',LAY,&
      !      ' FMATRIX_DN(LAY,K) = ',VQS_FixIn%Optical%FMATRIX_DN(LAY,K)
      !  END DO
      !END DO

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFR)  'LAMBERTIAN_ALBEDO = ',VQS_FixIn%Optical%LAMBERTIAN_ALBEDO

      WRITE(OUTUNIT,*)
      !WRITE(OUTUNIT,DWFR)  'ATMOS_WAVELENGTH  = ',VQS_FixIn%Optical%ATMOS_WAVELENGTH
      WRITE(OUTUNIT,DWFI)  'ATMOS_INDEX       = ',VQS_FixIn%Optical%ATMOS_INDEX

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,'(A)') '--------------------------'
      WRITE(OUTUNIT,'(A)') 'Standard Inputs - Variable'
      WRITE(OUTUNIT,'(A)') '--------------------------'

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFL)  'DO_RAYLEIGH_ONLY        = ',VQS_ModIn%MBool%DO_RAYLEIGH_ONLY
      WRITE(OUTUNIT,DWFL)  'DO_DELTAM_SCALING       = ',VQS_ModIn%MBool%DO_DELTAM_SCALING

      WRITE(OUTUNIT,*)
      WRITE(OUTUNIT,DWFI)  'NGREEK_MOMENTS_INPUT    = ',VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT

      WRITE(OUTUNIT,*)
      DO LAY=1,NLAYERS
        WRITE(OUTUNIT,DWFR1)  'LAY = ',LAY,' OMEGA_TOTAL_INPUT(LAY) = ',VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(LAY)
      END DO

!  Close output file

      CLOSE (OUTUNIT)

      END SUBROUTINE VQS_WRITE_STD_INPUT

      END MODULE vqs_writemodules_m
