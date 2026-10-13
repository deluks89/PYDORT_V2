
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

      MODULE VQS_Outputs_def_m

!  This module contains the following structures:

!                 VQS_Main_Outputs  nested in VQS_Outputs
!           VQS_Exception_Handling  nested in VQS_Outputs
!     VQS_Input_Exception_Handling  Intent(Out) from Input settings
!                      VQS_Outputs  Intent(Out)

!  4/24/24. Necessary to introduce the MAXSEGMENTS parameter from FOMS_pars.f90
!     -- Reason: want this for the complete set of Fluxes (all segement boundaries)

      USE FOMS_pars_m , only : Maxsegments

      USE MSSeg_PARS_m, only : fpk, MAX_USER_LEVELS, MAX_GEOMETRIES,  &
         MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAXLAYERS,&
         MAXBEAMS, MAXSTOKES, MAXSTREAMS, MAXMOMENTS, MAX_DIRECTIONS, &
         MAX_MESSAGES

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VQS_Main_Outputs

!  Stokes vector values: Total, First-Order (FO), and Multiple Scatter (MS)

      REAL(fpk), DIMENSION ( MAXSTOKES ) :: StokesTot
      REAL(fpk), DIMENSION ( MAXSTOKES ) :: StokesFO
      REAL(fpk), DIMENSION ( MAXSTOKES ) :: StokesMS

!  Actinic and Regular Fluxes: Diffuse (Diff) & Direct (Dir)

      REAL(fpk), DIMENSION ( MAXSTOKES ) :: ActFluxDiff
      REAL(fpk), DIMENSION ( MAXSTOKES ) :: RegFluxDiff

      REAL(fpk), DIMENSION ( MAXSTOKES ) :: ActFluxDir
      REAL(fpk), DIMENSION ( MAXSTOKES ) :: RegFluxDir

!  4/24/24. Fluxes defined at segment boundaries along the light path

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS ) :: ActFluxDiff_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS ) :: RegFluxDiff_AllSB

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS ) :: ActFluxDir_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS ) :: RegFluxDir_AllSB

!  Diagnostic information

      INTEGER :: Regime_idx
      INTEGER :: NSegments

      END TYPE VQS_Main_Outputs

! #####################################################################
! #####################################################################

      TYPE VQS_Exception_Handling

!  Exception handling for input checking.
!     Message Length should be at least 120 Characters

      INTEGER      :: STATUS_INPUTCHECK
      INTEGER      :: NCHECKMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: CHECKMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: ACTIONS

!  Exception handling for model calculation.
!mick note 6/6/2023 - this arrangement is a little different than VLIDORT

      INTEGER      :: STATUS_CALCULATION
      INTEGER      :: NMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: MESSAGES

      END TYPE VQS_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VQS_Input_Exception_Handling


!  Exception handling for Input Checking settings. New code, 18 May 2010
!     Message Length should be at least 120 Characters

      INTEGER      :: STATUS_INPUTREAD
      INTEGER      :: NINPUTMESSAGES

      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: INPUTMESSAGES
      CHARACTER (Len=120), dimension(0:MAX_MESSAGES)  :: INPUTACTIONS


      END TYPE VQS_Input_Exception_Handling

! #####################################################################
! #####################################################################

      TYPE VQS_Outputs

      TYPE(VQS_Main_Outputs)       :: Main
      TYPE(VQS_Exception_Handling) :: Status

      END TYPE VQS_Outputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VQS_Main_Outputs, &
                VQS_Exception_Handling, &
                VQS_Input_Exception_Handling, &
                VQS_Outputs

      END MODULE VQS_Outputs_def_m
