
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

      MODULE VQS_Inputs_def_m

!  This module contains the following VQS input structures,
!  with intents :

!        VQS_Fixed_Boolean    nested in VQS_Fixed_Inputs
!        VQS_Fixed_Control    nested in VQS_Fixed_Inputs
!        VQS_Fixed_Sunrays    nested in VQS_Fixed_Inputs
!     VQS_Fixed_UserValues    nested in VQS_Fixed_Inputs
!        VQS_Fixed_Geophys    nested in VQS_Fixed_Inputs
!        VQS_Fixed_Optical    nested in VQS_Fixed_Inputs
!          VQS_Fixed_Write    nested in VQS_Fixed_Inputs
!         VQS_Fixed_Inputs    Intent(In)

!     VQS_Modified_Boolean    nested in VQS_Modified_Inputs
!     VQS_Modified_Control    nested in VQS_Modified_Inputs
!     VQS_Modified_Optical    nested in VQS_Modified_Inputs
!      VQS_Modified_Inputs    Intent(InOut)

      USE MSSeg_PARS_m, Only : fpk, &
         MAXMOMENTS_INPUT, MAXSTOKES, MAXLAYERS, MAXSTREAMS, MAX_GEOMETRIES, &
         MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_USER_OBSGEOMS, &
         MAX_USER_LEVELS, MAXSTOKES_SQ

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Boolean

!  Compute multiple-scatter RT solution

      LOGICAL     :: DO_MSRT

!  Report progress of RT Calculation. Introduced 10/15/23.
!    Mk2: FO finished, MS Segments finished
!    Mk3: FO finished, MS Fouriers finished

      LOGICAL     :: DO_RT_PROGRESS

!  Compute regular & actinic Ffluxes

      LOGICAL     :: DO_CALCULATE_FLUXES

!  Flux input is non-polarized

      LOGICAL     :: DO_SUNLIGHT

!  Treat atmosphere layers as plane parallel (MS solution only)

      LOGICAL     :: DO_PLANE_PARALLEL

!  Surface reflectance

      LOGICAL     :: DO_LAMBERTIAN_SURFACE

!  Debugging (RT Solutions Use Only)

      LOGICAL     :: DO_DEBUG_FO_GEOMETRY
      LOGICAL     :: DO_DEBUG_FO_INPUT
      LOGICAL     :: DO_DEBUG_MS_INPUT
      LOGICAL     :: DO_DEBUG_FLUXES

      END TYPE VQS_Fixed_Boolean

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Control

!  Computational method for MS scattering

      INTEGER     :: METHOD_INDEX

!  Number of computational layers

      INTEGER     :: NLAYERS

!  Number of Stokes parameters

      INTEGER     :: NSTOKES

!  Number of discrete ordinate streams (in computational stream half-space)

      INTEGER     :: NSTREAMS

!  User-defined Number of Computational fine layers per atmospheric layer
!  (used if criticality below off)

      INTEGER     :: NFINE_INPUT

!  Criticality flag and threshold for automatic defining computational
!  fine layers per atmospheric layer

      LOGICAL     :: DOCRIT
      REAL(FPK)   :: ACRIT

      END TYPE VQS_Fixed_Control

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Sunrays

!  Stokes vector elements of input flux

      REAL(FPK) :: FLUX_FACTOR

      END TYPE VQS_Fixed_Sunrays

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_UserValues

!  Height of viewing camera

      REAL(FPK)   :: HEIGHT_CAMERA

!  Solar angle at viewing camera

      REAL(FPK)   :: THETA_CAMERA

!  Viewing angle of viewing camera

      REAL(FPK)   :: ALPHA_CAMERA

!  Azimuth angle of viewing camera

      REAL(FPK)   :: RELAZ_CAMERA

!  File name of file to store scene geometries

      CHARACTER(LEN=100):: GEO_FILENAME

      END TYPE VQS_Fixed_UserValues

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Geophys

!  Earth radius (in km)

      REAL(FPK)   :: PLANET_RADIUS

!  Height of atmospheric layer boundaries (in km)

      REAL(FPK), dimension (0:MAXLAYERS) :: HEIGHT_GRID

      END TYPE VQS_Fixed_Geophys

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Optical

!  Multilayer optical property (bulk)

      REAL(FPK), dimension ( MAXLAYERS ) :: DELTAU_VERT_INPUT

!  Phase function Legendre-polynomial expansion coefficients
!    Include all that you require for exact single scatter calculations

      REAL(FPK), dimension ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ ) :: GREEKMAT_TOTAL_INPUT

!  F-matrix input as an alternative to the use of expansion coefficients in
!   the single-scatter (FO) code.

      !REAL(FPK), dimension ( MAXLAYERS, 6 ) :: FMATRIX_UP
      !REAL(FPK), dimension ( MAXLAYERS, 6 ) :: FMATRIX_DN

!  Lambertian Surface control

      REAL(FPK) :: LAMBERTIAN_ALBEDO

!  Wavelength Index (RT Solutions Use Only)

      INTEGER   :: ATMOS_INDEX

      END TYPE VQS_Fixed_Optical

! #####################################################################
! #####################################################################

      TYPE VQS_Fixed_Inputs

      TYPE(VQS_Fixed_Boolean)    :: Bool
      TYPE(VQS_Fixed_Control)    :: Cont
      TYPE(VQS_Fixed_Sunrays)    :: Sunrays
      TYPE(VQS_Fixed_UserValues) :: UserVal
      TYPE(VQS_Fixed_Geophys)    :: Geophys
      TYPE(VQS_Fixed_Optical)    :: Optical

      END TYPE VQS_Fixed_Inputs

! #####################################################################
! #####################################################################

      TYPE VQS_Modified_Boolean

!  Rayleigh scattering flag

      LOGICAL     :: DO_RAYLEIGH_ONLY

!  Delta-M scaling flag

      LOGICAL     :: DO_DELTAM_SCALING

      END TYPE VQS_Modified_Boolean

! #####################################################################
! #####################################################################

      TYPE VQS_Modified_Control

!  Number of Legendre phase function expansion moments
!       May be re-set after Checking

      INTEGER   :: NGREEK_MOMENTS_INPUT

      END TYPE VQS_Modified_Control

! #####################################################################
! #####################################################################

      TYPE VQS_Modified_Optical

!  Multilayer optical property (bulk)

      REAL(FPK), dimension ( MAXLAYERS ) :: OMEGA_TOTAL_INPUT

      END TYPE VQS_Modified_Optical

! #####################################################################
! #####################################################################

      TYPE VQS_Modified_Timing

!  Timing for various VQS RT solution computations

      REAL :: GEO_TIME
      REAL :: SS_TIME
      REAL :: MS_TIME

      END TYPE VQS_Modified_Timing

! #####################################################################
! #####################################################################

      TYPE VQS_Modified_Inputs

      TYPE(VQS_Modified_Boolean) :: MBool
      TYPE(VQS_Modified_Control) :: MCont
      TYPE(VQS_Modified_Optical) :: MOptical
      TYPE(VQS_Modified_Timing)  :: MTiming

      END TYPE VQS_Modified_Inputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VQS_Fixed_Boolean, &
                VQS_Fixed_Control, &
                VQS_Fixed_Sunrays, &
                VQS_Fixed_UserValues, &
                VQS_Fixed_Geophys, &
                VQS_Fixed_Optical, &
                VQS_Fixed_Inputs,  &
                VQS_Modified_Boolean, &
                VQS_Modified_Control, &
                VQS_Modified_Optical, &
                VQS_Modified_Inputs

      END MODULE VQS_Inputs_def_m
