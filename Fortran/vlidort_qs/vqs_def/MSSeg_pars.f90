
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

!  Mk2, 6/21/21. .
!   -- module name now has the "_m" extension, Only use quantities you need
!   -- Here, Upgrade the Surface indices to VLIDORT 2.8.3.

!  File of constants for MSSEG model.

      MODULE MSSeg_pars_m
      IMPLICIT NONE

!  Real number type definitions

      INTEGER, PARAMETER :: VQS_SPKIND = SELECTED_REAL_KIND(6)
      INTEGER, PARAMETER :: VQS_DPKIND = SELECTED_REAL_KIND(15)
      INTEGER, PARAMETER :: FPK = VQS_DPKIND

!  RTS use

      INTEGER :: SDU,LDU

!  Version number
!  ==============

      CHARACTER (LEN=3), PARAMETER :: VQS_VERSION_NUMBER = '1.3'

!  File i/o unit numbers
!  ======================

      INTEGER, PARAMETER ::       VQS_INUNIT   = 21
      INTEGER, PARAMETER ::       VQS_SCENUNIT = 22
      INTEGER, PARAMETER ::       VQS_FUNIT    = 23
      INTEGER, PARAMETER ::       VQS_RESUNIT  = 24
      INTEGER, PARAMETER ::       VQS_ERRUNIT  = 25

!  Basic dimensions
!  ================

!  Computational dimensioning
!  --------------------------

!  Number of computational streams in the half-space

      INTEGER, PARAMETER :: MAXSTREAMS = 8

!  Maximum number of computational layers

      INTEGER, PARAMETER :: MAXLAYERS = 112 !51 !49

!  Maximum number of fine layers used in single scattering corrections

      INTEGER, PARAMETER :: MAXFINELAYERS = 20

!  Maximum number of input moments.
!    (Use full range for exact single scatter calculations)

      INTEGER, PARAMETER :: MAXMOMENTS_INPUT = 500

!  Geometrical and output parameters
!  ---------------------------------

!  Maximum number of solar zenith angles
!    -- Should equal MAXFINELAYERS

      INTEGER, PARAMETER :: MAX_SZANGLES = MAXFINELAYERS
			    
!  Maximum number of user-defined viewing zenith angles
!    -- Should equal MAXFINELAYERS

      INTEGER, PARAMETER :: MAX_USER_VZANGLES = MAXFINELAYERS

!  Maximum number of user-defined output relative azimuth angles
!    -- Should equal MAXFINELAYERS

      INTEGER, PARAMETER :: MAX_USER_RELAZMS = MAXFINELAYERS

!  Maximum number of Observational Geometries
!    -- Should equal MAXFINELAYERS

      INTEGER, PARAMETER :: MAX_USER_OBSGEOMS = MAXFINELAYERS

!  Maximum number of output optical depths

!      INTEGER, PARAMETER :: MAX_USER_LEVELS = 55  ! For multi-level Actinic Flux output
      INTEGER, PARAMETER :: MAX_USER_LEVELS = 1    ! For single-level A/C output

!  Version 2p7. Maximum number of Terms for Taylor series expansions
!    If you are retaining contributions of order EPS^n, 
!    then you need at least n+2 Taylor terms

   INTEGER, PARAMETER :: MAX_TAYLOR_TERMS = 7

!  Fixed parameters
!  ----------------

!  Number of Stokes parameters

      INTEGER, PARAMETER :: MAXSTOKES = 4

!  Two directions (Up and Down)

      INTEGER, PARAMETER :: MAX_DIRECTIONS = 2

!  Weighting functions
!  -------------------

!  Maximum number of profile/column weighting functions

      INTEGER, PARAMETER :: MAX_ATMOSWFS = 6 !5

!  Maximum number of error messages

      INTEGER, PARAMETER :: MAX_MESSAGES = 25

!  Derived dimensions
!  ==================

!  Copy Beam dimensioning

      INTEGER, PARAMETER :: MAXBEAMS = MAX_SZANGLES

!  Copy viewing zenith angle dimension

      INTEGER, PARAMETER :: MAX_USER_STREAMS = MAX_USER_VZANGLES

!  Maximum possible geometries. Observational Geometries only.

      INTEGER, PARAMETER :: MAX_GEOMETRIES = MAX_USER_OBSGEOMS

!  All streams

      INTEGER, PARAMETER :: MAX_ALLSTRMS = MAX_USER_STREAMS + MAXSTREAMS

!  All streams for the Legendre PI-matrix setup.
!   Refractive Geometry setting: Watch out for Kill. Memory Hog
!      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = &
!                            MAX_ALLSTRMS + MAXBEAMS*MAXLAYERS

!  All streams for the Legendre PI-matrix setup.
!   Straightline setting: This setting should avoid dimensioning error

      INTEGER, PARAMETER :: MAX_ALLSTRMS_P1 = MAX_ALLSTRMS + MAXBEAMS

!  Maximum number of moments in the diffuse field calculation
!   This is always 2*MAXSTREAMS, in case we need DELTA-M

      INTEGER, PARAMETER :: MAXMOMENTS = 2*MAXSTREAMS

!  Maximum number of Fourier components = 2*MAXSTREAMS - 1

      INTEGER, PARAMETER :: MAXFOURIER = 2*MAXSTREAMS - 1

!  Number of Stokes streams squared

      INTEGER, PARAMETER :: MAXSTOKES_SQ = 16

!  Other derived dimensions

      INTEGER, PARAMETER :: MAXSTREAMS_2   = 2*MAXSTREAMS
      INTEGER, PARAMETER :: MAXSTREAMS_21  = 2*MAXSTREAMS - 1
      INTEGER, PARAMETER :: MAXSTREAMS_P1  = MAXSTREAMS + 1
      INTEGER, PARAMETER :: MAXSTREAMS_P2  = MAXSTREAMS + 2

      INTEGER, PARAMETER :: MAXSTRMSTKS    = MAXSTREAMS * MAXSTOKES
      INTEGER, PARAMETER :: MAXSTRMSTKS_2  = 2*MAXSTRMSTKS
      INTEGER, PARAMETER :: MAXSTRMSTKS_21 = 2*MAXSTRMSTKS - 1

      INTEGER, PARAMETER :: MAXSTRMSTKS_P1 = MAXSTRMSTKS + 1
      INTEGER, PARAMETER :: MAXSTRMSTKS_P2 = MAXSTRMSTKS + 2
      INTEGER, PARAMETER :: MAXSTRMSTKS_P4 = MAXSTRMSTKS + 4

      INTEGER, PARAMETER :: MAX_USTRMSTKS = MAX_USER_STREAMS * MAXSTOKES

!  Maximum number of eigenvalues

      INTEGER, PARAMETER :: MAXEVALUES = MAXSTRMSTKS

!  For the BVP problem

      INTEGER, PARAMETER :: MAXTOTAL = MAXLAYERS*MAXSTRMSTKS_2
      INTEGER, PARAMETER :: MAXBANDTOTAL = 9*MAXSTRMSTKS - 2

!  Not so far used

      INTEGER, PARAMETER :: MAX_PSOLS = 2
      INTEGER, PARAMETER :: MAX_SCATPSOLS = MAX_PSOLS

!  Format constants
!  ================

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_HEADING = '( / T6, ''-----> '', A, /)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_INTEGER = '(T6, A, T58, I10)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_REAL    = '(T6, A, T58, 1PG14.6)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_CHAR    = '(T6, A, T48, A20)'

      CHARACTER (LEN=*), PARAMETER :: &
        FMT_SECTION = '( / T6, ''-----> '', A, /)'
!     $                  '( // T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                      / T6 A,
!     $                      / T6
!     $                  ''----------------------------------------'',
!     $                  ''-----------------------------------'',
!     $                          / )' )

!  Debug write format (DWF) constants
!  ==================================

      CHARACTER (LEN=*), PARAMETER :: DWFL  = '(A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1 = '(A,I3,A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL2 = '(2(A,I3),A,L1)'
      CHARACTER (LEN=*), PARAMETER :: DWFL1_2 = '(A,I3,2(A,L1))'

      CHARACTER (LEN=*), PARAMETER :: DWFI  = '(A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI1 = '(A,I3,A,I5)'
      CHARACTER (LEN=*), PARAMETER :: DWFI2 = '(2(A,I3),A,I5)'

      CHARACTER (LEN=*), PARAMETER :: DWFR  = '(A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR1 = '(A,I3,A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR2 = '(2(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR3 = '(3(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR4 = '(4(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR5 = '(5(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR6 = '(6(A,I3),A,ES13.6E2)'
      CHARACTER (LEN=*), PARAMETER :: DWFR7 = '(7(A,I3),A,ES13.6E2)'

      CHARACTER (LEN=*), PARAMETER :: DWFR1_2 = '(A,I3,2(A,ES13.6E2))'
      CHARACTER (LEN=*), PARAMETER :: DWFR1_3 = '(A,I3,3(A,ES13.6E2))'
      CHARACTER (LEN=*), PARAMETER :: DWFR1_4 = '(A,I3,4(A,ES13.6E2))'

      CHARACTER (LEN=*), PARAMETER :: DWFC  = '(2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC1 = '(A,I3,2A)'
      CHARACTER (LEN=*), PARAMETER :: DWFC2 = '(2(A,I3),2A)'

!  Numbers
!  =======

      DOUBLE PRECISION, PARAMETER :: &
        ONE = 1.0D0, ZERO = 0.0D0,  ONEP5 = 1.5D0
      DOUBLE PRECISION, PARAMETER :: &
        TWO = 2.0D0, THREE = 3.0D0, FOUR = 4.0D0
      DOUBLE PRECISION, PARAMETER :: &
        QUARTER = 0.25D0, HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_ONE = -ONE
      DOUBLE PRECISION, PARAMETER :: &
        MINUS_TWO = -TWO
      DOUBLE PRECISION, PARAMETER :: &
        DEG_TO_RAD = 1.7453292519943D-02
      DOUBLE PRECISION, PARAMETER :: &
        PIE = 180.0D0*DEG_TO_RAD
      DOUBLE PRECISION, PARAMETER :: &
        PI2 = 2.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PI4 = 4.0D0 * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO2 = HALF * PIE
      DOUBLE PRECISION, PARAMETER :: &
        PIO4 = QUARTER * PIE
      DOUBLE PRECISION, PARAMETER :: &
        EPS3 = 0.001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS4 = 0.0001D0
      DOUBLE PRECISION, PARAMETER :: &
        EPS5 = 0.00001D0
      DOUBLE PRECISION, PARAMETER :: &
        SMALLNUM = 1.0D-9 !1.0D-15
      DOUBLE PRECISION, PARAMETER :: &
        BIGEXP = 32.0D0

!  Rob fix 5/6/13 - Taylor series limiting values

      DOUBLE PRECISION, PARAMETER :: TAYLOR_SMALL = 1.0D-04
      DOUBLE PRECISION, PARAMETER :: TAYLOR_LARGE = 1.0d+04

!  Control for Using L'Hopital's Rule
!   Changed, January 2009 for Version 2.4..........

      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS3
!      DOUBLE PRECISION, PARAMETER :: HOPITAL_TOLERANCE = EPS5

!  Control for limits of single scatter albedo

      DOUBLE PRECISION, PARAMETER :: OMEGA_SMALLNUM = 1.0D-15

!  Control for limits of extinction optical depth along solar path

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_SPATH = 32.0D0

!  Control for limits of extinction optical depth along USER paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_UPATH = 32.0D0

!  Control for limits of extinction optical depth along QUADRATURE paths

      DOUBLE PRECISION, PARAMETER :: MAX_TAU_QPATH = 32.0D0

!  Error indices
!  =============

      INTEGER, PARAMETER :: VQS_SERIOUS  = 4
      INTEGER, PARAMETER :: VQS_WARNING  = 3
      INTEGER, PARAMETER :: VQS_INFO     = 2
      INTEGER, PARAMETER :: VQS_DEBUG    = 1
      INTEGER, PARAMETER :: VQS_SUCCESS  = 0

!  Directional indices

      INTEGER, PARAMETER :: UPIDX  = 1
      INTEGER, PARAMETER :: DNIDX  = 2

!  End of file.

      END MODULE MSSeg_pars_m
