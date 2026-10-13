
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

      MODULE MSSeg_Homog_Fourier_m

!  05 April 2023. Revision

      public  MSSEG_HOMOG_FOURIER

      contains

      SUBROUTINE MSSEG_HOMOG_FOURIER ( &
           DO_LAMBERTIAN_SURFACE, DO_FOURIER0_NSTOKES2, ASYMTX_TOLERANCE, &
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NMOMENTS, NSTREAMS_2,     &
           NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,    &
           MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
           QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,                      &
           LAMBERTIAN_ALBEDO, BRDF_F, DELTAU_VERT, OMEGA_GREEK,           &           
           DO_INCLUDE_SURFACE, HGSF, HSOL, HBVP,                          &
           STATUS, MESSAGE, TRACE_1, TRACE_2 )

!  Mk2, 6/21/21.
!   -- Use module names now have the "_m" extension, Only use quantities you need
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- Flux output (LOCAL_FLUXES) controlled by DO_LOCAL_FLUXES input flag
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.
!   -- Call subroutine specially written to get the discrete-ordinate field for Fluxes

!  Mk3. Experimental. 8/27/21.
!    -- Use internal Type structures. This is the homogeneous Fourier.

!  Mk3. 4/11/23. Modified for MK3 use

      USE MSSeg_PARS_m, Only : ZERO, ONE, TWO, MAXLAYERS, MAXMOMENTS, MAXSTOKES, MAXSTREAMS, &
                               MAXSTRMSTKS, MAXEVALUES, MAXSTOKES_SQ, VQS_SERIOUS, VQS_SUCCESS

      use MSSeg_HPTypes_Def_m
      use MSSeg_Homog_Solutions_m
      use MSSeg_Homog_BVProblem_m

      IMPLICIT NONE

!  INPUTS
!  ======

!  Input flags
!  -----------

      LOGICAL, INTENT (IN) ::          DO_LAMBERTIAN_SURFACE

!  Mk2, 6/21/21. Add the following three inputs
!   -- ASYMTX_TOLERANCE now introduced and passed into the QHOM call to ASYMTX
!   -- DO_FOURIER0_NSTOKES2 performance flag introduced.

      LOGICAL, INTENT (IN) ::          DO_FOURIER0_NSTOKES2
      DOUBLE PRECISION, INTENT (IN) :: ASYMTX_TOLERANCE

!  Integers
!  --------

!  Control integers (basic)

      INTEGER, INTENT (IN) ::          FOURIER
      INTEGER, INTENT (IN) ::          NSTOKES
      INTEGER, INTENT (IN) ::          NSTREAMS
      INTEGER, INTENT (IN) ::          NLAYERS
      INTEGER, INTENT (IN) ::          NMOMENTS ! derived

!  Control integers (derived). These are local valaues

      INTEGER, INTENT (IN) ::          NSTREAMS_2
      INTEGER, INTENT (IN) ::          NTOTAL
      INTEGER, INTENT (IN) ::          N_SUBDIAG
      INTEGER, INTENT (IN) ::          N_SUPDIAG
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS
      INTEGER, INTENT (IN) ::          NSTKS_NSTRMS_2

!  Derived book-keeping quantities
!  -------------------------------

      INTEGER, INTENT (IN) ::          MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION, INTENT (IN) :: DMAT  ( MAXSTOKES, MAXSTOKES )

      LOGICAL, INTENT (IN) ::          DO_REAL_EIGENSOLVER  ( 0:MAXMOMENTS, MAXLAYERS )
      LOGICAL, INTENT (IN) ::          DO_LAYER_SCATTERING  ( 0:MAXMOMENTS, MAXLAYERS )

!  Quadrature
!  ----------

      DOUBLE PRECISION, INTENT (IN) :: QUAD_STREAMS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_STRMWTS ( MAXSTREAMS )
      DOUBLE PRECISION, INTENT (IN) :: QUAD_HALFWTS ( MAXSTREAMS )

!  Optical inputs
!  --------------

!  surface reflection
!   -- Mk3. Reduce dimensioning, drop 0:MAXMOMENTS

      DOUBLE PRECISION, INTENT (IN) :: LAMBERTIAN_ALBEDO
      DOUBLE PRECISION, INTENT (IN) :: BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )

!  From the Deltam-scaling setup

      DOUBLE PRECISION, INTENT (INOUT) :: DELTAU_VERT ( MAXLAYERS )

!  From the SSALB_INIT setup

      DOUBLE PRECISION, INTENT (INOUT) :: OMEGA_GREEK  ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )

!  OUTPUTS
!  =======

!  Type structure output
!  ---------------------

      type(Homog_GSF)      , intent(inout) :: HGSF
      type(Homog_Solutions), intent(inout) :: HSOL
      type(Homog_BVPSetup) , intent(inout) :: HBVP

!  Bookkeeping and exception handling
!  ----------------------------------

      LOGICAL, INTENT (OUT) ::           DO_INCLUDE_SURFACE
      INTEGER, INTENT (OUT) ::           STATUS
      CHARACTER (LEN=*), INTENT (OUT) :: MESSAGE
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_1
      CHARACTER (LEN=*), INTENT (OUT) :: TRACE_2

!  Local variables
!  ===============

!  For the homogeneous solutions

      DOUBLE PRECISION :: REAL_KSQ    ( MAXSTRMSTKS )
      DOUBLE PRECISION :: IMAG_KSQ    ( MAXSTRMSTKS )
      DOUBLE PRECISION :: LEFT_EVEC   ( MAXSTRMSTKS, MAXSTRMSTKS )
      DOUBLE PRECISION :: RITE_EVEC   ( MAXSTRMSTKS, MAXSTRMSTKS )
      LOGICAL ::          EIGENDEGEN  ( MAXSTRMSTKS, MAXLAYERS )
      INTEGER ::          EIGENMASK_R ( MAXEVALUES )
      INTEGER ::          EIGENMASK_C ( MAXEVALUES )
      DOUBLE PRECISION :: KEIGEN_CSQ  ( MAXEVALUES )

!  others

      INTEGER          :: LAYER

!  Mk2, 6/21/21. Based on VLIDORT Version 2.8.3. Local NSTOKES and related control integers

      INTEGER :: LC_NSTK, LC_NSTKNSTRM, LC_NSTKNSTRM_2, LC_NTT, LC_NSUBD, LC_NSUPD

!  Flux multiplier and Fourier component numbers

      DOUBLE PRECISION ::  DELTA_FACTOR
      DOUBLE PRECISION ::  SURFACE_FACTOR

!  error tracing

      INTEGER ::           STATUS_SUB
      CHARACTER (LEN=2) :: CF

!  debug variables
!      integer :: i,j

!  ##############
!  initialization
!  ##############

!  module status and message initialization

      STATUS  = VQS_SUCCESS
      MESSAGE = ' '
      TRACE_1 = ' '
      TRACE_2 = ' '

!  Set local flags
!  ---------------

!  Albedo flag (for inclusion of some kind of reflecting boundary)

      DO_INCLUDE_SURFACE = .TRUE.
      IF ( DO_LAMBERTIAN_SURFACE ) THEN
        IF ( FOURIER .NE. 0 ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
        ELSE
          IF ( LAMBERTIAN_ALBEDO .EQ. ZERO ) THEN
            DO_INCLUDE_SURFACE = .FALSE.
          ENDIF
        ENDIF
      ENDIF

!  surface reflectance factors

      IF ( FOURIER .EQ. 0 ) THEN
        SURFACE_FACTOR = TWO
        DELTA_FACTOR   = ONE
      ELSE
        SURFACE_FACTOR = ONE
        DELTA_FACTOR   = TWO
      ENDIF

!  ###################
!  Spherical functions
!  ###################

!  Get Pi Matrices for this Fourier component

      CALL MSSEG_HOMOG_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NSTREAMS, NMOMENTS, QUAD_STREAMS, DMAT, HGSF )

!  ###################
!  LOCAL NSTOKES USAGE
!  ###################

!  Mk2, 6/21/21. Based on VLIDORT Version 2.8.3. 
!    -- Introduce control for Fast execution of Fourier-zero with LOCAL_NSTOKES = 2

      IF ( DO_Fourier0_NSTOKES2 .and. NSTOKES.eq.3 .and. FOURIER.eq.0 .and. DO_LAMBERTIAN_SURFACE ) THEN
        LC_NSTK = 2 ; LC_NSTKNSTRM = NSTREAMS * LC_NSTK ; LC_NSTKNSTRM_2 = LC_NSTKNSTRM * 2
        LC_NTT  = LC_NSTKNSTRM_2 * NLAYERS
        IF ( NLAYERS .EQ. 1 ) THEN
          LC_NSUBD = 2*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ELSE
          LC_NSUBD = 3*LC_NSTKNSTRM - 1 ; LC_NSUPD = LC_NSUBD
        ENDIF
      ELSE
        LC_NSTK = NSTOKES ; LC_NSTKNSTRM   = NSTKS_NSTRMS ; LC_NSTKNSTRM_2 = NSTKS_NSTRMS_2
        LC_NTT  = NTOTAL  ; LC_NSUBD       = N_SUBDIAG    ; LC_NSUPD       = N_SUPDIAG
      ENDIF

!  ########################################
!  RT differential equation Eigensolutions
!  ########################################

!  Start layer loop

      DO LAYER = 1, NLAYERS

!  Get Discrete ordinate solutions for this layer
!    MK2, 6/21/21. Add ASYMTX_TOLERANCE, Use local variables LC_NSTK, LC_NSTKNSTRM

         CALL MSSEG_QHOM_SOLUTION ( &
           LAYER, FOURIER, LC_NSTK, NSTREAMS, NMOMENTS, LC_NSTKNSTRM,            &
           QUAD_STREAMS, QUAD_HALFWTS, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
           ASYMTX_TOLERANCE, DELTAU_VERT, OMEGA_GREEK, HGSF, &
           HSOL, REAL_KSQ, IMAG_KSQ, LEFT_EVEC, RITE_EVEC,   &
           EIGENDEGEN, EIGENMASK_R, EIGENMASK_C, KEIGEN_CSQ, &
           STATUS_SUB, MESSAGE, TRACE_1 )

!if ( layer.gt.70)write(*,*)FOURIER,LAYER,HSOL%SOLA_XPOS(1,1,1:4,LAYER)
!if ( layer.eq.40)write(*,*)FOURIER,LAYER,HSOL%SAB(1:4,1,1,1,LAYER)
!  .. error tracing

         IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
           write(CF,'(I2)')FOURIER
           TRACE_2 =  'Called in MSSEG_HOMOG_FOURIER, Fourier component '//CF
           STATUS = VQS_SERIOUS ; RETURN
         ENDIF

!  end layer loop

      ENDDO

!  ############################################
!   boundary value problem - MATRIX PREPARATION
!  ############################################

!  standard case using compression of band matrices, etc..
!    MK2, 6/21/21.  Use local bookkeeping variables LC_NSTK etc.....

      CALL MSSEG_MATRIXSETUP_MASTER ( &
        DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, SURFACE_FACTOR,    &
        FOURIER, QUAD_STRMWTS, LAMBERTIAN_ALBEDO, BRDF_F,             &
        LC_NSTK, NSTREAMS, NLAYERS, LC_NSTKNSTRM, MUELLER_INDEX,      &
        NSTREAMS_2, LC_NTT, LC_NSUBD, LC_NSUPD, LC_NSTKNSTRM_2,       &
        HSOL, HBVP, STATUS_SUB, MESSAGE, TRACE_1 )

      IF ( STATUS_SUB .NE. VQS_SUCCESS ) THEN
        write(CF,'(I2)')FOURIER
        TRACE_2 = 'Error from MSSEG_MATRIXSETUP_MASTER, '// &
          'Called in MSSEG_HOMOG_FOURIER, Fourier # '//CF
        STATUS = VQS_SERIOUS ; RETURN
      ENDIF

!  finish

RETURN
END SUBROUTINE MSSEG_HOMOG_FOURIER

!  End module

END MODULE msseg_homog_fourier_m

