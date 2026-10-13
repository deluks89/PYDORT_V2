
! ###############################################################
! #                                                             #
! #                     VQS_QS_1p4                          #
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
! #      MS code upgraded to VQS 2.7 standards         (1.2)    #
! #      FO/MS codes upgraded to VQS 2.8 standards     (1.3)    #
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
! #    These routines are called by the Main VQS module         #
! #                                                             #
! #            VQS_CHECK_INPUT_DIMS                             #
! #            VQS_CHECK_INPUT                                  #
! #                                                             #
! ###############################################################

      MODULE vqs_inputs_m

      !PRIVATE
      PUBLIC

      CONTAINS

      SUBROUTINE VQS_CHECK_INPUT_DIMS &
        ( VQS_FixIn, VQS_ModIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  Module, dimensions and numbers

      USE MSSeg_PARS_m, only : MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT, MAXSTOKES, &
                               MAXFINELAYERS, MAX_MESSAGES, &
                               VQS_SUCCESS, VQS_SERIOUS

      USE VQS_Inputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VQS input structures

      TYPE(VQS_Fixed_Inputs)   , INTENT (IN) :: VQS_FixIn
      TYPE(VQS_Modified_Inputs), INTENT (IN) :: VQS_ModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , INTENT (OUT)   :: STATUS
      INTEGER      , INTENT (INOUT) :: NMESSAGES
      CHARACTER*(*), INTENT (INOUT) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), INTENT (INOUT) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VQS_SUCCESS
      NM = NMESSAGES

!  Check VQS input dimensions against maximum dimensions
!  =====================================================

!  1a. Basic dimensions - always checked

      IF ( VQS_FixIn%Cont%NSTOKES .GT. MAXSTOKES ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number Stokes parameters NSTOKES > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value to 4 or less'
        STATUS       = VQS_SERIOUS
      ENDIF

      IF ( VQS_FixIn%Cont%NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams NSTREAMS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension in MSSeg_PARS'
        STATUS       = VQS_SERIOUS
      ENDIF

      IF ( VQS_FixIn%Cont%NLAYERS .GT. MAXLAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of layers NLAYERS > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXLAYERS dimension in MSSeg_PARS'
        STATUS       = VQS_SERIOUS
      ENDIF

      IF ( VQS_FixIn%Cont%NFINE_INPUT .GT. MAXFINELAYERS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of fine layers NFINE_INPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXFINELAYERS dimension in MSSeg_PARS'
        STATUS       = VQS_SERIOUS
      ENDIF

      IF ( VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT .GT. MAXMOMENTS_INPUT ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of Expansion coefficients NGREEK_MOMENTS_INPUT > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXMOMENTS_INPUT dimension in MSSeg_PARS'
        STATUS       = VQS_SERIOUS
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      RETURN
      END SUBROUTINE VQS_CHECK_INPUT_DIMS

!

      SUBROUTINE VQS_CHECK_INPUT &
        ( VQS_FixIn, VQS_ModIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

      USE MSSeg_pars_m, Only : MAX_MESSAGES, &
                               VQS_WARNING, VQS_SUCCESS, VQS_SERIOUS, &
                               ZERO, ONE, OMEGA_SMALLNUM

      USE VQS_IO_DEFS_m

      IMPLICIT NONE

!  VQS Inputs
!  ==========

!  VQS input structures

      TYPE(VQS_Fixed_Inputs), INTENT (IN)       :: VQS_FixIn
      TYPE(VQS_Modified_Inputs), INTENT (INOUT) :: VQS_ModIn

!  Exception handling

      INTEGER, INTENT (OUT) ::              STATUS
      INTEGER, INTENT (INOUT) ::            NMESSAGES
      CHARACTER (LEN=*), INTENT (INOUT) ::  MESSAGES( 0:MAX_MESSAGES )
      CHARACTER (LEN=*), INTENT (INOUT) ::  ACTIONS ( 0:MAX_MESSAGES )

!  Local variables

      INTEGER :: I, L, NM
      LOGICAL :: LOOP
      CHARACTER (LEN=3) :: C3

!  Initialize output status and number of messages

      STATUS = VQS_SUCCESS
      NM = NMESSAGES

!  Check Delta-M scaling
!  ---------------------

!  Just FO alone - diffuse-field Delta-M scaling should be turned off

      IF ( .NOT.VQS_FixIn%Bool%DO_MSRT ) THEN
        IF ( VQS_ModIn%MBool%DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: First-Order ALONE --> diffuse-field delta-M on'
          ACTIONS(NM)  = 'First-Order ALONE: internal default to deltam_scaling = .false.'
          STATUS = VQS_WARNING
          VQS_ModIn%MBool%DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Just Rayleigh scattering - diffuse-field Delta-M scaling should be turned off

      IF ( VQS_ModIn%MBool%DO_RAYLEIGH_ONLY ) THEN
        IF ( VQS_ModIn%MBool%DO_DELTAM_SCALING ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: No delta-M scaling with Rayleigh-only'
          ACTIONS(NM)  = 'Warning: DO_DELTAM_SCALING turned off internally'
          STATUS = VQS_WARNING
          VQS_ModIn%MBool%DO_DELTAM_SCALING = .FALSE.
        ENDIF
      ENDIF

!  Check viewing geometry input
!  ----------------------------

!  Check Camera solar zenith angle

      IF ( VQS_FixIn%UserVal%theta_Camera.LT.ZERO .OR. &
           VQS_FixIn%UserVal%theta_Camera.GE.90.0D0 ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: out-of-range Camera solar zenith angle.'
        ACTIONS(NM)  = 'Look at theta_Camera input, should be [0,90)'
        STATUS = VQS_SERIOUS
      ENDIF

!  Check Camera relative azimuth angle

      IF ( VQS_FixIn%UserVal%relaz_Camera.LT.ZERO .OR. &
           VQS_FixIn%UserVal%relaz_Camera.GT.360.0D0 ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: out-of-range Camera azimuth angle'
        ACTIONS(NM)  = 'Look at relaz_Camera input, should be in [0,360]'
        STATUS = VQS_SERIOUS
      ENDIF

!  Check Camera view zenith angle

      IF ( VQS_FixIn%UserVal%alpha_Camera.LT.ZERO .OR. &
           VQS_FixIn%UserVal%alpha_Camera.GT.180.0D0 ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: out-of-range Camera view zenith angle'
        ACTIONS(NM)  = 'Look at alpha_Camera input, should be in [0,180]'
        STATUS = VQS_SERIOUS
      ENDIF

!  Check geophysical input
!  -----------------------

!  Check planet radius

      IF ( .NOT. VQS_FixIn%Bool%DO_PLANE_PARALLEL ) THEN
        IF ( VQS_FixIn%Geophys%PLANET_RADIUS .EQ. 0.0D0 ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: planet radius = zero'
          ACTIONS(NM)  = 'Look at planet radius'
          STATUS = VQS_SERIOUS
        ENDIF
      ENDIF

!  Check atmospheric height grid

      LOOP = .TRUE.
      I = 0
      DO WHILE (LOOP .AND. I.LT.VQS_FixIn%Cont%NLAYERS)
        I = I + 1
        IF ( VQS_FixIn%Geophys%HEIGHT_GRID(I-1).LE.VQS_FixIn%Geophys%HEIGHT_GRID(I) ) THEN
          WRITE(C3,'(I2)')I
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Height-grid not monotonically decreasing; Layer '//C3
          ACTIONS(NM)  = 'Look at Height-grid input'
          LOOP = .FALSE.
          STATUS = VQS_SERIOUS
        ENDIF
      ENDDO

!  Check Camera height (should always be within the atmosphere!)
!    --10/10/23. This is a check on CAMERA LEVEL. COde revised.
!      IF ( VQS_FixIn%UserVal%height_Camera.GT.DBLE(VQS_FixIn%Cont%NLAYERS) .OR. &
!           VQS_FixIn%UserVal%height_Camera.LT.ZERO )  THEN

      IF ( VQS_FixIn%UserVal%height_Camera.LT.ZERO )  THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Out of range for Camera height'
        ACTIONS(NM)  = 'Re-set Camera height'
        STATUS = VQS_SERIOUS
      ENDIF

      IF ( VQS_FixIn%UserVal%height_Camera.GT.VQS_FixIn%Geophys%HEIGHT_GRID(0) ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Out of range for Camera height'
        ACTIONS(NM)  = 'Re-set Camera height'
        STATUS = VQS_SERIOUS
      ENDIF

!  Check atmospheric optical properties
!  ------------------------------------

!  Check layer optical thicknesses

      DO L = 1, VQS_FixIn%Cont%NLAYERS
        IF ( VQS_FixIn%Optical%DELTAU_VERT_INPUT(L).LE.ZERO ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: optical thickness <= 0, layer '//C3
          ACTIONS(NM)  = 'Check optical thickness input'
          STATUS       = VQS_SERIOUS
        ENDIF
      ENDDO

!  Check layer single scatter albedos

      DO L = 1, VQS_FixIn%Cont%NLAYERS
        IF ( VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(L).GT.ONE-OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 1, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VQS_SERIOUS
        ELSE IF ( VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(L).LT.OMEGA_SMALLNUM ) THEN
          NM = NM + 1
          WRITE(C3,'(I3)')L
          MESSAGES(NM) = 'Bad input: SS-albedo too close to 0, layer '//C3
          ACTIONS(NM)  = 'Check SS-albedo input'
          STATUS = VQS_SERIOUS
        ENDIF
      ENDDO

!  Check number of input expansion coefficient moments

      IF ( .NOT.VQS_ModIn%MBool%DO_RAYLEIGH_ONLY ) THEN

        !Checks for general scattering case

        IF ( VQS_ModIn%MBool%DO_DELTAM_SCALING ) THEN
          IF ( VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT .LT. 2*VQS_FixIn%Cont%NSTREAMS ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Fewer than 2N expansion moments with delta-M'
            ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N internally'
            STATUS = VQS_WARNING
            VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT = 2*VQS_FixIn%Cont%NSTREAMS
          ENDIF
        ELSE
          IF ( VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT .LT. 2*VQS_FixIn%Cont%NSTREAMS-1 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Fewer than 2N-1 expansion moments without delta-M'
            ACTIONS(NM)  = 'Warning: Re-set NGREEK_MOMENTS_INPUT to 2N-1 internally'
            STATUS = VQS_WARNING
            VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT = 2*VQS_FixIn%Cont%NSTREAMS - 1
          ENDIF
        ENDIF

      ELSE

        !Checks for Rayleigh-only option

        IF ( VQS_ModIn%MBool%DO_RAYLEIGH_ONLY ) THEN
          IF ( VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT.NE.2 ) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: Rayleigh-only, expansion momemts NOT = 2'
            ACTIONS(NM)  = 'Warning: Set NGREEK_MOMENTS_INPUT = 2 internally'
            STATUS       = VQS_WARNING
            VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT = 2
          ENDIF
        ENDIF

      ENDIF

!  Check first moment of layer phase functions

      DO L = 1, VQS_FixIn%Cont%NLAYERS
        IF ( VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT(0,L,1).NE.ONE ) THEN
          WRITE(C3,'(I3)')L
          NM = NM + 1
          MESSAGES(NM) = 'First phase moment (GREEK_11) not 1 for layer '//C3
          ACTIONS(NM)  = 'Check First phase function moment'
          STATUS = VQS_SERIOUS
        ENDIF
      ENDDO

!  Check surface properties
!  ------------------------

!  Make sure the Lambertian surface is in range

      IF ( VQS_FixIn%Optical%LAMBERTIAN_ALBEDO.LT.ZERO .OR. VQS_FixIn%Optical%LAMBERTIAN_ALBEDO.GT.ONE ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Bad input: Lambertian albedo not in range [0,1]'
        ACTIONS(NM)  = 'Check albedo input'
        STATUS       = VQS_SERIOUS
      ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VQS_CHECK_INPUT

      END MODULE vqs_inputs_m

