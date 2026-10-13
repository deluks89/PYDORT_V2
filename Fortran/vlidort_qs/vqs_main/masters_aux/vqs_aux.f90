
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

! ###################################################################
! #                                                                 #
! #  This is VQS_AUX.f90. This is a collection of utility           #
! #  modules compiled from a variety of sources. This is the        #
! #  only part of the VQS package not written by the Author.        #
! #  The subroutines in VQS_AUX are listed below with their         #
! #  source of origin (order of appearance).                        #
! #                                                                 #
! #      vqs_read_error:        from M. Christi, 2017               #
! #      vqs_write_status:      from R.J.D. Spurr, RTS, 2005, 2010  #
! #                                                                 #
! ###################################################################

      MODULE vqs_aux_m

!  Parameter types

      USE MSSeg_PARS_m

      CONTAINS

      SUBROUTINE VQS_READ_ERROR ( ERRORFILE, VQS_InputStatus )

!  Module, dimensions and numbers

      USE MSSeg_PARS_m, only : VQS_ERRUNIT
      USE VQS_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file

      CHARACTER (LEN=*), intent(in) :: ERRORFILE

!  VQS file inputs status

      TYPE(VQS_Input_Exception_Handling) :: VQS_InputStatus

!  Local variables

      INTEGER :: N, W

!  Define some local variables

      W = VQS_ERRUNIT

!  Write VQS configuration file read errors to VQS error file 

      OPEN (UNIT = W, FILE = TRIM(ERRORFILE), STATUS = 'REPLACE')
      WRITE(W,*)' FATAL:   Wrong input from VQS input file-read'
      WRITE(W,*)'  ------ Here are the messages and actions '
      WRITE(W,'(A,I3)')'    ** Number of messages = ',VQS_InputStatus%NINPUTMESSAGES
      DO N = 1, VQS_InputStatus%NINPUTMESSAGES
        WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',&
          ADJUSTL(TRIM(VQS_InputStatus%INPUTMESSAGES(N)))
        WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',&
          ADJUSTL(TRIM(VQS_InputStatus%INPUTACTIONS(N)))
      ENDDO
      CLOSE(W)

      WRITE(*,'(/1X,A)') 'Read-input fail: Look at file ' // TRIM(ERRORFILE)
      STOP 'Program stopped in call to VQS_READ_ERROR'

      END SUBROUTINE VQS_READ_ERROR

!

      SUBROUTINE VQS_WRITE_STATUS &
         (ERRORUNIT, OPENFILEFLAG, VQS_Status )

      USE MSSeg_PARS_m, Only : MAX_MESSAGES, VQS_SUCCESS, &
                               VQS_SERIOUS,  VQS_WARNING
      USE VQS_Outputs_def_m

      IMPLICIT NONE

!  Subroutine arguments
!  --------------------

!  Error file and unit (will only opened if overall status not successfu

      INTEGER, INTENT(IN)           :: ERRORUNIT

      LOGICAL, INTENT(INOUT)        :: OPENFILEFLAG

!  VQS status

      TYPE(VQS_Exception_Handling), INTENT(IN) :: VQS_Status

!  Local variables

!  Exception handling for Input Checking.
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS_INPUTCHECK
      INTEGER ::             NCHECKMESSAGES
      CHARACTER (LEN=120) :: CKMESSAGES (0:MAX_MESSAGES)
      CHARACTER (LEN=120) :: ACTIONS (0:MAX_MESSAGES)

!  Exception handling for Model Calculation.

      INTEGER ::             STATUS_CALCULATION
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES(0:MAX_MESSAGES)

!  Help

      INTEGER ::           N, W
      LOGICAL ::           STATUS_OVERALL, VQS_ABORT

!  Copy inputs to local variables
      STATUS_INPUTCHECK            = VQS_Status%STATUS_INPUTCHECK
      NCHECKMESSAGES               = VQS_Status%NCHECKMESSAGES
      CKMESSAGES(0:NCHECKMESSAGES) = VQS_Status%CHECKMESSAGES(0:NCHECKMESSAGES)
      ACTIONS(0:NCHECKMESSAGES)    = VQS_Status%ACTIONS(0:NCHECKMESSAGES)

      IF (STATUS_INPUTCHECK .NE. VQS_SERIOUS) THEN
        STATUS_CALCULATION    = VQS_Status%STATUS_CALCULATION

        NMESSAGES             = VQS_Status%NMESSAGES
        MESSAGES(0:NMESSAGES) = VQS_Status%MESSAGES(0:NMESSAGES)

        STATUS_OVERALL = ( STATUS_INPUTCHECK  .NE. VQS_SUCCESS ) .or. &
                         ( STATUS_CALCULATION .NE. VQS_SUCCESS )
      ELSE
        STATUS_OVERALL = .true.
      ENDIF

!  Open if flagged
      IF ( .not.OPENFILEFLAG .and. STATUS_OVERALL ) THEN
        OPENFILEFLAG = .true.
        OPEN(ERRORUNIT, FILE = 'Debug/Execution.log', STATUS = 'unknown' )
      ENDIF

!  Define some local variables

      W = ERRORUNIT
      VQS_ABORT = .false.

!  Check status and write messages

      IF ( STATUS_INPUTCHECK .EQ. VQS_WARNING ) THEN
        WRITE(W,*)' WARNING: Suspicious input from VQS input check'
        WRITE(W,*)'  ------ VQS has executed with internal defaults'
        WRITE(W,*)'  ------ Here are the messages and actions'
        WRITE(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(Trim(CKMESSAGES(N)))
          WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',Adjustl(Trim(ACTIONS(N)))
        ENDDO
      ELSE IF ( STATUS_INPUTCHECK .EQ. VQS_SERIOUS ) THEN
        WRITE(W,*)' FATAL:   Wrong input from VQS input check'
        WRITE(W,*)'  ------ VQS will not execute'
        WRITE(W,*)'  ------ Here are the messages and actions'
        WRITE(W,'(A,I3)')'    ** Number of messages = ',NCHECKMESSAGES
        DO N = 1, NCHECKMESSAGES
          WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(Trim(CKMESSAGES(N)))
          WRITE(W,'(A,I3,A,A)')'Action  # ',N,' : ',Adjustl(Trim(ACTIONS(N)))
        ENDDO

        !Set abort flag and alert the user at the screen
        VQS_ABORT = .true.
        WRITE(*,*)
        WRITE(*,*)'VQS input abort'
      ENDIF

!  Model calculation status

      IF ( STATUS_INPUTCHECK .NE. VQS_SERIOUS ) THEN
        IF ( STATUS_CALCULATION .NE. VQS_SUCCESS ) THEN
          WRITE(W,*)' FATAL: VQS execution failure '
          WRITE(W,*)'  ------ Here are the messages'
          WRITE(W,'(A,I3)')'    ** Number of messages = ',NMESSAGES
          DO N = 1, NMESSAGES
            WRITE(W,'(A,I3,A,A)')'Message # ',N,' : ',Adjustl(Trim(MESSAGES(N)))
          ENDDO

          !Set abort flag and alert the user at the screen
          VQS_ABORT = .true.
          WRITE(*,*)
          WRITE(*,*)'VQS calculation abort'
        ENDIF
      ENDIF

!  If fatal error, abort

      IF ( VQS_ABORT ) THEN
        CLOSE(ERRORUNIT)
        WRITE(*,*)
        WRITE(*,'(/1x,a)')'Main program failed to execute, ' // &
          'error messages in Debug/Execution.log'
        STOP 'Program stopped in call to VQS_WRITE_STATUS'
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VQS_WRITE_STATUS

!  End module

      END MODULE vqs_aux_m
