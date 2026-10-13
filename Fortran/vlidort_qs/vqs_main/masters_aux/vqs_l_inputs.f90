
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
! #    These routines called by Main VQS LP module              #
! #                                                             #
! #            VQS_L_CHECK_INPUT_DIMS                           #
! #                                                             #
! ###############################################################

      MODULE vqs_l_inputs_m

      !PRIVATE
      PUBLIC

      CONTAINS

      SUBROUTINE VQS_L_CHECK_INPUT_DIMS &
        ( VQS_LinFixIn, &
          STATUS, NMESSAGES, MESSAGES, ACTIONS )

!  Check input dimensions

!  Module, dimensions and numbers

      USE MSSeg_PARS_m, only : MAX_ATMOSWFS, MAX_MESSAGES, &
                               VQS_SUCCESS, VQS_SERIOUS

      USE VQS_LinInputs_def_m

      IMPLICIT NONE

!  Subroutine inputs
!  -----------------

!  VQS input structures

      TYPE(VQS_Fixed_LinInputs)   , INTENT (IN) :: VQS_LinFixIn
      !TYPE(VQS_Modified_LinInputs), INTENT (IN) :: VQS_LinModIn

!  Exception handling.  Message Length should be at least 120 Characters

      INTEGER      , intent(out)   :: STATUS
      INTEGER      , intent(inout) :: NMESSAGES
      CHARACTER*(*), intent(inout) :: MESSAGES(0:MAX_MESSAGES)
      CHARACTER*(*), intent(inout) :: ACTIONS (0:MAX_MESSAGES)

!  Local variables

      INTEGER :: NM

!  Initialize Exception handling

      STATUS = VQS_SUCCESS
      NM = NMESSAGES

!  Check VQS linearized input dimensions against maximum dimensions
!  ================================================================

      !IF ( VQS_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION ) THEN
        IF ( VQS_LinFixIn%Cont%N_PROFILE_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Number of profile WFs N_PROFILE_WFS > maximum dimension'
         ACTIONS(NM)  = 'Re-set input value or increase MAX_ATMOSWFS in MSSeg_PARS'
         STATUS = VQS_SERIOUS
        ENDIF
      !ENDIF

!  9/21/23. Add check on N_COLUMN_WFS

       IF ( VQS_LinFixIn%Cont%N_COLUMN_WFS .GT. MAX_ATMOSWFS ) THEN
         NM = NM + 1
         MESSAGES(NM) = 'Number of column WFs N_COLUMN_WFS > maximum dimension'
         ACTIONS(NM)  = 'Re-set input value or increase MAX_ATMOSWFS in MSSeg_PARS'
         STATUS = VQS_SERIOUS
       ENDIF

!  Update NMESSAGES

      NMESSAGES = NM

!  Finish

      END SUBROUTINE VQS_L_CHECK_INPUT_DIMS

!  End Module

      END MODULE vqs_l_inputs_m

