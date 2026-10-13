
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

!  March 2020 update

module FOMS_pars_m

!mick mod 12/1/2021 - added debug write format (DWF) constants
!mick mod 6/6/2023  - removed (DWF) constants here
!                   - made all MSSEG parameters availabe to FOMS
!                   - trimmed real number constants defined here to those unique to FOMS

   use MSSeg_PARS_m

   implicit none

!  Precision
!  =========

   integer, parameter :: i4b = selected_int_kind(9) !4-byte integer
   integer, parameter :: ffp = fpk                  !8-byte real

!  Dimensioning for the DRIVERS
!  ============================
!@@! RTS 6/16/15. Add MaxViewAngles
!@@! RTS 6/26/15. Add MaxSunAngles, MaxAzimuths. MaxAllGeoms

   integer, parameter :: MaxViewAngles      = 1
   integer, parameter :: MaxSunAngles       = 1
   integer, parameter :: MaxAzimuths        = 1
   integer, parameter :: MaxAllGeoms        = MaxSunAngles * MaxViewAngles * MaxAzimuths
   integer, parameter :: maxmessages        = 10

!  External Dimensions with 'EXT_' prefix, distinguish from VLIDORT variables
!mick note 6/6/2023 - EXT_maxlayers, EXT_maxmoments, EXT_maxatmoswfs should equal
!                     "MSSeg_pars.f90" values MAXLAYERS, MAXMOMENTS_INPUT, MAX_ATMOSWFS.
!                     Also, affiliated FOMS parameters now dependent on MSSEG parameters
!                     in "MSSeg_pars.f90"

   integer, parameter :: EXT_maxlayers      = 100 !51 !53
   integer, parameter :: EXT_maxmoments     = 200
   integer, parameter :: EXT_maxatmoswfs    = 3  !5
   integer, parameter :: EXT_maxgeometries  = MaxAllGeoms

!  Dimensioning for the data
!  =========================
!  Removed for the general Version (March 2020)
!    4/5/16. Number of moments must be maxval(ic_maxmoms,aer_maxmoms,wc_maxmoms)

!  Gases, wavelengths

   integer, parameter :: maxgases       = 7   !3
   integer, parameter :: maxlambdas     = 101 !10001 !50001  !5001  !1001

!  Aerosol data Dimensioning

   integer, parameter :: aer_maxlambdas = 57
   integer, parameter :: aer_maxaggs    = 5
   integer, parameter :: aer_maxmoms    = 200 !800

!  Icecloud data Dimensioning. Updated 4/5/16  to 73

   integer, parameter :: ic_maxlambdas = 73
!   integer, parameter :: ic_maxlambdas = 3700 !73
   integer, parameter :: ic_maxaggs    = 20
   integer, parameter :: ic_maxmoms    = 100
!   integer, parameter :: ic_maxmoms    = 10

!  Water cloud data Dimensioning.

   integer, parameter :: wc_maxlambdas = 24
   integer, parameter :: wc_maxaggs    = 1
   integer, parameter :: wc_maxmoms    = 100

!  Dimensioning for the FOMS tools/routines
!  ========================================
!mick note 6/6/2023 - (1) MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT, MAX_ATMOSWFS, &
!                         MAXFINELAYERS defined in "MSSeg_pars.f90"
!                     (2) "MaxAllGeoms" defined above.

   integer, parameter   :: maxdisos    = MAXSTREAMS
   integer, parameter   :: maxlays     = MAXLAYERS
   integer, parameter   :: maxmoms     = MAXMOMENTS_INPUT
   integer, parameter   :: maxJacobs   = MAX_ATMOSWFS

   integer, parameter   :: maxfine     = MAXFINELAYERS
   integer, parameter   :: maxsegments = 200 !150        ! Caution, should = at least twice EXT_maxlayers
   integer, parameter   :: maxgeoms    = MaxAllGeoms

!  Numbers
!  =======
!mick mod 6/6/2023 - most previous parameters here eliminated (now using underlying
!                    parameters in "MSSeg_pars.f90").  Only unique FOMS parameters
!                    now active here.

   real(kind=fpk), parameter :: small = 1.0e-12_fpk
   real(kind=fpk), parameter :: dtr  = Pie / 180.0_fpk

end module FOMS_pars_m
