
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

      MODULE MSSeg_L_HPTypes_Def_m

      USE MSSEG_pars_m, only : MAXSTREAMS, MAXSTOKES, MAX_ATMOSWFS, MAXLAYERS, MAXBEAMS, MAX_USER_STREAMS,  &
                               MAXEVALUES, MAXSTRMSTKS, MAXSTREAMS_2, MAXSTRMSTKS_2, MAXTOTAL

      IMPLICIT NONE

!  Type structures

      TYPE L_Homog_Solutions

      DOUBLE PRECISION :: L_T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLA_XPOS ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_SOLB_XNEG ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_EIGENMAT  ( MAXEVALUES, MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_KEIGEN    ( MAXEVALUES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES,  MAXLAYERS,  MAX_ATMOSWFS )

      END TYPE L_Homog_Solutions

! #####################################################################
! #####################################################################

      TYPE L_Partic_Solutions

      DOUBLE PRECISION :: L_WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      DOUBLE PRECISION :: L_UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: L_UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LC_BVEC  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UPAR_DN_2 ( MAXSTOKES, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LC_UPAR_UP_2 ( MAXSTOKES, MAX_ATMOSWFS )

      DOUBLE PRECISION :: LP_BVEC  ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UPAR_DN_2 ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: LP_UPAR_UP_2 ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS )

      END TYPE L_Partic_Solutions

! #####################################################################
! #####################################################################

      TYPE L_Partic_BVPSolve

      DOUBLE PRECISION :: COL2_WF  ( MAXTOTAL, MAX_ATMOSWFS )
      DOUBLE PRECISION :: SCOL2_WF ( MAXSTRMSTKS_2, MAX_ATMOSWFS )

      DOUBLE PRECISION :: NCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )
      DOUBLE PRECISION :: PCON ( MAXSTRMSTKS, MAXLAYERS, MAX_ATMOSWFS )

      END TYPE L_Partic_BVPSolve

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: L_Homog_Solutions, L_Partic_BVPSolve, L_Partic_Solutions

      END MODULE MSSeg_L_HPTypes_Def_m

