
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

      MODULE MSSeg_HPTypes_Def_m

      USE MSSEG_pars_m, only : MAXSTREAMS, MAXSTOKES, MAXMOMENTS, MAXLAYERS, MAXEVALUES, MAXBEAMS, MAX_USER_STREAMS,  &
                               MAX_ALLSTRMS_P1, MAXSTRMSTKS, MAXSTREAMS_2, MAXSTRMSTKS_2, MAXSTOKES_SQ, MAXBANDTOTAL, MAXTOTAL

      IMPLICIT NONE

!  Type structures

! #####################################################################
! #####################################################################

      TYPE Homog_GSF

      DOUBLE PRECISION :: PIMM_11 ( MAXSTREAMS )
      DOUBLE PRECISION :: PIMM_KM ( MAXSTREAMS )

      DOUBLE PRECISION :: PI_XQP ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_POST ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQM_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XQP_PRE  ( 0:MAXMOMENTS, MAXSTREAMS, MAXSTOKES, MAXSTOKES )

!  Local matrices to be saved

      DOUBLE PRECISION :: XLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: YLM_DIAG ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION :: ZLM    ( 0:MAXMOMENTS, MAXSTOKES, MAXSTOKES )

      DOUBLE PRECISION :: PHI    ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PINORM ( 0:MAXMOMENTS, 0:MAXMOMENTS )
      DOUBLE PRECISION :: PISIGN ( 0:MAXMOMENTS, 0:MAXMOMENTS )

      DOUBLE PRECISION :: DF_L          ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: DF_LP1        ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: DF_2LP1       ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: DF_LSQM4      ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: DF_RT_LP3XLM1 ( 0:MAXMOMENTS )
      DOUBLE PRECISION :: ZHELP         ( 0:MAXMOMENTS )

      DOUBLE PRECISION :: UMXSQ         ( MAXSTREAMS )
      DOUBLE PRECISION :: RT_UMXSQ      ( MAXSTREAMS )
      DOUBLE PRECISION :: UPXSQ_D_UMXSQ ( MAXSTREAMS )
      DOUBLE PRECISION :: PLEG20        ( MAXSTREAMS )
      DOUBLE PRECISION :: M2X_D_UMXSQ   ( MAXSTREAMS )
      DOUBLE PRECISION :: XUMXSQ        ( MAXSTREAMS )
      DOUBLE PRECISION :: X_RT_UMXSQ    ( MAXSTREAMS )
      DOUBLE PRECISION :: X2LP1         ( 0:MAXMOMENTS, MAXSTREAMS )

      END TYPE Homog_GSF

! #####################################################################
! #####################################################################

!  04/07/23. Add Fourier dimension to PIMM_11 and PIMM_KM
 
      TYPE Partic_GSF

      DOUBLE PRECISION :: PIMM_11 ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )
      DOUBLE PRECISION :: PIMM_KM ( 0:MAXMOMENTS, MAX_ALLSTRMS_P1 )

      DOUBLE PRECISION :: PI_XUP ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_X0P ( 0:MAXMOMENTS, MAXBEAMS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUM_POST ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )
      DOUBLE PRECISION :: PI_XUP_PRE ( 0:MAXMOMENTS, MAX_USER_STREAMS, MAXSTOKES, MAXSTOKES )

      END TYPE Partic_GSF

! #####################################################################
! #####################################################################

      TYPE Homog_Solutions

      INTEGER ::          K_REAL       ( MAXLAYERS )
      INTEGER ::          K_COMPLEX    ( MAXLAYERS )
      DOUBLE PRECISION :: KEIGEN       ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: T_DELT_EIGEN ( MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: SOLA_XPOS    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: SOLB_XNEG    ( MAXSTREAMS_2, MAXSTOKES, MAXEVALUES, MAXLAYERS )

      DOUBLE PRECISION :: SAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: DAB ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION :: EIGENMAT_SAVE ( MAXEVALUES, MAXEVALUES, MAXLAYERS )
      DOUBLE PRECISION :: FWD_SUMVEC    ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: FWD_DIFVEC    ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )

      END TYPE Homog_Solutions

! #####################################################################
! #####################################################################

      TYPE Partic_Solutions

!mick fix 7/13/2015 - promoted HELPSTOKES back to an output variable
!                     for the linearized case (as in standard VLIDORT)

      DOUBLE PRECISION :: UHOM_DNDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: UHOM_DNUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: HELPSTOKES ( 0:MAXMOMENTS, MAXEVALUES, MAXSTOKES )

      DOUBLE PRECISION ::  QSUMVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION ::  QDIFVEC_SAVE ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION ::  QVEC_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION ::  QDIF_SAVE ( MAXSTRMSTKS )
      DOUBLE PRECISION ::  QMAT_SAVE ( MAXSTRMSTKS, MAXSTRMSTKS )

      INTEGER ::           QPIVOT ( MAXSTRMSTKS )
      DOUBLE PRECISION ::  BVEC   ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION ::  WUPPER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )
      DOUBLE PRECISION ::  WLOWER ( MAXSTREAMS_2, MAXSTOKES, MAXLAYERS )

      DOUBLE PRECISION ::  HELPSTOKES_BEAM ( 0:MAXMOMENTS, MAXSTOKES )
      DOUBLE PRECISION ::  UPAR_DN_2 ( MAXSTOKES )
      DOUBLE PRECISION ::  UPAR_UP_2 ( MAXSTOKES )

      END TYPE Partic_Solutions

! #####################################################################
! #####################################################################

      TYPE Homog_BVPSetup

      DOUBLE PRECISION :: BANDMAT2 ( MAXBANDTOTAL, MAXTOTAL )
      INTEGER ::          IPIVOT   ( MAXTOTAL )
      DOUBLE PRECISION :: SMAT2    ( MAXSTRMSTKS_2, MAXSTRMSTKS_2 )
      INTEGER ::          SIPIVOT  ( MAXSTRMSTKS_2 )

      DOUBLE PRECISION :: R2_HOMP ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: R2_HOMM ( MAXSTREAMS, MAXSTOKES, MAXEVALUES )
      DOUBLE PRECISION :: AXBID_F ( MAXSTREAMS, MAXSTREAMS, MAXSTOKES_SQ )

      END TYPE Homog_BVPSetup

! #####################################################################
! #####################################################################

      TYPE Partic_BVPSolve

      DOUBLE PRECISION ::  COL2  ( MAXTOTAL, MAXBEAMS )
      DOUBLE PRECISION ::  SCOL2 ( MAXSTRMSTKS_2, MAXBEAMS )

      DOUBLE PRECISION ::  R2_BEAM ( MAXSTREAMS, MAXSTOKES )
      DOUBLE PRECISION ::  LCON ( MAXSTRMSTKS, MAXLAYERS )
      DOUBLE PRECISION ::  MCON ( MAXSTRMSTKS, MAXLAYERS )

      END TYPE Partic_BVPSolve

! #####################################################################
! #####################################################################

      TYPE Homogeneous_RTE

      TYPE(Homog_GSF)      :: HGSF
      TYPE(Homog_BVPSetup) :: HBVP
      TYPE(Homog_Solutions):: HSOL

      END TYPE Homogeneous_RTE

! #####################################################################
! #####################################################################

! #####################################################################
! #####################################################################

      TYPE Particular_RTE

      TYPE(Partic_GSF)       :: PGSF
      TYPE(Partic_BVPSolve)  :: PBVP
      TYPE(Partic_Solutions) :: PSOL

      END TYPE Particular_RTE

! #####################################################################
! #####################################################################

      PRIVATE
      PUBLIC :: Homog_GSF,  Homog_BVPSetup,  Homog_Solutions,  Homogeneous_RTE, &
                Partic_GSF, Partic_BVPSolve, Partic_Solutions, Particular_RTE

      END MODULE MSSeg_HPTypes_Def_m

