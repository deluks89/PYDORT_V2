
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

      MODULE VQS_LinOutputs_def_m

!  This module contains the following VQS output structures,
!  with intents :

!         VQS_LinProf    nested in VQS_LinOutputs
!         VQS_LinCol     nested in VQS_LinOutputs (Added 9/21/23)
!      VQS_LinOutputs    Intent(Out)

!  1/31/21. Version 2.8.3. Add MAXLAYERS to this list (for MSST output)

      USE MSSeg_PARS_m, Only : fpk, MAX_USER_LEVELS, MAX_GEOMETRIES, &
        MAX_SZANGLES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAXLAYERS, &
        MAXLAYERS, MAXSTOKES, MAX_ATMOSWFS, MAX_DIRECTIONS
        
      USE FOMS_PARS_m, Only : MAXSEGMENTS

      IMPLICIT NONE

! #####################################################################
! #####################################################################

      TYPE VQS_LinProf

!  Stokes vector values: Total, First-Order (FO), and Multiple Scatter (MS)

      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_StokesTot
      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_StokesFO
      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_StokesMS

!  Actinic and Regular Fluxes: Diffuse (Diff) & Direct (Dir)

      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_ActFluxDiff
      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_RegFluxDiff

      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_ActFluxDir
      REAL(fpk), DIMENSION ( MAXSTOKES, MAXLAYERS, MAX_ATMOSWFS ) :: LP_RegFluxDir

!  4/24/24. Fluxes defined at segment boundaries along the light path

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAXLAYERS, MAX_ATMOSWFS ) :: LP_ActFluxDiff_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAXLAYERS, MAX_ATMOSWFS ) :: LP_RegFluxDiff_AllSB

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAXLAYERS, MAX_ATMOSWFS ) :: LP_ActFluxDir_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAXLAYERS, MAX_ATMOSWFS ) :: LP_RegFluxDir_AllSB

      END TYPE VQS_LinProf

! #####################################################################
! #####################################################################

!  9/21/23. Added type structure

      TYPE VQS_LinCol

!  Stokes vector values: Total, First-Order (FO), and Multiple Scatter (MS)

      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_StokesTot
      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_StokesFO
      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_StokesMS

!  Actinic and Regular Fluxes: Diffuse (Diff) & Direct (Dir)

      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_ActFluxDiff
      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_RegFluxDiff

      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_ActFluxDir
      REAL(fpk), DIMENSION ( MAXSTOKES, MAX_ATMOSWFS ) :: LC_RegFluxDir

!  4/24/24. Fluxes defined at segment boundaries along the light path

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAX_ATMOSWFS ) :: LC_ActFluxDiff_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAX_ATMOSWFS ) :: LC_RegFluxDiff_AllSB

      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAX_ATMOSWFS ) :: LC_ActFluxDir_AllSB
      REAL(fpk), DIMENSION ( MAXSTOKES, 0:MAXSEGMENTS, MAX_ATMOSWFS ) :: LC_RegFluxDir_AllSB

      END TYPE VQS_LinCol

! #####################################################################
! #####################################################################

      TYPE VQS_LinOutputs

      TYPE(VQS_LinProf)  :: Prof
      TYPE(VQS_LinCol)   :: Col

      END TYPE VQS_LinOutputs

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

      PRIVATE
      PUBLIC :: VQS_LinProf, VQS_LinCol, &
                VQS_LinOutputs

      END MODULE VQS_LinOutputs_def_m
