
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2024.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

!  Two subroutines here.
!    (1) use the Mie expansion coefficients to develop Fmatrix entries
!    (2) Use the Greekmatrix coefficients to develop Tmatrix entries
!  In both cases, F-matrix is expressed as a vector of the 6 non-zero entries

!  12/12/23. Linearization of the Use-Coefficient subroutines (New)
!    -- When the expansion coefficients have derivatives w.r.t various optical parameters,
!       then we use these coefficient derivatives to get the linearized Fmatrix.
!    -- Similarly with the Greek-matrix coefficients.
!    -- These subroutines do both tasks (regular and linearized coefficients/Fmatrices).

module vfzmat_L_UseCoeffs_m

public

contains

SUBROUTINE vfzmat_L_UseCoeffs_Fmatrices &
      ( max_InAngles, Max_Coeffs, Max_wfs, n_InAngles, ncoeffs, nstokes, n_wfs, & ! Input numbers
        InCosines, expcoeffs, L_expcoeffs, Fmatrices, L_Fmatrices )               ! Input/Outputs

!  Original Notes for the regular Fmatrices subroutine
!  ---------------------------------------------------

!  Use the expansion coefficients of the scattering matrix in 
!  generalized spherical functions to expand F matrix

!  Stand-alone routine to expand Fmatrices from expansion coefficients
!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER           , INTENT (IN) :: max_InAngles, Max_Coeffs, Max_wfs
   INTEGER           , INTENT (IN) :: n_InAngles, ncoeffs, nstokes, n_wfs
   REAL    (KIND=dpk), INTENT (IN) :: InCosines(Max_InAngles)
   REAL    (KIND=dpk), INTENT (IN) :: expcoeffs  (0:Max_Coeffs,6)
   REAL    (KIND=dpk), INTENT (IN) :: L_expcoeffs(Max_wfs,0:Max_Coeffs,6)

!  output, initialized here, so intent(out)

   REAL    (KIND=dpk), INTENT (OUT) :: FMatrices  (Max_InAngles,6)
   REAL    (KIND=dpk), INTENT (OUT) :: L_FMatrices(Max_Wfs,Max_InAngles,6)

!  local variables

   REAL    (KIND=dpk) :: P00(2), P02(2), P2p2(2), P2m2(2)
   REAL    (KIND=dpk) :: fmat(6), L_fmat(max_wfs,6)

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   real(dpk), parameter :: d_three = 3.0_dpk, d_four = 4.0_dpk

   INTEGER            :: l, k, lnew, lold, itmp, q
   INTEGER            :: index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk) :: dl, dl1, qroot6, fac1, fac2, uuu, FL2, FLL1, PERLL4, QFAC, WFACT, &
                         sql4, sql41, tmp1, tmp2, GK11, GK12, GK34, GK44, GK22, GK33, SUM23, DIF23, &
                         L_GK11(max_wfs), L_GK12(max_wfs), L_GK34(max_wfs), L_GK44(max_wfs), L_GK22(max_wfs), L_GK33(max_wfs) 

!  Initialization

   qroot6 = -0.25_dpk*SQRT(6.0_dpk)
   FMatrices   = d_zero
   L_FMatrices = d_zero

!  Indices

   index_11 = 1
   index_22 = 2
   index_33 = 3
   index_44 = 4
   index_12 = 5
   index_34 = 6

!  START LOOP OVER IN COSINES

   DO K = 1, N_InAngles

!  Cosine of the scattering angle
      
      FMAT = D_ZERO ; L_FMAT = D_ZERO
      UUU = InCosines(N_InAngles+1-k)

!  START LOOP OVER THE COEFFICIENT INDEX L
!  FIRST UPDATE GENERALIZED SPHERICAL FUNCTIONS, THEN CALCULATE COEFS.
!  LOLD AND LNEW ARE POINTER-LIKE INDICES USED IN RECURRENCE

      LNEW = 1
      LOLD = 2

      DO L = 0, NCOEFFS

        DL   = REAL(L,dpk)
        DL1  = DL - d_one

!  SET THE LOCAL COEFFICIENTS
!   44 AND 34 ARE NOT REQUIRED WITH NATURAL SUNLIGHT (DEFAULT HERE)
!   22 AND 33 REQUIRED FOR NON-MIE SPHEROIDAL PARTICLES

        GK11 = EXPCOEFFS(L,1) ; L_GK11(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,1)
        IF ( NSTOKES .GT. 1 ) THEN
          GK22 = + EXPCOEFFS(L,2) ; L_GK22(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,2)
          GK33 = + EXPCOEFFS(L,3) ; L_GK33(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,3)
          GK44 = + EXPCOEFFS(L,4) ; L_GK44(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,4)
          GK12 = + EXPCOEFFS(L,5) ; L_GK12(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,5)
          GK34 = - EXPCOEFFS(L,6) ; L_GK34(1:n_wfs) = L_EXPCOEFFS(1:n_wfs,L,6)
        ENDIF

!  FIRST MOMENT

        IF ( L .EQ. 0 ) THEN

!  ADDING PAPER EQS. (76) AND (77) WITH M=0
!   ADDITIONAL FUNCTIONS P2M2 AND P2P2 ZERO FOR M = 0

          P00(LOLD) = d_one
          P00(LNEW) = d_zero
          P02(LOLD) = d_zero
          P02(LNEW) = d_zero
          P2P2(LOLD) = d_zero
          P2P2(LNEW) = d_zero
          P2M2(LOLD) = d_zero
          P2M2(LNEW) = d_zero

        ELSE

          FAC1 = (d_two*DL-d_one)/DL
          FAC2 = DL1/DL

! ADDING PAPER EQ. (81) WITH M=0

          P00(LOLD) = FAC1*UUU*P00(LNEW) - FAC2*P00(LOLD)

        END IF

        IF ( L .EQ. 2 ) THEN

! ADDING PAPER EQ. (78)
! SQL4 CONTAINS THE FACTOR DSQRT((L+1)*(L+1)-4) NEEDED IN
! THE RECURRENCE EQS. (81) AND (82)

          P02(LOLD) = QROOT6*(d_one-UUU*UUU)
          P02(LNEW) = d_zero
          SQL41     = d_zero

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L = 2

          P2P2(LOLD)= 0.25_dpk*(d_one+UUU)*(d_one+UUU)
          P2M2(LOLD)= 0.25_dpk*(d_one-UUU)*(d_one-UUU)

        ELSE IF ( L .GT. 2) THEN

! ADDING PAPER EQ. (82) WITH M=0

          SQL4  = SQL41
          SQL41 = SQRT(DL*DL-d_four)
          TMP1  = (d_two*DL-d_one)/SQL41
          TMP2  = SQL4/SQL41
          P02(LOLD) = TMP1*UUU*P02(LNEW) - TMP2*P02(LOLD)

!  INTRODUCE THE P2P2 AND P2M2 FUNCTIONS FOR L > 2

          FL2 = d_two * DL - d_one
          FLL1 = DL * DL1
          PERLL4=d_one/(DL1*SQL41**d_two)
          QFAC  = DL  * ( DL1*DL1 - d_four)
          WFACT = FL2 * ( FLL1 * UUU - d_four )
          P2P2(LOLD) = (WFACT*P2P2(LNEW) - QFAC*P2P2(LOLD)) * PERLL4
          WFACT = FL2 * ( FLL1 * UUU + d_four )
          P2M2(LOLD) = (WFACT*P2M2(LNEW) - QFAC*P2M2(LOLD)) * PERLL4

        END IF

! SWITCH INDICES SO THAT LNEW INDICATES THE FUNCTION WITH
! THE PRESENT INDEX VALUE L, THIS MECHANISM PREVENTS SWAPPING
! OF ENTIRE ARRAYS.

        ITMP = LNEW
        LNEW = LOLD
        LOLD = ITMP

! NOW ADD THE L-TH TERM TO THE SCATTERING MATRIX.
! SEE DE HAAN ET AL. (1987) EQS. (68)-(73).

! SECTION FOR RANDOMLY-ORIENTED SPHEROIDS, ADDED 05 OCTOBER 2010
!  R. SPURR AND V. NATRAJ

        IF ( L.LE.NCOEFFS ) THEN

!  regular

          FMAT(INDEX_11) = FMAT(INDEX_11) + GK11 * P00(LNEW)
          FMAT(INDEX_12) = FMAT(INDEX_12) + GK12 * P02(LNEW)
          SUM23 = GK22 + GK33
          DIF23 = GK22 - GK33
          FMAT(INDEX_22) = FMAT(INDEX_22) + SUM23 * P2P2(LNEW)
          FMAT(INDEX_33) = FMAT(INDEX_33) + DIF23 * P2M2(LNEW)
          FMAT(INDEX_44) = FMAT(INDEX_44) + GK44 * P00(LNEW)
          FMAT(INDEX_34) = FMAT(INDEX_34) + GK34 * P02(LNEW)

!  Linearized

          do q = 1, n_wfs
            L_FMAT(q,INDEX_11) = L_FMAT(q,INDEX_11) + L_GK11(q) * P00(LNEW)
            L_FMAT(q,INDEX_12) = L_FMAT(q,INDEX_12) + L_GK12(q) * P02(LNEW)
            SUM23 = L_GK22(q) + L_GK33(q)
            DIF23 = L_GK22(q) - L_GK33(q)
            L_FMAT(q,INDEX_22) = L_FMAT(q,INDEX_22) + SUM23 * P2P2(LNEW)
            L_FMAT(q,INDEX_33) = L_FMAT(q,INDEX_33) + DIF23 * P2M2(LNEW)
            L_FMAT(q,INDEX_44) = L_FMAT(q,INDEX_44) + L_GK44(q) * P00(LNEW)
            L_FMAT(q,INDEX_34) = L_FMAT(q,INDEX_34) + L_GK34(q) * P02(LNEW)
          enddo

        ENDIF

!  END COEFFICIENT LOOP

      END DO

!   THIS MUST BE DONE AFTER THE MOMENT LOOP.

      FMAT(INDEX_22) = d_half * ( FMAT(INDEX_22) + FMAT(INDEX_33) )
      FMAT(INDEX_33) = FMAT(INDEX_22) - FMAT(INDEX_33)
      do q = 1, n_wfs
        L_FMAT(q,INDEX_22) = d_half * ( L_FMAT(q,INDEX_22) + L_FMAT(q,INDEX_33) )
        L_FMAT(q,INDEX_33) = L_FMAT(q,INDEX_22) - L_FMAT(q,INDEX_33)
      enddo

! REMEMBER FOR MIE SCATTERING : F11 = F22 AND F33 = F44
!  THIS CODE IS NO LONGER REQUIRED, AS WE HAVE INTRODUCED CODE NOW
!   FOR RANDOMLY ORIENTED SPHEROIDS. THE SYMMETRY SHOULD STILL OF
!   COURSE BE PRESENT FOR THE MIE PARTICLES, SO THIS WILL BE A
!   CHECK ON THE NEW CODE. R. SPURR AND V. NATRAJ,, 20 MARCH 2006
!        FMAT(INDEX_22) = FMAT(INDEX_11)
!        FMAT(INDEX_33) = FMAT(INDEX_44)

!  Assign output

      FMatrices(K,1:6) = FMAT(1:6)
      do q = 1, n_wfs
        L_FMatrices(q,K,1:6) = L_FMAT(q,1:6)
      enddo

!  End geometry loop

   ENDDO

!  Done

  RETURN
END SUBROUTINE vfzmat_L_UseCoeffs_Fmatrices

!

SUBROUTINE vfzmat_L_UseGkMats_Fmatrices &
   ( maxcoeffs, maxlayers, maxgeoms, max_wfs, nstokes, ncoeffs, nlayers, ngeoms, n_wfs, & ! Input Numbers
     do_sunlight, greekmats, L_greekmats, genspher, FMatrices, L_FMatrices )              ! Inputs/Outputs                   

!  This subroutine allows layer Greek Matrices to have their own derivatives. 
!   -- Thus, this subroutine could be used for LC or LP (column or profile) linearizations

!  Original Notes for the regular use-Greekmat subroutine
!  ------------------------------------------------------

!  Stand-alone routine to define the nonzero entries of the F-Matrix using
!  the generalized spherical function expansion coefficients of the nonzero
!  entries of the Greek matrix

!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!  Programmed 3/2/20 by M. Christi
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

   IMPLICIT NONE

!  Parameters
!  ----------

!  Precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Input
!  -----

!  Control

   INTEGER,   INTENT(IN) :: maxcoeffs, maxlayers, maxgeoms, max_wfs
   INTEGER,   INTENT(IN) :: nstokes, ncoeffs, nlayers, ngeoms, n_wfs(maxlayers)

   LOGICAL,   INTENT(IN) :: do_sunlight

!  Coefficients of the nonzero entries of the 4x4 Greek Matrix (vector version)

   REAL(dpk), INTENT(IN) :: Greekmats   ( 0:maxcoeffs, maxlayers, 16 )
   REAL(dpk), INTENT(IN) :: L_Greekmats ( max_wfs, 0:maxcoeffs, maxlayers, 16 )

!  Generalized spherical functions.
!    Rotations(1-4) = C1, S1, C2, S2

   REAL(dpk), INTENT(IN) :: genspher ( 0:maxcoeffs, 4, maxgeoms )

!  Output
!  ------

!  Nonzero entries of the F-matrix (vector version)

   REAL(dpk), INTENT(OUT) :: Fmatrices   ( maxlayers, maxgeoms, 6 )
   REAL(dpk), INTENT(OUT) :: L_Fmatrices ( max_wfs, maxlayers, maxgeoms, 6 )

!  Local variables
!  ---------------

   REAL(dpk), PARAMETER :: d_zero  = 0.0_dpk, d_half  = 0.5_dpk

   INTEGER   :: g1, g2, g3, g4, g5, g6, gmask(6), n, v, q
   REAL(dpk) :: fmat3, fmat4, gmsum(0:maxcoeffs), gmdif(0:maxcoeffs)

!  Initialization

   Fmatrices = d_zero ; L_Fmatrices = d_zero

!  Define some local variables

   !Entries in this version of gmask correspond to the following positions
   !in the 4x4 Greek Matrix:
   !          11, 12, 22, 33, 34, 44         
   gmask = (/  1,  2,  6, 11, 12, 16 /)

!  Get F-matrix

   if ( nstokes .eq. 1 ) then

!  Scalar case

     g1 = gmask(1)
     do v = 1, ngeoms
       do n = 1, nlayers
         Fmatrices(n,v,1) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
         do q = 1, n_wfs(n)
           L_Fmatrices(q,n,v,1) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
         enddo
       enddo
     enddo

   else

!  Vector cases (nstokes > 1)

     if ( do_sunlight ) then

!    Vector - with sunlight

       g1 = gmask(1); g2 = gmask(2)
       do v = 1, ngeoms
         do n = 1, nlayers
           Fmatrices(n,v,1) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
           Fmatrices(n,v,2) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g2),genspher(0:ncoeffs,2,v))
           do q = 1, n_wfs(n)
             L_Fmatrices(q,n,v,1) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
             L_Fmatrices(q,n,v,2) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g2),genspher(0:ncoeffs,1,v))
           enddo
         enddo
       enddo

     else

!    Vector - general case with possible linear and circular polarization
!      prepare to use full 4x4 matrix; code introduced but not tested, 05 October 2010

       g1 = gmask(1); g2 = gmask(2); g3 = gmask(3); g4 = gmask(4)
       do v = 1, ngeoms
         do n = 1, nlayers

           Fmatrices(n,v,1) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
           Fmatrices(n,v,2) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g2),genspher(0:ncoeffs,2,v))
           gmsum(0:ncoeffs) = Greekmats(0:ncoeffs,n,g3) + Greekmats(0:ncoeffs,n,g4)
           fmat3            = DOT_PRODUCT(gmsum(0:ncoeffs),genspher(0:ncoeffs,3,v))
           gmdif(0:ncoeffs) = Greekmats(0:ncoeffs,n,g3) - Greekmats(0:ncoeffs,n,g4)
           fmat4            = DOT_PRODUCT(gmdif(0:ncoeffs),genspher(0:ncoeffs,4,v))
           Fmatrices(n,v,3)   = ( fmat3 + fmat4 ) * d_half
           Fmatrices(n,v,4)   = ( fmat3 - fmat4 )
           
           do q = 1, n_wfs(n)
             L_Fmatrices(q,n,v,1) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g1),genspher(0:ncoeffs,1,v))
             L_Fmatrices(q,n,v,2) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g2),genspher(0:ncoeffs,1,v))
             gmsum(0:ncoeffs) = L_Greekmats(q,0:ncoeffs,n,g3) + L_Greekmats(q,0:ncoeffs,n,g4)
             fmat3            = DOT_PRODUCT(gmsum(0:ncoeffs),genspher(0:ncoeffs,3,v))
             gmdif(0:ncoeffs) = L_Greekmats(q,0:ncoeffs,n,g3) - L_Greekmats(q,0:ncoeffs,n,g4)
             fmat4            = DOT_PRODUCT(gmdif(0:ncoeffs),genspher(0:ncoeffs,4,v))
             L_Fmatrices(q,n,v,3)   = ( fmat3 + fmat4 ) * d_half
             L_Fmatrices(q,n,v,4)   = ( fmat3 - fmat4 )
           enddo

           if ( nstokes.eq.4 ) then
             g5 = gmask(5); g6 = gmask(6)
             Fmatrices(n,v,5) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g5),genspher(0:ncoeffs,2,v))
             Fmatrices(n,v,6) = DOT_PRODUCT(Greekmats(0:ncoeffs,n,g6),genspher(0:ncoeffs,1,v))
             do q = 1, n_wfs(n)
               L_Fmatrices(q,n,v,5) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g5),genspher(0:ncoeffs,2,v))
               L_Fmatrices(q,n,v,6) = DOT_PRODUCT(L_Greekmats(q,0:ncoeffs,n,g6),genspher(0:ncoeffs,1,v))
             enddo
           endif
         enddo
       enddo

     endif  !End sunlight IF block
   endif  !End stokes IF block

END SUBROUTINE vfzmat_L_UseGkMats_Fmatrices

!  End module

End Module vfzmat_L_UseCoeffs_m

