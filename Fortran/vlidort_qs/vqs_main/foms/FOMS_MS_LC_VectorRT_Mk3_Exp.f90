
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
! # MS (Multiple-scatter) Vector RT for FOMS                    #
! #                                                             #
! #        FOMS_MS_LC_VectorRT_Mk3_Exp                          #
! #                                                             #
! ###############################################################

!  Based on "Titan-Probe" FOMS code developed by R. Spurr in May 2012
!  This is a recapitulation and upgrade (VLIDORT Version 2.7)
!  All references to "Probe" have been replaced by "Attrex", now FOMS

!  Code uses the specialized MSSEG model
!  4/7/23. Mk3 Upgrade to VLIDORT Version 2.8.3 standards

module FOMS_MS_LC_VectorRT_Mk3_Exp_m

!mick fix 7/13/2015 - commented out since these already defined below in subroutine
   !use FOMS_pars_m
   !use vlidort_MSSEG_masters_m

!  For a given wavelength, this routine will calculate MS-only Stokes-Vector at Camera
!  This is based on Precalculated Geometrical quantities and appropriate Optical properties.
!       This will perform Enhanced-PS (LOS-path sphericity) calculations, using multiplier
!       results (SFI) from the main model FO calculation.

!  7/21/21. Mk2. Set additional MSSEG control
!   -- Flux output ontrolled by DO_NCHOICE_FLUXES input flag
!   -- ASYMTX_TOLERANCE flag fixed by hand
!   -- DO_FOURIER0_NSTOKES2 performance flag, can be turned on for Lambertian surfaces

!  9/21/21. Mk2a. New module for the Column Jacobians

!  4/24/24. Add "OutputDiffuse-fluxes_SB" and "LC_OutputDiffuse-fluxes_SB" 
!      -- (SEGMENT BOUNDARY) arrays as Output variables (2 extra lines)
!      -- This output was formerly written to file (fort.81) as part of a flux debug

public  :: FOMS_MS_LC_VectorRT_Mk3_Exp
private :: FOMS_WRITE_MS_LC_INPUTS, FOMS_MS_LC_Convergence_Mk3

contains

subroutine FOMS_MS_LC_VectorRT_Mk3_Exp ( do_debug_MS_input, RT_Progress,   & ! Inputs monitoring
       do_plane_parallel, do_Rayleigh_only, do_deltam_scaling,             & ! Inputs RT control flags
       do_lambertian_surface, do_sleave, do_sl_isotropic,                  & ! Inputs (surface flags)
       do_calculate_fluxes, do_debug_fluxes, Regime_idx,                   & ! Inputs (flux-flags/Regime)
       nstreams, nlayers, nstokes, ngreek_moments_input, eradius, heights, & ! Inputs (integers/heights)
       nsegments, seglayer_map, segupdown_map, alpha, theta, Relaz,        & ! Inputs (segments/geometry)
       nfinedivs, alphafine, thetafine, relazfine,                         & ! Inputs (fine-division control)
       do_Column_Wfs, n_Column_wfs, Lvaryflags, Lvarynumbers, npert, eps,  & ! Inputs Lin Control
       lostrans, sfi_factors, sfi_opdeps,                                  & ! Inputs Optical
       L_lostrans, L_sfi_factors, L_sfi_opdeps,                            & ! Inputs (SS Linearized)
       reflec, deltaus, omega, greekmat, solarflux,                        & ! Inputs Optical
       L_omega, L_deltaus, L_greekmat,                                     & ! Inputs (Optical Linearized)
       stokes_MS,    stokes_TOT,    Output_Regular_Fluxes,    Output_Actinic_Fluxes,    & ! Outputs    Stokes/fluxes
       LC_stokes_MS, LC_stokes_TOT, LC_Output_Regular_Fluxes, LC_Output_Actinic_Fluxes, & ! Outputs LC Stokes/fluxes
       OutputDiffuse_Actinic_fluxes_SB,    OutputDiffuse_Regular_fluxes_SB,             & ! Fluxes_SB    Outputs. Added, 4/24/24
       LC_OutputDiffuse_Actinic_fluxes_SB, LC_OutputDiffuse_Regular_fluxes_SB,          & ! LC_Fluxes_SB Outputs. Added, 4/24/24
       fail, status, nmessages, messages, actions )                                       ! Exception handling

!  Stand alone routine for MS field with Solar sources alone

!  7/21/21. Mk2. Some changes:
!    -- Introduce control for calculating fluxes. I/O list cleaned up and annotated
!    -- Output regular and actinic fluxes if control set
!    -- Internal performance control for VQS Call (ASYMTX_TOLERANCE, DO_FOURIER0_NSTOKES2) 

!  1/24/22. Use homogeneous and particular-integral modules and type structures

!  4/24/24. Add "OutputDiffuse-fluxes_SB" and "LC_OutputDiffuse-fluxes_SB" 
!      -- (SEGMENT BOUNDARY) arrays as Output variables (2 extra lines)
!      -- This output was formerly written to file (fort.81) as part of a flux debug

   use FOMS_pars_m
   USE MSSeg_PARS_m, Only : MAXSTREAMS, MAXLAYERS, MAXMOMENTS_INPUT, MAXSTOKES, MAX_MESSAGES, MAXBEAMS, &
                            MAX_USER_STREAMS, MAX_GEOMETRIES, MAX_ALLSTRMS_P1, MAXMOMENTS, MAXFOURIER, &
                            MAXSTOKES_SQ, MAX_ATMOSWFS, DEG_TO_RAD, VQS_SUCCESS, VQS_SERIOUS

   use MSSeg_HPTypes_Def_m
   use MSSeg_L_HPTypes_Def_m

   use msseg_homog_miscsetups_m
   use msseg_homog_Fourier_m
   use msseg_partic_miscsetups_m
   use msseg_partic_Fourier_m

   use MSSeg_LA_miscsetups_m
   use MSSeg_LA_solutions_m

   use MSSeg_LC_miscsetups_m
   use MSSeg_LC_bvproblem_m
   use MSSeg_LC_solutions_m
   use MSSeg_LC_PostProcess_m

   use MSSeg_LA_Homog_Fourier_m
   use MSSeg_LC_Fourier_m

   implicit none         

!  Inputs
!  ======

!  debug input and progress flags

   logical  , intent(in)     :: do_debug_MS_input
   logical  , intent(in)     :: RT_Progress

!  RT control flags

   logical  , intent(in)     :: do_plane_parallel
   logical  , intent(inout)  :: do_Rayleigh_only
   logical  , intent(inout)  :: do_deltam_scaling

!  surface flags

   logical  , intent(in)  :: do_lambertian_surface
   logical  , intent(in)  :: do_sleave
   logical  , intent(in)  :: do_sl_isotropic

!  7/21/21.Mk2. Introduce flux flags

   logical  , intent(in)  :: do_calculate_fluxes
   logical  , intent(in)  :: do_debug_fluxes

!  number of  layers, streams, Stokes, regime index

   integer  , intent(in)  :: nlayers, nstreams, nstokes, ngreek_moments_input, Regime_idx

!  Jacobian control

   LOGICAL, Intent(in) ::  Do_Column_WFs
   INTEGER, Intent(in) ::  N_Column_WFs
   LOGICAL, Intent(in) ::  Lvaryflags   (maxlays)
   INTEGER, Intent(in) ::  Lvarynumbers (maxlays)

!  Rob add 5/30/18, Perturbation layer and eps.(for debug only)

   INTEGER  , Intent(in) :: npert
   REAL(ffp), Intent(in) :: eps

!  Grid

   real(ffp), intent(in)  :: eradius, heights (0:maxlays)

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   integer  , intent(in)  :: nsegments
   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)
!   logical  , intent(in)  :: layermask_map (maxlays,2)

!  Number of fine divisions

   integer  , intent(in)  :: nfinedivs (maxsegments)

!  Angles 

   real(ffp), intent(in)  :: alpha  (0:maxsegments)
   real(ffp), intent(in)  :: theta  (0:maxsegments)
   real(ffp), intent(in)  :: relaz  (0:maxsegments)

!  Quadrature values

   real(ffp), intent(in)  :: alphafine (maxsegments,maxfine)
   real(ffp), intent(in)  :: thetafine (maxsegments,maxfine)
   real(ffp), intent(in)  :: Relazfine (maxsegments,maxfine)

!  Source function integration stuff

   real(ffp), Intent(in)  :: lostrans    ( maxsegments )
   real(ffp), Intent(in)  :: sfi_factors ( maxsegments, maxfine )
   real(ffp), Intent(in)  :: sfi_opdeps  ( maxsegments, maxfine )

!  Linearizations of source function integration material

   real(ffp), Intent(in)  :: L_lostrans    ( maxsegments, maxJacobs )
   real(ffp), Intent(in)  :: L_sfi_factors ( maxsegments, maxfine, maxJacobs )
   real(ffp), Intent(in)  :: L_sfi_opdeps  ( maxsegments, maxfine, maxJacobs )

!  Atmosphere optical inputs

   REAL(ffp), Intent(in) :: DELTAUS     ( maxlays )
   REAL(ffp), Intent(in) :: OMEGA       ( maxlays )
!   REAL(ffp), Intent(in) :: TRUNCFAC    ( maxlays )
   REAL(ffp), Intent(in) :: GREEKMAT    ( maxlays, 0:maxmoms, 4, 4 )

!  Linearized optical inputs

   REAL(ffp), Intent(in) :: L_DELTAUS     ( maxlays, maxJacobs )
   REAL(ffp), Intent(in) :: L_OMEGA       ( maxlays, maxJacobs )
   REAL(ffp), Intent(in) :: L_GREEKMAT    ( maxlays, 0:maxmoms, 4, 4, maxJacobs )

!  Solar Flux and Surface reflectivity (Regime 1, Lambertian-only)

   REAL(ffp), Intent(in) :: REFLEC(4,4), SOLARFLUX

!  BRDF supplemental material (not active currently). Will be pre-zeroed
!   -- Defined for each Fourier component locally

   real(ffp)  :: brdf_f      ( 16, maxdisos, maxdisos )
   real(ffp)  :: brdf_f_0    ( 16, maxdisos, maxfine  )
   real(ffp)  :: user_brdf_f ( 16, maxfine,  maxdisos )
  
!  SLEAVE supplemental material (not active currently). Will be pre-zeroed
!   -- Defined for each Fourier component locally

   real(ffp)  :: slterm_ISOTROPIC ( 4, maxfine )
   real(ffp)  :: slterm_f_0       ( 4, maxdisos, maxfine  )

!  outputs
!  -------

!  Stokes MS output 

   real(ffp), Intent(inout) :: stokes_MS(4) , LC_Stokes_MS  ( 4, maxJacobs )
   real(ffp), Intent(inout) :: stokes_TOT(4), LC_Stokes_TOT ( 4, maxJacobs )

!  7/21/21. MSSEG fluxes from the last segment calculation
!    -- 5/1/24. Must be Intent(InOut), otherwise you lose it

   real(ffp), Intent(inout) :: Output_Regular_Fluxes(4)
   real(ffp), Intent(inout) :: Output_Actinic_Fluxes(4)

!  8/19/21. Mk2, MSSEG Linearized fluxes from the last segment calculation
!    -- 9/21/23. Reversed ordering.
!    -- 5/1/24. Must be Intent(InOut), otherwise you lose it

   real(ffp), Intent(inout) :: LC_Output_Regular_Fluxes(4,maxJacobs)
   real(ffp), Intent(inout) :: LC_Output_Actinic_Fluxes(4,maxJacobs)

!  4/24/24.  Add Diffuse fluxes_SB and LC fluxes_SB arrays as Output.
!    -- Formerly written to dump file fort.81
!    -- 5/1/24. Must be Intent(InOut), otherwise you lose it

   real(ffp), Intent(inout) :: OutputDiffuse_Actinic_fluxes_SB(4,0:maxsegments)
   real(ffp), Intent(inout) :: OutputDiffuse_Regular_fluxes_SB(4,0:maxsegments)
   real(ffp), Intent(inout) :: LC_OutputDiffuse_Actinic_fluxes_SB(4,0:maxsegments,maxJacobs)
   real(ffp), Intent(inout) :: LC_OutputDiffuse_Regular_fluxes_SB(4,0:maxsegments,maxJacobs)

!  Exception handling
!mick fix 7/13/2015 - changed messages & actions from 100 to 25 elements to match the driver
!mick fix 6/6/2023  - changed messages & actions from 25 to 0:MAX_MESSAGES to make more general

   logical      , intent(inout) :: fail
   integer      , intent(inout) :: status, nmessages
   character*(*), intent(inout) :: messages(0:MAX_MESSAGES)
   character*(*), intent(inout) :: actions (0:MAX_MESSAGES)

!  LOCAL. USE VQS DIMENSIONING
!  ---------------------------

!  Additional Segmaster inputs

   logical   :: do_upwelling, do_surface_msst, do_double_convtest
   real(ffp) :: earth_radius, height_grid(0:maxlayers)
   real(ffp) :: flux_factor
   integer   :: N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES
   real(ffp) :: SZANGLES     (Max_geometries)
   real(ffp) :: USER_VZANGLES(Max_geometries)
   real(ffp) :: USER_RELAZMS (Max_geometries)

!  7/21/21. Mk2. Local flux control, internal VQS settings

   logical   :: do_nchoice_fluxes
   logical   :: DO_FOURIER0_NSTOKES2
   real(ffp) :: ASYMTX_TOLERANCE

!  Dummies

   real(ffp) :: sfi_opdeps_dum          ( Max_geometries )
   real(ffp) :: sfi_opdeps_dum_scaled   ( Max_geometries )
   real(ffp) :: L_sfi_opdeps_dum        ( MaxJacobs, Max_geometries )
   real(ffp) :: L_sfi_opdeps_dum_scaled ( MaxJacobs, Max_geometries )

!  VQS Optical inputs. Use VQS dimensioning

   real(ffp) :: deltau_vert_input (maxlayers)
   real(ffp) :: omega_total_input (maxlayers)
   real(ffp) :: greekmat_total_input (0:maxmoments, maxlayers, 16)
   real(ffp) :: ALBEDO

   real(ffp) :: L_DELTAU_VERT_INPUT ( maxJacobs, maxlays )
   real(ffp) :: L_OMEGA_TOTAL_INPUT ( maxJacobs, maxlays )
   real(ffp) :: L_GREEKMAT_TOTAL_INPUT ( maxJacobs, 0:maxmoments, maxlays, 16 )

!   --7/21/21. Mk2. Introduce local flux output from MSSEG

   real(ffp)  :: LOCAL_FLUXES    ( Max_geometries, 4, 2, 2 )
   real(ffp)  :: LOCAL_LC_FLUXES ( maxJacobs, Max_geometries, 4, 2, 2 )

!  Fourier components and Fourier-summed values

   real(ffp)  :: MSST_F(Max_geometries,4), SURFACE_MSST_F(Max_geometries,4)
   real(ffp)  :: LC_MSST_F(MaxJacobs, Max_geometries,4), LC_SURFACE_MSST_F(MaxJacobs, Max_geometries,4)
   real(ffp)  :: SAVED_MSST   ( Maxsegments, Max_geometries, 4)
   real(ffp)  :: SURFACE_MSST ( Max_geometries,4)
   real(ffp)  :: SAVED_LC_MSST   ( MaxJacobs, Maxsegments, Max_geometries, 4 )
   real(ffp)  :: LC_SURFACE_MSST ( MaxJacobs, Max_geometries, 4 )

!  sources output, recurrence relation variables

   real(ffp)  :: cumsource_ms(4), cumsource_db(4), sum1(4)
   real(ffp)  :: cumsource_ms_saved(maxsegments,4), cumsource_db_saved(0:maxsegments,4)
!   real(ffp)  :: temposs ( maxsegments, 4 )

   real(ffp)  :: L_cumsource_ms (maxJacobs,4), L_cumsource_db (maxJacobs,4), L_sum1(4)
   real(ffp)  :: sources   ( maxsegments, 4 ), L_sources ( maxsegments, 4, maxJacobs )

!  fourier index, accuracy

   Integer    :: FOURIER, Segment_converge_Count
   Real(ffp)  :: ACCURACY, DJ

!  convergence issues

   INTEGER    :: FOURIER_SAVED  ( MAXSEGMENTS, Max_geometries )
   LOGICAL    :: GEOM_ITERATION ( MAXSEGMENTS, Max_geometries )
   INTEGER    :: GEOM_TESTCONV  ( MAXSEGMENTS, Max_geometries )
   LOGICAL    :: DO_MGCONVERGE  ( Max_geometries, 0:MAXFOURIER )
   LOGICAL    :: SAVED_DO_MGCONVERGE  ( MAXSEGMENTS, Max_geometries, 0:MAXFOURIER )
   LOGICAL    :: do_Overall_Converge (MAXSEGMENTS)
   LOGICAL    :: MAIN_Converge

   INTEGER    :: SFOURIER_SAVED  ( Max_geometries )
   LOGICAL    :: SGEOM_ITERATION ( Max_geometries )
   INTEGER    :: SGEOM_TESTCONV  ( Max_geometries )
   LOGICAL    :: DO_MGSCONVERGE  ( Max_geometries, 0:MAXFOURIER )
   LOGICAL    :: SURF_Converge

   LOGICAL    :: GRAND_Converge

!  Local internal Type structure output

   type(Homog_GSF)       :: HGSF
   type(Homog_Solutions) :: HSOL
   type(Homog_BVPSetup)  :: HBVP

   type(Partic_GSF)       :: PGSF
   type(Partic_Solutions) :: PSOL
   type(Partic_BVPSolve)  :: PBVP

   type(L_Homog_Solutions)  :: L_HSOL
   type(L_Partic_Solutions) :: L_PSOL
   type(L_Partic_BVPSolve)  :: L_PBVP

!  Local Proxies or derived integers

   integer    :: ngeoms, nchoice, nflevel, nbeams, nstokes_SQ, nmoments, nsa, n_fourier_components, n_fourier_scomponents
   integer    :: nstreams_2, nstks_nstrms, nstks_nstrms_2, Ntotal, N_Subdiag, N_Supdiag

!  Local help variables

   CHARACTER(Len=2) :: CM
   CHARACTER(Len=3) :: CSEG
   integer          :: s, n, o1, o2, L, gk, j, igeom, M, M1, q, qv

!  Local Error handling

   LOGICAL             :: DO_INCLUDE_SURFACE
   CHARACTER (LEN=120) :: MESSAGE, TRACE_1, TRACE_2

!  1/24/22. Local for the Mk3 Calls
!  ================================

!  Bookkeeping and deltam-output from HOMOG Miscsetups
!  ---------------------------------------------------

!  deltam-scaling optical output

   Real(ffp) :: OMEGA_TOTAL ( MAXLAYERS )
   Real(ffp) :: DELTAU_VERT ( MAXLAYERS )
   Real(ffp) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
   Real(ffp) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
   Real(ffp) :: TRUNC_FACTOR ( MAXLAYERS ), FAC1 ( MAXLAYERS )
   Real(ffp) :: TAUGRID ( 0:MAXLAYERS ), TAUGRID_INPUT ( 0:MAXLAYERS )
   Real(ffp) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )

!  Quadrature and bookkeeping

   Real(ffp) ::   QUAD_STREAMS ( MAXSTREAMS )
   Real(ffp) ::   QUAD_WEIGHTS ( MAXSTREAMS )
   Real(ffp) ::   QUAD_STRMWTS ( MAXSTREAMS )
   Real(ffp) ::   QUAD_HALFWTS ( MAXSTREAMS )
   Integer   ::   MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
   Real(ffp) ::   DMAT ( MAXSTOKES, MAXSTOKES ), DFLUX ( MAXSTOKES)
   Real(ffp) ::   FLUXVEC ( MAXSTOKES )
   Logical   ::   DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
   Logical   ::   DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )

!  Particular-integral setups
!  --------------------------

!  general

   Real(ffp) :: COS_SZANGLES ( MAXBEAMS )
   Real(ffp) :: USER_STREAMS ( MAX_USER_STREAMS )
   LOGICAL   :: DO_DIRECT_BEAM, DO_NO_AZIMUTH

!  deltam scaling and Chapman factor output

   Real(ffp) :: DELTAU_SLANT    ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
   Real(ffp) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )

!  Beam parameterization output

   INTEGER   :: BEAM_CUTOFF        ( MAXBEAMS )
   Real(ffp) :: TRANS_SOLAR_BEAM   ( MAXBEAMS )
   LOGICAL   :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
   Real(ffp) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
   Real(ffp) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
   Real(ffp) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )

!  From the Linearized Deltam-scaling setup

   Real(ffp)  :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
   Real(ffp)  :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )
   Real(ffp)  :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
   Real(ffp)  :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS,    MAXLAYERS, MAXBEAMS )
   Real(ffp)  :: L_TRUNC_FACTOR   ( MAX_ATMOSWFS, MAXLAYERS )
   LOGICAL    :: DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
   Real(ffp)  :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )

!  From the Linearized SSALB_INIT setup

   Real(ffp)  :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )

!  From the Linearized QSPREP (average secant) setup.
!    -- 9/21/21. LC output from LAC Miscsetups

   Real(ffp)  :: LC_AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
   Real(ffp)  :: LC_INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
   Real(ffp)  :: LC_T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )

!  Saved PIMM_11 and PIMM_KM arrays

   Real(ffp) :: SAVED_PIMM_11 ( maxsegments, MAX_ALLSTRMS_P1 )
   Real(ffp) :: SAVED_PIMM_KM ( maxsegments, MAX_ALLSTRMS_P1 )

!  Zero the output and exception handling
!  --------------------------------------

!    -- 7/21/21. Mk2. Zero the fluxes

   Stokes_MS             = zero
   LC_Stokes_MS          = zero
   Output_Regular_Fluxes = zero
   Output_Actinic_Fluxes = zero
   LC_Output_Regular_Fluxes = zero
   LC_Output_Actinic_Fluxes = zero

!  4/24/24.  Initialize Diffuse Output fluxes_SB arrays

   OutputDiffuse_Regular_Fluxes_SB    = zero
   OutputDiffuse_Actinic_Fluxes_SB    = zero
   LC_OutputDiffuse_Regular_Fluxes_SB = zero
   LC_OutputDiffuse_Actinic_Fluxes_SB = zero

   fail      = .false.
   status    = 0
   nmessages = 0
   messages  = ' '

!  Check inputs if desired
!  -----------------------

!mick mod 8/3/2021 - add this section

   if ( do_debug_MS_input ) then
     CALL FOMS_WRITE_MS_LC_INPUTS ( do_debug_fluxes, &
       do_deltam_scaling, do_lambertian_surface, do_sleave, do_sl_isotropic, & ! Input flags
       nlayers, nstreams, nstokes, NGREEK_MOMENTS_INPUT, Regime_idx,         & ! Input control
       eradius, heights, nsegments, seglayer_map, segupdown_map,             & ! Input geophysical, segments
       nfinedivs, alpha, theta, relaz, alphafine, thetafine, relazfine,      & ! Input geometry
       lostrans, sfi_factors, sfi_opdeps,                                    & ! Input source func integ
       solarflux, deltaus, omega, greekmat, reflec,                          & ! Input flux, atmos & surf optical
       do_Column_Wfs, n_Column_Wfs, npert, eps,                              & ! Lin Control
       L_lostrans, L_sfi_factors, L_sfi_opdeps,                              & ! Inputs (SS Linearized)
       L_deltaus, L_omega, L_greekmat )                                        ! Inputs (Optical Linearized)
   endif

!  Step 1. Fix VQS inputs
!  ==========================

!  Initialize Saved PIMM arrays

   SAVED_PIMM_11(1:nsegments,:) = 0.0d0
   SAVED_PIMM_KM(1:nsegments,:) = 0.0d0

!   Accuracy, solar flux, grids

   ACCURACY         = 0.00001d0
   DO_DOUBLE_CONVTEST = .not. DO_RAYLEIGH_ONLY
   FLUX_FACTOR      = SOLARFLUX
   EARTH_RADIUS     = eradius
   HEIGHT_GRID(0:nlayers) = heights(0:nlayers)

!  Optical properties. Copy Greekmat up to 2*NSTREAMS

   GREEKMAT_TOTAL_INPUT = zero ; nstreams_2 = 2 * nstreams
   DO N = 1, NLAYERS
     DELTAU_VERT_INPUT(N) = deltaus(n)
     OMEGA_TOTAL_INPUT(N) = MIN(omega(n),0.9999999d0)
     DO L = 0, NSTREAMS_2
       DO O1 = 1, 4
         DO O2 = 1, 4 
           gk = 4*(o1-1)+O2
           GREEKMAT_TOTAL_INPUT(L,N,gk) = greekmat(N,L,O1,O2)
         ENDDO
       ENDDO
     ENDDO
   ENDDO

!  Linearized optical properties, double normalized

   L_GREEKMAT_TOTAL_INPUT = zero ; L_DELTAU_VERT_INPUT = zero ;L_OMEGA_TOTAL_INPUT    = zero
   if ( do_Column_WFs ) then
     DO N = 1, NLAYERS
       qv = lvarynumbers(n)
       if ( lvaryflags(n) ) then
         L_DELTAU_VERT_INPUT(1:qv,N) = L_deltaus(n,1:qv)/deltaus(n)
         L_OMEGA_TOTAL_INPUT(1:qv,N) = L_omega(n,1:qv)/omega(n)
         DO L = 0, NSTREAMS_2
           DO O1 = 1, 4
             DO O2 = 1, 4 
               gk = 4*(o1-1)+O2 ; if ( greekmat(N,L,O1,O2) .ne. zero ) &
               L_GREEKMAT_TOTAL_INPUT(1:qv,L,N,gk) = &
                         L_greekmat(N,L,O1,O2,1:qv) / greekmat(N,L,O1,O2)
             ENDDO
           ENDDO
         ENDDO
       ENDIF
     ENDDO
   ENDIF

!   surface reflectance. BRDF currently disabled

   IF ( DO_LAMBERTIAN_SURFACE ) then
      ALBEDO    = Reflec(1,1)
      brdf_f = zero ;  brdf_f_0 = zero ;  user_brdf_f = zero
   ELSE
      ALBEDO    = zero
      Stop ' Cannot have BRDFs for MS calculation, Placeholder'
!        PLaceholder - need a BRDF supplement to specify BRDF_F inputs
   ENDIF

!   Surface leaving. Not enabled at present

   IF ( .not. DO_SLEAVE ) then
      slterm_ISOTROPIC = zero ; slterm_f_0  = zero 
   ELSE
      Stop ' Cannot have SLEAVE for MS calculation, Placeholder'
!        PLaceholder - need a SLEAVE supplement to specify SLTERM inputs
   ENDIF

!  7/21/21. Mk2. Set additional MSSEG flags
!   -- Flux output (REGULAR_FLUX, ACTINIC_FLUX) controlled by DO_NCHOICE_FLUXES input flag
!   -- ASYMTX_TOLERANCE flag fixed by hand
!   -- DO_FOURIER0_NSTOKES2 performance flag, can be turned on for Lambertian surfaces

   ASYMTX_TOLERANCE     = 1.0d-20
   DO_FOURIER0_NSTOKES2 = do_Lambertian_Surface

!  Call the Homogeneous setups

   CALL MSSEG_HOMOG_MISCSETUPS ( & 
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING,                             & ! Input
        NSTOKES, NLAYERS, NSTREAMS, NGREEK_MOMENTS_INPUT, HEIGHT_GRID,   & ! Input
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,      & ! Input
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,              & ! Output
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTOKES_SQ,                        & ! Output
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,          & ! Output
        FLUXVEC, MUELLER_INDEX, DMAT, DFLUX, DO_REAL_EIGENSOLVER,        & ! Output
        DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,   & ! Output
        TAUGRID_INPUT, TAUGRID, TRUNC_FACTOR, FAC1, OMEGA_GREEK,         & ! Output
        STATUS, NMESSAGES, MESSAGES, ACTIONS )                             ! Output

!  Exception handling

    if ( STATUS.eq.VQS_SERIOUS ) RETURN
    
!  Now the linearized setups

    IF ( DO_COLUMN_WFS ) THEN
      CALL MSSEG_LA_MISCSETUPS ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS,        & ! Input
        MUELLER_INDEX, LVARYFLAGS, LVARYNUMBERS,              & ! Input
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, OMEGA_GREEK, & ! Input Optical/Deltam
        TRUNC_FACTOR, FAC1, L_DELTAU_VERT_INPUT,              & ! Input Optical/Deltam
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,          & ! Input Linearized Optical
        L_OMEGA_TOTAL, L_DELTAU_VERT, L_GREEKMAT_TOTAL,       & ! Output
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK )   ! Output
   ENDIF

!  number of Fourier components (general case)

   IF ( DO_RAYLEIGH_ONLY  ) THEN
      N_FOURIER_COMPONENTS = 2
   ELSE
      N_FOURIER_COMPONENTS = NMOMENTS
   ENDIF
   MAIN_CONVERGE = .true.

!  number of Fourier components (surface case)

   if ( Regime_idx .eq.1 ) then
      SURF_CONVERGE = .True. 
      if ( do_Lambertian_Surface ) then
         N_FOURIER_SCOMPONENTS = 0
      else
         IF ( DO_RAYLEIGH_ONLY  ) THEN
            N_FOURIER_SCOMPONENTS = 2
         ELSE
            N_FOURIER_SCOMPONENTS = NMOMENTS
         ENDIF
      ENDIF
   ELSE
      SURF_CONVERGE = .false.
   endif

!  set up Convergence flags (Surface case). Only require 1 geometry

   if ( Regime_idx .eq.1 ) then 
     SGEOM_TESTCONV  ( 1 )  = 0
     SGEOM_ITERATION ( 1 ) = .TRUE.
     DO L = 0, MAXFOURIER
       DO_MGSCONVERGE ( 1, L )                = .TRUE.
       DO_MGSCONVERGE ( 2:Max_geometries, L ) = .FALSE.
     ENDDO
   ENDIF

!  Have to set ngeoms here

   ngeoms = MAXVAL(nfinedivs(1:nsegments))

!  set up Convergence flags (general multi-segment case)

   GEOM_TESTCONV  ( 1:nsegments, 1:ngeoms )  = 0
   GEOM_ITERATION ( 1:nsegments, 1:ngeoms ) = .TRUE.
   DO L = 0, MAXFOURIER
      SAVED_DO_MGCONVERGE ( 1:nsegments, 1:ngeoms, L ) = .TRUE.
   ENDDO
   DO_OVERALL_CONVERGE(1:nsegments) = .true.

!  Start FOURIER Loop
!  ==================

   GRAND_CONVERGE = .true. ; FOURIER = -1
   DO WHILE ( GRAND_CONVERGE .and. FOURIER.lt. N_FOURIER_COMPONENTS )
      FOURIER = FOURIER + 1 ; M = FOURIER ; M1 = M - 1

!  10/15/23. Introduce progress control

      if ( RT_PROGRESS ) then
        write(*,'(A,I1,A,I3)')' -- Doing VQS Mk3 LC-Linearized MS calculation, Regime ',&
                                Regime_idx,', Homogeneous RT Solution, Fourier # ', FOURIER
      endif

!  Perform the Homogeneous solution calculation, with linearization if flagged

      IF ( do_Column_WFs ) THEN
        CALL MSSEG_LA_HOMOG_FOURIER ( &
           DO_LAMBERTIAN_SURFACE, DO_FOURIER0_NSTOKES2, ASYMTX_TOLERANCE,           & ! Input Control
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NMOMENTS, LVARYFLAGS, LVARYNUMBERS, & ! Input Control
           NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,  & ! Input Bookkeeping
           MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING,           & ! Input Bookkeeping
           QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS, ALBEDO, BRDF_F,                & ! Input Quads/Surface
           DELTAU_VERT, OMEGA_GREEK, L_DELTAU_VERT, L_OMEGA_GREEK,                  & ! Input Optical        
           DO_INCLUDE_SURFACE, HGSF, HSOL, HBVP, L_HSOL,                            & ! Output TS
           STATUS, MESSAGE, TRACE_1, TRACE_2 )                                        ! Exception handling
      ELSE
        CALL MSSEG_HOMOG_FOURIER ( &
           DO_LAMBERTIAN_SURFACE, DO_FOURIER0_NSTOKES2, ASYMTX_TOLERANCE, & ! Input
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NMOMENTS, NSTREAMS_2,     & ! Input
           NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,    & ! Input
           MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, & ! Input
           QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,                      & ! Input
           ALBEDO, BRDF_F, DELTAU_VERT, OMEGA_GREEK,                      & ! Input        
           DO_INCLUDE_SURFACE, HGSF, HSOL, HBVP,                          & ! Output
           STATUS, MESSAGE, TRACE_1, TRACE_2 )
      ENDIF

!  Exception handling

      if ( STATUS.eq.VQS_SERIOUS ) THEN
         write(CM,'(I2.2)')FOURIER
         NMESSAGES = NMESSAGES + 1
         if ( DO_Column_WFs ) then
            MESSAGES(NMESSAGES) = &
                  'Failure from MSSEG_LA_HOMOG_FOURIER, fourier # = '//CM//', Here are the errors --'
         else
           MESSAGES(NMESSAGES) = &
                  'Failure from MSSEG_HOMOG_FOURIER, fourier # = '//CM//', Here are the errors --'
         endif
         RETURN
      ENDIF

!  10/15/23. Introduce progress control

      if ( RT_PROGRESS ) then
        write(*,'(A,I1,A,I3)')' -- Doing VQS Mk3 LC-Linearized MS calculation, Regime ',&
                                Regime_idx,', Atmosphere contribution, Fourier # ', FOURIER
      endif

!  Start the Segment loop
!  ======================

!  NOTE: this segment loop will not be done if NSEGMENTS = 0
!        -- occurs for example is camera is at the ground, upwelling


!  initialize the segment convergence count

      Segment_converge_count = 0

!  start loop over segments

      do s = 1, nsegments

!  Only perform if segment has not yet converged

        if ( Do_overall_converge(s) ) then
          Segment_converge_count = Segment_converge_count + 1

!  VQS Setups, Use fine layer angles

          DO_SURFACE_MSST = .false.        
          DO_UPWELLING    = segupdown_map(s)         
          nchoice         = seglayer_map(s)
          ngeoms          = nfinedivs(s) ; nsa = 2 * ngeoms
          N_SZANGLES      = ngeoms ; SZANGLES     (1:ngeoms) = thetafine(s,1:ngeoms)/dtr
          N_USER_VZANGLES = ngeoms ; USER_VZANGLES(1:ngeoms) = alphafine(s,1:ngeoms)/dtr
          N_USER_RELAZMS  = ngeoms ; USER_RELAZMS (1:ngeoms) = relazfine(s,1:ngeoms)/dtr

!  7/21/21. Mk2. set the local flux control

          do_NChoice_fluxes = ( Do_Calculate_Fluxes .and.s.eq.nsegments )
          if ( do_debug_fluxes ) do_NChoice_fluxes = Do_Calculate_Fluxes

!  now the partic setups (has to be repeated)
!    10/4/23. Must use segment-saved SFI_OPDEPS values

          CALL MSSEG_PARTIC_MISCSETUPS ( &
            DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NSTOKES, NLAYERS,         & ! Input
            EARTH_RADIUS, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES,      & ! Input
            SZANGLES, USER_RELAZMS, USER_VZANGLES, HEIGHT_GRID,             & ! Input
            DELTAU_VERT_INPUT, DELTAU_VERT, TAUGRID, SFI_OPDEPS(s,:), FAC1, & ! Input
            DO_NO_AZIMUTH, DO_DIRECT_BEAM, COS_SZANGLES, USER_STREAMS,      & ! output
            NBEAMS, CHAPMAN_FACTORS, DELTAU_SLANT, SFI_OPDEPS_SCALED,       & ! output
            BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,         & ! output
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! output
            STATUS, NMESSAGES, MESSAGES, ACTIONS )                            ! output

!  Exception handling

          if ( STATUS.eq.VQS_SERIOUS ) RETURN

!  Now the linearized Partic setups
!    10/4/23. Must use segment-saved L_SFI_OPDEPS values

          IF ( DO_COLUMN_WFS ) THEN
            CALL MSSEG_LC_MISCSETUPS ( &
              DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NLAYERS, NBEAMS, & ! Input
              NCHOICE, LVARYFLAGS, LVARYNUMBERS, OMEGA_TOTAL_INPUT,  & ! Input
              DELTAU_VERT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,         & ! Input
              BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR,             & ! Input
              L_SFI_OPDEPS(s,:,:), L_TRUNC_FACTOR, L_DELTAU_VERT,    & ! Input
              L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,              & ! Input
              L_DELTAU_SLANT, L_SFI_OPDEPS_SCALED,                   & ! Output
              LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR )   ! Output
          ENDIF

!  Set the Local PIMM arrays to be used as input

          if ( FOURIER.gt.2 ) then
            PGSF%PIMM_11 (M1,1:nsa ) = SAVED_PIMM_11 (s,1:NSA)
            PGSF%PIMM_KM (M1,1:nsa ) = SAVED_PIMM_KM (s,1:NSA)
          endif

!  Set the local MGCONVERGE flag

          do igeom = 1, ngeoms
            DO_MGCONVERGE(igeom,0:MAXFOURIER) = SAVED_DO_MGCONVERGE(s,igeom,0:MAXFOURIER)
          enddo

!  call the Regular or Linearized Partic Fourier
!   -- 9/21/23. Add BRDF_F, DELTAU_SLANT and LC_T_DELT_MUBAR (for linearized BVP)

          IF ( DO_COLUMN_WFS ) THEN
            CALL MSSEG_LC_FOURIER ( &
              DO_SURFACE_MSST, DO_UPWELLING, DO_PLANE_PARALLEL, DO_DIRECT_BEAM,     &
              DO_MGCONVERGE, DO_NCHOICE_FLUXES, DO_FOURIER0_NSTOKES2,               &
              DO_INCLUDE_SURFACE,DO_LAMBERTIAN_SURFACE, DO_SLEAVE, DO_SL_ISOTROPIC, &
              DO_COLUMN_WFS, N_COLUMN_WFS, LVARYFLAGS, LVARYNUMBERS,                & ! Linearization control
              FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,                &
              NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,            &
              N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,                  &
              DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,             &
              USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,              &
              DELTAU_SLANT, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK,           & ! Input Optical        
              ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0, & ! Input surface 
              BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS, AVERAGE_SECANT,        & ! Input Beam param
              T_DELT_MUBAR, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,  & ! Input Beam param
              L_DELTAU_VERT, L_SFI_OPDEPS_SCALED, L_OMEGA_GREEK,                   & ! Input Lin. Optical
              HGSF, HSOL, HBVP, L_HSOL,                                            & ! Input TS
              MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,              & ! regular outputs and TS
              LC_MSST_F, LC_SURFACE_MSST_F, LOCAL_LC_FLUXES, L_PSOL, L_PBVP,       & ! Linearized outputs and TS
              STATUS, MESSAGE, TRACE_1, TRACE_2 )                                    ! Exception handling
          ELSE
            CALL MSSEG_PARTIC_FOURIER ( &
              DO_SURFACE_MSST, DO_UPWELLING, DO_DIRECT_BEAM, DO_MGCONVERGE,   &
              DO_NCHOICE_FLUXES, DO_FOURIER0_NSTOKES2, DO_INCLUDE_SURFACE,    &
              DO_LAMBERTIAN_SURFACE, DO_SLEAVE, DO_SL_ISOTROPIC,              &
              FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,           &
              NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,       &
              N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,             &
              DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,        &
              USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,         &
              ALBEDO, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK,            &           
              BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,            & 
              BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS,                   &
              AVERAGE_SECANT, T_DELT_MUBAR, HGSF, HSOL, HBVP,                 &
              MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,         &
              STATUS, MESSAGE, TRACE_1, TRACE_2 )
          ENDIF

!  10/04/23. Checks out now.
!if ( FOURIER.eq.0)write(701,*)FOURIER,s,MSST_F(1:6,1)
!if ( FOURIER.eq.0)write(700,*)FOURIER,s,LC_MSST_F(1,1:6,1)

!  Exception handling

          if ( STATUS.eq.VQS_SERIOUS ) THEN
            write(CM,'(I2.2)')FOURIER
            write(CSEG,'(I3.3)') S
            if ( DO_Column_WFs ) then
               MESSAGES(NMESSAGES) = &
               'Failure from MSSEG_LC_FOURIER, fourier # = '//CM//', in Segement # '//CSEG//'; Here are the errors --'
            else
              MESSAGES(NMESSAGES) = &
               'Failure from MSSEG_PARTIC_FOURIER, fourier # = '//CM//', in Segement # '//CSEG//'; Here are the errors --'
            endif
            NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(MESSAGE)
            NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(TRACE_1)
            NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(TRACE_2)
            RETURN
          ENDIF

!  Set the Local PIMM arrays to be used as input for the next Fourier term

          if ( FOURIER.ge.2 ) then
            SAVED_PIMM_11 (s,1:NSA) = PGSF%PIMM_11 (M,1:nsa)
            SAVED_PIMM_KM (s,1:NSA) = PGSF%PIMM_KM (M,1:nsa)
          endif

!!  7/21/21. Mk2. Output the Fluxes for the last segment only (averaged over fine-division outputs)
!   -- Output Top-of-layer fluxes (Regimes 1 and 2), Bottom-of-layer fluxes (Regime 3)

          if ( FOURIER.eq.0.and.DO_NCHOICE_FLUXES .and. s.eq.nsegments ) then
            J = nfinedivs(s) ; DJ = dble(J)
            nflevel = 1 ; if ( Regime_idx .eq. 3 ) nflevel = 2
            do o1 = 1, nstokes 
              OUTPUT_REGULAR_FLUXES(o1) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,2))/DJ
              OUTPUT_ACTINIC_FLUXES(o1) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,1))/DJ
            enddo
            if ( do_Column_WFs) then
              do q = 1, n_column_wfs
                do o1 = 1, nstokes 
                  LC_Output_REGULAR_FLUXES(o1,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,2))/ DJ
                  LC_Output_ACTINIC_FLUXES(o1,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,1))/ DJ
                enddo
              enddo
            endif
          endif

!  4/24/24. Diffuse atmospheric fluxes at segment boundaries
!  ---------------------------------------------------------

!   -- replaces the debug write-up to fort.81
!   -- Have to be careful with the segment boundary selection (regime dependent)

          if ( FOURIER.eq.0.and.DO_NCHOICE_FLUXES .and. do_debug_fluxes ) then
             J = nfinedivs(s) ; DJ = dble(J)

!  Starting value at fisrt segment needs special treatment 

             if ( s.eq.1 ) then
               nflevel = 1 ; if ( Regime_idx .eq. 1 ) nflevel = 2
               do o1 = 1, nstokes 
                 OUTPUTDIFFUSE_REGULAR_FLUXES_SB(o1,0) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,2))/DJ
                 OUTPUTDIFFUSE_ACTINIC_FLUXES_SB(o1,0) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,1))/DJ
               enddo
               if ( do_Column_WFs) then
                 do q = 1, n_column_wfs
                   do o1 = 1, nstokes 
                     LC_OUTPUTDIFFUSE_REGULAR_FLUXES_SB(o1,0,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,2))/ DJ
                     LC_OUTPUTDIFFUSE_ACTINIC_FLUXES_SB(o1,0,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,1))/ DJ
                   enddo
                 enddo
               endif
             endif

!  All other segments. Take care of regime dependence
!    -- for Regime 2, level = 1 on near side of tangent height, level 2 for far side

             nflevel = 2 ; if ( Regime_idx .eq. 1 ) nflevel = 1
             if ( Regime_idx .eq. 2 .and. segupdown_map(s) ) nflevel = 1
             do o1 = 1, nstokes 
                OUTPUTDIFFUSE_REGULAR_FLUXES_SB(o1,s) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,2))/DJ
                OUTPUTDIFFUSE_ACTINIC_FLUXES_SB(o1,s) = SUM(LOCAL_FLUXES(1:J,o1,nflevel,1))/DJ
             enddo
             if ( do_Column_WFs) then
               do q = 1, n_column_wfs
                 do o1 = 1, nstokes 
                   LC_OUTPUTDIFFUSE_REGULAR_FLUXES_SB(o1,s,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,2))/ DJ
                   LC_OUTPUTDIFFUSE_ACTINIC_FLUXES_SB(o1,s,q) = SUM(LOCAL_LC_FLUXES(q,1:J,o1,nflevel,1))/ DJ
                 enddo
               enddo
             endif

!  End flux output setting

          endif

!  7/21/21. Debug atmospheric fluxes. results to fort.81. Rob debug
!          if ( FOURIER.eq.0.and.DO_NCHOICE_FLUXES .and. do_debug_fluxes ) then
!            J = nfinedivs(s) ; DJ = dble(J)
!            if ( s.eq.1)write(81,'(/A/)')' segm. layer    ActFlux @TOP   ActFlux @BOT    RegFlux @TOP   RegFlux @BOT'
!            write(81,455)s,nchoice,SUM(LOCAL_FLUXES(1:J,1,1,1))/DJ,SUM(LOCAL_FLUXES(1:J,1,2,1))/DJ,&
!                                   SUM(LOCAL_FLUXES(1:J,1,1,1))/DJ,SUM(LOCAL_FLUXES(1:J,1,2,1))/DJ
!          endif 
!455       format(i4,3x,i4,2(1x,1p2e15.5))

!  Examine Convergence

          CALL FOMS_MS_LC_Convergence_Mk3 &
           ( Max_Geometries, MaxFourier, Max_Atmoswfs, do_Column_wfs,                  &
             do_Rayleigh_only, do_double_convtest, Fourier, nstokes, ngeoms, nmoments, &
             n_column_wfs, deg_to_rad, accuracy, USER_RELAZMS, msst_f, LC_msst_f,      &
             do_MGConverge, Geom_Testconv(s,:), Geom_Iteration(s,:), do_Overall_Converge(s), &
             Fourier_Saved(s,:), saved_msst(s,:,:), saved_LC_msst(:,s,:,:) )

!  Reset the global MGCONVERGE flag

           do igeom = 1, ngeoms
             SAVED_DO_MGCONVERGE(s,igeom,0:MAXFOURIER)= DO_MGCONVERGE(igeom,0:MAXFOURIER) 
           enddo

!  End segment loop

         endif
      enddo

!  Check Segment Converge count

      if ( Segment_converge_count .eq. 0 ) MAIN_CONVERGE = .false.

!  DONT FORGET THE SURFACE MS contribution ( Regime 1 )
!  ====================================================

      if ( regime_idx .eq. 1 .and. SURF_CONVERGE ) then

!  MSTerm Setups, Use single BOA geometry only, turn on surface MSST flag, no fluxes
!   10/10/23. Proper settings for geometrical variables; set SFI_OPDEPS dummies

        DO_SURFACE_MSST = .true.        
        DO_UPWELLING    = .true.        
        nchoice         = nlayers
        ngeoms          = 1
        N_SZANGLES      = ngeoms ; SZANGLES     (1) = theta(0)/dtr
        N_USER_VZANGLES = ngeoms ; USER_VZANGLES(1) = alpha(0)/dtr
        N_USER_RELAZMS  = ngeoms ; USER_RELAZMS (1) = relaz(0)/dtr
        DO_NCHOICE_FLUXES = .false.
        sfi_opdeps_dum   = 0.0_ffp
        L_sfi_opdeps_dum = 0.0_ffp

!  10/15/23. Introduce progress control

        if ( RT_PROGRESS ) then
          write(*,'(A,I1,A,I3)')' -- Doing VQS Mk3 LC-Linearized MS calculation, Regime ',&
                                Regime_idx,',   Surface  contribution, Fourier # ', FOURIER
        endif

!  now the partic setups. Has to be repeated
!    10/4/23. Must use segment-saved SFI_OPDEPS values (actually, dummy)

        CALL MSSEG_PARTIC_MISCSETUPS ( &
          DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NSTOKES, NLAYERS,         & ! Input
          EARTH_RADIUS, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES,      & ! Input
          SZANGLES, USER_RELAZMS, USER_VZANGLES, HEIGHT_GRID,             & ! Input
          DELTAU_VERT_INPUT, DELTAU_VERT, TAUGRID, sfi_opdeps_dum, FAC1, & ! Input
          DO_NO_AZIMUTH, DO_DIRECT_BEAM, COS_SZANGLES, USER_STREAMS,      & ! output
          NBEAMS, CHAPMAN_FACTORS, DELTAU_SLANT, sfi_opdeps_dum_SCALED,   & ! output
          BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,         & ! output
          INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                    & ! output
          STATUS, NMESSAGES, MESSAGES, ACTIONS )                            ! output

!  Exception handling

        if ( STATUS.eq.VQS_SERIOUS ) RETURN

!  Now the linearized Partic setups
!    10/4/23. Must use segment-saved L_SFI_OPDEPS values (actually, dummy)

        IF ( DO_COLUMN_WFS ) THEN
          CALL MSSEG_LC_MISCSETUPS ( &
            DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NLAYERS, NBEAMS, & ! Input
            NCHOICE, LVARYFLAGS, LVARYNUMBERS, OMEGA_TOTAL_INPUT,  & ! Input
            DELTAU_VERT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,         & ! Input
            BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR,             & ! Input
            L_sfi_opdeps_dum, L_TRUNC_FACTOR, L_DELTAU_VERT,       & ! Input
            L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,              & ! Input
            L_DELTAU_SLANT, L_SFI_OPDEPS_Dum_SCALED,               & ! Output
            LC_AVERAGE_SECANT, LC_INITIAL_TRANS, LC_T_DELT_MUBAR )   ! Output
        ENDIF

!  call the Regular or Linearized Partic Fourier
!   -- 9/21/23. Add BRDF_F, DELTAU_SLANT and LC_T_DELT_MUBAR (for linearized BVP)
!   -- 10/10/23. Use variable DO_MGS_CONVERGE

        IF ( DO_COLUMN_WFS ) THEN
          CALL MSSEG_LC_FOURIER ( &
            DO_SURFACE_MSST, DO_UPWELLING, DO_PLANE_PARALLEL, DO_DIRECT_BEAM,     &
            DO_MGSCONVERGE, DO_NCHOICE_FLUXES, DO_FOURIER0_NSTOKES2,              &
            DO_INCLUDE_SURFACE,DO_LAMBERTIAN_SURFACE, DO_SLEAVE, DO_SL_ISOTROPIC, &
            DO_COLUMN_WFS, N_COLUMN_WFS, LVARYFLAGS, LVARYNUMBERS,                & ! Linearization control
            FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,                &
            NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,            &
            N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,                  &
            DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,             &
            USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,              &
            DELTAU_SLANT, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK,           & ! Input Optical        
            ALBEDO, BRDF_F, BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0, & ! Input surface 
            BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS, AVERAGE_SECANT,        & ! Input Beam param
            T_DELT_MUBAR, LC_INITIAL_TRANS, LC_AVERAGE_SECANT, LC_T_DELT_MUBAR,  & ! Input Beam param
            L_DELTAU_VERT, L_SFI_OPDEPS_SCALED, L_OMEGA_GREEK,                   & ! Input Lin. Optical
            HGSF, HSOL, HBVP, L_HSOL,                                            & ! Input TS
            MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,              & ! regular outputs and TS
            LC_MSST_F, LC_SURFACE_MSST_F, LOCAL_LC_FLUXES, L_PSOL, L_PBVP,       & ! Linearized outputs and TS
            STATUS, MESSAGE, TRACE_1, TRACE_2 )                                    ! Exception handling
        ELSE
          CALL MSSEG_PARTIC_FOURIER ( &
            DO_SURFACE_MSST, DO_UPWELLING, DO_DIRECT_BEAM, DO_MGSCONVERGE,  &
            DO_NCHOICE_FLUXES, DO_FOURIER0_NSTOKES2, DO_INCLUDE_SURFACE,    &
            DO_LAMBERTIAN_SURFACE, DO_SLEAVE, DO_SL_ISOTROPIC,              &
            FOURIER, NSTOKES, NSTREAMS, NLAYERS, NGEOMS, NCHOICE,           &
            NMOMENTS, NSTREAMS_2, NTOTAL, MUELLER_INDEX, DMAT, DFLUX,       &
            N_SUBDIAG, N_SUPDIAG, NSTKS_NSTRMS, NSTKS_NSTRMS_2,             &
            DO_LAYER_SCATTERING, FLUX_FACTOR, FLUXVEC, COS_SZANGLES,        &
            USER_STREAMS, QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,         &
            ALBEDO, DELTAU_VERT, SFI_OPDEPS_SCALED, OMEGA_GREEK,            &           
            BRDF_F_0, USER_BRDF_F, SLTERM_ISOTROPIC, SLTERM_F_0,            & 
            BEAM_CUTOFF, TRANS_SOLAR_BEAM, INITIAL_TRANS,                   &
            AVERAGE_SECANT, T_DELT_MUBAR, HGSF, HSOL, HBVP,                 &
            MSST_F, SURFACE_MSST_F, LOCAL_FLUXES, PGSF, PSOL, PBVP,         &
            STATUS, MESSAGE, TRACE_1, TRACE_2 )
        ENDIF

!  Checks out, 10/10/23
!write(*,*)FOURIER, DO_SURFACE_MSST, SURFACE_MSST_F(1,1),LC_SURFACE_MSST_F(1:2,1,1)

!  Exception handling

        if ( STATUS.eq.VQS_SERIOUS ) THEN
          write(CM,'(I2.2)')FOURIER
          if ( DO_Column_WFs ) then
            MESSAGES(NMESSAGES) = &
              'Failure from MSSEG_LC_FOURIER, fourier # = '//CM//', for Surface regime 1 term; Here are the errors --'
          else
            MESSAGES(NMESSAGES) = &
              'Failure from MSSEG_PARTIC_FOURIER, fourier # = '//CM//', for Surface regime 1 term; Here are the errors --'
          endif
          NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(MESSAGE)
          NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(TRACE_1)
          NMESSAGES = NMESSAGES + 1 ; MESSAGES(NMESSAGES) = '  ** '//Trim(TRACE_2)
          RETURN
        ENDIF

!  8/26/24. Output the Fluxes for the ground segment only 
!   -- Output Bottom of layer (surface) fluxes (Regime 1)

        if ( FOURIER.eq.0.and.DO_NCHOICE_FLUXES ) then
          nflevel = 2
          do o1 = 1, nstokes 
            OUTPUT_REGULAR_FLUXES(o1) = LOCAL_FLUXES(1,o1,nflevel,2)
            OUTPUT_ACTINIC_FLUXES(o1) = LOCAL_FLUXES(1,o1,nflevel,1)
          enddo
          if ( do_Column_WFs) then
            do q = 1, n_column_wfs
              do o1 = 1, nstokes 
                LC_Output_REGULAR_FLUXES(o1,q) = LOCAL_LC_FLUXES(q,1,o1,nflevel,2)
                LC_Output_ACTINIC_FLUXES(o1,q) = LOCAL_LC_FLUXES(q,1,o1,nflevel,1)
              enddo
            enddo
          endif
        endif

!  Convergence (not needed for Lambertian_surface)

        if ( do_Lambertian_Surface .and. FOURIER.eq.0 ) then
          SURFACE_msst(1,1:nstokes) = SURFACE_msst_F(1,1:nstokes)
          SURF_CONVERGE = .false. ; SFourier_Saved(1) = 0
          IF ( DO_COLUMN_WFS ) THEN
            DO Q = 1, N_COLUMN_WFS
              LC_SURFACE_msst(q,1,1:nstokes) = LC_SURFACE_msst_F(q,1,1:nstokes)
            ENDDO
          ENDIF
        else
          CALL FOMS_MS_LC_Convergence_Mk3 &
           ( Max_Geometries, MaxFourier, Max_Atmoswfs, do_Column_wfs,                  &
             do_Rayleigh_only, do_double_convtest, Fourier, nstokes, ngeoms, nmoments, &
             n_column_wfs, deg_to_rad, accuracy, USER_RELAZMS, SURFACE_MSST_F, LC_SURFACE_MSST_F, &
             do_MGSConverge, SGeom_Testconv, SGeom_Iteration, Surf_Converge, SFourier_Saved,      &
             SURFACE_msst, LC_SURFACE_msst )
        endif

!  End Regime 1, surface term

      endif

!  Grand convergence

      GRAND_CONVERGE = ( MAIN_CONVERGE .or. SURF_CONVERGE ) 

!  End Fourier Loop

   enddo

!  Checks out, 1010/23
!write(*,*)'Converged Regime 1',LC_SURFACE_MSST(1:2,1,1)

!  Source function path integration for Fourier components
!  =======================================================

!  Fourier-summed Source Functions for all segments

   do s = 1, nsegments
     sum1(1:nstokes) = zero ; j = nfinedivs(s)
     do o1 = 1, nstokes 
       sum1(o1) = DOT_PRODUCT(sfi_factors(s,1:j),SAVED_MSST(s,1:j,o1))
     enddo
     sources(s,1:nstokes) = sum1(1:nstokes)
     if ( do_Column_WFs) then
       do q = 1, n_column_wfs
         do o1 = 1, nstokes 
           L_sum1(o1) = DOT_PRODUCT(sfi_factors(s,1:j),SAVED_LC_MSST(q,s,1:j,o1))
           L_sum1(o1) = L_sum1(o1) + DOT_PRODUCT(L_sfi_factors(s,1:j,q),SAVED_MSST(s,1:j,o1))
         enddo
         L_sources(s,1:nstokes,q) = L_sum1(1:nstokes)
       enddo
     endif
   enddo

!  Rob debug
!    fort.53. Validates sourceterm output and Fourier convergence
!   do s = 1, nsegments
!       write(53,*)S,sources(s,1), Fourier_Saved(s,1:6)
!   enddo

!  RT recursion. Have to save the cumulative sources at all stages

   CUMSOURCE_MS = zero
   do s = 1, nsegments
      CUMSOURCE_MS(1:nstokes) = LOSTRANS(S) * CUMSOURCE_MS(1:nstokes) + sources(S,1:nstokes)
      CUMSOURCE_MS_SAVED(s,:) = CUMSOURCE_MS(:)
    enddo
   STOKES_MS(1:nstokes) = SOLARFLUX * CUMSOURCE_MS(1:nstokes)

!  Add surface term for Regime 1
!   10/10/23.  Need to properly specify starting point for CUMSOURCE_DB_SAVED

   CUMSOURCE_DB = zero ; CUMSOURCE_DB_SAVED = zero
   if ( Regime_idx .eq. 1 ) then
      CUMSOURCE_DB(1:nstokes)         = SURFACE_msst(1,1:nstokes)
      CUMSOURCE_DB_SAVED(0,1:nstokes) = SURFACE_msst(1,1:nstokes)
      do s = 1, nsegments
         CUMSOURCE_DB(1:nstokes) = LOSTRANS(S) * CUMSOURCE_DB(1:nstokes)
         CUMSOURCE_DB_SAVED(s,:) = CUMSOURCE_DB(:)
      enddo
      STOKES_MS(1:nstokes) = STOKES_MS(1:nstokes) + SOLARFLUX * CUMSOURCE_DB(1:nstokes)
   endif

!  Column Wfs (atmospheric term). All Regimes
!mick fix 6/6/2023 - added "s.eq.1" porton to IF block to handle that case

   if ( do_Column_wfs ) then
     L_CUMSOURCE_MS = zero
     do s = 1, nsegments
       do q = 1, n_column_wfs
         if ( s.eq.1 ) then
           L_CUMSOURCE_MS(q,1:nstokes) = L_SOURCES(s,1:nstokes,q) + LOSTRANS(s) * L_CUMSOURCE_MS(q,1:nstokes)
         else
           L_CUMSOURCE_MS(q,1:nstokes) = L_SOURCES(s,1:nstokes,q) + LOSTRANS(s) * L_CUMSOURCE_MS(q,1:nstokes) &
                                                                    + L_LOSTRANS(s,q) * CUMSOURCE_MS_SAVED(s-1,1:nstokes)
         endif
       enddo
     enddo
     do q = 1, n_column_wfs
       LC_STOKES_MS(1:nstokes,q) = solarflux * L_CUMSOURCE_MS(q,1:nstokes)
     enddo
   endif

!  Add Column Wfs (Surface direct beam term, Regime 1 only)
!  BRDF capability added for FOMS.
!mick fix 7/13/2015 - added maxfine dimension to LC_SURFACE_MSST
!mick fix 6/6/2023  - added "s.eq.1" porton to IF block to handle that case

   if ( Regime_idx .eq. 1 .and. do_column_wfs ) then
     do q = 1, n_column_wfs
       L_CUMSOURCE_DB(q,1:nstokes) = LC_SURFACE_MSST(q,1,1:nstokes)
     enddo
     do s = 1, nsegments
       do q = 1, n_column_wfs
           L_CUMSOURCE_DB(q,1:nstokes) =   LOSTRANS(s)   * L_CUMSOURCE_DB(q,1:nstokes) &
                                       + L_LOSTRANS(s,q) *   CUMSOURCE_DB_SAVED(s-1,1:nstokes)
       enddo
     enddo
     do q = 1, n_column_wfs
       LC_STOKES_MS(1:nstokes,q) = LC_STOKES_MS(1:nstokes,q) + solarflux * L_CUMSOURCE_DB(q,1:nstokes)
     enddo
   endif

!  Update STOKES_TOT

   STOKES_TOT(1:nstokes) = STOKES_TOT(1:nstokes) + STOKES_MS(1:nstokes)
   if ( do_column_wfs ) then
     do q = 1, n_column_wfs
       LC_STOKES_TOT(1:nstokes,q) = LC_STOKES_TOT(1:nstokes,q) + LC_STOKES_MS(1:nstokes,q)
     enddo
   endif

!  Finish

   return
end subroutine FOMS_MS_LC_VectorRT_Mk3_Exp

!  Private convergence routine

subroutine FOMS_MS_LC_Convergence_Mk3 &
     ( Maxgeoms, MaxFourier, MaxJacs, do_Colwfs, do_Rayleigh_only, do_double_convtest, Fourier, &
       nstokes, ngeoms, nmoments, n_colwfs, deg_to_rad, accuracy, azimuths, msst_f, LC_msst_f,  &
       do_MGConverge, Geom_Testconv, Geom_Iteration, do_Overall_Converge, Fourier_Saved, msst, LC_msst )

   use FOMS_pars_m, only : ffp

   implicit none
!   integer, parameter :: ffp = selected_real_kind(15)

!  Inputs
!  ------

!  dimensions

   INTEGER, intent(in) :: Maxgeoms, MaxFourier, MaxJacs

!  Flags

   LOGICAL, intent(in) :: do_double_convtest, do_Colwfs, do_Rayleigh_only

!  fourier number

   INTEGER, intent(in) :: Fourier

!  Numbers

   INTEGER, intent(in) :: nstokes, ngeoms, nmoments, n_colwfs

!  constants, includeing convergence accuracy criterion

   real(ffp), intent(in) :: deg_to_rad, accuracy

!  Azimuth angles

   real(ffp), intent(in) :: azimuths(Maxgeoms)

!  source term fourier components

   real(ffp), intent(in) :: msst_f    (Maxgeoms,4)
   real(ffp), intent(in) :: LC_msst_f (MaxJacs,Maxgeoms,4)

!  Modified outputs
!  ----------------

!  Flags and bookkeeping

   LOGICAL, intent(inout) :: do_MGConverge  ( Maxgeoms, 0:MaxFourier )
   INTEGER, intent(inout) :: Geom_Testconv  ( Maxgeoms )
   LOGICAL, intent(inout) :: Geom_Iteration ( Maxgeoms )
   LOGICAL, intent(inout) :: do_Overall_Converge
   INTEGER, intent(inout) :: Fourier_Saved  ( Maxgeoms )

!  Fourier-integrated sources

   real(ffp), intent(inout) :: msst    (Maxgeoms,4)
   real(ffp), intent(inout) :: LC_msst (MaxJacs,Maxgeoms,4)
!  Local
!  -----

   integer   :: IGEOM, IGEOM_COUNT, COUNT, O, L, Q
   real(ffp) :: AZMFAC (Maxgeoms,4), Azm_argument(MaxGeoms)
   real(ffp) :: ACCUR, DFC, TOLD(4), TNEW(4), TAZM(4)

!  azimuth cosine/sine factors, using adjust geometries.
!      MISTAKE with the factor TWO. Corrected 17 November 2015.

      IF ( Fourier .GT. 0 ) THEN
        DFC = DBLE(Fourier)
        DO IGEOM = 1, NGEOMS
          Azm_argument(IGEOM) = azimuths(IGEOM) * DFC
          AZMFAC(IGEOM,1)     = 2.0d0*COS(DEG_TO_RAD*Azm_argument(IGEOM))
        ENDDO
        if ( NSTOKES.ge.3 ) then
          DO IGEOM = 1, NGEOMS
            AZMFAC(IGEOM,2) = AZMFAC(IGEOM,1)
            AZMFAC(IGEOM,3) = 2.0d0 * SIN(DEG_TO_RAD*Azm_argument(IGEOM))
            AZMFAC(IGEOM,4) = AZMFAC(IGEOM,3)
          ENDDO
        endif
      endif

!  Fourier summation and Convergence examination
!  ---------------------------------------------

!   -- only done for Solutions which are still not converged
!      This is controlled by flag DO_MULTIGEOM

      IGEOM_COUNT = 0
      DO IGEOM = 1, NGEOMS
         IF ( do_MGConverge ( IGEOM, FOURIER ) ) THEN

!  Fourier 0 component

            IF ( FOURIER.EQ.0 ) THEN
              MSST(IGEOM,1:nstokes) = MSST_F(IGEOM,1:nstokes)
              if ( do_ColWFs ) then
                DO Q = 1, n_ColWFs
                 LC_MSST(Q,IGEOM,1:nstokes) = LC_MSST_F(Q,IGEOM,1:nstokes)
                ENDDO
              endif
            ELSE

!  For Rayleigh atmosphere, skip convergence test on intensity

              IF ( do_Rayleigh_Only ) THEN
                DO O = 1, NSTOKES
                  MSST(IGEOM,O) = MSST(IGEOM,O) + AZMFAC(IGEOM,O) * MSST_F(IGEOM,O)
                ENDDO
                if ( do_ColWFs ) then
                  DO Q = 1, n_ColWFs ; DO O = 1, NSTOKES
                    LC_MSST(Q,IGEOM,O) = LC_MSST(Q,IGEOM,O) + AZMFAC(IGEOM,O) * LC_MSST_F(Q,IGEOM,O)
                  ENDDO ; ENDDO
                endif
              ENDIF

!  For the general case, Fourier > 0

              IF ( .not. do_Rayleigh_Only ) then
                COUNT = 0
                DO O = 1, NSTOKES
                  TOLD(O) = MSST(IGEOM,O)
                  TAZM(O) = AZMFAC(IGEOM,O) * MSST_F(IGEOM,O)
                  TNEW(O) = TOLD(O) + TAZM(O)
                ENDDO
                IF ( TAZM(1) .NE. 0.0d0 ) THEN
                  ACCUR     = ABS(TAZM(1)/TNEW(1))
                  IF ( ACCUR .LT. ACCURACY ) COUNT = COUNT + 1
                ELSE
                  COUNT   = COUNT + 1
                ENDIF
                MSST(IGEOM,1:nstokes)   = TNEW(1:nstokes)
                if ( do_ColWFs ) then
                  DO Q = 1, n_ColWFs ; DO O = 1, NSTOKES
                    LC_MSST(Q,IGEOM,O) = LC_MSST(Q,IGEOM,O) + AZMFAC(IGEOM,O) * LC_MSST_F(Q,IGEOM,O)
                  ENDDO ; ENDDO
                endif

!  set convergence counter TESTCONV

                IF ( COUNT .EQ. 1 ) THEN
                  Geom_Testconv(IGEOM) = Geom_Testconv(IGEOM) + 1
                  IF ( Do_Double_Convtest ) THEN
                    IF ( Geom_Testconv(IGEOM).EQ. 2 ) Geom_Iteration(IGEOM) = .FALSE.
                  ELSE
                    Geom_Iteration(IGEOM) = .FALSE.
                  ENDIF
                  IF ( .NOT. Geom_Iteration(IGEOM) ) THEN
                    Fourier_Saved(IGEOM) = Fourier
                  ENDIF
                ELSE
                  Geom_Testconv(IGEOM) = 0
                  Fourier_Saved(IGEOM) = NMoments
                ENDIF

!  End convergence test IF block

              ENDIF

!  For Rayleigh scattering alone, stop iteration after Fourier = 2

              IF ( do_Rayleigh_Only ) THEN
                IF ( Fourier .EQ. 2 ) THEN
                  Geom_Iteration(IGEOM) = .FALSE.
                  Fourier_Saved(IGEOM) = 2
                ENDIF
              ENDIF

!  Finish Fourier testing

            ENDIF
       
!  Check number of beams already converged

            IF ( Geom_Iteration(IGEOM) ) THEN
              IGEOM_COUNT = IGEOM_COUNT + 1
            ELSE
              DO L = Fourier + 1, MaxFourier
                do_MGConverge (IGEOM,L) = .FALSE.
              ENDDO
            ENDIF

!  End geometry count loop

         ENDIF
      END DO

!  If all beams have converged, stop iteration

      IF ( IGEOM_COUNT .EQ. 0 ) Do_Overall_Converge = .FALSE.

!  finish

      RETURN
END subroutine FOMS_MS_LC_Convergence_Mk3

! PRIVATE ROUTINE

SUBROUTINE FOMS_WRITE_MS_LC_INPUTS ( do_debug_fluxes, &
   do_deltam_scaling, do_lambertian_surface, do_sleave,do_sl_isotropic, & ! Input flags
   nlayers, nstreams, nstokes, nmoments_input, Regime_idx,              & ! Input control
   eradius, heights, nsegments, seglayer_map, segupdown_map,            & ! Input geophysical, segments
   nfinedivs, alpha, theta, relaz, alphafine, thetafine, relazfine,     & ! Input geometry
   lostrans, sfi_factors, sfi_opdeps,                                   & ! Input source func integ
   flux, deltaus, omega, greekmat, reflec,                              & ! Input flux, atmos & surf optical
   do_Column_Wfs, n_Column_Wfs, npert, eps,                             & ! Lin Control
   L_lostrans, L_sfi_factors, L_sfi_opdeps,                             & ! Inputs (SS Linearized)
   L_deltaus, L_omega, L_greekmat )                                       ! Inputs (Optical Linearized)

!  Check inputs if desired
!mick mod 8/3/2021 - added this subroutine, should be public.
!Rob mod 8/18/21. Add do_debug_fluxes

   use FOMS_pars_m, only : ffp, maxlays, maxsegments, maxfine, maxmoms, maxJacobs
   use MSSeg_pars_m, only: DWFL, DWFI, DWFL1, DWFI1, DWFR, DWFR1, DWFR2, DWFR3, DWFR5, DWFR6

   implicit none

!  Inputs
!  ======

!  Flags

   !logical  , intent(in)  :: do_sunlight
   logical  , intent(in)  :: do_debug_fluxes
   logical  , intent(in)  :: do_deltam_scaling
   logical  , intent(in)  :: do_lambertian_surface
   logical  , intent(in)  :: do_sleave
   logical  , intent(in)  :: do_sl_isotropic

!  Number of layers, streams, regime index

   integer  , intent(in)  :: nlayers, nstreams, Regime_idx

!  Number of Stokes components, moments

   integer  , intent(in)  :: nstokes
   integer  , intent(in)  :: nmoments_input

!  Geophysical

   real(ffp), intent(in)  :: eradius, heights (0:maxlays)

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   integer  , intent(in)  :: nsegments
   integer  , intent(in)  :: seglayer_map  (maxsegments)
   logical  , intent(in)  :: segupdown_map (maxsegments)
   !logical  , intent(in)  :: layermask_map (maxlays,2)

!  Number of fine divisions

   integer  , intent(in)  :: nfinedivs (maxsegments)

!  Angles - segment values

   real(ffp), intent(in)  :: alpha (0:maxsegments)
   real(ffp), intent(in)  :: theta (0:maxsegments)
   real(ffp), intent(in)  :: relaz (0:maxsegments)

!  Angles - quadrature values

   real(ffp), intent(in)  :: alphafine (maxsegments,maxfine)
   real(ffp), intent(in)  :: thetafine (maxsegments,maxfine)
   real(ffp), intent(in)  :: relazfine (maxsegments,maxfine)

!  Source function integration stuff

   real(ffp), intent(in)  :: lostrans    ( maxsegments )
   real(ffp), intent(in)  :: sfi_factors ( maxsegments, maxfine )
   real(ffp), intent(in)  :: sfi_opdeps  ( maxsegments, maxfine )

!  Atmosphere optical inputs

   real(ffp), intent(in) :: DELTAUS     ( maxlays )
   real(ffp), intent(in) :: OMEGA       ( maxlays )
!   real(ffp), intent(in) :: TRUNCFAC    ( maxlays )
   real(ffp), intent(in) :: GREEKMAT    ( maxlays,0:maxmoms, 4, 4 )

!  Solar flux and surface reflectivity (Regime 1, Lambertian-only)

   real(ffp), intent(in) :: REFLEC(4,4), FLUX

!  Linearized
!  ----------

!  Jacobian control

   logical, intent(in) ::  Do_Column_WFs
   integer, intent(in) ::  n_Column_WFs

!  Rob add 5/30/18, Perturbation layer and eps.(for debug only)

   integer  , intent(in) :: npert
   real(ffp), intent(in) :: eps

   real(ffp), intent(in)  :: L_lostrans    ( maxsegments, maxJacobs )
   real(ffp), intent(in)  :: L_sfi_factors ( maxsegments, maxfine, maxJacobs )
   real(ffp), intent(in)  :: L_sfi_opdeps  ( maxsegments, maxfine, maxJacobs )

!  Linearized optical inputs

   real(ffp), intent(in) :: L_deltaus     ( maxlays, maxJacobs )
   real(ffp), intent(in) :: L_omega       ( maxlays, maxJacobs )
   real(ffp), intent(in) :: L_greekmat    ( maxlays,0:maxmoms, 4, 4, maxJacobs )

!  Local variables

   integer :: outunit, n_wfs
   integer :: lay,mom,nf,ngreek,o1,o2,s,ss,wf
   integer :: gmask2(8,2)

!  Open output file

   OUTUNIT = 101
   OPEN (OUTUNIT,file = 'Debug/FOMS_MS_LC_Write_Input.dbg',status = 'replace')

!  Only output non-zero entries of GreekMat

   IF ( nstokes.eq.1 ) THEN
     ngreek = 1
   ELSE IF ( nstokes.eq.3 ) THEN
     ngreek = 5
   ELSE IF ( nstokes.eq.4 ) THEN
     ngreek = 8
   END IF

   !Non-zero GreekMat components (a 2-D 4x4 array):
   !             O1                 O2
   gmask2(1,1) = 1 ; gmask2(1,2) = 1
   gmask2(2,1) = 1 ; gmask2(2,2) = 2
   gmask2(3,1) = 2 ; gmask2(3,2) = 1
   gmask2(4,1) = 2 ; gmask2(4,2) = 2
   gmask2(5,1) = 3 ; gmask2(5,2) = 3
   gmask2(6,1) = 3 ; gmask2(6,2) = 4
   gmask2(7,1) = 4 ; gmask2(7,2) = 3
   gmask2(8,1) = 4 ; gmask2(8,2) = 4

   n_wfs = 0
   if (Do_Column_WFs) n_wfs = n_Column_Wfs

!  Write all VFO input to file

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,'(A)') '---------------'
   WRITE(OUTUNIT,'(A)') 'Standard Inputs'
   WRITE(OUTUNIT,'(A)') '---------------'

!  Control

   WRITE(OUTUNIT,*)
   !WRITE(OUTUNIT,DWFL)  'do_sunlight           = ',do_sunlight
   WRITE(OUTUNIT,DWFL)  'do_debug_fluxes       = ',do_debug_fluxes
   WRITE(OUTUNIT,DWFL)  'do_deltam_scaling     = ',do_deltam_scaling
   WRITE(OUTUNIT,DWFL)  'do_lambertian_surface = ',do_lambertian_surface
   WRITE(OUTUNIT,DWFL)  'do_sleave             = ',do_sleave
   WRITE(OUTUNIT,DWFL)  'do_sl_isotropic       = ',do_sl_isotropic

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFI)  'nlayers               = ',nlayers
   WRITE(OUTUNIT,DWFI)  'nstreams              = ',nstreams
   WRITE(OUTUNIT,DWFI)  'nstokes               = ',nstokes
   WRITE(OUTUNIT,DWFI)  'nmoments_input        = ',nmoments_input

!  Regime index

   WRITE(OUTUNIT,DWFI)  'Regime_idx            = ',Regime_idx

!  Geophysical

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFR)  'eradius = ',eradius

   WRITE(OUTUNIT,*)
   DO lay=0,nlayers
     WRITE(OUTUNIT,DWFR1)  'lay = ',lay,' heights(lay) = ',heights(lay)
   END DO

!  Segment to layer map, segment up/down switch, layer mask (1=up, 2=dn)

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFI)  'nsegments = ',nsegments

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFI1)  's = ',s,' seglayer_map(s) = ',seglayer_map(s)
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFL1)  's = ',s,' segupdown_map(s) = ',segupdown_map(s)
   END DO

   !WRITE(OUTUNIT,*)
   !DO dir=1,2
   !  DO n=1,nlayers
   !    WRITE(OUTUNIT,DWFL2)  'dir = ',dir,' n = ',n,' layermask_map(n,dir) = ',layermask_map(n,dir)
   !  END DO
   !END DO

!  Number of fine divisions

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFI1)  's = ',s,' nfinedivs(s) = ',nfinedivs(s)
   END DO

!  Angles - segment values (in degrees)

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFR1)  's = ',s,' alpha(s) = ',alpha(s)
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFR1)  's = ',s,' theta(s) = ',theta(s)
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFR1)  's = ',s,' relaz(s) = ',relaz(s)
   END DO

!  Angles - quadrature values (in degrees)

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     DO nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2)  's = ',s,' nf = ',nf,' alphafine(s,nf) = ',alphafine(s,nf)
     END DO
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     DO nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2)  's = ',s,' nf = ',nf,' thetafine(s,nf) = ',thetafine(s,nf)
     END DO
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     DO nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2)  's = ',s,' nf = ',nf,' relazfine(s,nf) = ',relazfine(s,nf)
     END DO
   END DO

!  Source function integration stuff

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     WRITE(OUTUNIT,DWFR1)  's = ',s,' lostrans(s) = ',lostrans(s)
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     DO nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2)  's = ',s,' nf = ',nf,' sfi_factors(s,nf) = ',sfi_factors(s,nf)
     END DO
   END DO

   WRITE(OUTUNIT,*)
   DO s=1,nsegments
     DO nf=1,nfinedivs(s)
       WRITE(OUTUNIT,DWFR2)  's = ',s,' nf = ',nf,' sfi_opdeps(s,nf) = ',sfi_opdeps(s,nf)
     END DO
   END DO

!  Atmosphere

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFR)  'flux = ',flux

   WRITE(OUTUNIT,*)
   DO lay=1,nlayers
     WRITE(OUTUNIT,DWFR1)  'lay = ',lay,' deltaus(lay) = ',deltaus(lay)
   END DO

   WRITE(OUTUNIT,*)
   DO lay=1,nlayers
     WRITE(OUTUNIT,DWFR1)  'lay = ',lay,' omega(lay) = ',omega(lay)
   END DO

   DO ss=1,ngreek
     o1 = gmask2(ss,1) ; o2 = gmask2(ss,2)
     WRITE(OUTUNIT,*)
     DO lay=1,nlayers
       DO mom=0,nmoments_input
         WRITE(OUTUNIT,DWFR5)  'ss = ',ss,' o1 = ',o1,' o2 = ',o2,' lay = ',lay,' mom = ',mom,&
           ' greekmat(lay,mom,o1,o2) = ',greekmat(lay,mom,o1,o2)
       END DO
     END DO
   END DO

!  Surface

   WRITE(OUTUNIT,*)
   DO ss=1,ngreek
     o1 = gmask2(ss,1) ; o2 = gmask2(ss,2)
     WRITE(OUTUNIT,DWFR3)  'ss = ',ss,' o1 = ',o1,' o2 = ',o2,&
       ' reflec(o1,o2) = ',reflec(o1,o2)
   END DO


!  Linearized Inputs
!  -----------------

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,'(A)') '-----------------'
   WRITE(OUTUNIT,'(A)') 'Linearized Inputs'
   WRITE(OUTUNIT,'(A)') '-----------------'

   WRITE(OUTUNIT,*)
   WRITE(OUTUNIT,DWFL) 'Do_Column_WFs = ',Do_Column_WFs

   if (Do_Column_WFs) then

!  Jacobian control

     WRITE(OUTUNIT,*)
     WRITE(OUTUNIT,DWFI) 'N_Column_Wfs = ',N_Column_Wfs

     WRITE(OUTUNIT,*)
     WRITE(OUTUNIT,DWFI) 'npert = ',npert
     WRITE(OUTUNIT,DWFR) 'eps   = ',eps

!  Source function integration stuff

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do s=1,nsegments
         WRITE(OUTUNIT,DWFR2)  'wf = ',wf,' s = ',s,' L_lostrans(s,wf) = ',L_lostrans(s,wf)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do s=1,nsegments
         do nf=1,nfinedivs(s)
           WRITE(OUTUNIT,DWFR3)  'wf = ',wf,' s = ',s,' nf = ',nf,' L_sfi_factors(s,nf,wf) = ',L_sfi_factors(s,nf,wf)
         enddo
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do s=1,nsegments
         do nf=1,nfinedivs(s)
           WRITE(OUTUNIT,DWFR3)  'wf = ',wf,' s = ',s,' nf = ',nf,' L_sfi_opdeps(s,nf,wf) = ',L_sfi_opdeps(s,nf,wf)
         enddo
       enddo
     enddo

!  Optical inputs

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFR2) 'wf = ',wf,' lay = ',lay,&
                             ' L_deltaus(lay,wf) = ',L_deltaus(lay,wf)
       enddo
     enddo

     WRITE(OUTUNIT,*)
     do wf=1,n_wfs
       do lay=1,nlayers
         WRITE(OUTUNIT,DWFR2) 'wf = ',wf,' lay = ',lay,&
                             ' L_omega(lay,wf)   = ',L_omega(lay,wf)
       enddo
     enddo

     do wf=1,n_wfs
       do ss=1,ngreek
         o1 = gmask2(ss,1); o2 = gmask2(ss,2); WRITE(OUTUNIT,*)
         do lay=1,nlayers
           do mom=0,nmoments_input
             WRITE(OUTUNIT,DWFR6) 'wf = ',wf,' ss = ',ss,' o1 = ',o1,' o2 = ',o2,' lay = ',lay,' mom = ',mom,&
               ' L_greekmat(lay,mom,o1,o2,wf) = ',L_greekmat(lay,mom,o1,o2,wf)
           enddo
         enddo
       enddo
     enddo
   endif

!  Close output file

   CLOSE (OUTUNIT)

   END SUBROUTINE FOMS_WRITE_MS_LC_INPUTS

! End module

end module FOMS_MS_LC_VectorRT_Mk3_Exp_m
