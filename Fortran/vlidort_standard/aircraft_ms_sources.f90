module aircraft_ms_sources_m
! Sparse AC uses the supplied QS homogeneous/particular MS kernels.
! Kernel call sequence adapted from R. Spurr's FOMS_MS_LP_VectorRT_Mk3_Exp;
! the supplied modules retain their original GPL headers.
! One homogeneous factorization per Fourier order serves every sample on both legs.
! Existing diffuse solutions are integrated on the actual curved chord, with
! analytic optical derivatives. No additional MS solve or QS master is called.
use aircraft_chord_quadrature_m
use MSSeg_PARS_m
use MSSeg_HPTypes_Def_m
use MSSeg_L_HPTypes_Def_m
use MSSEG_HOMOG_MISCSETUPS_m
use MSSEG_HOMOG_Fourier_m
use MSSEG_PARTIC_MISCSETUPS_m
use MSSEG_PARTIC_Solutions_m
use MSSEG_PARTIC_BVProblem_m
use MSSEG_PARTIC_PostProcess_m
use MSSEG_LA_MISCSETUPS_m
use MSSEG_LA_SOLUTIONS_m
use MSSEG_LA_HOMOG_Fourier_m
use MSSEG_LP_MISCSETUPS_m
use MSSEG_LP_SOLUTIONS_m
use MSSEG_LP_BVProblem_m
use MSSEG_LP_PostProcess_m
use VLIDORT_Inputs_def_m, only: VLIDORT_Fixed_Inputs, VLIDORT_Modified_Inputs
use VLIDORT_LinInputs_def_m, only: VLIDORT_Fixed_LinInputs
implicit none
private
integer,parameter::ffp=fpk
public :: aircraft_ms_sources
contains
subroutine aircraft_ms_sources(fix,modin,linfix,linearized,ngeoms,angles,up,needed, &
                              layer_ms,layer_jac,surface_ms,surface_jac, &
                              impact,node_mu,first_node,node_count,exit_t,chord_length)
type(VLIDORT_Fixed_Inputs),intent(in)::fix
type(VLIDORT_Modified_Inputs),intent(in)::modin
type(VLIDORT_Fixed_LinInputs),intent(in)::linfix
logical,intent(in)::linearized,up(:),needed(:,:)
integer,intent(in)::ngeoms
real(fpk),intent(in)::angles(:,:),impact,node_mu(:),exit_t(:,:),chord_length(:,:)
integer,intent(in)::first_node(:),node_count(:)
type(chord_rule),allocatable::chord(:,:)
real(fpk),intent(out)::layer_ms(:,:,:),layer_jac(:,:,:,:,:),surface_ms(:,:),surface_jac(:,:,:,:)
Real(ffp) :: OMEGA_TOTAL ( MAXLAYERS )
Real(ffp) :: DELTAU_VERT ( MAXLAYERS )
Real(ffp) :: GREEKMAT_TOTAL ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
Real(ffp) :: OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES )
Real(ffp) :: TRUNC_FACTOR ( MAXLAYERS ), FAC1 ( MAXLAYERS )
Real(ffp) :: TAUGRID ( 0:MAXLAYERS ), TAUGRID_INPUT ( 0:MAXLAYERS )
Real(ffp) :: SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES )
Real(ffp) ::   QUAD_STREAMS ( MAXSTREAMS )
Real(ffp) ::   QUAD_WEIGHTS ( MAXSTREAMS )
Real(ffp) ::   QUAD_STRMWTS ( MAXSTREAMS )
Real(ffp) ::   QUAD_HALFWTS ( MAXSTREAMS )
Integer   ::   MUELLER_INDEX ( MAXSTOKES, MAXSTOKES )
Real(ffp) ::   DMAT ( MAXSTOKES, MAXSTOKES ), DFLUX ( MAXSTOKES)
Real(ffp) ::   FLUXVEC ( MAXSTOKES )
Logical   ::   DO_REAL_EIGENSOLVER ( 0:MAXMOMENTS, MAXLAYERS )
Logical   ::   DO_LAYER_SCATTERING ( 0:MAXMOMENTS, MAXLAYERS )
Real(ffp) :: COS_SZANGLES ( MAXBEAMS )
Real(ffp) :: USER_STREAMS ( MAX_USER_STREAMS )
LOGICAL   :: DO_DIRECT_BEAM, DO_NO_AZIMUTH
Real(ffp) :: DELTAU_SLANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
Real(ffp) :: CHAPMAN_FACTORS ( MAXLAYERS, MAXLAYERS, MAXBEAMS )
INTEGER   :: BEAM_CUTOFF        ( MAXBEAMS )
Real(ffp) :: TRANS_SOLAR_BEAM   ( MAXBEAMS )
LOGICAL   :: DO_REFLECTED_DIRECTBEAM ( MAXBEAMS )
Real(ffp) :: INITIAL_TRANS  ( MAXLAYERS, MAXBEAMS )
Real(ffp) :: AVERAGE_SECANT ( MAXLAYERS, MAXBEAMS )
Real(ffp) :: T_DELT_MUBAR   ( MAXLAYERS, MAXBEAMS )
Real(ffp)  :: L_OMEGA_TOTAL    ( MAX_ATMOSWFS, MAXLAYERS )
Real(ffp)  :: L_DELTAU_VERT    ( MAX_ATMOSWFS, MAXLAYERS )
Real(ffp)  :: L_GREEKMAT_TOTAL ( MAX_ATMOSWFS, 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES_SQ )
Real(ffp)  :: L_DELTAU_SLANT   ( MAX_ATMOSWFS, MAXLAYERS,    MAXLAYERS, MAXBEAMS )
Real(ffp)  :: L_TRUNC_FACTOR   ( MAX_ATMOSWFS, MAXLAYERS )
LOGICAL    :: DO_SCATMAT_VARIATION ( MAXLAYERS, MAX_ATMOSWFS )
Real(ffp)  :: L_SFI_OPDEPS_SCALED  ( MAX_GEOMETRIES, MAX_ATMOSWFS )
Real(ffp)  :: L_OMEGA_GREEK ( 0:MAXMOMENTS, MAXLAYERS, MAXSTOKES, MAXSTOKES, MAX_ATMOSWFS )
Real(ffp)  :: LP_AVERAGE_SECANT ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
Real(ffp)  :: LP_INITIAL_TRANS  ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
Real(ffp)  :: LP_T_DELT_MUBAR   ( MAXLAYERS, MAXLAYERS, MAXBEAMS, MAX_ATMOSWFS )
type(Homog_GSF)::HGSF
type(Homog_Solutions)::HSOL
type(Homog_BVPSetup)::HBVP
type(Partic_GSF)::PGSF
type(Partic_Solutions)::PSOL
type(Partic_BVPSolve)::PBVP
type(L_Homog_Solutions)::L_HSOL
type(L_Partic_Solutions)::L_PSOL
type(L_Partic_BVPSolve)::L_PBVP
real(fpk)::HEIGHT_GRID(0:MAXLAYERS), EARTH_RADIUS,FLUX_FACTOR,ALBEDO,ASYMTX_TOLERANCE
real(fpk)::DELTAU_VERT_INPUT(MAXLAYERS),OMEGA_TOTAL_INPUT(MAXLAYERS)
real(fpk)::GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,MAXLAYERS,16)
real(fpk)::L_DELTAU_VERT_INPUT(MAX_ATMOSWFS,MAXLAYERS),L_OMEGA_TOTAL_INPUT(MAX_ATMOSWFS,MAXLAYERS)
real(fpk)::L_GREEKMAT_TOTAL_INPUT(MAX_ATMOSWFS,0:MAXMOMENTS,MAXLAYERS,16)
real(fpk)::BRDF_F(16,MAXSTREAMS,MAXSTREAMS),BRDF_F_0(16,MAXSTREAMS,MAX_GEOMETRIES)
real(fpk)::USER_BRDF_F(16,MAX_GEOMETRIES,MAXSTREAMS)
real(fpk)::SLTERM_ISOTROPIC(4,MAX_GEOMETRIES),SLTERM_F_0(4,MAXSTREAMS,MAX_GEOMETRIES)
real(fpk)::SZANGLES(MAX_GEOMETRIES),USER_VZANGLES(MAX_GEOMETRIES),USER_RELAZMS(MAX_GEOMETRIES)
real(fpk)::DUMMY_SFI(MAX_GEOMETRIES),DUMMY_L_SFI(MAX_GEOMETRIES,MAX_ATMOSWFS)
real(fpk)::DIRECT_BEAM(MAXSTREAMS,MAXBEAMS,MAXSTOKES),ATMOS_ATTN(MAXBEAMS)
real(fpk)::SURFACE_FACTOR,DELTA_FACTOR,STOKES_DOWNSURF(MAXSTREAMS,MAXSTOKES)
real(fpk)::SURFACE_MSST_F(MAX_GEOMETRIES,MAXSTOKES)
real(fpk)::LP_SURFACE_MSST_F(MAXLAYERS,MAX_ATMOSWFS,MAX_GEOMETRIES,MAXSTOKES)
real(fpk),allocatable::uc(:,:,:,:,:),luc(:,:,:,:,:,:),nc(:,:,:,:),pc(:,:,:,:)
real(fpk)::az(4),value(4),deriv(MAX_ATMOSWFS,MAXLAYERS,4)
integer::NLAYERS,NSTREAMS,NSTOKES,NGREEK_MOMENTS_INPUT,NMOMENTS,NSTREAMS_2,NSTKS_NSTRMS,NSTKS_NSTRMS_2
integer::NTOTAL,N_SUBDIAG,N_SUPDIAG,NSTOKES_SQ,N_SZANGLES,N_USER_RELAZMS,N_USER_VZANGLES,NBEAMS
integer::FOURIER,NCHOICE,LAYER,KV,QV,Q,IGEOM,N,NPARAMS_VARY,NQ,LL_NSTK
integer::LL_NSTKNSTRM,LL_NSTKNSTRM_2,LL_NTT,LL_NSUBD,LL_NSUPD,STATUS,STATUS_SUB,NMESSAGES
integer(kind=8)::prof_mark,prof_now,prof_rate,prof_ticks(6)
logical::profile
character(8)::profile_setting
logical::LVARYFLAGS(MAXLAYERS),DO_RTSOL_VARY,DO_INCLUDE_SURFACE,DO_FOURIER0_NSTOKES2
logical::DO_RAYLEIGH_ONLY,DO_DELTAM_SCALING,DO_PLANE_PARALLEL,DO_LAMBERTIAN_SURFACE
integer::LVARYNUMBERS(MAXLAYERS)
character(120)::MESSAGE,TRACE_1,TRACE_2,MESSAGES(0:MAX_MESSAGES),ACTIONS(0:MAX_MESSAGES)

call get_environment_variable('PYDORT_AC_SHARED_PROFILE',profile_setting)
profile=trim(profile_setting)=='1';prof_ticks=0
if(profile)call system_clock(prof_mark,prof_rate)
NLAYERS=fix%Cont%TS_NLAYERS;NSTREAMS=fix%Cont%TS_NSTREAMS;NSTOKES=fix%Cont%TS_NSTOKES
if(ngeoms>MAX_GEOMETRIES)error stop 'AC shared MS: too many samples'
if(.not.fix%Bool%TS_DO_LAMBERTIAN_SURFACE)error stop 'AC shared MS requires Lambertian surface'
NQ=0
if(linearized)NQ=linfix%Cont%TS_N_TOTALPROFILE_WFS
DO_LAMBERTIAN_SURFACE=.true.;DO_FOURIER0_NSTOKES2=.true.;DO_PLANE_PARALLEL=.false.
DO_RAYLEIGH_ONLY=modin%MBool%TS_DO_RAYLEIGH_ONLY
DO_DELTAM_SCALING=modin%MBool%TS_DO_DELTAM_SCALING
NGREEK_MOMENTS_INPUT=modin%MCont%TS_NGREEK_MOMENTS_INPUT
HEIGHT_GRID=fix%Chapman%TS_HEIGHT_GRID;EARTH_RADIUS=modin%MChapman%TS_EARTH_RADIUS
FLUX_FACTOR=fix%Sunrays%TS_FLUX_FACTOR;ALBEDO=fix%Optical%TS_LAMBERTIAN_ALBEDO
ASYMTX_TOLERANCE=1.0d-28
DELTAU_VERT_INPUT=fix%Optical%TS_DELTAU_VERT_INPUT
OMEGA_TOTAL_INPUT=modin%MOptical%TS_OMEGA_TOTAL_INPUT
GREEKMAT_TOTAL_INPUT=fix%Optical%TS_GREEKMAT_TOTAL_INPUT(0:MAXMOMENTS,:,:)
L_DELTAU_VERT_INPUT=0;L_OMEGA_TOTAL_INPUT=0;L_GREEKMAT_TOTAL_INPUT=0
LVARYFLAGS=.false.;LVARYNUMBERS=0
if(linearized)then
 LVARYFLAGS=linfix%Cont%TS_LAYER_VARY_FLAG;LVARYNUMBERS=linfix%Cont%TS_LAYER_VARY_NUMBER
 L_DELTAU_VERT_INPUT=linfix%Optical%TS_L_DELTAU_VERT_INPUT
 L_OMEGA_TOTAL_INPUT=linfix%Optical%TS_L_OMEGA_TOTAL_INPUT
 L_GREEKMAT_TOTAL_INPUT=linfix%Optical%TS_L_GREEKMAT_TOTAL_INPUT(:,0:MAXMOMENTS,:,:)
endif
BRDF_F=0;BRDF_F_0=0;USER_BRDF_F=0;SLTERM_ISOTROPIC=0;SLTERM_F_0=0
DUMMY_SFI=0;DUMMY_L_SFI=0
layer_ms=0;layer_jac=0;surface_ms=0;surface_jac=0
SURFACE_MSST_F(1:ngeoms,1:NSTOKES)=0
LP_SURFACE_MSST_F(1:NLAYERS,1:NQ,1:ngeoms,1:NSTOKES)=0
CALL MSSEG_HOMOG_MISCSETUPS ( &
        DO_RAYLEIGH_ONLY, DO_DELTAM_SCALING,                             &
        NSTOKES, NLAYERS, NSTREAMS, NGREEK_MOMENTS_INPUT, HEIGHT_GRID,   &
        DELTAU_VERT_INPUT, OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT,      &
        NMOMENTS, NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2,              &
        NTOTAL, N_SUBDIAG, N_SUPDIAG, NSTOKES_SQ,                        &
        QUAD_STREAMS, QUAD_WEIGHTS, QUAD_STRMWTS, QUAD_HALFWTS,          &
        FLUXVEC, MUELLER_INDEX, DMAT, DFLUX, DO_REAL_EIGENSOLVER,        &
        DO_LAYER_SCATTERING, DELTAU_VERT, OMEGA_TOTAL, GREEKMAT_TOTAL,   &
        TAUGRID_INPUT, TAUGRID, TRUNC_FACTOR, FAC1, OMEGA_GREEK,         &
        STATUS, NMESSAGES, MESSAGES, ACTIONS )
call checked(STATUS)
if(linearized)then
CALL MSSEG_LA_MISCSETUPS ( &
        DO_DELTAM_SCALING, NSTOKES, NLAYERS, NMOMENTS,        &
        MUELLER_INDEX, LVARYFLAGS, LVARYNUMBERS,              &
        OMEGA_TOTAL_INPUT, GREEKMAT_TOTAL_INPUT, OMEGA_GREEK, &
        TRUNC_FACTOR, FAC1, L_DELTAU_VERT_INPUT,              &
        L_OMEGA_TOTAL_INPUT, L_GREEKMAT_TOTAL_INPUT,          &
        L_OMEGA_TOTAL, L_DELTAU_VERT, L_GREEKMAT_TOTAL,       &
        DO_SCATMAT_VARIATION, L_TRUNC_FACTOR, L_OMEGA_GREEK )
endif
N_SZANGLES=ngeoms;N_USER_RELAZMS=ngeoms;N_USER_VZANGLES=ngeoms
SZANGLES(1:ngeoms)=angles(:,1);USER_VZANGLES(1:ngeoms)=angles(:,2);USER_RELAZMS(1:ngeoms)=angles(:,3)
NCHOICE=1
CALL MSSEG_PARTIC_MISCSETUPS ( &
            DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NSTOKES, NLAYERS,      &
            EARTH_RADIUS, N_SZANGLES, N_USER_RELAZMS, N_USER_VZANGLES,   &
            SZANGLES, USER_RELAZMS, USER_VZANGLES, HEIGHT_GRID,          &
            DELTAU_VERT_INPUT, DELTAU_VERT, TAUGRID, DUMMY_SFI, FAC1, &
            DO_NO_AZIMUTH, DO_DIRECT_BEAM, COS_SZANGLES, USER_STREAMS,   &
            NBEAMS, CHAPMAN_FACTORS, DELTAU_SLANT, SFI_OPDEPS_SCALED,    &
            BEAM_CUTOFF, TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM,      &
            INITIAL_TRANS, AVERAGE_SECANT, T_DELT_MUBAR,                 &
            STATUS, NMESSAGES, MESSAGES, ACTIONS )
call checked(STATUS)
if(linearized)then
CALL MSSEG_LP_MISCSETUPS ( &
              DO_DELTAM_SCALING, DO_PLANE_PARALLEL, NLAYERS, NBEAMS, &
              NCHOICE, LVARYFLAGS, LVARYNUMBERS, OMEGA_TOTAL_INPUT,  &
              DELTAU_VERT, DELTAU_SLANT, TRUNC_FACTOR, FAC1,         &
              BEAM_CUTOFF, AVERAGE_SECANT, T_DELT_MUBAR,             &
              DUMMY_L_SFI, L_TRUNC_FACTOR, L_DELTAU_VERT,    &
              L_DELTAU_VERT_INPUT, L_OMEGA_TOTAL_INPUT,              &
              L_DELTAU_SLANT, L_SFI_OPDEPS_SCALED,                   &
              LP_AVERAGE_SECANT, LP_INITIAL_TRANS, LP_T_DELT_MUBAR )
endif
allocate(chord(NLAYERS,ngeoms))
do IGEOM=1,ngeoms
 do N=1,NLAYERS
  if(.not.needed(IGEOM,N))cycle
  Q=first_node(IGEOM)
  call make_chord_rule(DELTAU_VERT(N),HEIGHT_GRID(N-1),HEIGHT_GRID(N),EARTH_RADIUS, &
       impact,exit_t(IGEOM,N),chord_length(IGEOM,N), &
       node_mu(Q:Q+node_count(IGEOM)-1),IGEOM-Q+1, &
       max(ONE/minval(QUAD_STREAMS(1:NSTREAMS)),AVERAGE_SECANT(N,IGEOM)),chord(N,IGEOM))
 enddo
enddo
allocate(uc(ngeoms,NSTOKES,NSTKS_NSTRMS,NLAYERS,4))
allocate(luc(ngeoms,NSTOKES,NSTKS_NSTRMS,NQ,NLAYERS,4))
allocate(nc(NSTKS_NSTRMS,NLAYERS,NQ,NLAYERS),pc(NSTKS_NSTRMS,NLAYERS,NQ,NLAYERS))
! All diffuse Fourier orders are retained; no geometry-dependent early exit can
! omit weak terms significant for Jacobians. PGSF recursion persists across orders.
if(profile)then
 call system_clock(prof_now);prof_ticks(1)=prof_ticks(1)+prof_now-prof_mark;prof_mark=prof_now
endif
do FOURIER=0,NMOMENTS
 if(linearized)then
CALL MSSEG_LA_HOMOG_FOURIER ( &
           DO_LAMBERTIAN_SURFACE, DO_FOURIER0_NSTOKES2, ASYMTX_TOLERANCE,           &
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NMOMENTS, LVARYFLAGS, LVARYNUMBERS, &
           NSTREAMS_2, NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,  &
           MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING,           &
           QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS, ALBEDO, BRDF_F,                &
           DELTAU_VERT, OMEGA_GREEK, L_DELTAU_VERT, L_OMEGA_GREEK,                  &
           DO_INCLUDE_SURFACE, HGSF, HSOL, HBVP, L_HSOL,                            &
           STATUS, MESSAGE, TRACE_1, TRACE_2 )
else
CALL MSSEG_HOMOG_FOURIER ( &
           DO_LAMBERTIAN_SURFACE, DO_FOURIER0_NSTOKES2, ASYMTX_TOLERANCE, &
           FOURIER, NSTOKES, NSTREAMS, NLAYERS, NMOMENTS, NSTREAMS_2,     &
           NSTKS_NSTRMS, NSTKS_NSTRMS_2, NTOTAL, N_SUBDIAG, N_SUPDIAG,    &
           MUELLER_INDEX, DMAT, DO_REAL_EIGENSOLVER, DO_LAYER_SCATTERING, &
           QUAD_STREAMS, QUAD_STRMWTS, QUAD_HALFWTS,                      &
           ALBEDO, BRDF_F, DELTAU_VERT, OMEGA_GREEK,           &
           DO_INCLUDE_SURFACE, HGSF, HSOL, HBVP,                          &
           STATUS, MESSAGE, TRACE_1, TRACE_2 )
endif
call checked(STATUS)
if(profile)then
 call system_clock(prof_now);prof_ticks(2)=prof_ticks(2)+prof_now-prof_mark;prof_mark=prof_now
endif
LL_NSTK=NSTOKES
if(DO_FOURIER0_NSTOKES2.and.NSTOKES==3.and.FOURIER==0)LL_NSTK=2
LL_NSTKNSTRM=NSTREAMS*LL_NSTK;LL_NSTKNSTRM_2=2*LL_NSTKNSTRM
LL_NTT=LL_NSTKNSTRM_2*NLAYERS;LL_NSUBD=3*LL_NSTKNSTRM-1
if(NLAYERS==1)LL_NSUBD=2*LL_NSTKNSTRM-1
LL_NSUPD=LL_NSUBD
SURFACE_FACTOR=ONE;DELTA_FACTOR=TWO
if(FOURIER==0)then
 SURFACE_FACTOR=TWO;DELTA_FACTOR=ONE
endif
DO_REFLECTED_DIRECTBEAM(1:ngeoms)=DO_DIRECT_BEAM.and.DO_INCLUDE_SURFACE
CALL MSSEG_DIRECTBEAM ( &
        DO_LAMBERTIAN_SURFACE, .false., .false., &
        FOURIER, NSTOKES, NSTREAMS, NGEOMS, DELTA_FACTOR,           &
        FLUX_FACTOR, FLUXVEC, MUELLER_INDEX, COS_SZANGLES,          &
        ALBEDO, BRDF_F_0, SLTERM_ISOTROPIC, SLTERM_F_0,             &
        TRANS_SOLAR_BEAM, DO_REFLECTED_DIRECTBEAM, &
        ATMOS_ATTN, DIRECT_BEAM )
CALL MSSEG_PARTIC_PIMATRIX_SETUP ( &
        FOURIER, NSTOKES, NLAYERS, NGEOMS, NMOMENTS, &
        COS_SZANGLES, USER_STREAMS, DMAT, HGSF, PGSF )
! User homogeneous coefficients are shared across every layer evaluation.
uc=0;luc=0
do NCHOICE=1,NLAYERS
if(.not.any(needed(:,NCHOICE)))cycle
PSOL%UHOM_DNDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)=0
PSOL%UHOM_DNUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)=0
PSOL%UHOM_UPDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)=0
PSOL%UHOM_UPUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)=0
L_PSOL%L_UHOM_DNDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=0
L_PSOL%L_UHOM_DNUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=0
L_PSOL%L_UHOM_UPDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=0
L_PSOL%L_UHOM_UPUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=0
CALL  MSSEG_UHOM_SOLUTION ( &
        NCHOICE, FOURIER, LL_NSTK, NSTREAMS, NMOMENTS, NGEOMS, &
        .true., .true., DO_LAYER_SCATTERING, &
        QUAD_HALFWTS, OMEGA_GREEK, HGSF, PGSF, HSOL, PSOL )
if(linearized)then
 DO_RTSOL_VARY=LVARYFLAGS(NCHOICE);NPARAMS_VARY=LVARYNUMBERS(NCHOICE)
CALL MSSEG_LAP_UHOM_SOLUTION ( &
            NCHOICE, FOURIER, DO_RTSOL_VARY, NPARAMS_VARY,   &
            LL_NSTK, NSTREAMS, NMOMENTS, NGEOMS,             &
            .true., .true., QUAD_HALFWTS,        &
            DO_LAYER_SCATTERING, OMEGA_GREEK, L_OMEGA_GREEK, &
            HGSF, PGSF, HSOL, PSOL, L_HSOL, L_PSOL )
endif
uc(:,:,:,NCHOICE,1)=PSOL%UHOM_DNDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)
if(linearized)luc(:,:,:,:,NCHOICE,1)=L_PSOL%L_UHOM_DNDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)
uc(:,:,:,NCHOICE,2)=PSOL%UHOM_DNUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)
if(linearized)luc(:,:,:,:,NCHOICE,2)=L_PSOL%L_UHOM_DNUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)
uc(:,:,:,NCHOICE,3)=PSOL%UHOM_UPDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)
if(linearized)luc(:,:,:,:,NCHOICE,3)=L_PSOL%L_UHOM_UPDN(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)
uc(:,:,:,NCHOICE,4)=PSOL%UHOM_UPUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS)
if(linearized)luc(:,:,:,:,NCHOICE,4)=L_PSOL%L_UHOM_UPUP(1:ngeoms,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)
enddo
! Each sample solves beam/BVP once; all source layers reuse those solutions.
if(profile)then
 call system_clock(prof_now);prof_ticks(3)=prof_ticks(3)+prof_now-prof_mark;prof_mark=prof_now
endif
do IGEOM=1,ngeoms
 do LAYER=1,NLAYERS
CALL MSSEG_QBEAM_SOLUTION ( &
              LAYER, FOURIER, IGEOM, LL_NSTK, NSTREAMS, NMOMENTS,            &
              NSTREAMS_2, LL_NSTKNSTRM, DO_LAYER_SCATTERING, DFLUX,          &
              FLUX_FACTOR, QUAD_STREAMS, OMEGA_GREEK, BEAM_CUTOFF,           &
              T_DELT_MUBAR, INITIAL_TRANS, AVERAGE_SECANT, PGSF, HGSF, HSOL, &
              PSOL, STATUS_SUB, MESSAGE, TRACE_1 )
call checked(STATUS_SUB)
 if(linearized)then
 DO_RTSOL_VARY=LVARYFLAGS(LAYER);NPARAMS_VARY=LVARYNUMBERS(LAYER)
CALL MSSEG_LP_QBEAM_SOLUTION ( &
                LAYER, FOURIER, IGEOM, DO_RTSOL_VARY, NPARAMS_VARY, DO_PLANE_PARALLEL, &
                LVARYFLAGS, LVARYNUMBERS, LL_NSTK, NSTREAMS, FLUX_FACTOR,              &
                NMOMENTS, NSTREAMS_2, LL_NSTKNSTRM, DFLUX,          &
                DO_LAYER_SCATTERING, QUAD_STREAMS, L_OMEGA_GREEK,   &
                BEAM_CUTOFF, AVERAGE_SECANT, LP_AVERAGE_SECANT,     &
                HGSF, PGSF, HSOL, PSOL, L_HSOL,                     &
                L_PSOL, STATUS_SUB, MESSAGE, TRACE_1 )
call checked(STATUS_SUB)
endif
enddo
CALL MSSEG_BVPSOLUTION_MASTER ( &
            DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IGEOM), DO_LAMBERTIAN_SURFACE, &
            IGEOM, FOURIER, NCHOICE, LL_NSTK, NSTREAMS, NLAYERS, SURFACE_FACTOR,       &
            NSTREAMS_2, LL_NTT, LL_NSUBD, LL_NSUPD, LL_NSTKNSTRM, LL_NSTKNSTRM_2,      &
            ALBEDO, QUAD_STRMWTS, MUELLER_INDEX, DIRECT_BEAM,                          &
            HSOL, PSOL, HBVP, PBVP, STATUS_SUB, MESSAGE, TRACE_1 )
call checked(STATUS_SUB)
if(profile)then
 call system_clock(prof_now);prof_ticks(4)=prof_ticks(4)+prof_now-prof_mark;prof_mark=prof_now
endif
nc=0;pc=0
if(linearized)then
 do KV=1,NLAYERS
  if(.not.LVARYFLAGS(KV))cycle
  QV=LVARYNUMBERS(KV)
CALL MSSEG_LP_BVPSOLUTION_MASTER ( &
                   DO_INCLUDE_SURFACE, DO_REFLECTED_DIRECTBEAM(IGEOM), DO_LAMBERTIAN_SURFACE,   &
                   DO_PLANE_PARALLEL, FOURIER, IGEOM, LL_NSTK, NSTREAMS, NLAYERS, KV, QV,       &
                   NSTREAMS_2, LL_NTT, LL_NSUBD, LL_NSUPD, LL_NSTKNSTRM, LL_NSTKNSTRM_2,        &
                   SURFACE_FACTOR, ALBEDO, BRDF_F, QUAD_STRMWTS, DIRECT_BEAM,                   &
                   DO_LAYER_SCATTERING, MUELLER_INDEX, DELTAU_SLANT, L_DELTAU_VERT,             &
                   BEAM_CUTOFF, T_DELT_MUBAR, INITIAL_TRANS, LP_INITIAL_TRANS, LP_T_DELT_MUBAR, &
                   HSOL, HBVP, PSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, STATUS_SUB, MESSAGE, TRACE_1 )
call checked(STATUS_SUB)
  nc(:,:,1:QV,KV)=L_PBVP%NCON(1:NSTKS_NSTRMS,1:NLAYERS,1:QV)
  pc(:,:,1:QV,KV)=L_PBVP%PCON(1:NSTKS_NSTRMS,1:NLAYERS,1:QV)
! L_WLOWER is specific to KV and is overwritten by the next linearized BVP.
  if(FOURIER==0.and.up(IGEOM))then
CALL MSSEG_LP_POSTPROCESS_SURFACEMSST ( &
                     DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE, IGEOM, FOURIER, &
                     LL_NSTK, NSTREAMS, NLAYERS, KV, QV, MUELLER_INDEX,         &
                     SURFACE_FACTOR, QUAD_STRMWTS, ALBEDO, USER_BRDF_F,         &
                     HSOL, PBVP, L_HSOL, L_PSOL, L_PBVP, LP_SURFACE_MSST_F )
   do Q=1,QV
    surface_jac(Q,KV,IGEOM,1:NSTOKES)=LP_SURFACE_MSST_F(KV,Q,IGEOM,1:NSTOKES)
   enddo
  endif
 enddo
endif
if(profile)then
 call system_clock(prof_now);prof_ticks(5)=prof_ticks(5)+prof_now-prof_mark;prof_mark=prof_now
endif
az=1
if(FOURIER>0)then
 az(1:2)=TWO*cos(DEG_TO_RAD*USER_RELAZMS(IGEOM)*FOURIER)
 az(3:4)=TWO*sin(DEG_TO_RAD*USER_RELAZMS(IGEOM)*FOURIER)
endif
if(FOURIER==0.and.up(IGEOM))then
CALL  MSSEG_POSTPROCESS_SURFACEMSST ( &
               DO_INCLUDE_SURFACE, DO_LAMBERTIAN_SURFACE,   &
               IGEOM, FOURIER, LL_NSTK, NSTREAMS, NLAYERS,  &
               MUELLER_INDEX, SURFACE_FACTOR, QUAD_STRMWTS, &
               ALBEDO, USER_BRDF_F, HSOL, PSOL, PBVP,       &
               STOKES_DOWNSURF, SURFACE_MSST_F )
 surface_ms(IGEOM,1:NSTOKES)=SURFACE_MSST_F(IGEOM,1:NSTOKES)

endif
 do NCHOICE=1,NLAYERS
 if(.not.needed(IGEOM,NCHOICE))cycle
PSOL%UHOM_DNDN(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS)=uc(IGEOM,:,:,NCHOICE,1)
if(linearized)L_PSOL%L_UHOM_DNDN(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=luc(IGEOM,:,:,:,NCHOICE,1)
PSOL%UHOM_DNUP(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS)=uc(IGEOM,:,:,NCHOICE,2)
if(linearized)L_PSOL%L_UHOM_DNUP(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=luc(IGEOM,:,:,:,NCHOICE,2)
PSOL%UHOM_UPDN(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS)=uc(IGEOM,:,:,NCHOICE,3)
if(linearized)L_PSOL%L_UHOM_UPDN(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=luc(IGEOM,:,:,:,NCHOICE,3)
PSOL%UHOM_UPUP(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS)=uc(IGEOM,:,:,NCHOICE,4)
if(linearized)L_PSOL%L_UHOM_UPUP(IGEOM,1:NSTOKES,1:NSTKS_NSTRMS,1:NQ)=luc(IGEOM,:,:,:,NCHOICE,4)
PSOL%UPAR_UP_2=0
PSOL%UPAR_DN_2=0
L_PSOL%LP_UPAR_UP_2=0
L_PSOL%LP_UPAR_DN_2=0
CALL MSSEG_UBEAM_SOLUTION ( &
            NCHOICE, FOURIER, IGEOM,                         &
            up(IGEOM), .not.up(IGEOM), DO_LAYER_SCATTERING, &
            LL_NSTK, NSTREAMS, NMOMENTS, FLUX_FACTOR,        &
            BEAM_CUTOFF, QUAD_HALFWTS, OMEGA_GREEK,          &
            HGSF, PGSF, PSOL )
if(linearized)then
 DO_RTSOL_VARY=LVARYFLAGS(NCHOICE);NPARAMS_VARY=LVARYNUMBERS(NCHOICE)
CALL MSSEG_LP_UBEAM_SOLUTION ( &
              NCHOICE, FOURIER, IGEOM, DO_RTSOL_VARY, NPARAMS_VARY, LL_NSTK, NSTREAMS, &
              DO_PLANE_PARALLEL, LVARYFLAGS, LVARYNUMBERS, up(IGEOM), .not.up(IGEOM), &
              FLUX_FACTOR, DO_LAYER_SCATTERING, NMOMENTS, BEAM_CUTOFF, QUAD_HALFWTS,   &
              OMEGA_GREEK, L_OMEGA_GREEK, HGSF, PGSF, PSOL, L_PSOL )
endif
call integrated_source(NCHOICE,IGEOM,NSTOKES,NLAYERS,NQ,up(IGEOM),chord(NCHOICE,IGEOM), &
 DELTAU_VERT,L_DELTAU_VERT,BEAM_CUTOFF,INITIAL_TRANS,AVERAGE_SECANT, &
 LP_INITIAL_TRANS,LP_AVERAGE_SECANT,HSOL,PSOL,PBVP,L_HSOL,L_PSOL,nc,pc,value,deriv)
layer_ms(IGEOM,1:NSTOKES,NCHOICE)=layer_ms(IGEOM,1:NSTOKES,NCHOICE)+az(1:NSTOKES)*value(1:NSTOKES)
do N=1,NLAYERS
 do Q=1,NQ
  layer_jac(Q,N,IGEOM,1:NSTOKES,NCHOICE)=layer_jac(Q,N,IGEOM,1:NSTOKES,NCHOICE)+ &
                                      az(1:NSTOKES)*deriv(Q,N,1:NSTOKES)
 enddo
enddo
 enddo
if(profile)then
 call system_clock(prof_now);prof_ticks(6)=prof_ticks(6)+prof_now-prof_mark;prof_mark=prof_now
endif
enddo
enddo
if(profile)write(*,*)'AC shared stages setup/hom/user/beam/lp/source:',real(prof_ticks,fpk)/real(prof_rate,fpk)
deallocate(uc,luc,nc,pc,chord)
end subroutine aircraft_ms_sources

subroutine checked(status)
integer,intent(in)::status
if(status/=VQS_SUCCESS)then
 write(*,*)'AC shared MS kernel failed, status=',status
 error stop 8
endif
end subroutine

subroutine integrated_source(n,g,ns,nl,nq,up,rule,tau,ltau,cutoff,initial,secant, &
                             linitial,lsecant,hs,ps,pb,lhs,lps,nc,pc,value,deriv)
integer,intent(in)::n,g,ns,nl,nq,cutoff(MAXBEAMS)
logical,intent(in)::up
type(chord_rule),intent(in)::rule
real(fpk),intent(in)::tau(MAXLAYERS),ltau(MAX_ATMOSWFS,MAXLAYERS)
real(fpk),intent(in)::initial(MAXLAYERS,MAXBEAMS),secant(MAXLAYERS,MAXBEAMS)
real(fpk),intent(in)::linitial(MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)
real(fpk),intent(in)::lsecant(MAXLAYERS,MAXLAYERS,MAXBEAMS,MAX_ATMOSWFS)
real(fpk),intent(in)::nc(:,:,:,:),pc(:,:,:,:)
type(Homog_Solutions),intent(in)::hs
type(Partic_Solutions),intent(in)::ps
type(Partic_BVPSolve),intent(in)::pb
type(L_Homog_Solutions),intent(in)::lhs
type(L_Partic_Solutions),intent(in)::lps
real(fpk),intent(out)::value(4),deriv(MAX_ATMOSWFS,MAXLAYERS,4)
real(fpk)::u(4,MAXEVALUES),v(4,MAXEVALUES),lu(4,MAXEVALUES,MAX_ATMOSWFS)
real(fpk)::lv(4,MAXEVALUES,MAX_ATMOSWFS),part(4),lpart(4,MAXLAYERS,MAX_ATMOSWFS),dt
complex(fpk)::k,lk,wd,wu,kd,ku,td,tu,lc,mc,uc,vc,dc,ec,dw,ew
integer::i,j,o,q,kv,last
value=ZERO;deriv=ZERO
u=ZERO;v=ZERO;lu=ZERO;lv=ZERO;part=ZERO;lpart=ZERO
if(up)then
 u(1:ns,:)=ps%UHOM_UPDN(g,1:ns,:);v(1:ns,:)=ps%UHOM_UPUP(g,1:ns,:)
 part(1:ns)=ps%UPAR_UP_2(1:ns)
 if(nq>0)then
  lu(1:ns,:,1:nq)=lps%L_UHOM_UPDN(g,1:ns,:,1:nq)
  lv(1:ns,:,1:nq)=lps%L_UHOM_UPUP(g,1:ns,:,1:nq)
  lpart(1:ns,1:nl,1:nq)=lps%LP_UPAR_UP_2(1:ns,1:nl,1:nq)
 endif
else
 u(1:ns,:)=ps%UHOM_DNDN(g,1:ns,:);v(1:ns,:)=ps%UHOM_DNUP(g,1:ns,:)
 part(1:ns)=ps%UPAR_DN_2(1:ns)
 if(nq>0)then
  lu(1:ns,:,1:nq)=lps%L_UHOM_DNDN(g,1:ns,:,1:nq)
  lv(1:ns,:,1:nq)=lps%L_UHOM_DNUP(g,1:ns,:,1:nq)
  lpart(1:ns,1:nl,1:nq)=lps%LP_UPAR_DN_2(1:ns,1:nl,1:nq)
 endif
endif
i=1;last=hs%K_REAL(n)+2*hs%K_COMPLEX(n)
do while(i<=last)
 j=0
 if(i>hs%K_REAL(n))j=i+1
 k=cmplx(hs%KEIGEN(i,n),ZERO,fpk)
 lc=cmplx(pb%LCON(i,n),ZERO,fpk);mc=cmplx(pb%MCON(i,n),ZERO,fpk)
 if(j>0)then
  k=cmplx(hs%KEIGEN(i,n),hs%KEIGEN(j,n),fpk)
  lc=cmplx(pb%LCON(i,n),pb%LCON(j,n),fpk);mc=cmplx(pb%MCON(i,n),pb%MCON(j,n),fpk)
 endif
 call chord_weights(k,tau(n),rule,wd,wu,kd,ku,td,tu,nq>0)
 do o=1,ns
  uc=cmplx(u(o,i),ZERO,fpk);vc=cmplx(v(o,i),ZERO,fpk)
  if(j>0)then
   uc=cmplx(u(o,i),u(o,j),fpk);vc=cmplx(v(o,i),v(o,j),fpk)
  endif
  value(o)=value(o)+real(lc*uc*wd+mc*vc*wu,fpk)
  do kv=1,nl
   do q=1,nq
    dc=cmplx(nc(i,n,q,kv),ZERO,fpk);ec=cmplx(pc(i,n,q,kv),ZERO,fpk)
    if(j>0)then
     dc=cmplx(nc(i,n,q,kv),nc(j,n,q,kv),fpk)
     ec=cmplx(pc(i,n,q,kv),pc(j,n,q,kv),fpk)
    endif
    dc=dc*uc;ec=ec*vc;dw=ZERO;ew=ZERO
    if(kv==n)then
     lk=cmplx(lhs%L_KEIGEN(i,n,q),ZERO,fpk)
     if(j>0)then
      lk=cmplx(lhs%L_KEIGEN(i,n,q),lhs%L_KEIGEN(j,n,q),fpk)
      dc=dc+lc*cmplx(lu(o,i,q),lu(o,j,q),fpk)
      ec=ec+mc*cmplx(lv(o,i,q),lv(o,j,q),fpk)
     else
      dc=dc+lc*lu(o,i,q);ec=ec+mc*lv(o,i,q)
     endif
     dt=tau(n)*ltau(q,n);dw=kd*lk+td*dt;ew=ku*lk+tu*dt
    endif
    deriv(q,kv,o)=deriv(q,kv,o)+real(dc*wd+lc*uc*dw+ec*wu+mc*vc*ew,fpk)
   enddo
  enddo
 enddo
 i=i+1
 if(j>0)i=i+1
enddo
if(n<=cutoff(g))then
 k=cmplx(secant(n,g),ZERO,fpk)
 call chord_weights(k,tau(n),rule,wd,wu,kd,ku,td,tu,nq>0)
 value(1:ns)=value(1:ns)+part(1:ns)*initial(n,g)*real(wd,fpk)
 do kv=1,n
  do q=1,nq
   dt=ZERO
   if(kv==n)dt=tau(n)*ltau(q,n)
   dw=kd*lsecant(n,kv,g,q)+td*dt+wd*linitial(n,kv,g,q)
   deriv(q,kv,1:ns)=deriv(q,kv,1:ns)+initial(n,g)* &
      (lpart(1:ns,kv,q)*real(wd,fpk)+part(1:ns)*real(dw,fpk))
  enddo
 enddo
endif
end subroutine integrated_source
end module aircraft_ms_sources_m
