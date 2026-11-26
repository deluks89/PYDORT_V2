! Direct in-process spherical FO adapter for the AC driver.
! Uses the bundled GPL VQS geometry and FO kernels by R. Spurr/M. Christi.
! Geometry setup and optical mapping follow vqs_lp_masters_V1p4.f90;
! no VQS master or MS code is linked or called.
module aircraft_fo_m
   use FOMS_pars_m
   use VLIDORT_PARS_m, only: AC_MAXLAYERS => MAXLAYERS, AC_MAXWFS => MAX_ATMOSWFS
   use VLIDORT_IO_DEFS_m, only: VLIDORT_Fixed_Inputs, VLIDORT_Modified_Inputs
   use VLIDORT_LIN_IO_DEFS_m, only: VLIDORT_Fixed_LinInputs
   use FOMS_Geometry_Masters_m, only: FOMS_Geometry_Master
   use FOMS_VectorSpherfuncs_m, only: FOMS_VectorSpherfuncs
   use FOMS_FO_VectorRT_m, only: FOMS_FO_VectorRT
   use FOMS_FO_LP_VectorRT_m, only: FOMS_FO_LP_VectorRT
   implicit none
   private
   public :: aircraft_first_order
contains
   subroutine aircraft_first_order(fix, modin, linfix, linearized, height, sza, vza, raa, fine, scaled, stokes, jac)
      type(VLIDORT_Fixed_Inputs), intent(in) :: fix
      type(VLIDORT_Modified_Inputs), intent(in) :: modin
      type(VLIDORT_Fixed_LinInputs), intent(in) :: linfix
      logical, intent(in) :: linearized, scaled
      real(fpk), intent(in) :: height,sza,vza,raa
      integer, intent(in) :: fine
      real(fpk), intent(out) :: stokes(4),jac(AC_MAXWFS,AC_MAXLAYERS,4)
      real(fpk) :: lp_stokes(4,maxlays,maxjacobs)
      real(fpk) :: actdir(4),regdir(4),actdir_sb(4,0:maxsegments),regdir_sb(4,0:maxsegments)
      real(fpk) :: lp_actdir(4,maxlays,maxjacobs),lp_regdir(4,maxlays,maxjacobs)
      real(fpk) :: lp_actdir_sb(4,0:maxsegments,maxlays,maxjacobs),lp_regdir_sb(4,0:maxsegments,maxlays,maxjacobs)
      integer :: nlayers,nstokes,FO_nmoments,nfine_input,n,q,l,o1,o2,g,m,nq
      real(fpk) :: FO_pradius,FO_heights(0:maxlays),height_camera,alpha_camera,theta_camera,relaz_camera
      real(fpk) :: FO_extinction(maxlays),extinction_m(maxlays),FO_omega(maxlays),FO_truncfac(maxlays)
      real(fpk) :: FO_greekmat(maxlays,0:maxmoms,4,4),FO_reflec(4,4),solarflux,fluxvec(4),tms
      real(fpk) :: L_extinction_m(maxlays,maxjacobs),FO_L_omega(maxlays,maxjacobs),FO_L_truncfac(maxlays,maxjacobs)
      real(fpk) :: FO_L_greekmat(maxlays,0:maxmoms,4,4,maxjacobs),Acrit
      logical :: FO_Lvaryflags(maxlays),FO_Lvarymoms(maxlays,maxjacobs),do_Profile_WFs,doCrit
      integer :: FO_Lvarynums(maxlays)
      logical :: do_sunlight,do_deltam_scaling,do_debug_FO_input,do_debug_fluxes,do_Calculate_fluxes
!  FO Geometry Output
!  ------------------

!  Regime index, number of segments

      integer    :: nsegments, Regime_idx
      logical    :: Exact_condition

!  Segment to layer map, segment up/down switch, Layer-mask (1=up,2=dn)

      integer    :: seglayer_map  ( maxsegments )
      logical    :: segupdown_map ( maxsegments )
      logical    :: layermask_map ( maxlays, 2 )

!  Radii, Ray constant, Lospaths

      real(fpk)  :: radii    ( 0:maxlays )
      real(fpk)  :: Raycon
      real(fpk)  :: Lospaths ( maxsegments )
      real(fpk)  :: Rsegs    ( 0:maxsegments )

!  LOS angles and trig ratios thereof
!    Convention, range always 0 to 90 degrees 

      real(fpk)  :: alpha    ( 0:maxsegments )
      real(fpk)  :: cosa     ( 0:maxsegments )
      real(fpk)  :: sina     ( 0:maxsegments )
      real(fpk)  :: cota     ( 0:maxsegments )

!  Critical values

!      logical    :: Crit_map    ( maxsegments )
!      real(fpk)  :: alpha_crits ( 0:maxsegments )
!      real(fpk)  :: Rad_crits   ( 0:maxsegments )

!  Quadrature values

      integer    :: nfinedivs ( maxsegments)
      real(fpk)  :: wfine     ( maxsegments, maxfine )
      real(fpk)  :: csqfine   ( maxsegments, maxfine )
      real(fpk)  :: cotfine   ( maxsegments, maxfine )
      real(fpk)  :: alphafine ( maxsegments, maxfine )
      real(fpk)  :: rsegsfine ( maxsegments, maxfine )

!  Solar Paths to segment boundaries

      integer    :: ntraverse ( 0:maxsegments )
      real(fpk)  :: sunpaths  ( 0:maxsegments, maxlays )

!  Solar Paths to fine levels

      integer    :: ntraversefine ( maxsegments, maxfine )
      real(fpk)  :: sunpathsfine  ( maxsegments, maxlays, maxfine )

!  Scattering angle cosine, segment-boundary and fine-grid angles

      real(fpk)  :: cosscat
      real(fpk)  :: theta     ( 0:maxsegments )
      real(fpk)  :: relaz     ( 0:maxsegments )
      real(fpk)  :: thetafine ( maxsegments, maxfine )
      real(fpk)  :: relazfine ( maxsegments, maxfine )

!  Exception handling

      logical            :: fail_geo
      integer            :: status_geo
      character(len=100) :: message_geo
      character(len=100) :: trace_geo

!  Spherical Functions
!  -------------------

!  Set starter flag before first call

      logical    :: starter

!  Solar zenith angle cosine

      real(fpk)  :: Mu0

!  Vsign

      real(fpk), parameter :: vsign_up =   one
      real(fpk), parameter :: vsign_dn = - one

!  Help

      real(fpk)  :: GSHELP(7,maxmoms)

!  Number of moments and spherical functions

      real(fpk)  :: Genspher_up ( 0:maxmoms, maxstokes )
      real(fpk)  :: Genspher_dn ( 0:maxmoms, maxstokes )
      real(fpk)  :: Rotations_up ( maxstokes )
      real(fpk)  :: Rotations_dn ( maxstokes )

!  FO code output
!  --------------

!  Stokes vector and integration factors

      real(fpk)  :: lostrans    ( maxsegments )
      real(fpk)  :: sfi_factors ( maxsegments, maxfine )
      real(fpk)  :: sfi_opdeps  ( maxsegments, maxfine )

      real(fpk)  :: L_lostrans    ( maxsegments, maxJacobs )
      real(fpk)  :: L_sfi_factors ( maxsegments, maxfine, maxJacobs )
      real(fpk)  :: L_sfi_opdeps  ( maxsegments, maxfine, maxJacobs )


      nlayers=fix%Cont%TS_NLAYERS; nstokes=fix%Cont%TS_NSTOKES
      FO_nmoments=modin%MCont%TS_NGREEK_MOMENTS_INPUT
      if (fine<1.or.fine>maxfine) error stop 'AC FO quadrature exceeds MAXFINELAYERS'
      nfine_input=fine; FO_pradius=modin%MChapman%TS_EARTH_RADIUS
      FO_heights=ZERO
      FO_heights(0:nlayers)=fix%Chapman%TS_HEIGHT_GRID(0:nlayers)
      height_camera=height; alpha_camera=max(1.0e-5_fpk,min(180.0_fpk-1.0e-5_fpk,180.0_fpk-vza))
      theta_camera=max(1.0e-5_fpk,sza); relaz_camera=raa
      doCrit=.true.; Acrit=1.0e-15_fpk
      do_sunlight=.true.; do_debug_FO_input=.false.; do_debug_fluxes=.false.; do_Calculate_fluxes=.false.
      do_deltam_scaling=scaled.and.modin%MBool%TS_DO_DELTAM_SCALING.and..not.modin%MBool%TS_DO_RAYLEIGH_ONLY
      nq=0
      if (linearized) nq=linfix%Cont%TS_N_TOTALPROFILE_WFS
      do_Profile_WFs=linearized
      FO_Lvaryflags=.false.; FO_Lvarynums=0; FO_Lvarymoms=.false.
      FO_Lvaryflags(1:nlayers)=linearized.and.linfix%Cont%TS_LAYER_VARY_FLAG(1:nlayers)
      FO_Lvarynums(1:nlayers)=nq
      FO_Lvarymoms(1:nlayers,1:nq)=.true.
      FO_extinction=ZERO; extinction_m=ZERO; FO_omega=ZERO; FO_truncfac=ZERO
      FO_greekmat=ZERO; FO_reflec=ZERO; FO_reflec(1,1)=fix%Optical%TS_LAMBERTIAN_ALBEDO
      L_extinction_m=ZERO; FO_L_omega=ZERO; FO_L_truncfac=ZERO
      if (linearized) FO_L_greekmat=ZERO
      fluxvec=ZERO; fluxvec(1)=ONE; solarflux=fix%Sunrays%TS_FLUX_FACTOR/(4.0_fpk*pie)
      m=2*fix%Cont%TS_NSTREAMS
      do n=1,nlayers
         FO_extinction(n)=fix%Optical%TS_DELTAU_VERT_INPUT(n)/(FO_heights(n-1)-FO_heights(n))
         FO_omega(n)=modin%MOptical%TS_OMEGA_TOTAL_INPUT(n)
         do l=0,FO_nmoments
            do o1=1,4
               do o2=1,4
                  g=4*(o1-1)+o2
                  FO_greekmat(n,l,o1,o2)=fix%Optical%TS_GREEKMAT_TOTAL_INPUT(l,n,g)
                  do q=1,nq
                     FO_L_greekmat(n,l,o1,o2,q)=FO_greekmat(n,l,o1,o2)* &
                        linfix%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,l,n,g)
                  enddo
               enddo
            enddo
         enddo
         if (do_deltam_scaling) FO_truncfac(n)=FO_greekmat(n,m,1,1)/real(2*m+1,fpk)
         tms=ONE-FO_omega(n)*FO_truncfac(n)
         extinction_m(n)=FO_extinction(n)*tms
         do q=1,nq
            FO_L_omega(n,q)=FO_omega(n)*linfix%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n)
            if (do_deltam_scaling) FO_L_truncfac(n,q)=FO_L_greekmat(n,m,1,1,q)/real(2*m+1,fpk)
            L_extinction_m(n,q)=extinction_m(n)*linfix%Optical%TS_L_DELTAU_VERT_INPUT(q,n) &
               -FO_extinction(n)*(FO_L_omega(n,q)*FO_truncfac(n)+FO_omega(n)*FO_L_truncfac(n,q))
         enddo
      enddo
      call FOMS_Geometry_Master &
       ( nlayers, FO_heights, FO_pradius, height_Camera,      & ! Input
         alpha_Camera, theta_Camera, relaz_Camera,            & ! Input
         nfine_input, doCrit, Acrit, FO_extinction,           & ! Input
         Regime_idx, Exact_condition, Raycon, nsegments,      & ! Output
         seglayer_map, segupdown_map, layermask_map, radii,   & ! Output
         Lospaths, Rsegs, alpha, sina, cosa, cota, nfinedivs, & ! Output
         wfine, rsegsfine, alphafine, csqfine, cotfine,       & ! Output
         sunpaths, ntraverse, sunpathsfine, ntraversefine,    & ! Output
         cosscat, theta, relaz, thetafine, relazfine,         & ! Output
         fail_geo, status_geo, message_geo, trace_geo )         ! Output

      if (fail_geo.or.status_geo==1) then
         write(*,*) 'AC first-order geometry: ',trim(message_geo),' ',trim(trace_geo)
         error stop 'Unsupported AC first-order geometry'
      endif
      if ( Regime_idx .eq. 1 ) then
        starter = .true.
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_UP, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_UP, GSHELP, GENSPHER_UP )                     ! Outputs
      else if ( Regime_idx .eq. 3 ) then
        starter = .true.
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_DN, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_DN, GSHELP, GENSPHER_DN )                     ! Outputs
      else
        starter = .true.
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, VSIGN_UP, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,    & ! Inputs
             ROTATIONS_UP, GSHELP, GENSPHER_UP )                     ! Outputs
        starter = .true.
        call FOMS_VectorSpherfuncs &
           ( STARTER, FO_NMOMENTS, NSTOKES, DO_SUNLIGHT, -VSIGN_DN, & ! Inputs
             THETA_Camera, ALPHA_Camera, RELAZ_Camera, COSSCAT,     & ! Inputs
             ROTATIONS_DN, GSHELP, GENSPHER_DN )                      ! Outputs
      endif


      Mu0=ZERO
      if (Regime_idx==1) Mu0=cos(theta(0))
      if (linearized) then
      call FOMS_FO_LP_VectorRT ( &
         do_debug_FO_input, do_debug_fluxes, do_Calculate_fluxes,                   & ! Inputs (5/30/22 flags added)
         do_sunlight, do_deltam_scaling, Regime_idx, nsegments, nlayers,            & ! Inputs
         nstokes, FO_nmoments, seglayer_map, segupdown_map, layermask_map,          & ! Inputs 
         nfinedivs, ntraverse, radii, Raycon, cota, theta, sunpaths,                & ! Inputs (5/30/22 theta added)
         Mu0, Genspher_up, Genspher_dn, Rotations_up, Rotations_dn,                 & ! Inputs
         ntraversefine, wfine, csqfine, cotfine, rsegsfine, sunpathsfine,           & ! Inputs
         do_Profile_Wfs, FO_Lvaryflags, FO_Lvarynums, FO_Lvarymoms,                 & ! Lin Control
         FO_reflec, extinction_m, FO_omega, FO_truncfac, FO_greekmat, solarflux,    & ! Inputs
         fluxvec, L_extinction_m, FO_L_omega, FO_L_truncfac, FO_L_greekmat,         & ! Inputs (Optical, linearized)
         stokes, lostrans, sfi_factors, sfi_opdeps,                  & ! Outputs
         lp_stokes, L_lostrans, L_sfi_factors, L_sfi_opdeps,      & ! Outputs (Linearized)
         actdir,             regdir,              & ! Flux    Outputs. Added, 5/30/22
         lp_actdir,       lp_regdir,        & ! LP_Flux Outputs. Added, 5/30/22
         actdir_sb,       regdir_sb,        & ! Flux_SB    Outputs. Added, 4/24/24
         lp_actdir_sb, lp_regdir_sb )   ! LP_Flux_SB Outputs. Added, 4/24/24


      else
         call FOMS_FO_VectorRT(.false.,.false.,.false., &
            do_sunlight,do_deltam_scaling,Regime_idx,nsegments,nlayers, &
            nstokes,FO_nmoments,seglayer_map,segupdown_map,layermask_map, &
            nfinedivs,ntraverse,radii,Raycon,cota,theta,sunpaths, &
            Mu0,Genspher_up,Genspher_dn,Rotations_up,Rotations_dn, &
            ntraversefine,wfine,csqfine,cotfine,rsegsfine,sunpathsfine, &
            FO_reflec,extinction_m,FO_omega,FO_truncfac,FO_greekmat,solarflux,fluxvec, &
            stokes,lostrans,sfi_factors,sfi_opdeps,actdir,regdir,actdir_sb,regdir_sb)
      endif
      jac=ZERO
      do q=1,nq
         do n=1,nlayers
            jac(q,n,:)=lp_stokes(:,n,q)
         enddo
      enddo
   end subroutine aircraft_first_order
end module aircraft_fo_m
