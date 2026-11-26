module sparse_aircraft_m
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use VLIDORT_PARS_m
   use VLIDORT_IO_DEFS_m
   use VLIDORT_LIN_IO_DEFS_m
   use aircraft_ms_sources_m, only: aircraft_ms_sources
   use VLIDORT_AC_SPHERCORR_ROUTINES_m
   implicit none
   private
   public :: sparse_aircraft_ms, limb_aircraft_ms
contains

   subroutine sparse_error(message)
      character(len=*), intent(in) :: message
      write(*,*) 'Sparse aircraft: ',message
      stop 8
   end subroutine sparse_error

   subroutine limb_aircraft_ms(fix, modin, sup, linfix, linmod, linsup, linearized, &
                              height, sza, vza, raa, points, stokes, jac)
      type(VLIDORT_Fixed_Inputs), intent(in) :: fix
      type(VLIDORT_Modified_Inputs), intent(in) :: modin
      type(VLIDORT_Sup_InOut), intent(in) :: sup
      type(VLIDORT_Fixed_LinInputs), intent(in) :: linfix
      type(VLIDORT_Modified_LinInputs), intent(in) :: linmod
      type(VLIDORT_LinSup_InOut), intent(in) :: linsup
      logical, intent(in) :: linearized
      real(fpk), intent(in) :: height, sza, vza, raa
      integer, intent(in) :: points
      real(fpk), intent(out) :: stokes(4), jac(MAX_ATMOSWFS,MAXLAYERS,4)
      ! Curved source quadrature evaluates partial layers directly. Keeping the
      ! original grid avoids vanishing tangent layers and their ill-conditioned BVPs.
      call sparse_aircraft_ms(fix,modin,sup,linfix,linmod,linsup,linearized, &
         height,sza,vza,raa,points,stokes,jac,.true.)
   end subroutine limb_aircraft_ms

   subroutine sparse_aircraft_ms(fix, modin, sup, linfix, linmod, linsup, linearized, &
                                height, sza, vza, raa, points, stokes, jac, limb)
      type(VLIDORT_Fixed_Inputs), intent(in) :: fix
      type(VLIDORT_Modified_Inputs), intent(in) :: modin
      type(VLIDORT_Sup_InOut), intent(in) :: sup
      type(VLIDORT_Fixed_LinInputs), intent(in) :: linfix
      type(VLIDORT_Modified_LinInputs), intent(in) :: linmod
      type(VLIDORT_LinSup_InOut), intent(in) :: linsup
      logical, intent(in) :: linearized
      ! Internal research mode: MS only on a ray whose impact parameter is
      ! above BOA. The caller must supply a spherical first-order solution.
      logical, intent(in), optional :: limb
      real(fpk), intent(in) :: height, sza, vza, raa ! km, degrees; VZA is 0..180
      integer, intent(in) :: points
      real(fpk), intent(out) :: stokes(4), jac(MAX_ATMOSWFS,MAXLAYERS,4)
      real(fpk), allocatable :: source(:,:,:), derivative(:,:,:,:,:)
      real(fpk) :: mus(33), tau(MAXLAYERS), dtau(MAX_ATMOSWFS,MAXLAYERS)
      real(fpk) :: path(MAXLAYERS), exit_path(MAXLAYERS), length_path(MAXLAYERS)
      real(fpk) :: radius, impact, alpha, dtr, top, bottom, rt, rb, distance
      real(fpk) :: low, high, mu, sample_height, theta, angle, phi, dummy(4)
      real(fpk) :: observer, solar_angle, azimuth, observer_alpha
      real(fpk) :: truncation, omega, tau0, denom, trans
      real(fpk) :: emissivity(4), previous(4), surface(4)
      real(fpk) :: surface_jac(MAX_ATMOSWFS,MAXLAYERS,4)
      integer :: nl, ns, nq, p, n, q, count, first, last, step, moment
      logical :: down
      logical :: turning_mode
      integer :: leg, nlegs, pass, ng, offsets(2), sample, nvisible
      integer :: all_first(66), all_count(66)
      real(fpk) :: all_angles(66,3), all_mu(66)
      real(fpk) :: all_exit(66,MAXLAYERS), all_length(66,MAXLAYERS)
      logical :: all_up(66), all_needed(66,MAXLAYERS)
      real(fpk), allocatable :: all_ms(:,:,:), all_jac(:,:,:,:,:), all_surface(:,:), all_surface_jac(:,:,:,:)
      real(fpk) :: tangent, clip_bottom, signed_t, observer_t, position(3), look(3), sun(3)
      real(fpk) :: costheta, cosscatter, cosphi, sign_mu

      nl = fix%Cont%TS_NLAYERS; ns = fix%Cont%TS_NSTOKES
      nq = linfix%Cont%TS_N_TOTALPROFILE_WFS
      if (.not. linearized) nq = 0
      if (points < 2 .or. points > 33) call sparse_error('points must be 2..33')
      down = vza < 90.0_fpk
      dtr = acos(-ONE)/180.0_fpk
      alpha = min(vza,180.0_fpk-vza)
      radius = modin%MChapman%TS_EARTH_RADIUS
      impact = (radius+height)*sin(alpha*dtr)
      turning_mode = .false.
      if (present(limb)) turning_mode = limb
      nlegs = 1
      look = (/-sin(vza*dtr),ZERO,cos(vza*dtr)/)
      sun = (/-sin(sza*dtr)*cos(raa*dtr),sin(sza*dtr)*sin(raa*dtr),cos(sza*dtr)/)
      observer_t = (radius+height)*cos(vza*dtr)
      cosscatter = dot_product(look,sun)
      if (turning_mode) then
         if (impact < radius+fix%Chapman%TS_HEIGHT_GRID(nl)) &
            call sparse_error('limb MS mode requires a ray that does not cross BOA')
         if (height > fix%Chapman%TS_HEIGHT_GRID(0)) &
            call sparse_error('limb MS observer must be at or below TOA')
         if (vza > 90.0_fpk) nlegs = 2
         tangent = impact-radius
         signed_t = sqrt(max(ZERO,(radius+fix%Chapman%TS_HEIGHT_GRID(0))**2-impact**2))
         position = (/ZERO,ZERO,radius+height/)+(signed_t-observer_t)*look
         ! The numerator of local solar cosine is affine along the ray.
         ! Positive endpoint values therefore exclude twilight on both legs.
         if (sza >= 90.0_fpk .or. dot_product(position,sun) <= ZERO) &
            call sparse_error('twilight limb paths are not supported')
      endif
      stokes = ZERO; jac = ZERO; surface = ZERO; surface_jac = ZERO

      ! Traverse from the far TOA to the observer. A turning ray crosses
      ! some shells twice; the near leg must retain the far-leg radiance
      ! and every nonlocal derivative, rather than resetting at the tangent.
      ng=0; all_needed=.false.
      ! First collect both legs, then solve their sparse geometries together.
      do pass=1,2
      do leg = 1, nlegs
      if (turning_mode) down = leg == 1
      path = ZERO; exit_path=ZERO; length_path=ZERO
      low=ONE;high=ZERO;nvisible=0

      ! Clip the observer's layer geometrically, retaining the full atmosphere
      ! for the local diffuse-field solves and their nonlocal Jacobians.
      do n = 1, nl
         top = fix%Chapman%TS_HEIGHT_GRID(n-1)
         bottom = fix%Chapman%TS_HEIGHT_GRID(n)
         if (turning_mode) then
            clip_bottom = tangent
            if (nlegs == 1) clip_bottom = height
            bottom = max(bottom,clip_bottom)
            if (.not. down) top = min(top,height)
         else
         if (down) then
            bottom = max(bottom,height)
         else
            top = min(top,height)
         endif
         endif
         if (top <= bottom) cycle
         rt = radius+top; rb = radius+bottom
         ! A tangent obtained by subtracting Earth radius can differ from
         ! an identical observer boundary by a few rounding bits. After
         ! adding the radius back, this is an empty shell: avoid 0/0 at
         ! its tangent instead of perturbing the physical viewing angle.
         if (rt <= rb) cycle
         if (.not. turning_mode .and. impact >= rb) call sparse_error('turning/horizontal rays require QS')
         ! Rationalized shell intersection avoids cancellation near zenith.
         distance = (rt-rb)*(rt+rb)/(sqrt(max(ZERO,rt*rt-impact*impact))+sqrt(max(ZERO,rb*rb-impact*impact)))
         path(n) = distance/(fix%Chapman%TS_HEIGHT_GRID(n-1)-fix%Chapman%TS_HEIGHT_GRID(n))
         length_path(n)=distance
         exit_path(n)=sqrt(max(ZERO,(rb-impact)*(rb+impact)))
         if(.not.down)exit_path(n)=-sqrt(max(ZERO,(rt-impact)*(rt+impact)))
         low=min(low,sqrt(max(ZERO,(rb-impact)*(rb+impact)))/rb)
         high=max(high,sqrt(max(ZERO,(rt-impact)*(rt+impact)))/rt)
         nvisible=nvisible+1
      enddo
      if (.not. any(path(1:nl)>ZERO)) then
         ! At BOA, retain diffuse surface reflection for a zero-length path.
         if (turning_mode.or.down.or.height/=fix%Chapman%TS_HEIGHT_GRID(nl)) cycle
         low=cos(alpha*dtr);high=low
      endif
      ! Sample actual local viewing cosines, including the tangent limit.
      ! Unlike whole-layer effective cosines these bounds remain continuous at 90 deg.
      count=points
      if(high-low<1.0e-12_fpk)count=1
      if(turning_mode.and.leg==2.and.nvisible==1)count=1
      do p=1,count
         if(count==1)then
            mus(p)=0.5_fpk*(low+high)
         else
            mus(p)=low+(high-low)*real(p-1,fpk)/real(count-1,fpk)
         endif
      enddo
      if(pass==1)then
         offsets(leg)=ng; ng=ng+count
      endif

      if(pass==1)then
         do p=1,count
         mu = mus(p)
         if (count == 1) then
            sample_height = height
         else
            sample_height = impact/sqrt(ONE-mu*mu)-radius
         endif
         if (mu >= ONE.or.alpha < 1.0e-10_fpk) then
            theta=sza;angle=alpha;phi=raa
         else if (turning_mode.or.down) then
            ! Local geometry from Cartesian vectors, with no extrapolation
            ! to a fictitious ground intersection. The scattering angle and
            ! LOS meridian plane are invariant along this straight ray.
            sign_mu = ONE
            if (.not. down) sign_mu = -ONE
            signed_t = sign_mu*impact*mu/sqrt(ONE-mu*mu)
            position = (/ZERO,ZERO,radius+height/)+(signed_t-observer_t)*look
            costheta = dot_product(position,sun)/sqrt(dot_product(position,position))
            costheta = max(-ONE,min(ONE,costheta))
            theta = acos(costheta)/dtr
            angle = acos(mu)/dtr
            denom = sqrt(max(ZERO,ONE-costheta*costheta))*sqrt(ONE-mu*mu)
            phi = ZERO
            if (denom > 1.0e-20_fpk) then
               cosphi = (cosscatter-sign_mu*mu*costheta)/denom
               phi = acos(max(-ONE,min(ONE,cosphi)))/dtr
               if (raa > 180.0_fpk) phi = 360.0_fpk-phi
            endif
         else if (alpha < 1.0e-10_fpk) then
            theta = sza; angle = alpha; phi = raa
         else if (down) then
            ! A virtual BOA at the observer avoids the supplier geometry
            ! helper's unnecessary extension of an upward-looking ray.
            call AC_to_BOATOA_geometry_conversion(radius,height,height,sample_height,sza,alpha,raa, &
               dummy(1),dummy(2),dummy(3),dummy(4),theta,angle,phi)
         else
            observer = height; solar_angle = sza; azimuth = raa; observer_alpha = alpha
            call AC_to_BOA_geometry_conversion(radius,sample_height,observer,solar_angle,observer_alpha,azimuth, &
               dummy(1),theta,angle,phi)
         endif
         if (theta >= 90.0_fpk) call sparse_error('twilight sources require a different solver')

            sample=offsets(leg)+p
            all_angles(sample,:)=(/theta,angle,phi/)
            all_mu(sample)=mu
            all_first(sample)=offsets(leg)+1;all_count(sample)=count
            all_exit(sample,:)=exit_path;all_length(sample,:)=length_path
            all_up(sample)=.not.down
            all_needed(sample,1:nl)=path(1:nl)>ZERO
         enddo
         cycle
      endif
      allocate(source(count,4,nl),derivative(count,nq,nl,4,nl))
      source=ZERO;derivative=ZERO;dtau=ZERO
      do n=1,nl
         tau0=fix%Optical%TS_DELTAU_VERT_INPUT(n)
         omega=modin%MOptical%TS_OMEGA_TOTAL_INPUT(n)
         truncation=ZERO;moment=2*fix%Cont%TS_NSTREAMS
         if(modin%MBool%TS_DO_DELTAM_SCALING) &
            truncation=fix%Optical%TS_GREEKMAT_TOTAL_INPUT(moment,n,1)/real(2*moment+1,fpk)
         tau(n)=tau0*(ONE-omega*truncation)
         if(tau(n)<=ZERO)call sparse_error('nonpositive scaled optical depth')
         do q=1,nq
            dtau(q,n)=tau(n)*linfix%Optical%TS_L_DELTAU_VERT_INPUT(q,n)-tau0*omega*truncation * &
               (linfix%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n)+linfix%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,moment,n,1))
            if(.not.linfix%Cont%TS_LAYER_VARY_FLAG(n))dtau(q,n)=ZERO
         enddo
      enddo
      do p=1,count
         sample=offsets(leg)+p;mu=mus(p)
         if(.not.down.and.p==1.and..not.turning_mode)then
            surface(1:ns)=all_surface(sample,1:ns)
            if(linearized)surface_jac(1:nq,1:nl,1:ns)=all_surface_jac(:,:,sample,1:ns)
         endif
         source(p,1:ns,:)=all_ms(sample,1:ns,:)
         if(linearized)derivative(p,1:nq,1:nl,1:ns,:)=all_jac(:,:,sample,1:ns,:)

      enddo

      if (.not. turning_mode) then
         stokes = surface; jac = surface_jac
      endif
      first = 1; last = nl; step = 1
      if (.not. down) then
         first = nl; last = 1; step = -1
      endif
      do n = first, last, step
         if (path(n) == ZERO) cycle
         emissivity=sum(source(:,:,n),dim=1)
         trans=exp(-tau(n)*path(n))
         previous=stokes
         stokes=trans*previous+emissivity
         if(linearized)then
            jac(1:nq,1:nl,:)=trans*jac(1:nq,1:nl,:)+sum(derivative(:,:,:,:,n),dim=1)
            do q=1,nq
               jac(q,n,:)=jac(q,n,:)-trans*path(n)*dtau(q,n)*previous
            enddo
         endif
      enddo
      deallocate(source,derivative)
      enddo
      if(pass==1)then
      allocate(all_ms(ng,4,nl),all_jac(nq,nl,ng,4,nl))
      allocate(all_surface(ng,4),all_surface_jac(nq,nl,ng,4))
      if(ng>0)call aircraft_ms_sources(fix,modin,linfix,linearized,ng, &
         all_angles(1:ng,:),all_up(1:ng),all_needed(1:ng,1:nl), &
         all_ms,all_jac,all_surface,all_surface_jac, &
         impact,all_mu(1:ng),all_first(1:ng),all_count(1:ng),all_exit(1:ng,:),all_length(1:ng,:))
      if(.not.all(ieee_is_finite(all_ms)))call sparse_error('nonfinite shared MS source')
      if(linearized.and..not.all(ieee_is_finite(all_jac)))call sparse_error('nonfinite shared MS derivative')
      endif
      enddo
      deallocate(all_ms,all_jac,all_surface,all_surface_jac)
   end subroutine sparse_aircraft_ms
end module sparse_aircraft_m
