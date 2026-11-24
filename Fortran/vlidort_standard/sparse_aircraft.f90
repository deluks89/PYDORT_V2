module sparse_aircraft_m
   use VLIDORT_PARS_m
   use VLIDORT_IO_DEFS_m
   use VLIDORT_LIN_IO_DEFS_m
   use VLIDORT_MASTERS_m
   use VLIDORT_LPS_MASTERS_m
   use VLIDORT_AC_SPHERCORR_ROUTINES_m
   implicit none
   private
   public :: sparse_aircraft_ms
contains

   subroutine sparse_error(message)
      character(len=*), intent(in) :: message
      write(*,*) 'Sparse aircraft: ',message
      stop 8
   end subroutine sparse_error

   ! Stable 1-exp(-x), including optically thin layers.
   pure function loss(x) result(value)
      real(fpk), intent(in) :: x
      real(fpk) :: value
      if (abs(x) < 1.0e-5_fpk) then
         value = x*(ONE-x*(0.5_fpk-x*(ONE/6.0_fpk-x/24.0_fpk)))
      else
         value = ONE-exp(-x)
      endif
   end function loss

   subroutine sparse_aircraft_ms(fix, modin, sup, linfix, linmod, linsup, linearized, &
                                height, sza, vza, raa, points, stokes, jac)
      type(VLIDORT_Fixed_Inputs), intent(in) :: fix
      type(VLIDORT_Modified_Inputs), intent(in) :: modin
      type(VLIDORT_Sup_InOut), intent(in) :: sup
      type(VLIDORT_Fixed_LinInputs), intent(in) :: linfix
      type(VLIDORT_Modified_LinInputs), intent(in) :: linmod
      type(VLIDORT_LinSup_InOut), intent(in) :: linsup
      logical, intent(in) :: linearized
      real(fpk), intent(in) :: height, sza, vza, raa ! km, degrees; VZA is 0..180
      integer, intent(in) :: points
      real(fpk), intent(out) :: stokes(4), jac(MAX_ATMOSWFS,MAXLAYERS,4)
      type(VLIDORT_Fixed_Inputs) :: sf
      type(VLIDORT_Modified_Inputs) :: sm
      type(VLIDORT_Sup_InOut) :: ss
      type(VLIDORT_Fixed_LinInputs) :: lf
      type(VLIDORT_Modified_LinInputs) :: lm
      type(VLIDORT_LinSup_InOut) :: ls
      type(VLIDORT_Outputs) :: out
      type(VLIDORT_LinOutputs) :: lo
      real(fpk), allocatable :: source(:,:,:), derivative(:,:,:,:,:)
      real(fpk) :: mus(33), tau(MAXLAYERS), dtau(MAX_ATMOSWFS,MAXLAYERS)
      real(fpk) :: path(MAXLAYERS), target(MAXLAYERS)
      real(fpk) :: radius, impact, alpha, dtr, top, bottom, rt, rb, distance
      real(fpk) :: low, high, mu, sample_height, theta, angle, phi, dummy(4)
      real(fpk) :: observer, solar_angle, azimuth, observer_alpha
      real(fpk) :: truncation, omega, tau0, denom, ddenom, trans, frac, slab, dslab
      real(fpk) :: emissivity(4), previous(4), surface(4)
      real(fpk) :: surface_jac(MAX_ATMOSWFS,MAXLAYERS,4)
      integer :: nl, ns, nq, p, n, k, q, j, count, first, last, step, moment
      logical :: down

      nl = fix%Cont%TS_NLAYERS; ns = fix%Cont%TS_NSTOKES
      nq = linfix%Cont%TS_N_TOTALPROFILE_WFS
      if (.not. linearized) nq = 0
      if (points < 2 .or. points > 33) call sparse_error('points must be 2..33')
      down = vza < 90.0_fpk
      dtr = acos(-ONE)/180.0_fpk
      alpha = min(vza,180.0_fpk-vza)
      radius = modin%MChapman%TS_EARTH_RADIUS
      impact = (radius+height)*sin(alpha*dtr)
      stokes = ZERO; jac = ZERO; surface = ZERO; surface_jac = ZERO
      path = ZERO; target = ONE

      ! Clip the observer's layer geometrically, retaining the full atmosphere
      ! for the local diffuse-field solves and their nonlocal Jacobians.
      do n = 1, nl
         top = fix%Chapman%TS_HEIGHT_GRID(n-1)
         bottom = fix%Chapman%TS_HEIGHT_GRID(n)
         if (down) then
            bottom = max(bottom,height)
         else
            top = min(top,height)
         endif
         if (top <= bottom) cycle
         rt = radius+top; rb = radius+bottom
         if (impact >= rb) call sparse_error('turning/horizontal rays require QS')
         ! Rationalized shell intersection avoids cancellation near zenith.
         distance = (rt-rb)*(rt+rb)/(sqrt(rt*rt-impact*impact)+sqrt(rb*rb-impact*impact))
         path(n) = distance/(fix%Chapman%TS_HEIGHT_GRID(n-1)-fix%Chapman%TS_HEIGHT_GRID(n))
         target(n) = (top-bottom)/distance
      enddo
      if (.not. any(path(1:nl)>ZERO)) return
      low = minval(target(1:nl),mask=path(1:nl)>ZERO)
      high = maxval(target(1:nl),mask=path(1:nl)>ZERO)
      if (.not. down) then
         rb = radius+fix%Chapman%TS_HEIGHT_GRID(nl)
         low = min(low,sqrt(ONE-(impact/rb)**2))
      endif
      count = points
      if (high-low < 1.0e-12_fpk) count = 1
      do p = 1, count
         mus(p) = exp(log(low)+(log(high)-log(low))*real(p-1,fpk)/real(max(1,count-1),fpk))
      enddo
      allocate(source(count,4,nl),derivative(count,nq,nl,4,nl))
      source = ZERO; derivative = ZERO; dtau = ZERO

      do p = 1, count
         sf = fix; sm = modin; ss = sup
         lf = linfix; lm = linmod; ls = linsup
         sf%Bool%TS_DO_SPHERICITY_CORRECTION = .true.
         sf%Bool%TS_DO_MSSTS = .true.
         sf%Bool%TS_DO_PLANE_PARALLEL = .false.
         sf%Bool%TS_DO_FULLRAD_MODE = .true.
         sf%Bool%TS_DO_UPWELLING = .not. down
         sf%Bool%TS_DO_DNWELLING = down
         sf%Cont%TS_SPHERCORR_TYPE = '2'
         ! Only native MS source outputs are used here. FO is separately
         ! computed along the actual observer ray by the caller.
         sm%MBool%TS_DO_FOCORR = .true.
         sm%MBool%TS_DO_FOCORR_EXTERNAL = .false.
         sm%MBool%TS_DO_FOCORR_OUTGOING = .true.
         sm%MUserVal%TS_USER_LEVELS(1) = ZERO
         if (down) sm%MUserVal%TS_USER_LEVELS(1) = real(nl,fpk)
         mu = mus(p)
         if (count == 1) then
            sample_height = height
         else
            sample_height = impact/sqrt(ONE-mu*mu)-radius
         endif
         if (alpha < 1.0e-10_fpk) then
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
         sm%MSunrays%TS_SZANGLES(1) = theta
         sm%MUserVal%TS_USER_VZANGLES_INPUT(1) = angle
         sm%MUserVal%TS_USER_RELAZMS(1) = phi
         sm%MUserVal%TS_USER_OBSGEOMS_INPUT(1,:) = (/theta,angle,phi/)
         if (linearized) then
            call VLIDORT_LPS_MASTER(.false.,sf,sm,ss,out,lf,lm,ls,lo)
         else
            call VLIDORT_MASTER(.false.,sf,sm,ss,out)
         endif
         if (out%Status%TS_STATUS_INPUTCHECK /= VLIDORT_SUCCESS .or. &
             out%Status%TS_STATUS_CALCULATION /= VLIDORT_SUCCESS) then
            write(*,*) 'Sparse aircraft local source failed:',p
            write(*,*) trim(out%Status%TS_CHECKMESSAGES(1)),trim(out%Status%TS_MESSAGE)
            stop 8
         endif
         if (.not. down .and. p == 1) then
            surface(1:ns) = out%Main%TS_SURF_MSSTS(1,1:ns)
            if (linearized) surface_jac(1:nq,1:nl,1:ns) = lo%Prof%TS_LP_SURF_MSSTS(1:nq,1:nl,1,1:ns)
         endif
         do n = 1, nl
            ! Match VLIDORT delta-M extinction and its normalized optical
            ! derivatives. Interpolate J=M/(1-exp(-tau/mu)), not M itself.
            tau0 = sf%Optical%TS_DELTAU_VERT_INPUT(n)
            omega = sm%MOptical%TS_OMEGA_TOTAL_INPUT(n)
            truncation = ZERO
            moment = 2*sf%Cont%TS_NSTREAMS
            if (sm%MBool%TS_DO_DELTAM_SCALING) &
               truncation = sf%Optical%TS_GREEKMAT_TOTAL_INPUT(moment,n,1)/real(2*moment+1,fpk)
            tau(n) = tau0*(ONE-omega*truncation)
            if (tau(n) <= ZERO) call sparse_error('nonpositive scaled optical depth')
            do q = 1, nq
               dtau(q,n) = tau(n)*lf%Optical%TS_L_DELTAU_VERT_INPUT(q,n) - tau0*omega*truncation * &
                  (lf%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n)+lf%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,moment,n,1))
               if (.not. lf%Cont%TS_LAYER_VARY_FLAG(n)) dtau(q,n) = ZERO
            enddo
            denom = loss(tau(n)/mu)
            source(p,1:ns,n) = out%Main%TS_LAYER_MSSTS(1,1:ns,n)/denom
            if (linearized) then
               derivative(p,1:nq,1:nl,1:ns,n) = lo%Prof%TS_LP_LAYER_MSSTS(1:nq,1:nl,1,1:ns,n)/denom
               do q = 1, nq
                  ddenom = exp(-tau(n)/mu)*dtau(q,n)/mu
                  derivative(p,q,n,1:ns,n) = derivative(p,q,n,1:ns,n)-source(p,1:ns,n)*ddenom/denom
               enddo
            endif
         enddo
      enddo

      stokes = surface; jac = surface_jac
      first = 1; last = nl; step = 1
      if (.not. down) then
         first = nl; last = 1; step = -1
      endif
      do n = first, last, step
         if (path(n) == ZERO) cycle
         j = 1; frac = ZERO
         if (count > 1) then
            do while (j < count-1 .and. mus(j+1) < target(n))
               j = j+1
            enddo
            frac = (log(target(n))-log(mus(j)))/(log(mus(j+1))-log(mus(j)))
            frac = max(ZERO,min(ONE,frac))
         endif
         k = min(j+1,count)
         emissivity = (ONE-frac)*source(j,:,n)+frac*source(k,:,n)
         trans = exp(-tau(n)*path(n)); slab = loss(tau(n)*path(n))
         previous = stokes
         stokes = trans*previous+slab*emissivity
         if (linearized) then
            jac(1:nq,1:nl,:) = trans*jac(1:nq,1:nl,:)+slab * &
               ((ONE-frac)*derivative(j,1:nq,1:nl,:,n)+frac*derivative(k,1:nq,1:nl,:,n))
            do q = 1, nq
               dslab = trans*path(n)*dtau(q,n)
               jac(q,n,:) = jac(q,n,:)+dslab*(emissivity-previous)
            enddo
         endif
      enddo
      deallocate(source,derivative)
   end subroutine sparse_aircraft_ms
end module sparse_aircraft_m
