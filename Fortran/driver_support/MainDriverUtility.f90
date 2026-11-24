MODULE MainDriverUtility_m
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: find_string_in_array
   PUBLIC :: external_observer_to_toa
   PUBLIC :: open_irradiance_files, close_irradiance_files
   PUBLIC :: write_irradiance_orders

CONTAINS

   subroutine open_irradiance_files(labels, do_jacobians, do_orders)
     character(len=*), intent(in) :: labels(:)
     logical, intent(in) :: do_jacobians
     logical, intent(in) :: do_orders
     integer :: q
     open(21,file='OUTPUT/irradiance_diffuse_upwelling.out',status='replace')
     open(22,file='OUTPUT/irradiance_diffuse_downwelling.out',status='replace')
     open(23,file='OUTPUT/irradiance_direct.out',status='replace')
     if (do_orders) then
        open(24,file='OUTPUT/irradiance_single_upwelling.out',status='replace')
        open(25,file='OUTPUT/irradiance_single_downwelling.out',status='replace')
        open(26,file='OUTPUT/irradiance_multiple_upwelling.out',status='replace')
        open(27,file='OUTPUT/irradiance_multiple_downwelling.out',status='replace')
     endif
     if (do_jacobians) then
        do q = 1,size(labels)
           open(200+q,file='OUTPUT/jacs_irradiance_diffuse_upwelling_'//trim(labels(q))//'.out',status='replace')
           open(300+q,file='OUTPUT/jacs_irradiance_diffuse_downwelling_'//trim(labels(q))//'.out',status='replace')
           open(400+q,file='OUTPUT/jacs_irradiance_direct_'//trim(labels(q))//'.out',status='replace')
           if (do_orders) then
              open(500+q,file='OUTPUT/jacs_irradiance_single_upwelling_'//trim(labels(q))//'.out',status='replace')
              open(600+q,file='OUTPUT/jacs_irradiance_single_downwelling_'//trim(labels(q))//'.out',status='replace')
              open(700+q,file='OUTPUT/jacs_irradiance_multiple_upwelling_'//trim(labels(q))//'.out',status='replace')
              open(800+q,file='OUTPUT/jacs_irradiance_multiple_downwelling_'//trim(labels(q))//'.out',status='replace')
           endif
        enddo
     endif
   end subroutine open_irradiance_files

   subroutine close_irradiance_files(n_jacobians, do_jacobians, do_orders)
     integer, intent(in) :: n_jacobians
     logical, intent(in) :: do_jacobians
     logical, intent(in) :: do_orders
     integer :: q
     close(21); close(22); close(23)
     if (do_orders) then
        close(24); close(25); close(26); close(27)
     endif
     if (do_jacobians) then
        do q = 1,n_jacobians
           close(200+q); close(300+q); close(400+q)
           if (do_orders) then
              close(500+q); close(600+q); close(700+q); close(800+q)
           endif
        enddo
     endif
   end subroutine close_irradiance_files

   subroutine write_irradiance_orders(samples, diffuse, step, nlayers, njacs, do_jacobians)
     ! At fixed extinction, scaling all SSA and surface albedo by x gives
     ! diffuse(x) = x*first_order + x**2*second_order + ... .
     ! Extrapolate diffuse(x)/x to zero, cancelling the x and x**2 terms.
     ! Samples are at x = h, 2h, 4h, h/2. Compare two cubic-error estimates.
     ! Apply the same linear combination to native profile Jacobians.
     real(kind=kind(1d0)), intent(in) :: samples(:,:), diffuse(:), step
     integer, intent(in) :: nlayers, njacs
     logical, intent(in) :: do_jacobians
     real(kind=kind(1d0)) :: single(size(diffuse)), coarse(size(diffuse)), multiple(size(diffuse))
     integer :: q, lo, hi, offset
     if (.not. all(ieee_is_finite(samples)) .or. .not. all(ieee_is_finite(diffuse))) &
        error stop 'Nonfinite irradiance scattering-order input.'
     coarse = (8d0*samples(:,1)-3d0*samples(:,2)+0.25d0*samples(:,3))/(3d0*step)
     single = (16d0*samples(:,4)-6d0*samples(:,1)+0.5d0*samples(:,2))/(3d0*step)
     if (.not. all(ieee_is_finite(single))) error stop 'Nonfinite irradiance scattering-order estimate.'
     if (any(abs(single-coarse) > 1d-9 + 5d-6*abs(single))) then
        q = maxloc(abs(single-coarse)/(1d-9+5d-6*abs(single)), dim=1)
        write(*,*) 'Scattering expansion index, fine, coarse: ', q, single(q), coarse(q)
        error stop 'Irradiance scattering expansion did not converge; adjust irradiance_scattering_step.'
     endif
     multiple = diffuse-single
     write(24,*) single(1)
     write(25,*) single(2)
     write(26,*) multiple(1)
     write(27,*) multiple(2)
     if (do_jacobians) then
        offset = njacs*nlayers
        do q = 1,njacs
           lo = 3+(q-1)*nlayers
           hi = lo+nlayers-1
           write(500+q,*) single(lo:hi)
           write(600+q,*) single(lo+offset:hi+offset)
           write(700+q,*) multiple(lo:hi)
           write(800+q,*) multiple(lo+offset:hi+offset)
        enddo
     endif
   end subroutine write_irradiance_orders

   subroutine external_observer_to_toa(radius, toa_m, observer_m, sza, vza, raa, do_flux)
     ! Angles are defined at the observer: VZA=0 zenith, 180 nadir.
     ! Trace the forward line of sight to the near TOA intersection and
     ! rotate the fixed solar and viewing directions into that local frame.
     ! real(8) follows each driver's compiler precision (including binary128).
     real(kind=kind(1d0)), intent(in) :: radius, toa_m
     real, intent(inout) :: observer_m
     real(kind=kind(1d0)), intent(inout) :: sza, vza, raa
     logical, intent(in) :: do_flux
     real(kind=kind(1d0)) :: ro, rt, dtr, impact, distance, root, mu, sunmu, normv, norms
     real(kind=kind(1d0)) :: los(3), sun(3), up(3), vt(3), st(3)

     if (real(observer_m,kind(1d0)) <= toa_m) return
     if (do_flux) error stop 'Above-TOA actinic flux or irradiance is not supported; request radiance only.'
     dtr = acos(-1d0) / 180d0
     ro = radius + real(observer_m,kind(1d0)) / 1000d0
     rt = radius + toa_m / 1000d0
     los = [sin(vza*dtr), 0d0, cos(vza*dtr)]
     sun = [sin(sza*dtr)*cos(raa*dtr), sin(sza*dtr)*sin(raa*dtr), cos(sza*dtr)]
     impact = ro * abs(los(1))
     if (los(3) >= 0d0 .or. impact >= rt) &
        error stop 'Above-TOA viewing ray does not enter the atmosphere.'
     root = sqrt((rt-impact)*(rt+impact))
     ! Rationalized near root avoids cancellation for observers close to TOA.
     distance = ((ro-rt)*(ro+rt)) / (-ro*los(3)+root)
     up = [distance*los(1), 0d0, ro+distance*los(3)] / rt
     up = up / sqrt(dot_product(up,up))
     mu = max(-1d0,min(1d0,dot_product(los,up)))
     sunmu = max(-1d0,min(1d0,dot_product(sun,up)))
     if (sunmu < 0d0) error stop 'Above-TOA conversion gives a nighttime TOA entry (SZA > 90 degrees).'
     vt = los - mu*up
     st = sun - sunmu*up
     normv = sqrt(dot_product(vt,vt))
     norms = sqrt(dot_product(st,st))
     if (normv > 1d-14 .and. norms > 1d-14) then
        raa = acos(max(-1d0,min(1d0,dot_product(vt,st)/(normv*norms)))) / dtr
     else
        raa = 0d0  ! Azimuth is immaterial for a vertical solar or viewing ray.
     endif
     sza = acos(sunmu) / dtr
     vza = acos(mu) / dtr
     observer_m = real(toa_m)
   end subroutine external_observer_to_toa

   FUNCTION find_string_in_array(target_string, string_array) RESULT(index)

     ! Arguments:
     !   target_string (CHARACTER, INTENT(IN)) :: The string to search for.
     !   string_array  (CHARACTER, DIMENSION(:), INTENT(IN)) :: The array of strings to search within.
     !
     ! Returns:
     !   index (INTEGER) :: The 1-based index of the first occurrence of the
     !                      target_string in string_array. Returns 0 if not found.

     CHARACTER(LEN=*), INTENT(IN) :: target_string
     CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: string_array
     INTEGER :: index
     INTEGER :: i

     index = 0 ! Initialize index to 0 (not found)

     DO i = 1, SIZE(string_array)
       IF (TRIM(string_array(i)) == TRIM(target_string)) THEN
         index = i
         EXIT ! Exit the loop once the string is found
       END IF
     END DO

   END FUNCTION find_string_in_array
   
END MODULE MainDriverUtility_m
