! Independent geometric identities, exercised through the public QS master.
program check_qs_limb_geometry
   use FOMS_pars_m
   use FOMS_Geometry_Masters_m, only : FOMS_Geometry_Master
   implicit none
   integer :: regime, ns, nl, status, i, j, s, failures, cases
   integer :: layers(maxsegments), nf(maxsegments), nt(0:maxsegments), ntf(maxsegments,maxfine)
   logical :: exact, up(maxsegments), mask(maxlays,2), fail
   real(ffp) :: heights(0:maxlays), extinction(maxlays), radii(0:maxlays), b, scatter
   real(ffp) :: paths(maxsegments), r(0:maxsegments), a(0:maxsegments)
   real(ffp) :: sina(0:maxsegments), cosa(0:maxsegments), cota(0:maxsegments)
   real(ffp) :: wf(maxsegments,maxfine), rf(maxsegments,maxfine), af(maxsegments,maxfine)
   real(ffp) :: csqf(maxsegments,maxfine), cotf(maxsegments,maxfine)
   real(ffp) :: sun(0:maxsegments,maxlays), sunf(maxsegments,maxlays,maxfine)
   real(ffp) :: theta(0:maxsegments), azimuth(0:maxsegments)
   real(ffp) :: thetaf(maxsegments,maxfine), azimuthf(maxsegments,maxfine)
   real(ffp) :: observer, vza, radial_error, chord_error, x0, x1, expected
   real(ffp), parameter :: observers(3) = (/0.55_ffp,15.0_ffp,16.0_ffp/)
   real(ffp), parameter :: views(9) = (/10.0_ffp,89.0_ffp,89.999_ffp,90.0_ffp, &
                                      90.001_ffp,90.5_ffp,91.0_ffp,93.0_ffp,170.0_ffp/)
   character(len=256) :: message, trace

   nl = 6
   heights = zero
   heights(0:nl) = (/70.0_ffp,30.0_ffp,16.0_ffp,14.0_ffp,2.0_ffp,1.0_ffp,0.0_ffp/)
   extinction = 0.001_ffp
   failures = 0; cases = 0
   do i = 1, size(observers)
      observer = observers(i)
      do j = 1, size(views)
         vza = views(j); cases = cases+1
         call FOMS_Geometry_Master(nl,heights,6371.0_ffp,observer,180.0_ffp-vza,60.0_ffp,90.0_ffp, &
            4,.false.,0.001_ffp,extinction,regime,exact,b,ns,layers,up,mask,radii,paths,r, &
            a,sina,cosa,cota,nf,wf,rf,af,csqf,cotf,sun,nt,sunf,ntf,scatter,theta,azimuth, &
            thetaf,azimuthf,fail,status,message,trace)
         if (fail) then
            write(*,*) 'Geometry master failed:',observer,vza,trim(message)
            failures = failures+1
            cycle
         endif
         ! Every boundary must lie on the SAME ray (constant impact parameter).
         radial_error = maxval(abs(r(0:ns)*sina(0:ns)-b))
         chord_error = zero
         do s = 1, ns
            ! Distances from the tangent, independently from shell radii.
            x0 = sqrt(max(zero,(r(s-1)-b)*(r(s-1)+b)))
            x1 = sqrt(max(zero,(r(s)-b)*(r(s)+b)))
            expected = abs(x1-x0)
            chord_error = max(chord_error,abs(paths(s)-expected))
         enddo
         if (radial_error > 1.0e-7_ffp .or. chord_error > 1.0e-5_ffp) then
            failures = failures+1
            write(*,*) 'FAIL: observer, VZA, regime, impact error km, segment error km:', &
               observer,vza,regime,radial_error,chord_error
         endif
      enddo
   enddo
   if (failures > 0) then
      write(*,*) failures,'of',cases,'geometry cases failed'
      stop 1
   endif
   write(*,*) 'PASS:',cases,'QS geometry cases; constant ray impact and exact shell chords'
end program check_qs_limb_geometry
