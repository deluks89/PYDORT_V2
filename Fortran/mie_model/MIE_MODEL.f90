program RTSMie_Trimodal_test

   use RTSMie_parameters_m
   use RTSMie_master_trimodal_m

   implicit none

!  List of Inputs
!  ==============


!      Do_Fmatrix        - 

   logical    :: Do_Expcoeffs                ! Boolean flag for computing Expansion Coefficients
   logical    :: Do_Fmatrix                  ! Boolean flag for computing F-matrix at equal-angles
   logical    :: do_refr_idx_linearization   ! Boolean flag for refractive index linearizations
   logical    :: do_psd_linearization        ! Boolean flag for PSD linearizations

!    psd_index      - Index for particle size distribution of spheres
!    psd_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    psd_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    psd_index = 3 : BIMODAL GAMMA with equal mode weights
!    psd_index = 4 : LOG-NORMAL with rg and sigma given
!    psd_index = 5 : LOG-NORMAL with reff and veff given
!    psd_index = 6 : POWER LAW
!    psd_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    psd_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

   integer        :: psd_indexes (3)
   real(KIND=dp)  :: psd_pars (3,3) ! first index is mode, second is parameter (r,sig,dummy)

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!  then Internal routine to calculate R1 and R2 (outputs)
!  If Not set, Use Input R1 and R2 for PSD limits.

   logical            :: FixR1R2(3)

! Minimum and maximum radii (Microns) of calculatio range, for each mode

   real    (KIND=dp)  :: R1(3), R2(3)

!  R1R2_cutoff particle size for setting R1 and R2 internally

   REAL    (KIND=dp)  :: R1R2_cutoff(3)

!  PSD quadrature control
!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.

   INTEGER            :: nblocks(3)
   INTEGER            :: nweights(3)

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (KIND=dp)  :: xparticle_limit

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER         :: n_Fmatrix_angles
   REAL (KIND=dp)  :: Fmatrix_angles(max_Mie_angles)

! PSD parameters. These will be filled from file before sroted into psd_pars

   real(kind=dp)   :: modal_fractions(3)
   real(kind=dp)   :: modal_radii(3)
   real(kind=dp)   :: modal_widths(3)

! Dummy variables to make main model function work

   logical, parameter       :: do_monodisperse = .false.
   real(kind=dp), parameter :: monoradius = 0.1
   logical, parameter       :: Do_Save_PSD = .false.


!  OUTPUT ARGUMENTS
!  ================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(kind=dp)  :: TriMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter
!  ----------------------------------------------

!  Regular quantities

   integer       :: TriMie_ncoeffs
   real(KIND=dp) :: TriMie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp) :: TriMie_asymm

!  F-matrix,  optional output
!  --------------------------

!  F-matrix

   real(KIND=dp) :: TriMie_Fmatrix(4,max_Mie_angles)
   real(kind=dp)  :: TriMie_dist (5,3)

!  Exception handling
   logical       :: fail
   integer       :: istatus
   character*90  :: Trimessages(3)
   character*90  :: trace_3

!  Local variables
!  ===============

   logical            :: apply
   integer            :: nwav, iwav, L, i, j
   real(kind=dp)      :: wavelengths(max_wavelengths)
   real(kind=dp)      :: n_real_arr(max_wavelengths), n_imag_arr(max_wavelengths)
   real(kind=dp)      :: n_real(3), n_imag(3), wl
   
   real(kind=dp)      :: div

   real(kind=dp) :: moment_cutoff = 1.0d-20

! Read wavelength grid and refractive indexes
! ==========================================

   open (1, file = 'INPUT/refractive_indexes.inp', status = 'old')
   read(1,*) nwav
   read(1,*)
   read(1,*) ! Read header away
   do iwav = 1, nwav
      read(1,*) wavelengths(iwav), n_real_arr(iwav), n_imag_arr(iwav)
   enddo
   
   ! Convert wavelengths to microns
   do iwav = 1, nwav
      wavelengths(iwav) = wavelengths(iwav)/1e3
   enddo


! Read static configuration file
! ==============================

   open(1,file = 'static_configuration.inp', status = 'old')
   read(1,*) ! Header
   read(1,*) ! Header
   read(1,*) ! Blank
   read(1,*) Do_expcoeffs
   read(1,*) Do_Fmatrix
   read(1,*) do_refr_idx_linearization
   read(1,*) do_psd_linearization
   read(1,*) psd_indexes
   read(1,*) FixR1R2
   read(1,*) R1R2_cutoff
   read(1,*) nblocks
   read(1,*) nweights
   read(1,*) xparticle_limit
   read(1,*) n_Fmatrix_angles
   close(1)
   !  Fmatrix angles are regular here between 0 and 180.
   if ( do_Fmatrix ) then
      div = 180.0d0 / dble( n_Fmatrix_angles - 1 )
      do j = 1, n_Fmatrix_angles
         Fmatrix_angles(j) = dble(j-1) * div
      enddo
   else
      Fmatrix_angles = 0.0d0
   endif


!  Read configuration
!  ==================

   open(1,file = 'INPUT/configuration.inp', status = 'old')
   read(1,*) ! Header
   read(1,*) ! Header
   read(1,*) ! Header
   read(1,*) 
   read(1,*) modal_radii
   read(1,*) modal_widths
   read(1,*) modal_fractions
   read(1,*) R1
   read(1,*) R2
   close(1)

! Open result files and write headers
! ===================================

   open(1,file='OUTPUT/bulk_properties.out',status = 'unknown')
   write(1,*) 'wavelength   extinction_xs   single_scattering_albedo   asymmetry_parameter'

   open(2,file='OUTPUT/phasefunction_moments.out',status = 'unknown')
   write(2,*)  'wavelength   moment_order   alpha1   alpha2   alpha3   alpha4   beta1   beta2'
   
   open(3,file='OUTPUT/phasefunction.out', status = 'unknown')
   write(3,*) 'wavelength   scattering_angle   F11   F21   F33   F34'


! Start wavelength loop
! =====================

   do iwav = 1, nwav

      write(*,*)'Doing wavelength #', iwav, wavelengths(iwav)

      psd_pars(1,:) = modal_radii(:)
      psd_pars(2,:) = modal_widths(:)
      psd_pars(3,:) = 0.0d0
      
      wl = wavelengths(iwav)
      ! Apply same refracive index to all modes
      do i = 1,3
         n_real(i) = n_real_arr(iwav)
         n_imag(i) = n_imag_arr(iwav)
      enddo
      
      call RTSMie_master_trimodal &
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         psd_indexes, psd_pars, monoradius, R1, R2, FixR1R2,   & ! I
         nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         wl, n_real, n_imag, modal_fractions, Do_Save_PSD,   & ! I
         TriMie_bulk, TriMie_asymm, TriMie_ncoeffs,          & ! O
         TriMie_expcoeffs, TriMie_Fmatrix, TriMie_dist,      & ! O
         fail, istatus, Trimessages, trace_3 )                 ! O
        
      if ( fail ) then
        write(*,'(/a,i3)')' ** FAILURE! Trimodal test, wav # ',iwav
        write(*,'(a,i4)') '      - Status output (2 = fatal): ',istatus
        write(*,'(a,a)')  '      - First Message : ',adjustl(trim(Trimessages(1)))
        write(*,'(a,a)')  '      - First  trace  : ',adjustl(trim(Trimessages(2)))
        write(*,'(a,a)')  '      - Second trace  : ',adjustl(trim(Trimessages(3)))
        write(*,'(a,a)')  '      - Third  trace  : ',adjustl(trim(trace_3))
        stop
      endif

      !  Apply moment_cutoff

      apply = .true. ; L = 0
      do while ( apply )
         L = L + 1 ; 
         if ( abs(TriMie_expcoeffs(1,L)).lt.moment_cutoff ) apply = .false.
      enddo
      TriMie_ncoeffs = L - 1
      
      !  Save results
      
      ! Write bulk properties to file
      write(1,*) wl*1d3, TriMie_bulk(1), TriMie_bulk(3), TriMie_asymm
      
      ! Write phasefunction moments to file
      do L = 0, TriMie_ncoeffs
         write(2,*) wl*1d3, L, TriMie_expcoeffs(1:6,L)
      enddo
      
      ! write phasefunction to file
      do L = 1, n_Fmatrix_angles
         write(3,*) wl*1d3, Fmatrix_angles(L), TriMie_Fmatrix(1:4,L)
      enddo

   enddo

! Close result files
close(1)
close(2)
close(3)
close(4)

stop

end program RTSMie_Trimodal_test

