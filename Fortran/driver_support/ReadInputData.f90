module ReadInputData_m

public  :: ReadConfigFile
public  :: ReadArrayFromFile
public  :: ReadIntegerArrayFromFile
public  :: ReadVerticalGrid
public  :: ReadSpeciesNames
public  :: ReadGasXS
public  :: ReadAerosolMoments

integer, parameter, private :: fpk = Selected_Real_Kind(15)

contains


subroutine ReadConfigfile (configfile,                                     & ! Input filename for reading
       do_ss_only, do_vector_calculation, do_calculate_radiance,           & ! Main RT flags
       do_calculate_flux_at_observer, do_calculate_flux_profiles,          & ! Main RT flags
       do_rayleigh_only, do_deltam_scaling,                                & ! Main RT flags
       n_stokes, n_streams,n_fine_layers,                                  & ! Main RT flags
       planet_radius,                                                      & ! Planet properties
       observer_height,                                                    & ! Viewing geometry
       solar_zenith_angle, viewing_zenith_angle, relative_azimuth_angle,   & ! Viewing geometry
       do_radiance_jacobians, do_actinic_flux_jacobians,                   & ! Jacobian control
       backend_mode, sphericity_interpolation_points)                        ! Optional backend control


   implicit none

!  Inputs only
!  ===========

!  Input filename

   character*(*), intent(in) :: Configfile

!  Outputs (to be filled from file)
!  ================================

!  RT options

   logical, intent(out) :: do_ss_only
   logical, intent(out) :: do_Vector_Calculation
   logical, intent(out) :: do_calculate_radiance
   logical, intent(out) :: do_calculate_flux_at_observer
   logical, intent(out) :: do_calculate_flux_profiles
   
   logical, intent(out) :: do_rayleigh_only
   logical, intent(out) :: do_deltam_scaling
   logical, intent(out) :: do_radiance_jacobians
   logical, intent(out) :: do_actinic_flux_jacobians

   integer, intent(out) :: n_streams
   integer, intent(out) :: n_stokes
   integer, intent(out) :: n_fine_layers
   integer, optional, intent(out) :: sphericity_interpolation_points

   real(kind=fpk), intent(out) :: planet_radius
   
   real          , intent(out) :: observer_height
   real(kind=fpk), intent(out) :: solar_zenith_angle
   real(kind=fpk), intent(out) :: viewing_zenith_angle
   real(kind=fpk), intent(out) :: relative_azimuth_angle
   character(len=*), optional, intent(out) :: backend_mode


!  Initialize everything
!  =====================

   do_ss_only         = .false.
   do_vector_Calculation = .false.
   do_calculate_radiance = .false.
   do_calculate_flux_at_observer   = .false.
   do_calculate_flux_profiles = .false.
   
   do_deltam_scaling     = .false.

   n_streams = 0
   n_stokes  = 0
	n_fine_layers = 0

   observer_height     = 0
   Planet_radius    = 0.0_fpk

   solar_zenith_angle = 0.0_fpk
   if (present(backend_mode)) backend_mode = 'vlidort_qs'
   if (present(sphericity_interpolation_points)) sphericity_interpolation_points = 3

!  Read the configuration file
!  ===========================

!  Open file

   Open(1, file=Trim(Configfile), form = 'formatted', status='old' )

!  Read Radiative transfer specifications

   read(1,*) do_ss_only
   read(1,*) do_vector_calculation

   if (do_vector_Calculation) then 
      n_stokes = 3
   else 
      n_stokes = 1
   endif

   read(1,*) do_calculate_radiance
   read(1,*) do_calculate_flux_at_observer
   read(1,*) do_calculate_flux_profiles
   
   read(1,*) do_rayleigh_only
   read(1,*) do_deltam_scaling
   read(1,*) n_streams
	read(1,*) n_fine_layers

   if ( do_ss_only ) then
      do_deltam_scaling = .false.
   endif

   read(1,*) observer_height
   read(1,*) planet_radius

   read(1,*) solar_zenith_angle
   read(1,*) viewing_zenith_angle
   read(1,*) relative_azimuth_angle

   read(1,*) do_radiance_jacobians
   read(1,*) do_actinic_flux_jacobians
   if (present(backend_mode)) read(1,*) backend_mode
   if (present(sphericity_interpolation_points)) read(1,*) sphericity_interpolation_points


   close(1)
   return

end subroutine ReadConfigfile



subroutine ReadArrayFromFile (filename,                 &
     MAXROWS, MAXCOLS,                             & ! Input     
     n_rows, n_cols, array)             ! Outputs

   implicit none

! Input
   character*(*), intent(in) :: filename
   integer, intent(in) :: MAXROWS
   integer, intent(in) :: MAXCOLS

! Output
   integer, intent(out) :: n_rows
   integer, intent(out) :: n_cols
   real(kind=fpk), intent(out) :: array (MAXROWS, MAXCOLS )

! Local variables
   integer        :: n

! Functionality

   array   = 0.0_fpk

   open(18,file=TRIM(filename), status='old')
   read (18,*) n_rows
   read (18,*) n_cols
   read (18,*) ! Ignore column header
   do n = 1, n_rows
      read(18,*) array(n,1:n_cols)
   enddo
   close(18)
   
   return
end subroutine ReadArrayFromFile



subroutine ReadIntegerArrayFromFile (filename,                 &
     MAXROWS, MAXCOLS,                             & ! Input     
     n_rows, n_cols, array)             ! Outputs

   implicit none

! Input
   character*(*), intent(in) :: filename
   integer, intent(in) :: MAXROWS
   integer, intent(in) :: MAXCOLS

! Output
   integer, intent(out) :: n_rows
   integer, intent(out) :: n_cols
   integer, intent(out) :: array (MAXROWS, MAXCOLS )

! Local variables
   integer        :: n

! Functionality

   array = 0

   open(18,file=TRIM(filename), status='old')
   read (18,*) n_rows
   read (18,*) n_cols
   read (18,*) ! Ignore column header
   do n = 1, n_rows
      read(18,*) array(n,1:n_cols)
   enddo
   close(18)

   return
end subroutine ReadIntegerArrayFromFile



subroutine ReadLogicalArrayFromFile (filename,                 &
     MAXROWS, MAXCOLS,                             & ! Input     
     n_rows, n_cols, array)             ! Outputs

   implicit none

! Input
   character*(*), intent(in) :: filename
   integer, intent(in) :: MAXROWS
   integer, intent(in) :: MAXCOLS

! Output
   integer, intent(out) :: n_rows
   integer, intent(out) :: n_cols
   logical, intent(out) :: array (MAXROWS, MAXCOLS )

! Local variables
   integer        :: n

! Functionality

   array   = .false.

   open(18,file=TRIM(filename), status='old')
   read (18,*) n_rows
	read (18,*) n_cols
   read (18,*) ! Ignore column header
   do n = 1, n_rows
      read(18,*) array(n,1:n_cols)
   enddo
   close(18)

   return
end subroutine ReadLogicalArrayFromFile


subroutine ReadVerticalGrid (grid_file,                 &
     MAXLAYERS,                                         & ! Input     
     n_layers, layer_boundary_altitudes)             ! Outputs

   implicit none

! Input

   character*(*), intent(in) :: grid_file
   integer, intent(in) :: MAXLAYERS

! Output

   integer, intent(out) :: n_layers
   real(kind=fpk), intent(out) :: layer_boundary_altitudes (0:MAXLAYERS)

! Local variables

   integer        :: n

! Function

   layer_boundary_altitudes   = 0.0_fpk

   open(18,file=TRIM(grid_file), status='old')
   read (18,*) n_layers
   read (18,*) ! Ignore header
   do n = 0, n_layers
      read(18,*) layer_boundary_altitudes(n)
   enddo
   close(18)

   return
   
end subroutine ReadVerticalGrid


subroutine ReadSpeciesNames(filename,            &
   MAX_SPECIES, MAX_SPECIESNAME_LENGTH,               & ! input
   species_names, n_species)                            ! output
   
   implicit none

! Input
   character*(*), intent(in) :: filename
   integer, intent(in) :: MAX_SPECIES
   integer, intent(in) :: MAX_SPECIESNAME_LENGTH

! Output
   integer, intent(out)   :: n_species
   character(len=MAX_SPECIESNAME_LENGTH), dimension(MAX_SPECIES), intent(out) :: species_names

! Local variables
   integer        :: n

   open(18,file=TRIM(filename), status='old')
   read (18,*) n_species
   do n = 1, n_species
      read(18,*) species_names(n)
   enddo
   close(18)

   return

end subroutine ReadSpeciesNames



subroutine ReadGasXS(pathname,                  &
   MAX_SPECIES, MAX_WAVELENGTHS, MAX_GASXS_TEMPCOEFFS, MAX_SPECIESNAME_LENGTH, &
   n_species, n_wavelengths, species_names,                 &
   gas_xs_data)
   
   implicit none

! Input
   character*(*), intent(in) :: pathname
   integer, intent(in) :: MAX_SPECIES, MAX_WAVELENGTHS, MAX_GASXS_TEMPCOEFFS, MAX_SPECIESNAME_LENGTH
   integer, intent(in) :: n_species
   character(len=*), dimension(*), intent(in) :: species_names

! Output
   real(fpk), intent(out) :: gas_xs_data(MAX_SPECIES,MAX_WAVELENGTHS,0:MAX_GASXS_TEMPCOEFFS)
   integer, intent(out)   :: n_wavelengths
   
! Local variables
   integer             :: w
   integer             :: g
   integer             :: i
   integer             :: n_columns
   character(len=MAX_SPECIESNAME_LENGTH) :: filename
   
   do g = 1,n_species ! Iterate over gases
      filename = trim(adjustl(pathname)) // trim(adjustl(species_names(g))) // ".inp"
      open(18,file=trim(filename), status='old')
      read (18,*) n_wavelengths
      read (18,*) n_columns
      read (18,*) ! Ignore column header
      do w = 1, n_wavelengths
         read(18,*) gas_xs_data(g,w,0:n_columns-1)
      enddo
      close(18)
   enddo 

   return

end subroutine ReadGasXS



subroutine ReadAerosolMoments(pathname,                            &
   MAX_AEROSOL_TYPES, MAX_WAVELENGTHS, MAXMOMENTS_INPUT, MAX_SPECIESNAME_LENGTH,  & ! Input
   n_aertypes, n_wavelengths, aertype_names,                 & ! Input
   aer_moments)                                                ! Output
   
   implicit none

! Input
   character*(*), intent(in) :: pathname
   integer, intent(in) :: MAX_AEROSOL_TYPES, MAX_WAVELENGTHS, MAXMOMENTS_INPUT, MAX_SPECIESNAME_LENGTH
   integer, intent(in) :: n_aertypes
   character(len=*), dimension(*), intent(in) :: aertype_names

! Output
   real(fpk), intent(out) :: aer_moments(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT,6)
   integer, intent(out)   :: n_wavelengths
   
! Local variables
   integer             :: a ! iterates over aerosol types
   integer             :: l ! iterates over all (all wavelengths) moments of current aerosol type
   integer             :: n_all_moments
   integer             :: w ! buffer for wavelength index as read from file
   integer             :: m ! buffer for moment index (current wavelength) as read from file
   character(len=MAX_SPECIESNAME_LENGTH) :: filename
   real(fpk)           :: mom_values(6)
   
   aer_moments = 0.0 ! Initialize with zero
   
   do a = 1,n_aertypes ! Iterate over aerosol types
      filename = trim(adjustl(pathname)) // trim(adjustl(aertype_names(a))) // ".inp"
      open(18,file=trim(filename), status='old')
      read (18,*) n_all_moments
      read (18,*) ! Ignore column header
      do l = 1, n_all_moments ! iterate over all lines of the file
         read(18,*) w, m, mom_values
         aer_moments(w, a, m, :) = mom_values
      enddo
      close(18)
   enddo 

   return

end subroutine ReadAerosolMoments


subroutine ReadAerosolMomentJacobians(filename,                            &
   MAX_WAVELENGTHS, MAXMOMENTS_INPUT,                        & ! Input
   aer_moment_jacobians)                                       ! Output
   
   implicit none

! Input
   character*(*), intent(in) :: filename
   integer, intent(in) :: MAX_WAVELENGTHS, MAXMOMENTS_INPUT

! Output
   real(fpk), intent(out) :: aer_moment_jacobians(MAX_WAVELENGTHS, 0:MAXMOMENTS_INPUT,6)
   
! Local variables
   integer             :: n_wavelengths
   integer             :: a ! iterates over aerosol types
   integer             :: l ! iterates over all (all wavelengths) moments of current aerosol type
   integer             :: n_all_moments
   integer             :: w ! buffer for wavelength index as read from file
   integer             :: m ! buffer for moment index (current wavelength) as read from file
   real(fpk)           :: mom_values(6)
   
   aer_moment_jacobians = 0.0 ! Initialize with zero
   
   open(18,file=trim(filename), status='old')
   read (18,*) n_all_moments
   read (18,*) ! Ignore column header
   do l = 1, n_all_moments ! iterate over all lines of the file
      read(18,*) w, m, mom_values
      aer_moment_jacobians(w, m, :) = mom_values
   enddo
   close(18)

   return

end subroutine ReadAerosolMomentJacobians


end module ReadInputData_m
