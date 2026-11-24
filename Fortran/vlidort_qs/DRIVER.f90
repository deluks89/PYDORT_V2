
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

program DRIVER

   use MSSeg_PARS_m, only : MAX_MESSAGES
   use FOMS_pars_m
   use DRIVER_pars_m

   use MainDriverUtility_m
   use ReadInputData_m
   use RayleighScattering_m

   use VQS_IO_DEFS_m
   use VQS_LIN_IO_DEFS_m
   use VQS_MASTERS_m
   use VQS_LC_MASTERS_m
   use VQS_LP_MASTERS_m
   use VQS_AUX_m

   implicit none


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  VLIDORT-QS I/O
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VQS input structures

   TYPE(VQS_Fixed_Inputs)    :: VQS_FixIn
   TYPE(VQS_Modified_Inputs) :: VQS_ModIn
   TYPE(VQS_Fixed_LinInputs) :: VQS_LinFixIn

!  VQS output structure

   TYPE(VQS_Outputs)         :: VQS_Out
   TYPE(VQS_LinOutputs)      :: VQS_LinOut

!  Stokes vector elements

   real(kind=ffp), dimension(1) :: Ivalues
   real(kind=ffp), dimension(1) :: Qvalues
   real(kind=ffp), dimension(1) :: Uvalues

!  Exception handling for Tool

   logical            :: fail
   integer            :: nmessages
   character(len=200) :: messages (maxmessages)
   logical            :: openfileflag = .false.
   logical            :: do_debug_master_input = .false.


! PARAMETERS:

   real(kind=fpk), parameter :: k_b = 1.3806488d-23    ! Boltzmann constant
   real(kind=fpk), parameter :: pi  = 3.14159265359    
   
   character (len = *), parameter :: output_path = "OUTPUT/"

!  Mask for correct sorting of Mie outputs into greek matrix

   integer        :: GMASK(8), CMASK(8), SMASK(8)
   integer        :: gk,ck,sk ! for iteration

! Data read from configuration and data input files

   Logical        :: do_ss_only
   Logical        :: do_vector_calculation
   Logical        :: do_calculate_radiance
   Logical        :: do_calculate_flux_at_observer
   Logical        :: do_calculate_flux_profiles
   Logical        :: do_deltam_scaling
   Logical        :: do_rayleigh_only
   real(kind=ffp) :: planet_radius
   integer        :: n_stokes, n_streams
   real(kind=ffp) :: solar_zenith_angle, viewing_zenith_angle, relative_azimuth_angle
   Logical        :: do_radiance_jacobians
   Logical        :: do_actinic_flux_jacobians
   logical :: do_irradiance, do_irradiance_jacobians
   logical :: do_irradiance_orders
   integer :: flux_pass, n_order_values, order_index
   real(fpk) :: scattering_step, scatter_scale, scatter_albedo, order_direct
   real(fpk) :: ground_up, ground_up_jac(MAXLAYERS)
   real(fpk) :: scatter_omega(MAXLAYERS)
   real(fpk) :: order_samples(2+2*MAX_ATMOSWFS*MAXLAYERS,4)
   real(fpk) :: order_diffuse(2+2*MAX_ATMOSWFS*MAXLAYERS)

   integer        :: n_fine_layers

   integer        :: n_nonzero_greekmat_entries

   integer        :: n_wavelengths
   real(kind=ffp) :: simulated_wavelengths (1:MAX_WAVELENGTHS)  ! [nm]

   integer        :: n_layers
   real(fpk)      :: layer_boundary_altitudes(0:MAXLAYERS)  ! [m]
   real(fpk)      :: layer_thicknesses(MAXLAYERS)  ! [m]
   
   integer        :: n_aertypes
   character(len=MAX_SPECIESNAME_LENGTH), dimension(MAX_AEROSOL_TYPES) :: aertype_names
   real(fpk)      :: aer_extinction_profiles(MAXLAYERS, MAX_AEROSOL_TYPES) ! [km-1]
   real(fpk)      :: aer_extinction_multipliers(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES)
   real(fpk)      :: aer_ssas(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES)
   real(fpk)      :: aer_moments(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT, 6)

   ! For reading and storage of aerosol model jacobians:
   real(fpk)      :: d_ssa_d_param

   character(len=MAX_SPECIESNAME_LENGTH), dimension(MAX_GASES) :: gas_names
   integer        :: n_gases
   real(fpk)      :: gas_concentration_profiles(MAXLAYERS, MAX_GASES) ! [molec/cm3]
   real(fpk)      :: gas_xs_data(MAX_GASES, MAX_WAVELENGTHS, 0:MAX_GASXS_TEMPCOEFFS)  ! [molec/cm2]
   real(fpk)      :: surface_albedos(MAX_WAVELENGTHS)

   real(fpk)      :: pressure_profile(MAXLAYERS) ! [hPa]
   real(fpk)      :: temperature_profile(MAXLAYERS) ! [K]
   real(fpk)      :: vapor_pressure_profile(MAXLAYERS) ! [hPa]

   integer        :: dummy_integer
	logical        :: dummy_logical

   real           :: observer_height ! [m]

   character(len=MAX_SPECIESNAME_LENGTH), dimension(MAX_ATMOSWFS) :: jacobian_list
   integer        :: n_jacobians
   logical        :: jacobian_layer_flags(MAXLAYERS)
   integer        :: aerosol_model_jacobians_aertype_indexes(MAX_ATMOSWFS)
   real(fpk)      :: aerosol_model_jacobians_extmul(MAX_ATMOSWFS,MAX_WAVELENGTHS)
   real(fpk)      :: aerosol_model_jacobians_ssa(MAX_ATMOSWFS,MAX_WAVELENGTHS)
   real(fpk)      :: aerosol_model_jacobians_moms(MAX_ATMOSWFS,MAX_WAVELENGTHS,0:MAXMOMENTS_INPUT,6)

! UPDATED FOR EACH WAVELENGTH:
   
   real(fpk)      :: wavelength
   real(fpk)      :: rayleigh_greekcoeffs(0:2,6)
   INTEGER        :: n_moments_input
   
! EVEN MORE LOCAL:

   real(fpk)      :: rayleigh_xs ! for current layer and wavelength [cm2]
   real(fpk)      :: air_density ! for current layer [m-3]
   real(fpk)      :: rayleigh_layer_od ! for current layer
   
   ! Local buffer for aerosol model jacobians:
   real(fpk)      :: d_extmul_d_p
   real(fpk)      :: d_ssa_d_p
   real(fpk)      :: d_moms_d_p(0:MAXMOMENTS_INPUT,16)
   
   real(fpk)      :: aer_layer_extinction_ods(MAX_AEROSOL_TYPES)
   real(fpk)      :: aer_layer_scattering_ods(MAX_AEROSOL_TYPES)
   real(fpk)      :: aer_greekcoeffs(MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT, 6) ! Updated for each wavelength
   real(fpk)      :: aer_greekcoeffs_matrix(MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT, 16) ! Will contain same information as aer_greekcoeffs but in the 16-entry matrix form required by VLIDORT master routine. Required for linearisations.

   ! variables for aggregated aerosol
   real(fpk)      :: aer_agg_layer_extinction_od
   real(fpk)      :: aer_agg_layer_scattering_od
   real(fpk)      :: aer_agg_layer_greekcoeffs(0:MAXMOMENTS_INPUT, 6)
   real(fpk)      :: gas_agg_layer_od

   ! Total optical properties of currently iterated layer:
   real(fpk)      :: layer_thickness
   real(fpk)      :: layer_greekmat(0:MAXMOMENTS_INPUT, 16) ! greekmat for layer (rayleigh + aerosol)
   real(fpk)      :: layer_extinction_od
   real(fpk)      :: layer_scattering_od
   real(fpk)      :: layer_absorption_od
   real(fpk)      :: layer_ssa


! VERY VERY LOCAL:

   real(fpk)         :: xs, oh, vna, sza, h, mom
   integer           :: n, l, g, w, i, a, q
   character(MAX_SPECIESNAME_LENGTH) :: jstr ! for jacobian labels
   character(MAX_SPECIESNAME_LENGTH) :: astr ! for aerosol type labels
   character(MAX_SPECIESNAME_LENGTH) :: pstr ! for aerosol type labels
   integer           :: idx

!  ##################################################################
!      S T A R T    of    C O D E
!  ##################################################################

   GMASK = (/  1,  2,  5, 6, 11, 12, 15, 16 /)     ! Greek Matrix indices
   CMASK = (/  1,  5,  5, 2,  3,  6,  6,  4 /)     ! Associated coefficients from Mie model output
   SMASK = (/  1, -1, -1, 1,  1, -1,  1,  1 /)     ! Signs to apply to Mie model output


!  Initialize output and exception handling
!  ----------------------------------------

!  Initialize (important)

   fail = .false.
   nmessages = 0
   messages = ' '
   OPENFILEFLAG = .false.



!  Read configuration file
!  -----------------------
   
   Call ReadConfigfile ('INPUT/configuration.inp',                         & ! Input filename for reading
       do_ss_only, do_vector_calculation, do_calculate_radiance,           & ! Main RT flags
       do_calculate_flux_at_observer, do_calculate_flux_profiles,          & ! Main RT flags
       do_rayleigh_only, do_deltam_scaling,                                & ! Main RT flags
       n_stokes, n_streams,n_fine_layers,                                  & ! Main RT flags
       planet_radius,                                                      & ! Planet properties
       observer_height,                                                    & ! Viewing geometry
       solar_zenith_angle, viewing_zenith_angle, relative_azimuth_angle,   & ! Viewing geometry
       do_radiance_jacobians, do_actinic_flux_jacobians, &
       do_irradiance=do_irradiance, do_irradiance_jacobians=do_irradiance_jacobians, &
       do_irradiance_orders=do_irradiance_orders, scattering_step=scattering_step)
   if (do_irradiance_jacobians .and. .not. do_irradiance) &
      error stop 'Enable irradiance when requesting its Jacobians.'
   ! Both angular moments share the native up/down flux solves.
   do_calculate_flux_at_observer = do_calculate_flux_at_observer .or. do_irradiance
   do_actinic_flux_jacobians = do_actinic_flux_jacobians .or. do_irradiance_jacobians
   if (do_irradiance_orders) then
      if (.not. do_irradiance .or. do_ss_only) error stop 'Orders require full-scattering irradiance.'
      if (do_deltam_scaling .and. .not. do_rayleigh_only) error stop 'Disable delta-M for scattering orders.'
      if (scattering_step < 1d-6 .or. scattering_step > 1d-2) error stop 'Invalid scattering-order step.'
   endif


!  Entries in Greekmat

   if ( n_stokes .eq. 1 ) then 
      n_nonzero_greekmat_entries = 1
   else if ( n_stokes .eq. 2 ) then 
      n_nonzero_greekmat_entries = 4
   else if ( n_stokes .eq. 3 ) then 
      n_nonzero_greekmat_entries = 5
   else if ( n_stokes .eq. 4 ) then 
      n_nonzero_greekmat_entries = 8
   endif
   

   ! Read wavelength grid:
   Call ReadArrayFromFile ('INPUT/wavelengths.inp',        &
      MAX_WAVELENGTHS, 1,                                            & ! Inputs
      n_wavelengths, dummy_integer, simulated_wavelengths)            ! Outputs
   !write(*,*) simulated_wavelengths(1:n_wavelengths)

   ! Read vertical grid:
   Call ReadVerticalGrid ('INPUT/vertical_grid.inp',    &
     MAXLAYERS,                                      & ! Input     
     n_layers, layer_boundary_altitudes) ! Outputs
   ! write(*,*) layer_boundary_altitudes(0:n_layers)
   call external_observer_to_toa(planet_radius, layer_boundary_altitudes(0), observer_height, &
      solar_zenith_angle, viewing_zenith_angle, relative_azimuth_angle, do_calculate_flux_at_observer)
   layer_thicknesses = layer_boundary_altitudes(0:n_layers) - layer_boundary_altitudes(1:)
	if (n_layers .gt. MAXLAYERS) then
		write(*,*) 'Maximum number of layers defined by MSEG_pars.f90 exceeded!'
		stop
	endif

   ! Read pressure profile
   Call ReadArrayFromFile('INPUT/pressure_profile.inp',    &
      MAXLAYERS, 1,                                           &
      n_layers, dummy_integer, pressure_profile)
   !write(*,*) pressure_profile(1:n_layers)
   
   ! Read temperature profile
   Call ReadArrayFromFile('INPUT/temperature_profile.inp',    &
      MAXLAYERS, 1,                                           &
      n_layers, dummy_integer, temperature_profile)
   ! write(*,*) temperature_profile(1:n_layers)
   
   ! Read water vapor partial pressure profile
   Call ReadArrayFromFile('INPUT/vapor_pressure_profile.inp',    &
      MAXLAYERS, 1,                                           &
      n_layers, dummy_integer, vapor_pressure_profile)
   ! write(*,*) vapor_pressure_profile(1:n_layers)

   ! Read names of considered aerosols
   Call ReadSpeciesNames('INPUT/aertype_names.inp',        &
      MAX_AEROSOL_TYPES, MAX_SPECIESNAME_LENGTH,               & ! input
      aertype_names, n_aertypes)                            ! output
   !do n = 1,n_aertypes
   !   write(*,*) trim(aertype_names(n))
   !enddo
	if (n_aertypes .gt. MAX_AEROSOL_TYPES) then
		write(*,*) 'Maximum number of aerosol types defined by DRIVER_pars.f90 exceeded!'
		stop
	endif

   ! Read aerosol extinction profiles for all aerosol types:
   Call ReadArrayFromFile ('INPUT/aerosol_extinction_profiles.inp',    &
        MAXLAYERS, MAX_AEROSOL_TYPES,                                  & ! Input     
        n_layers, n_aertypes, aer_extinction_profiles)                 ! Outputs
   !do n = 1,n_layers
   !   write(*,*) aer_extinction_profiles(n,:)
   !enddo

   ! Read aerosol extinction multipliers for all aerosol types:
   Call ReadArrayFromFile ('INPUT/aerosol_extinction_multipliers.inp',    &
        MAX_WAVELENGTHS, MAX_AEROSOL_TYPES,                                  & ! Input     
        n_wavelengths, n_aertypes, aer_extinction_multipliers)                 ! Outputs
   !do n = 1,n_wavelengths
   !   write(*,*) aer_extinction_multipliers(n,:)
   !enddo
   
   ! Read aerosol extinction multipliers for all aerosol types:
   Call ReadArrayFromFile ('INPUT/aerosol_single_scattering_albedos.inp',    &
        MAX_WAVELENGTHS, MAX_AEROSOL_TYPES,                                  & ! Input     
        n_wavelengths, n_aertypes, aer_ssas)                 ! Outputs
   !do n = 1,n_wavelengths
   !   write(*,*) aer_ssas(n,:)
   !enddo
   
   ! Read aerosol scattering moments:
   Call ReadAerosolMoments('INPUT/aerosol_moments/',                            &
      MAX_AEROSOL_TYPES, MAX_WAVELENGTHS, MAXMOMENTS_INPUT, MAX_SPECIESNAME_LENGTH,  & ! Input
      n_aertypes, n_wavelengths, aertype_names,                 & ! Input
      aer_moments)                                                ! Output
	! do l = 0,MAXMOMENTS_INPUT
   !    write (*,*) aer_moments(1,1,l,1) ! wavlength, aertype, moment_order, matrix_element
	! enddo
   
   ! Read names of considered trace gases:
   Call ReadSpeciesNames('INPUT/gas_names.inp',        &
      MAX_GASES, MAX_SPECIESNAME_LENGTH,               & ! input
      gas_names, n_gases)                              ! output
   ! do n = 1,n_gases
   !    write(*,*) trim(gas_names(n))
   ! enddo
	if (n_gases .gt. MAX_GASES) then
		write(*,*) 'Maximum number of gases defined by DRIVER_pars.f90 exceeded!'
		stop
	endif
   
	! Read names of considered aerosol types:
   Call ReadSpeciesNames('INPUT/aertype_names.inp',        &
      MAX_AEROSOL_TYPES, MAX_SPECIESNAME_LENGTH,               & ! input
      aertype_names, n_aertypes)                              ! output
   ! do n = 1,n_gases
   !    write(*,*) trim(gas_names(n))
   ! enddo
	
   ! Read trace gas absorption cross-sections including temperature coefficients:
   Call ReadGasXS('INPUT/xs/',                            &
      MAX_GASES, MAX_WAVELENGTHS, MAX_GASXS_TEMPCOEFFS, MAX_SPECIESNAME_LENGTH,  & ! Input
      n_gases, n_wavelengths, gas_names,                 & ! Input
      gas_xs_data)                                     ! Output
   ! do n = 1,n_wavelengths
   !   write(*,*) gas_xs_data(3,n,:)
   ! enddo

   ! Read trace gas concentration profiles for all gas species:
   Call ReadArrayFromFile ('INPUT/gas_concentration_profiles.inp',     &
      MAXLAYERS, MAX_GASES,                                               & ! Input     
      n_layers, n_gases, gas_concentration_profiles)                        ! Outputs
   ! do n = 1,n_layers
   !   write(*,*) gas_concentration_profiles(n,:)
   ! enddo
   
   ! Read lambertian surface albedos:
   Call ReadArrayFromFile ('INPUT/surface_albedos.inp',        &
      MAX_WAVELENGTHS, 1,                                            & ! Inputs
      n_wavelengths, dummy_integer, surface_albedos)           ! Outputs
   !write(*,*) surface_albedos(1:n_wavelengths)
      
   ! Read list of requested profile jacobians:
   Call ReadSpeciesNames ('INPUT/jacobian_list.inp',        &
      MAX_ATMOSWFS, MAX_SPECIESNAME_LENGTH,                                            & ! Inputs
      jacobian_list, n_jacobians)           ! Outputs
	! do q = 1,n_jacobians
   !    write(*,*) trim(jacobian_list(q))
   ! enddo
	
   ! Read Jacobian layer flags:
   Call ReadLogicalArrayFromFile('INPUT/jacobian_layer_flags.inp', &
								 MAXLAYERS, 1, &
								 n_layers, dummy_integer, jacobian_layer_flags)
   ! write(*,*) jacobian_layer_flags(1:n_layers)	
   
   ! Load aerosol model jacobians
   aerosol_model_jacobians_aertype_indexes = 0
   aerosol_model_jacobians_extmul = ZERO
   aerosol_model_jacobians_ssa = ZERO
   aerosol_model_jacobians_moms = ZERO
   if (do_radiance_jacobians .or. do_actinic_flux_jacobians) then
      do q = 1,n_jacobians
         jstr = jacobian_list(q)
         if (index(trim(jstr), "_") .ne. 0 .and. find_string_in_array(trim(jstr), aertype_names) <= 0) then 
            ! Determine and save aerosol type index for the respective jacobian:
            idx = index(trim(jstr),'_', back=.true.)
            astr = jstr(1:idx-1)
            aerosol_model_jacobians_aertype_indexes(q) = find_string_in_array(trim(astr), aertype_names)
            ! Read d_extmul_d_param
            Call ReadArrayFromFile ( &
               'INPUT/aerosol_model_jacobians/' // trim(jstr) // '_extmul.inp', &
               MAX_WAVELENGTHS, 1,                           & ! Inputs
               n_wavelengths, dummy_integer, aerosol_model_jacobians_extmul(q,:)) ! Outputs
            ! Read d_ssa_d_param
            Call ReadArrayFromFile ( &
               'INPUT/aerosol_model_jacobians/' // trim(jstr) // '_ssa.inp', &
               MAX_WAVELENGTHS, 1,                           & ! Inputs
               n_wavelengths, dummy_integer, aerosol_model_jacobians_ssa(q,:)) ! Outputs
            ! Read d_ssa_d_param
            Call ReadAerosolMomentJacobians( &
               'INPUT/aerosol_model_jacobians/' // trim(jstr) // '_moms.inp', &
               MAX_WAVELENGTHS, MAXMOMENTS_INPUT,                           & ! Inputs
               aerosol_model_jacobians_moms(q,:,:,:)) ! Outputs
         endif
      enddo
   endif

	
   VQS_FixIn%Bool%DO_SUNLIGHT           = .true.
   VQS_FixIn%Bool%DO_PLANE_PARALLEL     = .false.
   VQS_FixIn%Bool%DO_LAMBERTIAN_SURFACE = .true.
   VQS_FixIn%Cont%METHOD_INDEX          = 1     ! 1 = Fast (Mk3 implementation), 2 = Slow (Mk2 implementation)
   VQS_FixIn%Cont%DOCRIT                = .true.
   VQS_FixIn%Cont%ACRIT                 = 1.0d-15
   VQS_FixIn%Sunrays%FLUX_FACTOR        = 1d0
   VQS_FixIn%UserVal%GEO_FILENAME       = 'debug/geometry_debug.log'

   VQS_FixIn%Bool%DO_DEBUG_FO_GEOMETRY  = .false.
   VQS_FixIn%Bool%DO_DEBUG_FO_INPUT     = .false.
   VQS_FixIn%Bool%DO_DEBUG_MS_INPUT     = .false.
   VQS_FixIn%Bool%DO_RT_PROGRESS        = .false.


   VQS_FixIn%Bool%DO_MSRT               = .not. do_ss_only
   
   VQS_FixIn%Geophys%PLANET_RADIUS = planet_radius
   VQS_FixIn%Geophys%HEIGHT_GRID   = layer_boundary_altitudes / 1d3 ! vertical grid [km]

   VQS_FixIn%Cont%nlayers      = n_layers
   VQS_FixIn%Cont%nstokes      = n_stokes
   VQS_FixIn%Cont%nstreams     = n_streams
   VQS_FixIn%Cont%NFINE_INPUT  = n_fine_layers

   VQS_FixIn%UserVal%RELAZ_CAMERA  = relative_azimuth_angle

   ! Preapre viewing nadir angle. But assingment to VLIDORT happens later because sign has to be changed to obtain up-/downwelling flux.
   vna = 180d0 - viewing_zenith_angle ! VLIDORT requires viewing nadir angle
   if (vna < 1d-5) vna = 1d-5 ! 
   if (vna > 180d0 - 1d-5) vna = 180d0 - 1d-5

   VQS_FixIn%Optical%ATMOS_INDEX          = 0

   ! Similar for solar zenith angle
   sza = solar_zenith_angle
   if (sza < 1d-5) sza = 1d-5 ! 
   if (sza > 90d0 - 1d-5) sza = 90d0 - 1d-5
   VQS_FixIn%UserVal%THETA_CAMERA  = sza

   ! Open result files
   ! Open all files despite some will not be filled with data (depending on VQS flags) to ensure that old files are replaced.

   ! Files actinic fluxes at observer altitude
   OPEN(10, FILE=trim(output_path) // 'flux_diffuse_upwelling.out', STATUS='replace')
   OPEN(11, FILE=trim(output_path) // 'flux_diffuse_downwelling.out', STATUS='replace')
   OPEN(12, FILE=trim(output_path) // 'flux_direct.out', STATUS='replace')
   ! Files actinic flux profiles
   OPEN(13, FILE=trim(output_path) // 'flux_profile_diffuse_upwelling.out', STATUS='replace')
   OPEN(14, FILE=trim(output_path) // 'flux_profile_diffuse_downwelling.out', STATUS='replace')
   OPEN(15, FILE=trim(output_path) // 'flux_profile_direct.out', STATUS='replace')

   OPEN(20, FILE=trim(output_path) // 'radiances.out', STATUS='replace')
   if (do_irradiance) call open_irradiance_files(jacobian_list(1:n_jacobians), &
      do_irradiance_jacobians, do_irradiance_orders)
   
   do q = 1,n_jacobians
		if (do_radiance_jacobians) then
			OPEN(30+q, FILE=trim(output_path) & 
				// 'jacs_radiance_' // trim(jacobian_list(q)) // '.out', STATUS='replace')
		endif
		if (do_actinic_flux_jacobians) then
			OPEN(50+q, FILE=trim(output_path) & 
				// 'jacs_flux_diffuse_upwelling_' // trim(jacobian_list(q)) // '.out', STATUS='replace')
			OPEN(70+q, FILE=trim(output_path) & 
				// 'jacs_flux_diffuse_downwelling_' // trim(jacobian_list(q)) // '.out', STATUS='replace')
			OPEN(90+q, FILE=trim(output_path) & 
				// 'jacs_flux_direct_' // trim(jacobian_list(q)) // '.out', STATUS='replace')
		endif
	enddo
   
!  Timing inputs

   VQS_ModIn%MTiming%GEO_TIME = 0.0
   VQS_ModIn%MTiming%SS_TIME  = 0.0
   VQS_ModIn%MTiming%MS_TIME  = 0.0


   do w = 1, n_wavelengths
   
      wavelength = simulated_wavelengths(w)
      
      write (*,*)
      write (*,*) "Simulating wavelength ", w ," of ",n_wavelengths, &
                  ": ", wavelength, " nm"

      ! Find number of phase function moments to actually consider:
      n_moments_input = 0
      do a = 1,n_aertypes
         do l = 0,MAXMOMENTS_INPUT
            do i = 1,6
               if (aer_moments(w,a,l,i) .gt. 0d0) then
                  if (l>n_moments_input) n_moments_input = l
               endif
            enddo
         enddo
      enddo
      

      VQS_ModIn%MCont%NGREEK_MOMENTS_INPUT = n_moments_input


      rayleigh_greekcoeffs = RayleighGreekCoeffs(wavelength) ! Rayleigh greekmat is independent of altitude
      aer_greekcoeffs = aer_moments(w,:,:,:) ! no zero initilaisation required here. Invalid moments are zeroed from read-in routine
		
		! Sort aer_greekcoeffs into aer_greekcoeffs_matrix (in 16-entry matrix format):
		aer_greekcoeffs_matrix = ZERO
		
		do i = 1, n_nonzero_greekmat_entries
			gk = gmask(i); ck = cmask(i); sk = smask(i) ! position in the 16 element greek matrix; index of coefficient in mie model output; sign to apply for aerosol coefficients
			aer_greekcoeffs_matrix(1:n_aertypes, 0:n_moments_input, gk) & 
					= aer_greekcoeffs(1:n_aertypes, 0:n_moments_input, ck) * sk
		enddo
      
      ! Initialize surface albedo:
      VQS_FixIn%Optical%LAMBERTIAN_ALBEDO    = surface_albedos(w) !defined later

      ! Initialize layer properties:
      VQS_FixIn%Optical%DELTAU_VERT_INPUT = ZERO
      VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT = ZERO
      VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT = ZERO
      
      ! Initialize layer linearizarions:
	   VQS_LinFixIn%Cont%layer_vary_flag = jacobian_layer_flags
		VQS_LinFixIn%Cont%n_profile_wfs = n_jacobians
      VQS_LinFixIn%Cont%n_column_wfs    = ZERO
      VQS_LinFixIn%Optical%L_DELTAU_VERT_INPUT    = ZERO
      VQS_LinFixIn%Optical%L_OMEGA_TOTAL_INPUT    = ZERO
      VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT = ZERO

      do n = 1, n_layers
			
			layer_thickness = layer_thicknesses(n) ! in [m]
			
         ! Layer specific rayleigh quantities
         rayleigh_xs = RayleighXS(wavelength, &
                                  pressure_profile(n), &
                                  temperature_profile(n), &
                                  vapor_pressure_profile(n))
         air_density = (pressure_profile(n) * 1d2)/(k_b * temperature_profile(n)) ! [molec / m3]
         rayleigh_layer_od = rayleigh_xs * 1d-4 &  ! Convert to [m2]
                             * air_density * layer_thickness
			  
         ! Layer specific aerosol quantities:
         aer_layer_extinction_ods = 0d0
         aer_layer_extinction_ods(1:n_aertypes) = aer_extinction_profiles(n,1:n_aertypes) &
                                       * aer_extinction_multipliers(w,1:n_aertypes) &
                                       * layer_thickness / 1d3
         aer_agg_layer_extinction_od = sum(aer_layer_extinction_ods(1:n_aertypes))
			
			aer_layer_scattering_ods = 0d0
         aer_layer_scattering_ods(1:n_aertypes) = aer_layer_extinction_ods(1:n_aertypes) &
                                                  * aer_ssas(w,1:n_aertypes)
         aer_agg_layer_scattering_od = sum(aer_layer_scattering_ods(1:n_aertypes))
         
							
         aer_agg_layer_greekcoeffs = 0d0
         if (aer_agg_layer_scattering_od .gt. 0d0) then
            do a = 1,n_aertypes
               aer_agg_layer_greekcoeffs(0:n_moments_input,:) &
                  = aer_agg_layer_greekcoeffs(0:n_moments_input,:) &
                  + aer_layer_scattering_ods(a) * aer_greekcoeffs(a,0:n_moments_input,:)
            enddo
            aer_agg_layer_greekcoeffs = aer_agg_layer_greekcoeffs / aer_agg_layer_scattering_od
         endif
         
         ! Layer specific gas quantities:
         gas_agg_layer_od = 0d0
         do g = 1, n_gases
            xs = gas_xs_data(g,w,0) ! Read gas cross section value [cm2]
            do i = 1, MAX_GASXS_TEMPCOEFFS ! Apply temperature dependence
               xs = xs + gas_xs_data(g,w,i) * (temperature_profile(n)-273.15)**i
            enddo
            gas_agg_layer_od = gas_agg_layer_od & 
               + gas_concentration_profiles(n,g) * xs * layer_thickness * 1d2 ! Conversion to cm
         enddo
         
         ! Get total layer optical properties:
         layer_extinction_od = rayleigh_layer_od + aer_agg_layer_extinction_od + gas_agg_layer_od
         layer_scattering_od = rayleigh_layer_od + aer_agg_layer_scattering_od
         layer_absorption_od = layer_extinction_od - layer_scattering_od
         layer_ssa           = layer_scattering_od/layer_extinction_od
			
         ! Get total layer greekmat coefficients (same order as mie model output):
         layer_greekmat = 0d0
         do i = 1, n_nonzero_greekmat_entries
            gk = gmask(i); ck = cmask(i); sk = smask(i) ! position in the 16 element greek matrix; index of coefficient in mie model output; sign to apply for aerosol coefficients
            ! Insert weighted rayliegh greekmat coefficients first:
				layer_greekmat(0:2,gk) = rayleigh_layer_od * rayleigh_greekcoeffs(0:2,ck)
				! Add weighted aerosol greekmat coefficients: 
				layer_greekmat(0:n_moments_input,gk) = layer_greekmat(0:n_moments_input,gk) &
                  + aer_agg_layer_scattering_od &
                 * sk * aer_agg_layer_greekcoeffs(0:n_moments_input,ck)
         enddo
         layer_greekmat = layer_greekmat / (rayleigh_layer_od + aer_agg_layer_scattering_od) ! weight normalisation
   
         ! Correct minor deviations from one in first phase function moment. Seems to come from numerical issues in the Mie model itself.
         mom = layer_greekmat(0,1)
         if (mom .ne. 1d0) then
            if (abs(mom - 1d0) .lt. 1e-14) then
               layer_greekmat(:,:) = layer_greekmat(:,:) / mom
            endif
         endif
         
         
			VQS_FixIn%Optical%DELTAU_VERT_INPUT(n) = layer_extinction_od
         VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(n) = layer_ssa
         VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT(:,n,:) = layer_greekmat
			
      
         ! ####################
         ! LINEARISATION INPUTS
         ! ####################

		do q = 1, n_jacobians ! Index of Python Jacobian list (counts entries provided by Pydort)
            
         jstr = jacobian_list(q)
            
         if (trim(jstr) == 'absorption') then ! Layer absorption Jacobians (for bAMFs)
            VQS_LinFixIn%Optical%L_DELTAU_VERT_INPUT(q,n) &
               = 1 / layer_extinction_od
            VQS_LinFixIn%Optical%L_OMEGA_TOTAL_INPUT(q,n) &
               = - 1 / layer_extinction_od
            VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT(q,:,:,:) = ZERO
         endif
				
			do a = 1,n_aertypes
				if (trim(jstr) == trim(aertype_names(a))) then ! Layer aerosol extinction Jacobians
					VQS_LinFixIn%Optical%L_DELTAU_VERT_INPUT(q,n) &
						= aer_extinction_multipliers(w,a) &
						  * layer_thickness / 1d3 / layer_extinction_od
					VQS_LinFixIn%Optical%L_OMEGA_TOTAL_INPUT(q,n) &
						= aer_extinction_multipliers(w,a) &
						  * layer_thickness / 1d3 / layer_scattering_od & 
						  * (aer_ssas(w,a) - layer_ssa)
					do l = 0,n_moments_input
						do i = 1,16
							if (layer_greekmat(l,i) .ne. 0d0) then
								VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT(q,l,n,i) &
									= aer_extinction_multipliers(w,a) &
										* layer_thickness / 1d3 * aer_ssas(w,a) & 
										/ layer_greekmat(l,i) / layer_scattering_od &
										* (aer_greekcoeffs_matrix(a,l,i) - layer_greekmat(l,i))
							endif
						enddo
					enddo
				endif
			enddo

         ! Is the current weighting function q, an aerosol property weighting function?
         if (aerosol_model_jacobians_aertype_indexes(q) .gt. 0) then 
            a = aerosol_model_jacobians_aertype_indexes(q)
            
            ! Pick and arange aerosol model jacobians
            d_extmul_d_p = aerosol_model_jacobians_extmul(q,w)
            d_ssa_d_p = aerosol_model_jacobians_ssa(q,w)
            do i = 1, n_nonzero_greekmat_entries
               gk = gmask(i); ck = cmask(i); sk = smask(i) ! position in the 16 element greek matrix; index of coefficient in mie model output; sign to apply for aerosol coefficients
               d_moms_d_p(0:n_moments_input, gk) & 
                  = aerosol_model_jacobians_moms(q, w, 0:n_moments_input, ck) * sk
            enddo
            
            ! Helper variables
            h = layer_thickness/1d3
            
            VQS_LinFixIn%Optical%L_DELTAU_VERT_INPUT(q,n) &
               = h * aer_extinction_profiles(n,a) / layer_extinction_od &
               * d_extmul_d_p
            VQS_LinFixIn%Optical%L_OMEGA_TOTAL_INPUT(q,n) &
               = h * aer_extinction_profiles(n,a) / layer_scattering_od & 
                  * ( (aer_ssas(w,a) - layer_ssa) * d_extmul_d_p &
                     + aer_extinction_multipliers(w,a) * d_ssa_d_p )
            do l = 0,n_moments_input
               do i = 1,16
                  if (layer_greekmat(l,i) .ne. 0d0) then
                     VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT(q,l,n,i) &
                        = 1 / layer_greekmat(l,i) / layer_scattering_od &
                          * ( h * aer_extinction_profiles(n,a) &
                           * (aer_greekcoeffs_matrix(a,l,i) - layer_greekmat(l,i)) & 
                           * (aer_ssas(w,a) * d_extmul_d_p + aer_extinction_multipliers(w,a) * d_ssa_d_p) &
                           + aer_layer_scattering_ods(a) * d_moms_d_p(l,i))
                  endif
               enddo
            enddo
            
            ! VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT = ZERO
            ! VQS_LinFixIn%Optical%L_GREEKMAT_TOTAL_INPUT(q,1,n,1) &
            !             = 1 / layer_greekmat(1,1) / layer_scattering_od &
            !             * aer_layer_scattering_ods(a)
         endif
			

		enddo ! layer loop


      
	enddo ! wavelength loop

      ! ##############################
      ! Calculate radiance at observer
      ! ##############################
      
      if (observer_height .gt. layer_boundary_altitudes(0) - 1d-3) then 
         observer_height = layer_boundary_altitudes(0) - 1d-3
      elseif (observer_height .lt. 0d0) then
         observer_height = 0d0
      endif

      if (do_calculate_radiance) then

         write(*,*) '  Simulations for radiance:'

         VQS_FixIn%Bool%DO_CALCULATE_FLUXES   = .false.
         VQS_FixIn%UserVal%ALPHA_CAMERA  = vna
         VQS_ModIn%MBool%DO_RAYLEIGH_ONLY     = do_rayleigh_only
         VQS_ModIn%MBool%DO_DELTAM_SCALING    = do_deltam_scaling
         VQS_FixIn%UserVal%HEIGHT_CAMERA = observer_height / 1d3 ! Observer height in [km]
	      VQS_FixIn%Bool%DO_DEBUG_FLUXES       = .false. ! MUST BE TRUE FOR ACTINIC FLUX PROFILE CALCULATIONS!
         
         !do n = 1, n_layers
         !   write (*,*) VQS_FixIn%Geophys%HEIGHT_GRID(n), &
         !      VQS_FixIn%Optical%DELTAU_VERT_INPUT(n), &
         !      VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(n), &
         !      VQS_FixIn%Optical%GREEKMAT_TOTAL_INPUT(0:2,n,1)
         !enddo
         
         if (do_radiance_jacobians) then
            write (*,*) '    Calling linearized VQS master.'
            CALL VQS_LP_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out, &
               VQS_LinFixIn, VQS_LinOut)
            do q = 1, n_jacobians
               write(30+q, *) VQS_LinOut%Prof%LP_StokesTot(1,1:n_layers,q)
            enddo
         else
            write (*,*) '    Calling non-linearized VQS master.'
            CALL VQS_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out )
         endif
      
         CALL VQS_WRITE_STATUS ( &
           VQS_ERRUNIT, OPENFILEFLAG, VQS_Out%Status)

         ! write (20, '(4(2X ES22.15E2))') VQS_Out%Main%StokesTot(1)
         write (20, *) VQS_Out%Main%StokesTot(1)

         write (*,*) '      Radiance total:', VQS_Out%Main%StokesTot(1)
         write (*,*) '      Radiance first order:', VQS_Out%Main%StokesFO(1)
         write (*,*) '      Radiance multiple scattering:', VQS_Out%Main%StokesMS(1)

         
      endif


      ! ####################################
      ! Calculate actinic fluxes at observer
      ! ####################################
      
      if (do_calculate_flux_at_observer) then
         n_order_values = 2
         if (do_irradiance_jacobians) n_order_values = 2+2*n_jacobians*n_layers
         scatter_omega(1:n_layers) = VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(1:n_layers)
         scatter_albedo = VQS_FixIn%Optical%LAMBERTIAN_ALBEDO
         do flux_pass = 0,merge(4,0,do_irradiance_orders)
            if (flux_pass > 0) then
               scatter_scale = scattering_step * 2d0**(flux_pass-1)
               if (flux_pass == 4) scatter_scale = scattering_step/2d0
               VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(1:n_layers) = scatter_omega(1:n_layers)*scatter_scale
               VQS_FixIn%Optical%LAMBERTIAN_ALBEDO = scatter_albedo*scatter_scale
               ! For an unpolarized Sun, first-order Stokes I depends only on
               ! P11. Scalar probes avoid nearly degenerate weak-scattering
               ! vector eigenmodes; the baseline total remains polarized.
               VQS_FixIn%Cont%nstokes = 1
            endif

         if (flux_pass == 0) write(*,*) '  Simulations for actinic flux at observer:'

         ! UPWELLING (nadir view)
         VQS_FixIn%Bool%DO_CALCULATE_FLUXES = .true.
         VQS_FixIn%UserVal%ALPHA_CAMERA = 1d-5
         VQS_ModIn%MBool%DO_RAYLEIGH_ONLY     = do_rayleigh_only
         VQS_ModIn%MBool%DO_DELTAM_SCALING    = do_deltam_scaling
         VQS_FixIn%UserVal%HEIGHT_CAMERA = observer_height / 1d3 ! Observer height in [km]
		   VQS_FixIn%Bool%DO_DEBUG_FLUXES       = .true. ! MUST BE TRUE FOR ACTINIC FLUX PROFILE CALCULATIONS!

         if (do_actinic_flux_jacobians) then
            if (flux_pass == 0) write(*,*) '    Calling linearized VQS master for upwelling radiation.'
            CALL VQS_LP_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out, &
               VQS_LinFixIn, VQS_LinOut)
            do q = 1, n_jacobians
               if (flux_pass == 0 .and. observer_height > 0) &
                  write(50+q, *) VQS_LinOut%Prof%LP_ActFluxDiff(1,1:n_layers,q)*4d0*pi
               if (do_irradiance_jacobians .and. flux_pass == 0 .and. observer_height > 0) &
                  write(200+q,*) VQS_LinOut%Prof%LP_RegFluxDiff(1,1:n_layers,q)
            enddo
         else
            if (flux_pass == 0) write(*,*) '    Calling non-linearized VQS master for upwelling radiation.'
            CALL VQS_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out )
         endif
         
         ! write (10, '(2X ES22.15E2)') VQS_Out%Main%ActFluxDiff(1)*4d0*pi
         if (flux_pass == 0 .and. observer_height > 0) write(10, *) VQS_Out%Main%ActFluxDiff(1)*4d0*pi
         if (do_irradiance .and. flux_pass == 0 .and. observer_height > 0) write(21,*) VQS_Out%Main%RegFluxDiff(1)
         if (do_irradiance_orders) then
            if (flux_pass == 0) then
               order_diffuse(1) = VQS_Out%Main%RegFluxDiff(1)
            else
               order_samples(1,flux_pass) = VQS_Out%Main%RegFluxDiff(1)
            endif
            if (do_irradiance_jacobians) then
               do q = 1,n_jacobians
                  do n = 1,n_layers
                     order_index = 2 + (q-1)*n_layers+n
                     if (flux_pass == 0) then
                        order_diffuse(order_index) = VQS_LinOut%Prof%LP_RegFluxDiff(1,n,q)
                     else
                        order_samples(order_index,flux_pass) = VQS_LinOut%Prof%LP_RegFluxDiff(1,n,q)
                     endif
                  enddo
               enddo
            endif
         endif


         CALL VQS_WRITE_STATUS ( &
           VQS_ERRUNIT, OPENFILEFLAG, VQS_Out%Status)
         
         ! DOWNWELLING (zenith view)
         VQS_FixIn%UserVal%ALPHA_CAMERA  = 180d0 - 1d-5 ! Look upward
         VQS_ModIn%MBool%DO_RAYLEIGH_ONLY     = do_rayleigh_only
         VQS_ModIn%MBool%DO_DELTAM_SCALING    = do_deltam_scaling
         VQS_FixIn%UserVal%HEIGHT_CAMERA = observer_height / 1d3 ! Observer height in [km]
         
         if (do_actinic_flux_jacobians) then
            if (flux_pass == 0) write(*,*) '    Calling linearized VQS master for downwelling radiation.'
            CALL VQS_LP_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out, &
               VQS_LinFixIn, VQS_LinOut)
            do q = 1, n_jacobians
               if (flux_pass == 0) write(70+q, *) VQS_LinOut%Prof%LP_ActFluxDiff(1,1:n_layers,q)*4d0*pi
               if (flux_pass == 0) write(90+q, *) VQS_LinOut%Prof%LP_ActFluxDir(1,1:n_layers,q)*4d0*pi
               if (do_irradiance_jacobians) then
                  if (flux_pass == 0) write(300+q,*) VQS_LinOut%Prof%LP_RegFluxDiff(1,1:n_layers,q)
                  if (flux_pass == 0) write(400+q,*) VQS_LinOut%Prof%LP_RegFluxDir(1,1:n_layers,q)
               endif
            enddo
         else
            if (flux_pass == 0) write(*,*) '    Calling non-linearized VQS master for downwelling radiation.'
            CALL VQS_MASTER ( do_debug_master_input, &
               VQS_FixIn, VQS_ModIn, VQS_Out )
         endif
         
         CALL VQS_WRITE_STATUS ( &
           VQS_ERRUNIT, OPENFILEFLAG, VQS_Out%Status )
         
         !write (11, '(2X ES22.15E2)') VQS_Out%Main%ActFluxDiff(1)*4d0*p
         !write (12, '(2X ES22.15E2)') VQS_Out%Main%ActFluxDir(1)*4d0*pi

			if (flux_pass == 0) write(11, *) VQS_Out%Main%ActFluxDiff(1)*4d0*pi
			if (flux_pass == 0) write(12, *) VQS_Out%Main%ActFluxDir(1)*4d0*pi
         if (do_irradiance) then
            if (flux_pass == 0) write(22,*) VQS_Out%Main%RegFluxDiff(1)
            if (flux_pass == 0) write(23,*) VQS_Out%Main%RegFluxDir(1)
         endif
         
         !write (*,*) 'Radiance total:', VQS_Out%Main%StokesTot(1)
         !write (*,*) 'Radiance first order:', VQS_Out%Main%StokesFO(1)
         !write (*,*) 'Radiance multiple scattering:', VQS_Out%Main%StokesMS(1)
         ! At the surface the upward viewing ray has no atmospheric segments,
         ! so QS leaves its local upward flux empty. Use the exact Lambertian
         ! boundary condition from the downward flux instead. Atmospheric
         ! profile derivatives hold surface albedo fixed.
         if (observer_height == 0) then
            ground_up = VQS_FixIn%Optical%LAMBERTIAN_ALBEDO * &
               (VQS_Out%Main%RegFluxDiff(1)+VQS_Out%Main%RegFluxDir(1))
            if (flux_pass == 0) then
               write(10,*) 2d0*ground_up
               if (do_irradiance) write(21,*) ground_up
               if (do_irradiance_orders) order_diffuse(1) = ground_up
            else
               order_samples(1,flux_pass) = ground_up
            endif
            if (do_actinic_flux_jacobians) then
               do q = 1,n_jacobians
                  ground_up_jac(1:n_layers) = VQS_FixIn%Optical%LAMBERTIAN_ALBEDO * &
                     (VQS_LinOut%Prof%LP_RegFluxDiff(1,1:n_layers,q)+VQS_LinOut%Prof%LP_RegFluxDir(1,1:n_layers,q))
                  if (flux_pass == 0) then
                     write(50+q,*) 2d0*ground_up_jac(1:n_layers)
                     if (do_irradiance_jacobians) write(200+q,*) ground_up_jac(1:n_layers)
                  endif
                  if (do_irradiance_orders .and. do_irradiance_jacobians) then
                     do n = 1,n_layers
                        order_index = 2+(q-1)*n_layers+n
                        if (flux_pass == 0) then
                           order_diffuse(order_index) = ground_up_jac(n)
                        else
                           order_samples(order_index,flux_pass) = ground_up_jac(n)
                        endif
                     enddo
                  endif
               enddo
            endif
         endif
         

         if (do_irradiance_orders) then
            if (flux_pass == 0) then
               order_diffuse(2) = VQS_Out%Main%RegFluxDiff(1)
            else
               order_samples(2,flux_pass) = VQS_Out%Main%RegFluxDiff(1)
            endif
            if (do_irradiance_jacobians) then
               do q = 1,n_jacobians
                  do n = 1,n_layers
                     order_index = 2 + n_jacobians*n_layers + (q-1)*n_layers+n
                     if (flux_pass == 0) then
                        order_diffuse(order_index) = VQS_LinOut%Prof%LP_RegFluxDiff(1,n,q)
                     else
                        order_samples(order_index,flux_pass) = VQS_LinOut%Prof%LP_RegFluxDiff(1,n,q)
                     endif
                  enddo
               enddo
            endif
         endif

         if (do_irradiance_orders) then
            if (flux_pass == 0) then
               order_direct = VQS_Out%Main%RegFluxDir(1)
            else
               if (abs(VQS_Out%Main%RegFluxDir(1)-order_direct) > 1d-10+1d-8*abs(order_direct)) &
                  error stop 'Direct irradiance changed during fixed-extinction scattering expansion.'
            endif
         endif
         enddo
         VQS_ModIn%MOptical%OMEGA_TOTAL_INPUT(1:n_layers) = scatter_omega(1:n_layers)
         VQS_FixIn%Optical%LAMBERTIAN_ALBEDO = scatter_albedo
         VQS_FixIn%Cont%nstokes = n_stokes
         if (do_irradiance_orders) call write_irradiance_orders( &
            order_samples(1:n_order_values,:), order_diffuse(1:n_order_values), &
            scattering_step, n_layers, n_jacobians, do_irradiance_jacobians)
      endif
      ! ###############################
      ! Calculate actinic flux profiles
      ! ###############################
      
      if (do_calculate_flux_profiles) then
      
         write(*,*) '  Simulations for actinic flux profiles:'
      
         ! UPWELLING (nadir view)
         VQS_FixIn%Bool%DO_CALCULATE_FLUXES = .true.
         VQS_FixIn%UserVal%ALPHA_CAMERA    = 1d-5
         VQS_FixIn%UserVal%HEIGHT_CAMERA   = (layer_boundary_altitudes(0) - 1d-2)/1e3 ! Put camera at TOA (1 cm below actually)
         VQS_ModIn%MBool%DO_RAYLEIGH_ONLY  = do_rayleigh_only
         VQS_ModIn%MBool%DO_DELTAM_SCALING = do_deltam_scaling
         VQS_FixIn%Bool%DO_DEBUG_FLUXES       = .true. ! MUST BE TRUE FOR ACTINIC FLUX PROFILE CALCULATIONS!
         
         write (*,*) '    Calling non-linearized VQS master for upwelling radiation.'
         CALL VQS_MASTER ( do_debug_master_input, &
            VQS_FixIn, VQS_ModIn, VQS_Out )
         
         CALL VQS_WRITE_STATUS ( &
           VQS_ERRUNIT, OPENFILEFLAG, VQS_Out%Status )

         write (13,*) VQS_Out%Main%ActFluxDiff_AllSB(1,0:n_layers)*4d0*pi ! Actinic flux profile (all segment boundaries)

         ! DOWNWELLING (zenith view)
         VQS_FixIn%UserVal%ALPHA_CAMERA  = 180d0 - 1d-5 ! Look upward
         VQS_FixIn%UserVal%HEIGHT_CAMERA = 0d0 ! Put camera on ground
         VQS_ModIn%MBool%DO_RAYLEIGH_ONLY     = do_rayleigh_only
         VQS_ModIn%MBool%DO_DELTAM_SCALING    = do_deltam_scaling
         
         write (*,*) '    Calling non-linearized VQS master for downwelling radiation.'
         CALL VQS_MASTER ( do_debug_master_input, &
            VQS_FixIn, VQS_ModIn, VQS_Out )

         CALL VQS_WRITE_STATUS ( &
           VQS_ERRUNIT, OPENFILEFLAG, VQS_Out%Status )
         
         write (14, *) VQS_Out%Main%ActFluxDiff_AllSB(1,0:n_layers) * 4d0*pi
         write (15, *) VQS_Out%Main%ActFluxDir_AllSB(1,0:n_layers) * 4d0*pi
         
      endif
      
   enddo
   
   ! CLOSE FILES
   Close(10)
   Close(11)
   Close(12)
   Close(13)
   Close(14)
   Close(15)
   Close(20)
   if (do_irradiance) call close_irradiance_files(n_jacobians, do_irradiance_jacobians, do_irradiance_orders)
   do q = 1, n_jacobians
      if (do_radiance_jacobians) then
			Close(30+q)
      endif
		if (do_actinic_flux_jacobians) then
			Close(50+q)
			Close(70+q)
			Close(90+q)
		endif
   enddo
   
   write (*,*) ' '

   write (*,*) "Execution successful!"

   stop

end program DRIVER
