program DRIVER

   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use VLIDORT_PARS_m, only : fpk, MAXLAYERS, MAXMOMENTS_INPUT, MAX_ATMOSWFS, &
      MAXSTOKES, MAX_DIRECTIONS, MAX_MESSAGES, MAX_SURFACEWFS, VLIDORT_SUCCESS, &
      ZERO, ONE, PI4, UPIDX, DNIDX
   use DRIVER_pars_m
   use MainDriverUtility_m
   use ReadInputData_m
   use RayleighScattering_m
   use VLIDORT_IO_DEFS_m
   use VLIDORT_LIN_IO_DEFS_m
   use VLIDORT_MASTERS_m
   use VLIDORT_LPS_MASTERS_m
   use VLIDORT_SPHERCORR_ROUTINES_m
   use VLIDORT_LPS_SPHERCORR_ROUTINES_m
   use VLIDORT_AC_SPHERCORR_ROUTINES_m
   use VLIDORT_AC_LPS_SPHERCORR_ROUTINES_m
   use VLIDORT_L_INPUTS_m, only : VLIDORT_L_INIT_INPUTS, VLIDORT_LinSup_Init

   implicit none

   real(kind=fpk), parameter :: k_b = 1.3806488d-23
   character(len=*), parameter :: output_path = 'OUTPUT/'

   type(VLIDORT_Fixed_Inputs)    :: VLIDORT_FixIn
   type(VLIDORT_Modified_Inputs) :: VLIDORT_ModIn
   type(VLIDORT_Sup_InOut)       :: VLIDORT_Sup
   type(VLIDORT_Outputs)         :: VLIDORT_Out
   type(VLIDORT_Fixed_LinInputs)    :: VLIDORT_LinFixIn
   type(VLIDORT_Modified_LinInputs) :: VLIDORT_LinModIn
   type(VLIDORT_LinSup_InOut)       :: VLIDORT_LinSup
   type(VLIDORT_LinOutputs)         :: VLIDORT_LinOut

   logical :: do_ss_only, do_vector_calculation, do_calculate_radiance
   logical :: do_calculate_flux_at_observer, do_calculate_flux_profiles
   logical :: do_deltam_scaling, do_rayleigh_only
   logical :: do_radiance_jacobians, do_actinic_flux_jacobians
   logical :: do_sphericity_correction, do_aircraft_correction
   logical :: fail, twilight_flag
   integer :: n_stokes, n_streams, n_fine_layers, sphercorr_points
   integer :: n_wavelengths, n_layers, n_aertypes, n_gases, n_jacobians
   integer :: dummy_integer
   real(kind=fpk) :: planet_radius, solar_zenith_angle, viewing_zenith_angle, relative_azimuth_angle
   real :: observer_height

   real(kind=fpk) :: simulated_wavelengths(MAX_WAVELENGTHS)
   real(kind=fpk) :: layer_boundary_altitudes(0:MAXLAYERS)
   real(kind=fpk) :: layer_thicknesses(MAXLAYERS)
   real(kind=fpk) :: pressure_profile(MAXLAYERS), temperature_profile(MAXLAYERS)
   real(kind=fpk) :: vapor_pressure_profile(MAXLAYERS)
   character(len=MAX_SPECIESNAME_LENGTH) :: aertype_names(MAX_AEROSOL_TYPES)
   character(len=MAX_SPECIESNAME_LENGTH) :: gas_names(MAX_GASES)
   character(len=MAX_SPECIESNAME_LENGTH) :: jacobian_list(MAX_ATMOSWFS)
   real(kind=fpk) :: aer_extinction_profiles(MAXLAYERS, MAX_AEROSOL_TYPES)
   real(kind=fpk) :: aer_extinction_multipliers(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES)
   real(kind=fpk) :: aer_ssas(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES)
   real(kind=fpk) :: aer_moments(MAX_WAVELENGTHS, MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT, 6)
   real(kind=fpk) :: gas_concentration_profiles(MAXLAYERS, MAX_GASES)
   real(kind=fpk) :: gas_xs_data(MAX_GASES, MAX_WAVELENGTHS, 0:MAX_GASXS_TEMPCOEFFS)
   real(kind=fpk) :: surface_albedos(MAX_WAVELENGTHS)
   logical :: jacobian_layer_flags(MAXLAYERS)

   integer :: GMASK(8), CMASK(8), SMASK(8)
   integer :: w, n, a, g, i, l, q, idx, gk, ck, sk, output_diridx
   integer :: n_nonzero_greekmat_entries, n_moments_input
   real(kind=fpk) :: wavelength, layer_thickness, rayleigh_xs, air_density, rayleigh_layer_od
   real(kind=fpk) :: aer_layer_extinction_ods(MAX_AEROSOL_TYPES)
   real(kind=fpk) :: aer_layer_scattering_ods(MAX_AEROSOL_TYPES)
   real(kind=fpk) :: aer_greekcoeffs_matrix(MAX_AEROSOL_TYPES, 0:MAXMOMENTS_INPUT, 16)
   real(kind=fpk) :: aer_agg_layer_extinction_od, aer_agg_layer_scattering_od
   real(kind=fpk) :: aer_agg_layer_greekcoeffs(0:MAXMOMENTS_INPUT, 6)
   real(kind=fpk) :: rayleigh_greekcoeffs(0:2, 6)
   real(kind=fpk) :: layer_greekmat(0:MAXMOMENTS_INPUT, 16)
   real(kind=fpk) :: gas_agg_layer_od, layer_extinction_od, layer_scattering_od
   real(kind=fpk) :: layer_ssa, xs, mom, h
   real(kind=fpk) :: d_extmul_d_p, d_ssa_d_p
   real(kind=fpk) :: d_moms_d_p(0:MAXMOMENTS_INPUT, 16)
   real(kind=fpk) :: observer_level, obs_height_clamped
   real(kind=fpk) :: solver_vza
   real(kind=fpk) :: radiance_stokes(4), corrected_stokes(4), scratch1(4), scratch2(4), scratch3(4)
   real(kind=fpk) :: profilewf_stokes(MAX_ATMOSWFS, MAXLAYERS)
   integer :: aerosol_model_jacobians_aertype_indexes(MAX_ATMOSWFS)
   real(kind=fpk) :: aerosol_model_jacobians_extmul(MAX_ATMOSWFS, MAX_WAVELENGTHS)
   real(kind=fpk) :: aerosol_model_jacobians_ssa(MAX_ATMOSWFS, MAX_WAVELENGTHS)
   real(kind=fpk) :: aerosol_model_jacobians_moms(MAX_ATMOSWFS, MAX_WAVELENGTHS, 0:MAXMOMENTS_INPUT, 6)
   real(kind=fpk) :: lps_stokes_a(4), lps_stokes_b(4), lps_stokes_c(4)
   real(kind=fpk) :: lps_stokes_d(4), lps_stokes_e(4), lps_stokes_f(4)
   real(kind=fpk) :: lps_profile_a(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_profile_b(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_profile_c(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_profile_d(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_profile_e(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_profile_f(MAX_ATMOSWFS, MAXLAYERS, 4)
   real(kind=fpk) :: lps_surface_a(MAX_SURFACEWFS, 4), lps_surface_b(MAX_SURFACEWFS, 4)
   real(kind=fpk) :: lps_surface_c(MAX_SURFACEWFS, 4), lps_surface_d(MAX_SURFACEWFS, 4)
   real(kind=fpk) :: lps_surface_e(MAX_SURFACEWFS, 4), lps_surface_f(MAX_SURFACEWFS, 4)
   character(len=20) :: backend_mode
   character(len=MAX_SPECIESNAME_LENGTH) :: jstr, astr
   character(len=200) :: main_message, local_message, local_action

   GMASK = (/  1,  2,  5, 6, 11, 12, 15, 16 /)
   CMASK = (/  1,  5,  5, 2,  3,  6,  6,  4 /)
   SMASK = (/  1, -1, -1, 1,  1, -1,  1,  1 /)

   call ReadConfigfile('INPUT/configuration.inp', do_ss_only, do_vector_calculation, &
      do_calculate_radiance, do_calculate_flux_at_observer, do_calculate_flux_profiles, &
      do_rayleigh_only, do_deltam_scaling, n_stokes, n_streams, n_fine_layers, &
      planet_radius, observer_height, solar_zenith_angle, viewing_zenith_angle, &
      relative_azimuth_angle, do_radiance_jacobians, do_actinic_flux_jacobians, &
      backend_mode, sphercorr_points)
   call decode_backend_mode(backend_mode, do_sphericity_correction, do_aircraft_correction)
   solver_vza = min(viewing_zenith_angle, 180.0_fpk - viewing_zenith_angle)
   if (do_ss_only) then
      do_sphericity_correction = .false.
      do_aircraft_correction = .false.
   endif

   if (n_stokes == 1) then
      n_nonzero_greekmat_entries = 1
   else if (n_stokes == 2) then
      n_nonzero_greekmat_entries = 4
   else if (n_stokes == 3) then
      n_nonzero_greekmat_entries = 5
   else
      n_nonzero_greekmat_entries = 8
   endif

   call ReadArrayFromFile('INPUT/wavelengths.inp', MAX_WAVELENGTHS, 1, &
      n_wavelengths, dummy_integer, simulated_wavelengths)
   call ReadVerticalGrid('INPUT/vertical_grid.inp', MAXLAYERS, n_layers, layer_boundary_altitudes)
   layer_thicknesses(1:n_layers) = layer_boundary_altitudes(0:n_layers-1) - layer_boundary_altitudes(1:n_layers)

   call ReadArrayFromFile('INPUT/pressure_profile.inp', MAXLAYERS, 1, n_layers, dummy_integer, pressure_profile)
   call ReadArrayFromFile('INPUT/temperature_profile.inp', MAXLAYERS, 1, n_layers, dummy_integer, temperature_profile)
   call ReadArrayFromFile('INPUT/vapor_pressure_profile.inp', MAXLAYERS, 1, n_layers, dummy_integer, vapor_pressure_profile)
   call ReadSpeciesNames('INPUT/aertype_names.inp', MAX_AEROSOL_TYPES, MAX_SPECIESNAME_LENGTH, aertype_names, n_aertypes)
   call ReadArrayFromFile('INPUT/aerosol_extinction_profiles.inp', MAXLAYERS, MAX_AEROSOL_TYPES, &
      n_layers, n_aertypes, aer_extinction_profiles)
   call ReadArrayFromFile('INPUT/aerosol_extinction_multipliers.inp', MAX_WAVELENGTHS, MAX_AEROSOL_TYPES, &
      n_wavelengths, n_aertypes, aer_extinction_multipliers)
   call ReadArrayFromFile('INPUT/aerosol_single_scattering_albedos.inp', MAX_WAVELENGTHS, MAX_AEROSOL_TYPES, &
      n_wavelengths, n_aertypes, aer_ssas)
   call ReadAerosolMoments('INPUT/aerosol_moments/', MAX_AEROSOL_TYPES, MAX_WAVELENGTHS, MAXMOMENTS_INPUT, &
      MAX_SPECIESNAME_LENGTH, n_aertypes, n_wavelengths, aertype_names, aer_moments)

   call ReadSpeciesNames('INPUT/gas_names.inp', MAX_GASES, MAX_SPECIESNAME_LENGTH, gas_names, n_gases)
   call ReadGasXS('INPUT/xs/', MAX_GASES, MAX_WAVELENGTHS, MAX_GASXS_TEMPCOEFFS, &
      MAX_SPECIESNAME_LENGTH, n_gases, n_wavelengths, gas_names, gas_xs_data)
   call ReadArrayFromFile('INPUT/gas_concentration_profiles.inp', MAXLAYERS, MAX_GASES, &
      n_layers, n_gases, gas_concentration_profiles)
   call ReadArrayFromFile('INPUT/surface_albedos.inp', MAX_WAVELENGTHS, 1, n_wavelengths, dummy_integer, surface_albedos)
   call ReadSpeciesNames('INPUT/jacobian_list.inp', MAX_ATMOSWFS, MAX_SPECIESNAME_LENGTH, jacobian_list, n_jacobians)
   call ReadLogicalArrayFromFile('INPUT/jacobian_layer_flags.inp', MAXLAYERS, 1, n_layers, dummy_integer, jacobian_layer_flags)

   aerosol_model_jacobians_aertype_indexes = 0
   aerosol_model_jacobians_extmul = ZERO
   aerosol_model_jacobians_ssa = ZERO
   aerosol_model_jacobians_moms = ZERO
   if (do_radiance_jacobians .or. do_actinic_flux_jacobians) then
      do q = 1, n_jacobians
         jstr = jacobian_list(q)
         if (index(trim(jstr), '_') /= 0 .and. find_string_in_array(trim(jstr), aertype_names) <= 0) then
            idx = index(trim(jstr), '_', back=.true.)
            astr = jstr(1:idx-1)
            aerosol_model_jacobians_aertype_indexes(q) = find_string_in_array(trim(astr), aertype_names)
            if (aerosol_model_jacobians_aertype_indexes(q) <= 0) then
               write(*,*) 'Unknown aerosol type in Jacobian label: ', trim(jstr)
               stop 2
            endif
            call ReadArrayFromFile('INPUT/aerosol_model_jacobians/' // trim(jstr) // '_extmul.inp', &
               MAX_WAVELENGTHS, 1, n_wavelengths, dummy_integer, aerosol_model_jacobians_extmul(q,:))
            call ReadArrayFromFile('INPUT/aerosol_model_jacobians/' // trim(jstr) // '_ssa.inp', &
               MAX_WAVELENGTHS, 1, n_wavelengths, dummy_integer, aerosol_model_jacobians_ssa(q,:))
            call ReadAerosolMomentJacobians('INPUT/aerosol_model_jacobians/' // trim(jstr) // '_moms.inp', &
               MAX_WAVELENGTHS, MAXMOMENTS_INPUT, aerosol_model_jacobians_moms(q,:,:,:))
         endif
      enddo
   endif

   call initialize_vlidort_inputs()
   call initialize_linearized_inputs()

   open(10, file=trim(output_path)//'flux_diffuse_upwelling.out', status='replace')
   open(11, file=trim(output_path)//'flux_diffuse_downwelling.out', status='replace')
   open(12, file=trim(output_path)//'flux_direct.out', status='replace')
   open(20, file=trim(output_path)//'radiances.out', status='replace')

   do q = 1, n_jacobians
      if (do_radiance_jacobians) then
         open(30+q, file=trim(output_path)//'jacs_radiance_'//trim(jacobian_list(q))//'.out', status='replace')
      endif
      if (do_actinic_flux_jacobians) then
         open(50+q, file=trim(output_path)//'jacs_flux_diffuse_upwelling_'//trim(jacobian_list(q))//'.out', &
            status='replace')
         open(70+q, file=trim(output_path)//'jacs_flux_diffuse_downwelling_'//trim(jacobian_list(q))//'.out', &
            status='replace')
         open(90+q, file=trim(output_path)//'jacs_flux_direct_'//trim(jacobian_list(q))//'.out', status='replace')
      endif
   enddo

   do w = 1, n_wavelengths
      wavelength = simulated_wavelengths(w)
      write(*,*) 'Simulating wavelength ', w, ' of ', n_wavelengths, ': ', wavelength, ' nm'

      n_moments_input = 0
      do a = 1, n_aertypes
         do l = 0, MAXMOMENTS_INPUT
            do i = 1, 6
               if (aer_moments(w,a,l,i) > 0.0_fpk) n_moments_input = max(n_moments_input, l)
            enddo
         enddo
      enddo
      n_moments_input = min(n_moments_input, MAXMOMENTS_INPUT)
      if (do_rayleigh_only) n_moments_input = 2

      rayleigh_greekcoeffs = RayleighGreekCoeffs(wavelength)
      aer_greekcoeffs_matrix = ZERO
      do i = 1, n_nonzero_greekmat_entries
         gk = GMASK(i); ck = CMASK(i); sk = SMASK(i)
         aer_greekcoeffs_matrix(1:n_aertypes,0:n_moments_input,gk) = &
            aer_moments(w,1:n_aertypes,0:n_moments_input,ck) * sk
      enddo
      VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = surface_albedos(w)
      VLIDORT_FixIn%Optical%TS_ATMOS_WAVELENGTH = wavelength / 1000.0_fpk
      VLIDORT_FixIn%Optical%TS_ATMOS_INDEX = w
      VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = n_moments_input
      VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT = ZERO
      VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT = ZERO
      VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_UP = ZERO
      VLIDORT_LinFixIn%Optical%TS_L_FMATRIX_DN = ZERO

      do n = 1, n_layers
         layer_thickness = layer_thicknesses(n)
         rayleigh_xs = RayleighXS(wavelength, pressure_profile(n), temperature_profile(n), vapor_pressure_profile(n))
         air_density = (pressure_profile(n) * 1.0d2) / (k_b * temperature_profile(n))
         rayleigh_layer_od = rayleigh_xs * 1.0d-4 * air_density * layer_thickness

         aer_layer_extinction_ods = ZERO
         aer_layer_scattering_ods = ZERO
         aer_layer_extinction_ods(1:n_aertypes) = aer_extinction_profiles(n,1:n_aertypes) * &
            aer_extinction_multipliers(w,1:n_aertypes) * layer_thickness / 1.0d3
         aer_layer_scattering_ods(1:n_aertypes) = aer_layer_extinction_ods(1:n_aertypes) * aer_ssas(w,1:n_aertypes)
         aer_agg_layer_extinction_od = sum(aer_layer_extinction_ods(1:n_aertypes))
         aer_agg_layer_scattering_od = sum(aer_layer_scattering_ods(1:n_aertypes))

         aer_agg_layer_greekcoeffs = ZERO
         if (aer_agg_layer_scattering_od > ZERO) then
            do a = 1, n_aertypes
               aer_agg_layer_greekcoeffs(0:n_moments_input,:) = aer_agg_layer_greekcoeffs(0:n_moments_input,:) + &
                  aer_layer_scattering_ods(a) * aer_moments(w,a,0:n_moments_input,:)
            enddo
            aer_agg_layer_greekcoeffs = aer_agg_layer_greekcoeffs / aer_agg_layer_scattering_od
         endif

         gas_agg_layer_od = ZERO
         do g = 1, n_gases
            xs = gas_xs_data(g,w,0)
            do i = 1, MAX_GASXS_TEMPCOEFFS
               xs = xs + gas_xs_data(g,w,i) * (temperature_profile(n) - 273.15_fpk)**i
            enddo
            gas_agg_layer_od = gas_agg_layer_od + gas_concentration_profiles(n,g) * xs * layer_thickness * 1.0d2
         enddo

         layer_extinction_od = rayleigh_layer_od + aer_agg_layer_extinction_od + gas_agg_layer_od
         layer_scattering_od = rayleigh_layer_od + aer_agg_layer_scattering_od
         if (layer_extinction_od <= ZERO) layer_extinction_od = 1.0d-12
         layer_ssa = layer_scattering_od / layer_extinction_od
         if (.not. ieee_is_finite(layer_ssa) .or. .not. ieee_is_finite(layer_extinction_od) .or. &
             layer_ssa < ZERO .or. layer_ssa > ONE + 1.0d-12) then
            write(*,*) 'Invalid optical properties at wavelength/layer:', w, n, layer_extinction_od, layer_ssa
            stop 7
         endif

         layer_greekmat = ZERO
         if (layer_scattering_od > ZERO) then
            do i = 1, n_nonzero_greekmat_entries
               gk = GMASK(i); ck = CMASK(i); sk = SMASK(i)
               layer_greekmat(0:2,gk) = rayleigh_layer_od * rayleigh_greekcoeffs(0:2,ck)
               layer_greekmat(0:n_moments_input,gk) = layer_greekmat(0:n_moments_input,gk) + &
                  aer_agg_layer_scattering_od * sk * aer_agg_layer_greekcoeffs(0:n_moments_input,ck)
            enddo
            layer_greekmat = layer_greekmat / layer_scattering_od
         endif
         mom = layer_greekmat(0,1)
         if (abs(mom - ONE) < 1.0d-14 .and. mom /= ZERO) layer_greekmat = layer_greekmat / mom

         VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n) = layer_extinction_od
         VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = layer_ssa
         VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(:,n,:) = layer_greekmat

         do q = 1, n_jacobians
            jstr = jacobian_list(q)

             if (trim(jstr) == 'absorption') then
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = ONE / layer_extinction_od
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = -ONE / layer_extinction_od
             endif

            do a = 1, n_aertypes
               if (trim(jstr) == trim(aertype_names(a))) then
                  VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = &
                     aer_extinction_multipliers(w,a) * layer_thickness / 1.0d3 / layer_extinction_od
                  if (layer_scattering_od > ZERO) then
                     VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = &
                        aer_extinction_multipliers(w,a) * layer_thickness / 1.0d3 / layer_scattering_od * &
                        (aer_ssas(w,a) - layer_ssa)
                     do l = 0, n_moments_input
                        do i = 1, 16
                           if (layer_greekmat(l,i) /= ZERO) then
                              VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,l,n,i) = &
                                 aer_extinction_multipliers(w,a) * layer_thickness / 1.0d3 * aer_ssas(w,a) / &
                                 layer_greekmat(l,i) / layer_scattering_od * &
                                 (aer_greekcoeffs_matrix(a,l,i) - layer_greekmat(l,i))
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo

            if (aerosol_model_jacobians_aertype_indexes(q) > 0) then
               a = aerosol_model_jacobians_aertype_indexes(q)
               d_extmul_d_p = aerosol_model_jacobians_extmul(q,w)
               d_ssa_d_p = aerosol_model_jacobians_ssa(q,w)
               d_moms_d_p = ZERO
               do i = 1, n_nonzero_greekmat_entries
                  gk = GMASK(i); ck = CMASK(i); sk = SMASK(i)
                  d_moms_d_p(0:n_moments_input,gk) = &
                     aerosol_model_jacobians_moms(q,w,0:n_moments_input,ck) * sk
               enddo
               h = layer_thickness / 1.0d3

               VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = &
                  h * aer_extinction_profiles(n,a) / layer_extinction_od * d_extmul_d_p
               if (layer_scattering_od > ZERO) then
                  VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = &
                     h * aer_extinction_profiles(n,a) / layer_scattering_od * &
                     ((aer_ssas(w,a) - layer_ssa) * d_extmul_d_p + &
                      aer_extinction_multipliers(w,a) * d_ssa_d_p)
                  do l = 0, n_moments_input
                     do i = 1, 16
                        if (layer_greekmat(l,i) /= ZERO) then
                           VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,l,n,i) = &
                              ONE / layer_greekmat(l,i) / layer_scattering_od * &
                              (h * aer_extinction_profiles(n,a) * &
                               (aer_greekcoeffs_matrix(a,l,i) - layer_greekmat(l,i)) * &
                               (aer_ssas(w,a) * d_extmul_d_p + aer_extinction_multipliers(w,a) * d_ssa_d_p) + &
                               aer_layer_scattering_ods(a) * d_moms_d_p(l,i))
                        endif
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo

      if (do_calculate_radiance) then
         if (do_radiance_jacobians) then
            call run_vlidort_jacobians_for_geometry(.true., radiance_stokes, profilewf_stokes)
            if (.not. all(ieee_is_finite(profilewf_stokes(1:n_jacobians,1:n_layers)))) then
               write(*,*) 'Nonfinite native radiance Jacobian.'
               stop 7
            endif
            do q = 1, n_jacobians
               write(30+q,*) profilewf_stokes(q,1:n_layers)
            enddo
         else
            call run_vlidort_for_geometry(.true., radiance_stokes)
         endif
         if (.not. all(ieee_is_finite(radiance_stokes(1:n_stokes)))) then
            write(*,*) 'Nonfinite radiance; check geometry and optical properties.'
            stop 7
         endif
         write(20,*) radiance_stokes(1)
      endif

      if (do_calculate_flux_at_observer) then
         if (do_actinic_flux_jacobians) then
            call run_vlidort_jacobians_for_geometry(.false., radiance_stokes, profilewf_stokes)
            do q = 1, n_jacobians
               write(50+q,*) VLIDORT_LinOut%Prof%TS_MEANST_DIFFUSE_PROFWF(q,1:n_layers,1,1,1,UPIDX) * PI4
               write(70+q,*) VLIDORT_LinOut%Prof%TS_MEANST_DIFFUSE_PROFWF(q,1:n_layers,1,1,1,DNIDX) * PI4
               write(90+q,*) VLIDORT_LinOut%Prof%TS_DNMEANST_DIRECT_PROFWF(q,1:n_layers,1,1,1) * PI4
            enddo
         else
            call run_vlidort_for_geometry(.false., radiance_stokes)
         endif
         write(10,*) VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1,1,1,UPIDX) * PI4
         write(11,*) VLIDORT_Out%Main%TS_MEANST_DIFFUSE(1,1,1,DNIDX) * PI4
         write(12,*) VLIDORT_Out%Main%TS_DNMEANST_DIRECT(1,1,1) * PI4
      endif
   enddo

   do q = 1, n_jacobians
      if (do_radiance_jacobians) close(30+q)
      if (do_actinic_flux_jacobians) then
         close(50+q); close(70+q); close(90+q)
      endif
   enddo
   close(10); close(11); close(12); close(20)

contains

   subroutine initialize_vlidort_inputs()
      integer :: n

      VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE = .not. do_ss_only
      VLIDORT_FixIn%Bool%TS_DO_THERMAL_EMISSION = .false.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_EMISSION = .false.
      VLIDORT_FixIn%Bool%TS_DO_PLANE_PARALLEL = .not. do_sphericity_correction
      VLIDORT_FixIn%Bool%TS_DO_UPWELLING = .true.
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING = .true.
      VLIDORT_FixIn%Bool%TS_DO_TOA_CONTRIBS = .false.
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .true.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_1 = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_2 = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPECIALIST_OPTION_3 = .false.
      VLIDORT_FixIn%Bool%TS_DO_SURFACE_LEAVING = .false.
      VLIDORT_FixIn%Bool%TS_DO_SL_ISOTROPIC = .false.
      VLIDORT_FixIn%Bool%TS_DO_WATER_LEAVING = .false.
      VLIDORT_FixIn%Bool%TS_DO_FLUORESCENCE = .false.
      VLIDORT_FixIn%Bool%TS_DO_TF_ITERATION = .false.
      VLIDORT_FixIn%Bool%TS_DO_WLADJUSTED_OUTPUT = .false.
      VLIDORT_FixIn%Bool%TS_DO_TOA_ILLUMINATION = .false.
      VLIDORT_FixIn%Bool%TS_DO_BOA_ILLUMINATION = .false.
      VLIDORT_FixIn%Bool%TS_DO_ALBTRN_MEDIA = .false.
      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .false.
      VLIDORT_FixIn%Bool%TS_DO_SPHERICITY_CORRECTION = do_sphericity_correction
      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.
      ! The supplier uses the reduced I/Q system for the azimuth-independent mode.
      VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2 = (n_stokes > 1)

      VLIDORT_FixIn%Cont%TS_TAYLOR_ORDER = 2
      VLIDORT_FixIn%Cont%TS_NSTOKES = n_stokes
      VLIDORT_FixIn%Cont%TS_NSTREAMS = n_streams
      VLIDORT_FixIn%Cont%TS_NLAYERS = n_layers
      VLIDORT_FixIn%Cont%TS_NFINELAYERS = n_fine_layers
      VLIDORT_FixIn%Cont%TS_N_THERMAL_COEFFS = 2
      VLIDORT_FixIn%Cont%TS_VLIDORT_ACCURACY = 1.0d-6
      ! Near-conservative vector Jacobians need tightly converged eigenvectors
      ! as well as the binary128 arithmetic selected by the gfortran build.
      VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-28
      VLIDORT_FixIn%Cont%TS_NLAYERS_NOMS = n_layers
      VLIDORT_FixIn%Cont%TS_NLAYERS_CUTOFF = n_layers
      VLIDORT_FixIn%Cont%TS_TF_MAXITER = 0
      VLIDORT_FixIn%Cont%TS_TF_CRITERION = ZERO
      VLIDORT_FixIn%Cont%TS_TOA_ILLUMINATION = ZERO
      VLIDORT_FixIn%Cont%TS_BOA_ILLUMINATION = ZERO
      if (do_aircraft_correction) then
         VLIDORT_FixIn%Cont%TS_SPHERCORR_TYPE = '2'
      else if (sphercorr_points <= 2) then
         VLIDORT_FixIn%Cont%TS_SPHERCORR_TYPE = '2'
      else
         VLIDORT_FixIn%Cont%TS_SPHERCORR_TYPE = '3'
      endif
      VLIDORT_FixIn%Cont%TS_SPHERCORR_3PTFACTOR = 0.5d0
      VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = ONE
      VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = 1

      VLIDORT_FixIn%Chapman%TS_HEIGHT_GRID(0:n_layers) = layer_boundary_altitudes(0:n_layers) / 1000.0_fpk
      do n = 0, n_layers
         VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(n) = pressure_at_boundary(n)
         VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(n) = temperature_at_boundary(n)
      enddo
      VLIDORT_FixIn%Chapman%TS_FINEGRID(1:n_layers) = max(1, n_fine_layers)
      VLIDORT_FixIn%Chapman%TS_RFINDEX_PARAMETER = ZERO
      VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT = ZERO
      VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT = ZERO
      VLIDORT_FixIn%Write%TS_DO_DEBUG_WRITE = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_INPUT = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_SCENARIO = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_FOURIER = .false.
      VLIDORT_FixIn%Write%TS_DO_WRITE_RESULTS = .false.

      VLIDORT_ModIn%MBool%TS_DO_FOCORR = .true.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL = .false.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR = .false.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .true.
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_TRUNCATION = .false.
      VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT = .false.
      VLIDORT_ModIn%MBool%TS_DO_EXTERNAL_WLEAVE = .false.
      VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST = .false.
      VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES = .true.
      VLIDORT_ModIn%MBool%TS_DO_CLASSICAL_SOLUTION = .true.
      VLIDORT_ModIn%MBool%TS_DO_REFRACTIVE_GEOMETRY = .false.
      VLIDORT_ModIn%MBool%TS_DO_CHAPMAN_FUNCTION = .true.
      VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY = do_rayleigh_only
      VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING = do_deltam_scaling
      VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = .false.
      VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = .false.
      VLIDORT_ModIn%MBool%TS_DO_USER_VZANGLES = .true.
      VLIDORT_ModIn%MBool%TS_DO_ADDITIONAL_MVOUT = .true.
      VLIDORT_ModIn%MBool%TS_DO_MVOUT_ONLY = .false.
      VLIDORT_ModIn%MBool%TS_DO_THERMAL_TRANSONLY = .false.
      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY = .true.
      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY = .false.

      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES = 1
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1) = solar_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = 1
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1) = relative_azimuth_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 1
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1) = viewing_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 1
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = solar_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = viewing_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = relative_azimuth_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS = 0
      VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = ZERO
      obs_height_clamped = min(max(real(observer_height,kind=fpk), layer_boundary_altitudes(n_layers)), &
         layer_boundary_altitudes(0))
      observer_level = height_to_user_level(obs_height_clamped)
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = observer_level
      VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS = planet_radius
   end subroutine initialize_vlidort_inputs

   subroutine initialize_linearized_inputs()
      logical :: do_jacobians
      integer :: layer_index

      call VLIDORT_L_INIT_INPUTS(VLIDORT_LinFixIn, VLIDORT_LinModIn)
      call VLIDORT_LinSup_Init(VLIDORT_LinSup)

      do_jacobians = do_radiance_jacobians .or. do_actinic_flux_jacobians
      if (.not. do_jacobians) return

      VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = n_jacobians
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES = ' '
      VLIDORT_LinFixIn%Cont%TS_PROFILEWF_NAMES(1:n_jacobians) = jacobian_list(1:n_jacobians)
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG = .false.
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER = 0
      VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG(1:n_layers) = jacobian_layer_flags(1:n_layers)
      do layer_index = 1, n_layers
         if (jacobian_layer_flags(layer_index)) then
            VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(layer_index) = n_jacobians
         endif
      enddo

      VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION = .false.
      VLIDORT_LinModIn%MCont%TS_DO_PROFILE_LINEARIZATION = .true.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION = .true.
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .false.
      VLIDORT_LinModIn%MCont%TS_DO_LINEARIZATION = .true.
      VLIDORT_LinModIn%MCont%TS_DO_SIMULATION_ONLY = .false.
      VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LBBF = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LBBF = .false.
      VLIDORT_LinModIn%MCont%TS_DO_SLEAVE_WFS = .false.
   end subroutine initialize_linearized_inputs

   subroutine run_vlidort_for_geometry(for_radiance, stokes_out)
      logical, intent(in) :: for_radiance
      real(kind=fpk), intent(out) :: stokes_out(4)

      stokes_out = ZERO
      fail = .false.; twilight_flag = .false.
      main_message = ' '; local_message = ' '; local_action = ' '

      if (for_radiance) then
         call configure_radiance_geometry(output_diridx)
      else
         call configure_flux_geometry(output_diridx)
      endif

      if (do_sphericity_correction .and. for_radiance) then
         if (do_aircraft_correction) then
            call AC_SPHERICITY_CORRECTION_MASTER(.false., VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, &
               scratch1, scratch2, scratch3, lps_stokes_d, corrected_stokes, stokes_out, &
               VLIDORT_Out, twilight_flag, fail, main_message, local_message, local_action)
         else
            call SPHERICITY_CORRECTION_MASTER(.false., VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, &
               scratch1, scratch2, scratch3, VLIDORT_Out, lps_stokes_d, stokes_out, &
               twilight_flag, fail, main_message, local_message, local_action)
         endif
      else
         if (for_radiance) call prepare_standard_first_order(.false.)
         if (for_radiance .and. do_ss_only) then
            stokes_out = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,:,output_diridx)
            if (output_diridx == UPIDX) stokes_out = stokes_out + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,:)
            return
         endif
         call VLIDORT_MASTER(.false., VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_Out)
         stokes_out(1:MAXSTOKES) = VLIDORT_Out%Main%TS_STOKES(1,1,1:MAXSTOKES,output_diridx)
      endif

      if (fail) then
         write(*,*) trim(main_message)
         write(*,*) trim(local_message)
         write(*,*) trim(local_action)
         stop 3
      endif
      if (VLIDORT_Out%Status%TS_STATUS_INPUTCHECK /= VLIDORT_SUCCESS) then
         write(*,*) 'VLIDORT input check failed:'
         write(*,*) trim(VLIDORT_Out%Status%TS_CHECKMESSAGES(1))
         write(*,*) trim(VLIDORT_Out%Status%TS_ACTIONS(1))
         stop 4
      endif
      if (VLIDORT_Out%Status%TS_STATUS_CALCULATION /= VLIDORT_SUCCESS) then
         write(*,*) 'VLIDORT calculation failed: ', trim(VLIDORT_Out%Status%TS_MESSAGE)
         write(*,*) trim(VLIDORT_Out%Status%TS_TRACE_1)
         stop 5
      endif
      if (for_radiance) call write_rt_diagnostics(.false., stokes_out)
   end subroutine run_vlidort_for_geometry

   subroutine run_vlidort_jacobians_for_geometry(for_radiance, stokes_out, profilewf_out)
      logical, intent(in) :: for_radiance
      real(kind=fpk), intent(out) :: stokes_out(4)
      real(kind=fpk), intent(out) :: profilewf_out(MAX_ATMOSWFS, MAXLAYERS)

      stokes_out = ZERO
      profilewf_out = ZERO
      fail = .false.; twilight_flag = .false.
      main_message = ' '; local_message = ' '; local_action = ' '

      if (for_radiance) then
         call configure_radiance_geometry(output_diridx)
      else
         call configure_flux_geometry(output_diridx)
      endif

      if (do_sphericity_correction .and. for_radiance) then
         if (do_aircraft_correction) then
            call AC_LPS_SPHERICITY_CORRECTION_MASTER(.false., &
               VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, &
               lps_stokes_a, lps_stokes_b, lps_stokes_c, lps_stokes_d, lps_stokes_e, lps_stokes_f, &
               lps_profile_a, lps_profile_b, lps_profile_c, lps_profile_d, lps_profile_e, lps_profile_f, &
               lps_surface_a, lps_surface_b, lps_surface_c, lps_surface_d, lps_surface_e, lps_surface_f, &
               VLIDORT_Out, VLIDORT_LinOut, twilight_flag, fail, main_message, local_message, local_action)
            stokes_out = lps_stokes_f
            profilewf_out(1:n_jacobians,1:n_layers) = lps_profile_f(1:n_jacobians,1:n_layers,1)
         else
            call SPHERICITY_CORRECTION_LPS_MASTER(.false., &
               VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, &
               lps_stokes_a, lps_stokes_b, lps_stokes_c, lps_profile_a, lps_profile_b, lps_profile_c, &
               lps_surface_a, lps_surface_b, lps_surface_c, VLIDORT_Out, VLIDORT_LinOut, &
               lps_stokes_d, lps_stokes_e, lps_profile_d, lps_profile_e, lps_surface_d, lps_surface_e, &
               twilight_flag, fail, main_message, local_message, local_action)
            stokes_out = lps_stokes_e
            profilewf_out(1:n_jacobians,1:n_layers) = lps_profile_e(1:n_jacobians,1:n_layers,1)
         endif
      else
         if (for_radiance) call prepare_standard_first_order(.true.)
         if (for_radiance .and. do_ss_only) then
            stokes_out = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,:,output_diridx)
            profilewf_out(1:n_jacobians,1:n_layers) = &
               VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_SS(1:n_jacobians,1:n_layers,1,1,1,output_diridx)
            if (output_diridx == UPIDX) then
               stokes_out = stokes_out + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,:)
               profilewf_out(1:n_jacobians,1:n_layers) = profilewf_out(1:n_jacobians,1:n_layers) + &
                  VLIDORT_LinSup%SS%Prof%TS_PROFILEWF_DB(1:n_jacobians,1:n_layers,1,1,1)
            endif
            return
         endif
         call VLIDORT_LPS_MASTER(.false., VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, VLIDORT_Out, &
            VLIDORT_LinFixIn, VLIDORT_LinModIn, VLIDORT_LinSup, VLIDORT_LinOut)
         stokes_out(1:MAXSTOKES) = VLIDORT_Out%Main%TS_STOKES(1,1,1:MAXSTOKES,output_diridx)
         profilewf_out(1:n_jacobians,1:n_layers) = &
            VLIDORT_LinOut%Prof%TS_PROFILEWF(1:n_jacobians,1:n_layers,1,1,1,output_diridx)
      endif

      if (fail) then
         write(*,*) trim(main_message)
         write(*,*) trim(local_message)
         write(*,*) trim(local_action)
         stop 3
      endif
      if (VLIDORT_Out%Status%TS_STATUS_INPUTCHECK /= VLIDORT_SUCCESS) then
         write(*,*) 'VLIDORT linearized input check failed:'
         write(*,*) trim(VLIDORT_Out%Status%TS_CHECKMESSAGES(1))
         write(*,*) trim(VLIDORT_Out%Status%TS_ACTIONS(1))
         stop 4
      endif
      if (VLIDORT_Out%Status%TS_STATUS_CALCULATION /= VLIDORT_SUCCESS) then
         write(*,*) 'VLIDORT linearized calculation failed: ', trim(VLIDORT_Out%Status%TS_MESSAGE)
         write(*,*) trim(VLIDORT_Out%Status%TS_TRACE_1)
         stop 5
      endif
      if (for_radiance) call write_rt_diagnostics(.true., stokes_out)
   end subroutine run_vlidort_jacobians_for_geometry

   subroutine write_rt_diagnostics(linearized, stokes)
      logical, intent(in) :: linearized
      real(kind=fpk), intent(in) :: stokes(4)
      character(len=8) :: enabled
      integer :: unit, layer, status
      real(kind=fpk) :: fo, ms, components(6)

      call get_environment_variable('PYDORT_RT_DIAGNOSTICS', enabled, status=status)
      if (status /= 0 .or. trim(enabled) /= '1') return
      if (do_aircraft_correction) then
         if (linearized) then
            components = (/lps_stokes_a(1), lps_stokes_b(1), lps_stokes_c(1), &
                           lps_stokes_d(1), lps_stokes_e(1), lps_stokes_f(1)/)
         else
            components = (/scratch1(1), scratch2(1), scratch3(1), &
                           lps_stokes_d(1), corrected_stokes(1), stokes(1)/)
         endif
      else
         fo = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1,output_diridx)
         if (output_diridx == UPIDX) fo = fo + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1)
         ms = stokes(1) - fo
         components = (/fo, ms, stokes(1), fo, ms, stokes(1)/)
      endif
      open(newunit=unit, file='OUTPUT/rt_components.out', status='replace')
      write(unit,*) '# original_FO original_MS original_total corrected_FO corrected_MS corrected_total'
      write(unit,*) components
      close(unit)
      if (.not. do_aircraft_correction) return
      open(newunit=unit, file='OUTPUT/rt_layers.out', status='replace')
      write(unit,*) '# layer z_top_m z_bottom_m tau mu_top mu_bottom trans_boa MS_boa MS_ac MS_toa'
      do layer = 1, n_layers
         write(unit,*) layer, layer_boundary_altitudes(layer-1), layer_boundary_altitudes(layer), &
            VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(layer), &
            cos(VLIDORT_Out%Main%TS_PATHGEOMS(2,layer-1)), cos(VLIDORT_Out%Main%TS_PATHGEOMS(2,layer)), &
            VLIDORT_Out%Main%TS_LOSTRANS(1,layer), VLIDORT_Out%Main%TS_LAYER_MSSTS(1:3,1,layer)
      enddo
      close(unit)
      if (linearized) then
         open(newunit=unit, file='OUTPUT/rt_profile_components.out', status='replace')
         write(unit,*) '# layer dFO_dtau_abs dMS_dtau_abs dtotal_dtau_abs'
         do layer = 1, n_layers
            write(unit,*) layer, lps_profile_d(1,layer,1), lps_profile_e(1,layer,1), lps_profile_f(1,layer,1)
         enddo
         close(unit)
      endif
      call sample_ac_sources()
   end subroutine write_rt_diagnostics

   subroutine sample_ac_sources()
      type(VLIDORT_Fixed_Inputs) :: sample_fix
      type(VLIDORT_Modified_Inputs) :: sample_mod
      type(VLIDORT_Sup_InOut) :: sample_sup
      type(VLIDORT_Outputs) :: sample_out
      character(len=8) :: enabled
      integer :: status, unit, layer, first_layer, last_layer
      real(kind=fpk) :: radius, boa, altitude, height, theta, alpha, phi, scattering
      real(kind=fpk) :: theta_boa, alpha_boa, phi_boa, theta_local, alpha_local, phi_local

      call get_environment_variable('PYDORT_RT_SAMPLE_SOURCES', enabled, status=status)
      if (status /= 0 .or. trim(enabled) /= '1') return
      sample_fix = VLIDORT_FixIn
      sample_mod = VLIDORT_ModIn
      sample_sup = VLIDORT_Sup
      first_layer = 1; last_layer = nint(observer_level)
      if (output_diridx == UPIDX) then
         first_layer = nint(observer_level) + 1; last_layer = n_layers
      endif
      open(newunit=unit, file='OUTPUT/rt_sampled_sources.out', status='replace')
      write(unit,*) '# layer local_sza_deg local_vza_deg local_raa_deg directly_evaluated_MS'
      do layer = first_layer, last_layer
         radius = planet_radius
         boa = layer_boundary_altitudes(n_layers) / 1000.0_fpk
         altitude = obs_height_clamped / 1000.0_fpk
         height = 0.5_fpk * (layer_boundary_altitudes(layer-1) + layer_boundary_altitudes(layer)) / 1000.0_fpk
         theta = solar_zenith_angle; alpha = solver_vza; phi = relative_azimuth_angle
         if (output_diridx == DNIDX) then
            call AC_to_BOATOA_geometry_conversion(radius, boa, altitude, height, theta, alpha, phi, &
               scattering, theta_boa, alpha_boa, phi_boa, theta_local, alpha_local, phi_local)
         else
            call AC_to_BOA_geometry_conversion(radius, height, altitude, theta, alpha, phi, &
               scattering, theta_local, alpha_local, phi_local)
         endif
         sample_mod%MSunrays%TS_SZANGLES(1) = theta_local
         sample_mod%MUserVal%TS_USER_VZANGLES_INPUT(1) = alpha_local
         sample_mod%MUserVal%TS_USER_RELAZMS(1) = phi_local
         sample_mod%MUserVal%TS_USER_OBSGEOMS_INPUT(1,:) = (/theta_local, alpha_local, phi_local/)
         call VLIDORT_MASTER(.false., sample_fix, sample_mod, sample_sup, sample_out)
         if (sample_out%Status%TS_STATUS_INPUTCHECK /= VLIDORT_SUCCESS .or. &
             sample_out%Status%TS_STATUS_CALCULATION /= VLIDORT_SUCCESS) then
            write(*,*) 'Diagnostic source sampling failed at layer', layer
            stop 8
         endif
         write(unit,*) layer, theta_local, alpha_local, phi_local, sample_out%Main%TS_LAYER_MSSTS(1,1,layer)
      enddo
      close(unit)
   end subroutine sample_ac_sources

   subroutine prepare_standard_first_order(linearized)
      logical, intent(in) :: linearized
      type(VLIDORT_Fixed_Inputs) :: fo_fix
      type(VLIDORT_Modified_Inputs) :: fo_mod
      type(VLIDORT_Modified_LinInputs) :: fo_linmod
      integer :: bottom
      real(kind=fpk) :: fraction, radius, boa, altitude, theta, alpha, phi
      real(kind=fpk) :: cosscat, theta_boa, alpha_boa, phi_boa

      ! MS uses the observer's local geometry. Its spherical FO supplement
      ! must use that same ray, not the same numeric angle at the surface.
      fo_fix = VLIDORT_FixIn
      fo_mod = VLIDORT_ModIn
      fo_linmod = VLIDORT_LinModIn
      fo_fix%Bool%TS_DO_FULLRAD_MODE = .false.
      fo_fix%Bool%TS_DO_PLANE_PARALLEL = .false.
      fo_fix%Bool%TS_DO_SPHERICITY_CORRECTION = .false.
      fo_fix%Bool%TS_DO_MSSTS = .false.
      fo_mod%MBool%TS_DO_FOCORR_EXTERNAL = .false.

      if (viewing_zenith_angle <= 90.0_fpk) then
         bottom = ceiling(observer_level)
         if (bottom == 0) then
            VLIDORT_Sup%SS%TS_STOKES_SS = ZERO
            VLIDORT_Sup%SS%TS_STOKES_DB = ZERO
            call VLIDORT_LinSup_Init(VLIDORT_LinSup)
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL = .true.
            VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .false.
            return
         endif
         fraction = observer_level - real(bottom-1,kind=fpk)
         fo_fix%Cont%TS_NLAYERS = bottom
         fo_fix%Cont%TS_NLAYERS_NOMS = bottom
         fo_fix%Cont%TS_NLAYERS_CUTOFF = bottom
         fo_fix%Chapman%TS_HEIGHT_GRID(bottom) = obs_height_clamped / 1000.0_fpk
         fo_fix%Optical%TS_DELTAU_VERT_INPUT(bottom) = &
            fraction * fo_fix%Optical%TS_DELTAU_VERT_INPUT(bottom)
         ! Fractional optical derivatives retain the original layer's
         ! normalization, so the result is still dI/d(original layer OD).
         fo_mod%MUserVal%TS_USER_LEVELS(1) = real(bottom,kind=fpk)
         fo_mod%MUserVal%TS_GEOMETRY_SPECHEIGHT = fo_fix%Chapman%TS_HEIGHT_GRID(bottom)
      else
         radius = planet_radius
         boa = layer_boundary_altitudes(n_layers) / 1000.0_fpk
         altitude = obs_height_clamped / 1000.0_fpk
         theta = solar_zenith_angle; alpha = solver_vza; phi = relative_azimuth_angle
         if ((radius+altitude)*sin(alpha*acos(-ONE)/180.0_fpk) >= radius+boa) then
            write(*,*) 'Standard spherical FO: limb ray misses surface; use VLIDORT-QS.'
            stop 6
         endif
         call AC_to_BOA_geometry_conversion(radius, boa, altitude, theta, alpha, phi, &
            cosscat, theta_boa, alpha_boa, phi_boa)
         fo_mod%MSunrays%TS_SZANGLES(1) = theta_boa
         fo_mod%MUserVal%TS_USER_VZANGLES_INPUT(1) = alpha_boa
         fo_mod%MUserVal%TS_USER_RELAZMS(1) = phi_boa
         fo_mod%MUserVal%TS_USER_OBSGEOMS_INPUT(1,:) = (/theta_boa, alpha_boa, phi_boa/)
      endif

      if (linearized) then
         call VLIDORT_LPS_MASTER(.false., fo_fix, fo_mod, VLIDORT_Sup, VLIDORT_Out, &
            VLIDORT_LinFixIn, fo_linmod, VLIDORT_LinSup, VLIDORT_LinOut)
      else
         call VLIDORT_MASTER(.false., fo_fix, fo_mod, VLIDORT_Sup, VLIDORT_Out)
      endif
      if (VLIDORT_Out%Status%TS_STATUS_INPUTCHECK /= VLIDORT_SUCCESS .or. &
          VLIDORT_Out%Status%TS_STATUS_CALCULATION /= VLIDORT_SUCCESS) then
         write(*,*) 'Observer-referenced spherical first-order calculation failed.'
         write(*,*) trim(VLIDORT_Out%Status%TS_CHECKMESSAGES(1))
         write(*,*) trim(VLIDORT_Out%Status%TS_MESSAGE)
         stop 6
      endif
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL = .true.
      ! The spherical FO is already in the supplement; do not recompute it.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .false.
   end subroutine prepare_standard_first_order

   subroutine configure_radiance_geometry(diridx)
      integer, intent(out) :: diridx

      call reset_single_geometry()
      VLIDORT_FixIn%Bool%TS_DO_SPHERICITY_CORRECTION = do_sphericity_correction
      VLIDORT_FixIn%Bool%TS_DO_MSSTS = do_sphericity_correction
      if (do_aircraft_correction) then
         if (observer_level <= ZERO .or. observer_level >= real(n_layers,kind=fpk) .or. &
             abs(observer_level-real(nint(observer_level),kind=fpk)) > 1.0d-8) then
            write(*,*) 'Aircraft correction requires an interior layer boundary; use VLIDORT-QS.'
            stop 6
         endif
         if ((planet_radius+obs_height_clamped/1000.0_fpk)*sin(solver_vza*acos(-ONE)/180.0_fpk) &
             >= planet_radius+layer_boundary_altitudes(n_layers)/1000.0_fpk) then
            write(*,*) 'Aircraft correction ray misses the surface; use VLIDORT-QS.'
            stop 6
         endif
      endif
      if (viewing_zenith_angle <= 90.0_fpk) then
         VLIDORT_FixIn%Bool%TS_DO_UPWELLING = .false.
         VLIDORT_FixIn%Bool%TS_DO_DNWELLING = .true.
         diridx = DNIDX
         if (do_sphericity_correction) then
            if (do_aircraft_correction) then
               VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = &
                  real(nint(observer_level), kind=fpk)
            else
               VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = real(n_layers, kind=fpk)
            endif
         endif
      else
         VLIDORT_FixIn%Bool%TS_DO_UPWELLING = .true.
         VLIDORT_FixIn%Bool%TS_DO_DNWELLING = .false.
         diridx = UPIDX
         if (do_sphericity_correction) then
            if (do_aircraft_correction) then
               VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = &
                  real(nint(observer_level), kind=fpk)
            else
               VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = ZERO
            endif
         endif
      endif
   end subroutine configure_radiance_geometry

   subroutine configure_flux_geometry(diridx)
      integer, intent(out) :: diridx

      call reset_single_geometry()
      ! LOS interpolation applies to radiances, not angularly integrated fluxes.
      VLIDORT_FixIn%Bool%TS_DO_SPHERICITY_CORRECTION = .false.
      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .false.
      VLIDORT_FixIn%Bool%TS_DO_UPWELLING = .true.
      VLIDORT_FixIn%Bool%TS_DO_DNWELLING = .true.
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = observer_level
      diridx = UPIDX
   end subroutine configure_flux_geometry

   subroutine reset_single_geometry()
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_EXTERNAL = .false.
      VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = .true.
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES = 1
      VLIDORT_ModIn%MSunrays%TS_SZANGLES = ZERO
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(1) = solar_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = 1
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS = ZERO
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1) = relative_azimuth_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 1
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT = ZERO
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1) = solver_vza
      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 1
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT = ZERO
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1) = solar_zenith_angle
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,2) = solver_vza
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,3) = relative_azimuth_angle
      VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS = 0
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS = ZERO
      VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1) = observer_level
   end subroutine reset_single_geometry

   function height_to_user_level(height_m) result(user_level)
      real(kind=fpk), intent(in) :: height_m
      real(kind=fpk) :: user_level
      integer :: k

      user_level = real(n_layers, kind=fpk)
      do k = 1, n_layers
         if (height_m <= layer_boundary_altitudes(k-1) .and. height_m >= layer_boundary_altitudes(k)) then
            user_level = real(k-1, kind=fpk) + &
               (layer_boundary_altitudes(k-1) - height_m) / &
               (layer_boundary_altitudes(k-1) - layer_boundary_altitudes(k))
            return
         endif
      enddo
   end function height_to_user_level

   function pressure_at_boundary(idx) result(value)
      integer, intent(in) :: idx
      real(kind=fpk) :: value
      if (idx == 0) then
         value = pressure_profile(1)
      else if (idx >= n_layers) then
         value = pressure_profile(n_layers)
      else
         value = 0.5_fpk * (pressure_profile(idx) + pressure_profile(idx+1))
      endif
   end function pressure_at_boundary

   function temperature_at_boundary(idx) result(value)
      integer, intent(in) :: idx
      real(kind=fpk) :: value
      if (idx == 0) then
         value = temperature_profile(1)
      else if (idx >= n_layers) then
         value = temperature_profile(n_layers)
      else
         value = 0.5_fpk * (temperature_profile(idx) + temperature_profile(idx+1))
      endif
   end function temperature_at_boundary

   subroutine decode_backend_mode(mode, spherical, aircraft)
      character(len=*), intent(inout) :: mode
      logical, intent(out) :: spherical, aircraft

      spherical = .false.
      aircraft = .false.
      mode = adjustl(mode)
      aircraft = index(mode, 'aircraft_spherical') > 0
      spherical = index(mode, 'spherical') > 0
      sphercorr_points = max(2, sphercorr_points)
   end subroutine decode_backend_mode

end program DRIVER
