MODULE RayleighScattering_m

IMPLICIT NONE


PUBLIC :: RayleighXs          ! Calculate Rayleigh cross section
PUBLIC :: AirRefractiveIndex       ! Calculate refractive index of moist air
PUBLIC :: RayleighGreekCoeffs    ! Rayleigh phase matrix expansion coefficients

!####################
!  Some constants  #
!####################

integer, parameter, private :: fpk = Selected_Real_Kind(15)
real(kind=fpk), PARAMETER, private :: PI  = 3.141592653589793  
real(kind=fpk), PARAMETER, private :: k_B = 1.3806488d-23    ! Boltzmann constant
real(kind=fpk), PARAMETER, private :: h   = 6.62606957d-34   ! Planck constant
real(kind=fpk), PARAMETER, private :: c   = 299792458        ! Speed of light

CONTAINS

FUNCTION RayleighXS(wavelength, pressure, temperature, water_vapor_pressure)

   ! ######################################################################################################
   ! # Calculate Rayleigh cross section for given wavelength, pressure, temperature and relative humidity #
   ! #                                      Adapted from Bucholtz [1995]                                  #
   ! ######################################################################################################

    real(kind=fpk) RayleighXs ! RAYLEIGH CROSS SECTION [cm^2]

   ! # Input parameter 

    real(kind=fpk), INTENT(IN) :: wavelength   ! WAVELENGTH [nm]
    real(kind=fpk), INTENT(IN) :: pressure   ! AIR PRESSURE [hPa]
    real(kind=fpk), INTENT(IN) :: temperature   ! TEMPERATURE [K]
    real(kind=fpk), INTENT(IN) :: water_vapor_pressure   ! WATER VAPOUR PRESSURE [hPa]
    
   ! # Local variables
   
    real(kind=fpk) :: Rho ! Air density [molec/m^3]
    real(kind=fpk) :: K   ! King correction factor
    real(kind=fpk) :: N   ! Refractive index of air
    
    ! Air density
    Rho = (pressure * 1d2)/(k_B * temperature)
    
    ! King correction factor
    K = 1.046d0 + 5.036d-4 * (1.0d3/wavelength)**2 + 2.896d-5 * (1.0d3/wavelength)**4
   
    ! Air refractive index
    N = AirRefractiveIndex(wavelength, pressure, temperature, water_vapor_pressure)
    
    RayleighXS = 1d0/Rho**2 * 24d0 * PI**3                 &
                 / (wavelength*1d-9)**4 * ((N**2 - 1d0)    &
                 / (N**2 + 2d0))**2 * K * 1d4

END FUNCTION RayleighXS


! ########################################################################################


FUNCTION AirRefractiveIndex(wavelength, pressure, temperature, water_vapor_pressure)

   ! #####################################
   ! # Calculate refractive index of air #
   ! #    Adapted from Edlen [1966]      #
   ! #         and Ciddor [1996]         #
   ! #####################################

   real(kind=fpk) :: AirRefractiveIndex 

   ! Input parameter 

   real(kind=fpk), INTENT(IN) :: wavelength   ! WAVELENGTH [nm]
   real(kind=fpk), INTENT(IN) :: pressure   ! AIR PRESSURE [hPa]
   real(kind=fpk), INTENT(IN) :: temperature   ! TEMPERATURE [K]
   real(kind=fpk), INTENT(IN) :: water_vapor_pressure   ! WATER VAPOUR PRESSURE [hPa]
    
   ! Local variables
   real(kind=fpk) :: W2 

   W2 = (1d3/wavelength)**2 ! square of wavenumber in micron

   ! Calculate refractive index of dry air at 15°C  
   AirRefractiveIndex = 5792105d0/(238.0185d0 - W2) + 167917d0/(57.362d0 - W2)
   AirRefractiveIndex = AirRefractiveIndex*1d-8 + 1d0
                       
   ! Calculate refractive index at given pressure and temperature
   AirRefractiveIndex =  1 + (AirRefractiveIndex - 1d0)   &
      * 0.00104126d0 * pressure        & 
      / (1d0 + 0.003671d0 * (temperature - 273.15d0)) 

   ! Calculate refractive index of moist air   
   AirRefractiveIndex = AirRefractiveIndex &
      - water_vapor_pressure * (4.292d0 - 0.0341d0 * W2) * 1d-8

END FUNCTION AirRefractiveIndex


! ########################################################################################

FUNCTION RayleighGreekCoeffs(wavelength)

   ! ##############################################
   ! # Calculate Rayleigh scattering phase matrix #
   ! #    adapted from Landgraf et al. [2004]     #
   ! ##############################################

   real(kind=fpk) :: RayleighGreekCoeffs(0:2,6) ! Dimensions: moment order, matrix element

   ! Input parameter 

    real(kind=fpk), INTENT(IN) :: wavelength         ! Wavelength[nm]
  
   ! Local variables 

   real(kind=fpk) :: K      ! King correction factor
   real(kind=fpk) :: depol  ! Depolarisation factor
   real(kind=fpk) :: beta2
   real(kind=fpk) :: gamma

   K = 1.046d0 + 5.036d-4 * (1d3/wavelength)**2 + 2.896d-5 * (1d3/wavelength)**4
   depol = (6d0*K - 6d0) / (7d0*K + 3)
   
   RayleighGreekCoeffs = 0.0d0
   
   beta2 = ( 1D0 - depol ) / ( 2.0d0 + depol )
   gamma = -sqrt(6.0d0) * beta2
   
   ! Apply same order as for coefficients from mie model output: 1 to 4 is diagnoal, 4 & 5 are off-diagonal elements
   RayleighGreekCoeffs(0,1) =  1d0
   RayleighGreekCoeffs(2,1) =  beta2
   RayleighGreekCoeffs(2,2) =  6.0d0 * beta2
   RayleighGreekCoeffs(1,4) =  3.D0 * ( 1D0 - 2.D0*DEPOL ) / (2D0 + depol )
   RayleighGreekCoeffs(2,5) =  gamma

END FUNCTION RayleighGreekCoeffs


END MODULE RayleighScattering_m



