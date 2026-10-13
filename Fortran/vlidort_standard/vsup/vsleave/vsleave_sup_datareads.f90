
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
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
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2024.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            get_foQ_data                                     #
! #            get_fluorescence_data                            #
! #            get_solarspec_irradiance                         #
! #                                                             #
! ###############################################################

Module vsleave_sup_datareads_m

use VLIDORT_PARS_m, Only : fpk, Zero

!  Created 3/28/23 (VLIDORT 2.8.3 Upgrade). Three subroutines are PUBLIC

public

contains

subroutine get_foQ_data ( fOQ_File, fOQ_Full, fail, message )

!  3/28/23. This routine just obtains the data for later use.
!     -- Adapted for tool use by R. Spurr.
!     -- Output data array is double precision

   implicit none

!  Input

   character(LEN=*), intent(in) :: fOQ_File ! file containing foQ full dadta

!  Output (pre-initialized)

   real(fpk), intent(inout)   :: fOQ_Full(7,6,6,17,13)

!  Exception  handling

   logical      , intent(out) :: fail
   character*(*), intent(out) :: message

!  Local

   integer   :: i, j, j1, k, m, m1
   real(fpk) :: lam, sun, pigc, fval, FoQ_Table(17,13)

!  Initialize

   fOQ_Full = Zero
   fail = .false. ; message = ' '

!  Read table, store Full Table only.

   open(45,file=Trim(fOQ_File),err=88,status='old',action='read')
   do i = 1, 7 !lam
     do j = 1, 6 !sun
       do k = 1, 6 !pigc
         read(45,*)lam, sun, pigc, fval ; j1 = 7-j
         do m = 1, 17
           m1 = 18-m ; read(45,*)FoQ_Table(m,1:13)
           fOQ_Full(i,j1,k,m1,1:13) = FoQ_Table(m,1:13) 
         enddo
       enddo
     enddo
   enddo
   close(45)
   return

!  Error return

88 continue
   close(45) ; fail = .true.
   message = 'Openfile error in get_foQ_data ; file not found: ' // Trim(FoQ_File)
   return

end subroutine get_foQ_data

!

subroutine get_fluorescence_data & 
     ( NLAT_FLUO_FILE, NLON_FLUO_FILE, fluo_file, use_nag_compiler, &
       lat_data, lon_data, fluo_data, fail, message )

!  Function from Chris O'Dell, 12 July 12 ==> fluorescence data at 755 nm in W/m^2/um/sr
!    Minor reprogramming for the SL supplement by R. Spurr, 12 July 2012

!  3/28/23. This routine just obtains the data for later use.
!     -- Adapted for tool use by R. Spurr.
!     -- Output data and lat/lon grids are now both double precision (M. Christi)

      implicit none

!  Input

      integer         , intent(in) :: NLAT_FLUO_FILE, NLON_FLUO_FILE
      logical         , intent(in) :: use_nag_compiler
      character(LEN=*), intent(in) :: fluo_file ! file containing fluorescence climatology

!  Output (pre-initialized)

      real(fpk), intent(inout) :: lat_data (NLAT_FLUO_FILE)
      real(fpk), intent(inout) :: lon_data (NLON_FLUO_FILE)
      real(fpk), intent(inout) :: fluo_data(NLON_FLUO_FILE, NLAT_FLUO_FILE, 12)

!  Exception  handling

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Local variables

      integer, parameter, dimension(12) :: DAYS_IN_MONTH = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

      integer :: i, j, mon, k, kmax, FUNIT, ios
      real    :: fluo_data_in(NLON_FLUO_FILE, NLAT_FLUO_FILE, 12) !3-D
      logical :: verbose

!  data array for use with NAG compiler read

      real    :: fluo_data_temp(NLON_FLUO_FILE*NLAT_FLUO_FILE*12) !1-D

!  initialize

      fail = .false. ; message = ' '
      verbose = .false. !.true.
      FUNIT = 1001

!  Set Lat/Lon Grids

      do i = 1, NLAT_FLUO_FILE
        lat_data(i) = 90.0d0 - dble(i-1) * 0.5d0
      enddo
      do j = 1, NLON_FLUO_FILE
        lon_data(j) = -180.0d0 + 0.5d0 * dble(j-1)
      enddo

!  Open File of Fluorescence data

      if (verbose) write(*,*) 'Fluo File = ' // trim(fluo_file)
      open(FUNIT, file=trim(fluo_file), form='UNFORMATTED', status='OLD', IOSTAT=ios)
      if (ios /=0) then
         close(FUNIT) ; message = 'Error Opening Fluorescence file ' // trim(fluo_file)
         fail = .true. ; return
      endif

!  Read file using normal compiler

      if (.not. use_nag_compiler) then
         read(FUNIT, IOSTAT=ios) fluo_data_in
         fluo_data = real(fluo_data_in,fpk)
         if (ios /=0) then
           close(FUNIT) ; fail = .true.
           message = 'Error Reading Normal Fluorescence file ' // trim(fluo_file)
         else
           close(FUNIT)
         endif
      endif

!  Read file using NAG-compiler
!    - prepare to read from array "fluo_data_temp" starting at position 2 (not 1!) since NAG reads
!      the binary file record header as a data point

      if (use_nag_compiler) then
        k = 1 ; kmax = NLON_FLUO_FILE*NLAT_FLUO_FILE*12
        read(FUNIT, IOSTAT=ios) fluo_data_temp
        if (ios /=0) then
          close(FUNIT) ; fail = .true.
          message = 'Error NAG-compiler Reading Fluorescence file ' // trim(fluo_file)
        else
          do mon=1,12 ; do i=1,NLAT_FLUO_FILE ; do j=1,NLON_FLUO_FILE
            k = k + 1
            if (k <= kmax) then
              fluo_data(j,i,mon) = real(fluo_data_temp(k),fpk)
            else
              fluo_data(j,i,mon) = 0.0_fpk
            end if
          enddo ;  enddo ; enddo
          close(FUNIT)
        endif
      endif

!  done

      return
end subroutine get_fluorescence_data

!

subroutine get_solarspec_irradiance &
         ( maxsolarfiledata, solarfile_choice, obs_period, solarfile_1, solarfile_2, &
           numsolardata, wavdata, solardata, fail, message )  

!Reads a solar spectral irradiance file and returns the solar
!spectral irradiance for an air mass of zero in
!units of W m^-2 μm^-1 for the wavelength input

! Function by Mick Christi - 16 July 2012
! Add filenames externally, Rob Spurr 10/5/15

!  3/28/23. This routine just obtains the data for later use.
!     -- Adapted for tool use by R. Spurr.
!     -- Output data are double precision

!  solar_file = 1  1985 Wehrli Standard Extraterrestrial Spectral Solar Irradiance  (TSI = 1367 W/m2)
!  solar_file = 2  Solar Spectral Irradiance Reference Spectra for Whole Heliosphere Interval (WHI) 2008
!                     (TSI = ?; however, may be normalized by user)

      implicit none

!  Input

      integer         , intent(in) :: maxsolarfiledata
      integer         , intent(in) :: solarfile_choice, obs_period
      character(LEN=*), intent(in) :: solarfile_1, solarfile_2 ! files containing solar data

!  Output (pre-initialized)

      integer  , intent(inout) :: numsolardata
      real(fpk), intent(inout) :: wavdata  (maxsolarfiledata)
      real(fpk), intent(inout) :: solardata(maxsolarfiledata)

!  Exception handling

      logical      , intent(out) :: fail
      character*(*), intent(out) :: message

!  Local variables

      !Output variable

      integer           :: i, funit
      real(fpk)         :: wav,data_src,norm_factor,solar_integ_irad_in, solar_spec_irad_in(3)
      logical           :: normalize
      character*250     :: solarfile

!  Initialize.

      numsolardata = 0 ; wavdata = Zero ; solardata = Zero
      fail = .false. ; message = ' '

!  Obtain solar data npoints

      if (solarfile_choice == 1) then
        numsolardata = 920
      else if (solarfile_choice == 2) then
        numsolardata = 24000
      else
        message = 'Solar file choice not recognized - abort' ; fail = .true.
      endif

!  Unit number, normalization factor for second solar file (set)

      FUNIT = 50
      normalize = .true. !.false.

!  First choice ( 1985 Wehrli Standard Extraterrestrial Spectral Solar Irradiance file )
!     Read solar data and convert (1) wavelength from nm to um
!                                 (2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1

      if ( solarfile_choice == 1) then
         solarfile = Trim(solarfile_1)
         open(unit=Funit, file=Trim(solarfile),status='old',err=96,action='read')
         do i=1,5
            read(Funit,*)
         enddo
         do i = 1, numsolardata
            read(Funit,*) wav, solar_spec_irad_in(1), solar_integ_irad_in
            wavdata(i)   = wav * 1.0d-3
            solardata(i) = solar_spec_irad_in(1) * 1.0d3
         enddo
         close(Funit) ; return
      endif

!  Second choice ( Solar Spectral Irradiance Reference Spectra for Whole Heliosphere Interval (WHI) 2008 file)
!     Note: observation period for this data set -
          ! obs_period = 1  Moderately low solar activity with sunspot darkening.   TSI = 1360.696 W/m^2
          ! obs_period = 2  Moderately low solar activity with faculae brightening. TSI = 1360.944 W/m^2
          ! obs_period = 3  Very close to solar cycle minimum condition.            TSI = 1360.840 W/m^2

!     Read solar data and convert (1) wavelength from nm to um
!                                 (2) spectral irradiance data from W m^-2 nm^-1 to W m^-2 μm^-1
!     Create normalization factor to adjust radiances to correspond to a TSI of 1366.1 W/m^2

      if ( solarfile_choice == 2) then
         solarfile = Trim(solarfile_2)
         open(unit=funit,file=Trim(solarfile),status='old',err=96,action='read')
         do i=1,144
            read(Funit,*)
         enddo
         if (obs_period == 1) then
            norm_factor = 1366.1d0/1360.696d0
         else if (obs_period == 2) then
            norm_factor = 1366.1d0/1360.944d0
         else if (obs_period == 3) then
            norm_factor = 1366.1d0/1360.840d0
         end if
         do i=1,numsolardata
            read(funit,'(F8.2,3E12.4,F4.0)') wav,solar_spec_irad_in(1:3),data_src
            wavdata(i)   = wav * 1.0d-3
            solardata(i) = solar_spec_irad_in(obs_period) * 1.0d3
            if (normalize) then
              solardata(i) = solardata(i)*norm_factor
            end if
         enddo
         close(funit) ; return
      end if

!  error returns

96    continue
      close(funit)
      message = 'Following solar file not found ; '//Trim(solarfile)
      fail = .true. ; return

end subroutine get_solarspec_irradiance

!  End module

END MODULE vsleave_sup_datareads_m

