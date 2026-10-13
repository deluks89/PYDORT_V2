module RTSMie_master_trimodal_m

!  This is the Trimodal master #1 for RTS MIE code

!  history of trimodal and Trimodal
!    ** RT Solutions, Version 1.4, 30 June     2011 (trimodal control Tmatrix)
!    ** RT Solutions, Version 1.5, 17 August   2011 (trimodal control Mie)
!    ** RT Solutions, 03 March 2016. Trimodal development, UCLA for ATTREX project

  use RTSMie_parameters_m
  use RTSMie_sourcecode_m

!  Everything PUBLIC here
!  ----------------------

public

contains

subroutine RTSMie_master_trimodal &
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index, PSD_pars, MonoRadius, R1, R2, FixR1R2,   & ! I
         nblocks, nweights, xparticle_limit, R1R2_cutoff,    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real, n_imag, Trifractions, Do_Save_PSD,               & ! I
         TriMie_bulk, TriMie_asymm, TriMie_ncoeffs,          & ! O
         TriMie_expcoeffs, TriMie_Fmatrix, TriMie_dist,      & ! O
         fail, istatus, Trimessages, trace_3 )                 ! O

   USE RTSMie_parameters_m

!  implicit none statement

   IMPLICIT NONE

!  List of Inputs
!  ==============

!  Flag inputs
!  -----------

!      Do_Expcoeffs      - Boolean flag for computing Expansion Coefficients
!      Do_Fmatrix        - Boolean flag for computing F-matrix at equal-angles

   logical  , intent(in)  :: Do_Expcoeffs
   logical  , intent(in)  :: Do_Fmatrix

!      Do_Monodisperse   - Boolean flag for Doing a Monodisperse calculation
!                          If set, the PSD stuff will be turned off internally

   LOGICAL  , INTENT (IN) :: do_Monodisperse

!  PSD inputs
!  ----------

!  FixR1R2 : If  set, Use "R1R2_cutoff" for smallest particle size
!                     then Internal routine to calculate R1 and R2 (outputs)
!            If Not set, Use Input R1 and R2 for PSD limits.

   logical, intent(inout)  :: FixR1R2(3)

!  R1, R2         - Minimum and Maximum radii (Microns)

   real    (KIND=dp), intent(inout)  :: R1(3), R2(3)

!  Limiting particle size value. Set to 10000.0 default.
!   If you exceed this, program will tell you to increase dimensioning.

   REAL    (KIND=dp), INTENT (IN) :: xparticle_limit

!      Monoradius     - Monodisperse radius size (Microns)

   real    (KIND=dp), intent(in)  :: Monoradius

!      psd_Index      - Index for particle size distribution of spheres
!      psd_pars       - Parameters characterizing PSD (up to 3 allowed)

!  Mie inputs (distribution index, PSD parameters)
!    PSD_index = 1 : TWO-PARAMETER GAMMA with alpha and b given
!    PSD_index = 2 : TWO-PARAMETER GAMMA with par(1)= reff and par(2)= veff given
!    PSD_index = 3 : trimodal GAMMA with equal mode weights
!    PSD_index = 4 : LOG-NORMAL with rg and sigma given
!    PSD_index = 5 : LOG-NORMAL with reff and veff given
!    PSD_index = 6 : POWER LAW
!    PSD_index = 7 : MODIFIED GAMMA with alpha, rc and gamma given
!    PSD_index = 8 : MODIFIED GAMMA with alpha, b and gamma given

   integer          , intent(in)  :: PSD_Index(3)
   real    (KIND=dp), intent(in)  :: PSD_pars (3,3)

!  PSD quadrature control
!  ----------------------

!    PSD integration ranges is divided into so many blocks.
!    For each block, integrate using Gaussian quadrature, with so many weights.
!    R1R2_cutoff particle size for setting R1 and R2 internally

   INTEGER          , INTENT (IN) :: nblocks(3)
   INTEGER          , INTENT (IN) :: nweights(3)
   REAL    (KIND=dp), INTENT (IN) :: R1R2_cutoff(3)

!  Optical: Wavelength, refractive index
!  -------------------------------------

!      LAMBDA         - wavelength of light (microns)
!      N_REAL, N_IMAG - real and imaginary parts, refractive index (N-i.GE.0)

   real    (KIND=dp), intent(in)  :: lambda, n_real(3), n_imag(3)

!  Trimodal Fraction (by volume weight). Must add up to 1 (this is checked)

   real    (KIND=dp), intent(in)  :: Trifractions(3)
   
!  If saving PSD data to file

   logical    :: Do_Save_PSD

!  F-matrix Angular control input
!  ------------------------------

!  Calculate F-matrix at user-defined angles (do_Fmatrix flag MUST BE set)
!       n_Fmatrix_angles = number of user-defined angles.
!       Fmatrix_angles   = user-defined angles, in DEGREES between [0, 180]

   INTEGER          , INTENT (IN) :: n_Fmatrix_angles
   REAL    (KIND=dp), INTENT (IN) :: Fmatrix_angles(max_Mie_angles)

!  Output arguments
!  ================

!  Bulk distribution parameters
!    1 = Extinction coefficient
!    2 = Scattering coefficient
!    3 = Single scattering albedo

   real(KIND=dp), intent(out) :: TriMie_bulk (3)

!  Expansion coefficients and Asymmetry parameter, optional output

   integer      , intent(out) :: TriMie_ncoeffs
   real(KIND=dp), intent(out) :: TriMie_expcoeffs (6,0:max_Mie_angles)
   real(KIND=dp), intent(out) :: TriMie_asymm

!  F-matrix, optional output

   real(KIND=dp), intent(out) :: TriMie_Fmatrix(4,max_Mie_angles)

!  Distribution parameters. BOTH MODES !!!!!!!!!!
!    1 = Normalization
!    2 = Cross-section
!    3 = Volume
!    4 = REFF
!    5 = VEFF

   real(KIND=dp), intent(out) :: TriMie_dist (5,3)

!  Exception handling

   LOGICAL          , INTENT (OUT)   :: fail
   INTEGER          , INTENT (OUT)   :: istatus
   CHARACTER*(*)    , INTENT (OUT)   :: Trimessages(3), trace_3

!  Local Arrays
!  ------------

!  Bulk distribution parameters
!  Expansion coefficients and Asymmetry parameter
!  F-matrix,  optional output

   real(KIND=dp) :: Mie1_bulk (3,3)
   integer       :: Mie1_ncoeffs(3)
   real(KIND=dp) :: Mie1_expcoeffs (6,0:max_Mie_angles,3)
   real(KIND=dp) :: Mie1_asymm(3)
   real(KIND=dp) :: Mie1_Fmatrix (4,max_Mie_angles,3)

!  Other local variables
!  ---------------------

   integer       :: k, k12, k13, k2, L, nmommin3, nmommin2
   real(KIND=dp) :: FF1, FF2, FF3, WW1, WW2, WW3, WW(3)
   real(KIND=dp) :: Csca1_FF1, Csca2_FF2, Csca3_FF3, Csca_total
   
!  Output variable to save to file and plot size distribution
!  -------------------

   REAL    (KIND=dp)  :: Total_Radius(nblocks(1), max_Mie_sizes), Total_PSD(nblocks(1), max_Mie_sizes)
   character*100      :: filename, X1
   
!  Zero the output
!  ---------------

   TriMie_bulk      = d_zero
   TriMie_Fmatrix   = d_zero
   TriMie_expcoeffs = d_zero
   TriMie_asymm     = d_zero
   TriMie_ncoeffs   = 0
   TriMie_dist      = d_zero

   FF1 = Trifractions(1)
   FF2 = Trifractions(2)
   FF3 = Trifractions(3)

   fail           = .false.
   istatus        = 0
   Trimessages(1:3) = ' '
   trace_3        = ' '

!  Check fractions
!  ---------------

   if ( abs(sum(Trifractions(1:3)))- d_one .gt. 1.0d-12 ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_trimodal module: Three fractions do not add up to 1.0 !!'
      return
   endif

!  Check: No Monodisperse here !
!  -----------------------------

   if ( do_monodisperse ) then
      fail = .true.; istatus = 1
      trace_3 = 'RTSMie_master_trimodal module: Input error: MONODISPERSE FLAG must be Off!'
      return
   endif

!  Three Calls
!  -----------

   do k = 1, 3
     
      write(*,'(A,i2)')' ** Doing Trimodal Mie for PSD/Mode # ',k

      CALL RTSMie_main                                          & 
       ( Do_Expcoeffs, Do_Fmatrix, do_Monodisperse,          & ! I
         PSD_Index(k), PSD_pars(:,k), MonoRadius,            & ! I
         R1(k), R2(k), FixR1R2(k), nblocks(k), nweights(k),  & ! I
         xparticle_limit, R1R2_cutoff(k),                    & ! I
         n_Fmatrix_angles, Fmatrix_angles,                   & ! I
         lambda, n_real(k), n_imag(k),                       & ! I
         Mie1_bulk(:,k), Mie1_asymm(k), Mie1_ncoeffs(k),     & ! O
         Mie1_expcoeffs(:,0:,k), Mie1_Fmatrix(:,:,k), TriMie_dist(:,k),  & ! O
         fail, istatus, Trimessages(1), Trimessages(2), Trimessages(3) , & ! O
         Total_Radius, Total_PSD)   ! O
         
         if (Do_Save_PSD) then
             write (X1,'(i1)') k ! converting integer to string using a 'internal file'
             filename='SizeInfo/output_PSD'//trim(X1)//'.dat'
             open(1,file=Trim(filename),status = 'unknown')
             write(1,*)Total_Radius
             write(1,*)Total_PSD
         endif

!  Exception handling
        
      if ( fail ) then
         if (k.eq.1)trace_3 = 'RTSMie_master_trimodal module: First  PSD call, Error'
         if (k.eq.2)trace_3 = 'RTSMie_master_trimodal module: Second PSD call, Error'
         if (k.eq.3)trace_3 = 'RTSMie_master_trimodal module: Third  PSD call, Error'
         if ( Istatus .eq. 1 ) return
      endif

   enddo

!  trimodal determination
!  ---------------------

!  Revision 20 September 2011 (Bimodal)
!    Correct definition for Expcoeffs/Fmatrix: WW1/WW2 in place of FF1/FF2

   Csca1_FF1  = FF1 * Mie1_bulk(2,1)   
   Csca2_FF2  = FF2 * Mie1_bulk(2,2) 
   Csca3_FF3  = FF3 * Mie1_bulk(2,3)   

   Csca_total =  Csca1_FF1 + Csca2_FF2 + Csca3_FF3
   WW1   = Csca1_FF1 / Csca_total ; WW(1) = WW1
   WW2   = Csca2_FF2 / Csca_total ; WW(2) = WW2
   WW3   = Csca3_FF3 / Csca_total ; WW(3) = WW3

!  @@@ Rob Fix 21 Sep 12, combined Single-scatter-albedo was wrong (BIMODAL)
!  original code
!   TriMie_bulk(1:3) = FF1 * Mie1_bulk(1:3) + FF2 * Mie2_bulk(1:3)

   TriMie_bulk(1:2) = FF1 * Mie1_bulk(1:2,1) + FF2 * Mie1_bulk(1:2,2) + FF3 * Mie1_bulk(1:2,3)
   TriMie_bulk(3)   = TriMie_bulk(2) / TriMie_bulk(1)

!write(0,*) TriMie_bulk(1:3)

!  set filter for computing moments

   nmommin3 = minval(Mie1_ncoeffs)
   if ( Mie1_ncoeffs(1) .eq. nmommin3 ) then
      nmommin2 = min(Mie1_ncoeffs(2),Mie1_ncoeffs(3)) ;  k12 = 2 ; k13 = 3
      k2 = 2 ; if ( Mie1_ncoeffs(2).eq. nmommin2 ) k2 = 3
   else if ( Mie1_ncoeffs(2) .eq. nmommin3 ) then
      nmommin2 = min(Mie1_ncoeffs(1),Mie1_ncoeffs(3)) ;  k12 = 1 ; k13 = 3
      k2 = 1 ; if ( Mie1_ncoeffs(1).eq. nmommin2 ) k2 = 3
   else if ( Mie1_ncoeffs(3) .eq. nmommin3 ) then
      nmommin2 = min(Mie1_ncoeffs(1),Mie1_ncoeffs(2)) ;  k12 = 1 ; k13 = 2
      k2 = 1 ; if ( Mie1_ncoeffs(1).eq. nmommin2 ) k2 = 2
   endif

!  Moments; use the filter

   if ( Do_Expcoeffs ) then
      TriMie_asymm = WW1 * Mie1_asymm(1) + WW2 * Mie1_asymm(2) + WW3 * Mie1_asymm(3)
      TriMie_ncoeffs = maxval(Mie1_ncoeffs)      
      do L = 0, nmommin3
         TriMie_expcoeffs(1:6,L) = WW1 * Mie1_expcoeffs(1:6,L,1) + & 
                                   WW2 * Mie1_expcoeffs(1:6,L,2) + & 
                                   WW3 * Mie1_expcoeffs(1:6,L,3)
      enddo
      do L = nmommin3 + 1,nmommin2
         TriMie_expcoeffs(1:6,L) = WW(k12) * Mie1_expcoeffs(1:6,L,k12) + WW(k13) * Mie1_expcoeffs(1:6,L,k13)
      enddo
      do L = nmommin2 + 1,Mie1_ncoeffs(k2)
         TriMie_expcoeffs(1:6,L) = WW(k2) * Mie1_expcoeffs(1:6,L,k2)
      enddo
   endif

!  F-matrix

   if ( Do_Fmatrix ) then
      do L = 1, n_Fmatrix_angles
         TriMie_Fmatrix(1:4,L) = WW1 * Mie1_Fmatrix(1:4,L,1) + & 
                                 WW2 * Mie1_Fmatrix(1:4,L,2) + & 
                                 WW3 * Mie1_Fmatrix(1:4,L,3)
      enddo
   endif   

!  Finish

   return
end subroutine RTSMie_master_trimodal

!  End module

end module RTSMie_master_trimodal_m
