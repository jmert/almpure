
module s2hat_types_internal

  ! definitions needed for s2hat package
  ! Radek Stompor (APC) September, 2008

  use s2hat_defs

#ifdef HEALPIX_fft
     use healpix_types
     use healpix_fft
     use misc_utils
#endif  

  implicit none

#ifdef FFTW3_C2R
     include "fftw3.f"
#endif

#ifdef FFTW3_HC2R
     include "fftw3.f"
#endif

  include "mpif.h"

  ! HEALPIX-like type defs
  
#ifndef HEALPIX_fft
    INTEGER, PARAMETER :: i8b = s2hat_i8b
    INTEGER, PARAMETER :: i4b = s2hat_i4b
    INTEGER, PARAMETER :: i2b = s2hat_i2b
    INTEGER, PARAMETER :: i1b = s2hat_i1b
    INTEGER, PARAMETER :: sp  = s2hat_sp
    INTEGER, PARAMETER :: dp  = s2hat_dp
    INTEGER, PARAMETER :: lgt = s2hat_lgt
    INTEGER, PARAMETER :: spc = s2hat_spc
    INTEGER, PARAMETER :: dpc = s2hat_dpc

    ! and constants

    REAL(kind=dp), PARAMETER :: HALFPI= 1.570796326794896619231321691639751442099_dp
    REAL(kind=dp), PARAMETER :: PI    = 3.141592653589793238462643383279502884197_dp
    REAL(kind=dp), PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394_dp
    REAL(kind=dp), PARAMETER :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
    REAL(kind=dp), PARAMETER :: max_dp  = HUGE(1.0_dp)
    real(kind=dp), parameter ::  KvS = 1.0_dp ! 1.0 : CMBFAST (Healpix 1.2)
#endif
  
  ! define large and small numbers used to renormalise the recursion 
  ! on the Legendre Polynomials
  integer(I4B),      public, parameter :: LOG2LG   = 100
  real(KIND=DP),     public, parameter :: FL_LARGE = 2.0_dp **   LOG2LG
  real(KIND=DP),     public, parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
  ! declare array for dynamic rescaling of the Ylm
  integer(kind=i4b), public, parameter :: RSMAX = 20, RSMIN = -20
  real(dp),          public, dimension(RSMIN:RSMAX) :: rescale_tab
  real(DP),          public, parameter :: ALN2_INV = 1.4426950408889634073599246810_dp ! 1/log(2)

  integer(kind=i4b), public, parameter :: HPX_MXL0 = 40 ! minimum used, choose <=0 to do full loop
  real   (kind=dp),  public, parameter :: HPX_MXL1 = 1.35_dp

  type( s2hat_pixeltype), target :: globalPixelization
  type( s2hat_scandef), target :: globalScan

  integer(s2hat_i4b) :: fft_mc_flag=0   ! no MC by default

end module s2hat_types_internal
