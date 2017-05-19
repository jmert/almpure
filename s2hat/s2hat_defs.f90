
module s2hat_defs

  ! definitions needed for s2hat package
  ! Radek Stompor (APC) September, 2008

  ! HEALPIX-like type defs
  
  INTEGER, PARAMETER, public :: s2hat_i8b = SELECTED_INT_KIND(16)
  INTEGER, PARAMETER, public :: s2hat_i4b = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER, public :: s2hat_i2b = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER, public :: s2hat_i1b = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER, public :: s2hat_sp  = SELECTED_REAL_KIND(5,30)
  INTEGER, PARAMETER, public :: s2hat_dp  = SELECTED_REAL_KIND(12,200)
  INTEGER, PARAMETER, public :: s2hat_lgt = KIND(.TRUE.)
  INTEGER, PARAMETER, public :: s2hat_spc = KIND((1.0_s2hat_sp, 1.0_s2hat_sp))
  INTEGER, PARAMETER, public :: s2hat_dpc = KIND((1.0_s2hat_dp, 1.0_s2hat_dp))  

  ! define pixel choices here - HEALPIX the only one so far ...
  integer(s2hat_i4b), parameter :: PIXCHOICE_HEALPIX = 0  ! for HEALPIX pixelization
  integer(s2hat_i4b), parameter :: PIXCHOICE_GLESP = 1    ! for GLESP pixelization
  integer(s2hat_i4b), parameter :: PIXCHOICE_ECP = 2      ! for ECP pixelization
  integer(s2hat_i4b), parameter :: PIXCHOICE_GLCP = 3     ! for GLCP pixelization

  ! spin convention 
  integer(s2hat_i4b), parameter :: SPIN_CONV_SIGN = -1    ! for Healpix like choice

  type s2hat_pixparameters
    integer(s2hat_i4b) :: par1
    integer(s2hat_i4b) :: par2
  end type s2hat_pixparameters

  type s2hat_pixeltype                        ! pixelization is assumed to be iso-latitudal with pixels evenly spaced for each latitude
		     			      ! and a symmetry present between northern and southern hemispheres.
    ! pixel types defined as above - 8 byte ints to avoid padding problems in the C-interface
    integer(s2hat_i8b) :: type
    ! a total number of pixels
    integer(s2hat_i8b) :: npixsall
    ! a total number of iso-latitude rings in the north hemisphere (including equator)
    integer(s2hat_i8b) :: nringsall
    ! a total maximum number of pixels per iso-ring
    integer(s2hat_i8b) :: nphmx
    ! a number of the first pixel for each ring
    integer(s2hat_i8b), dimension(:), pointer :: fpix  ! ranges [1:nringsall]
    ! a number of pixel for each iso-ring
    integer(s2hat_i8b), dimension(:), pointer :: nph   ! ranges [1:nringsall]
    ! an offset of the 1st pixel for each ring with respect to a meridian 0 (radians)
    real(s2hat_dp), dimension(:), pointer :: kphi  ! ranges [1:nringsall]
    ! quadrature weights for each ring
    real(s2hat_dp), dimension(:), pointer :: qwght  ! ranges [1:nringsall]    
    ! pixel center separations for each iso-ring
    real(s2hat_dp), dimension(:), pointer :: pixphi    ! ranges [1:nringsall]
    ! pixel area (assumed to be constant) for each iso-ring
    real(s2hat_dp), dimension(:), pointer :: parea     ! ranges [1:nringsall]
    ! cosines of the polar angle for each iso-ring
    real(s2hat_dp), dimension(:), pointer :: cth       ! ranges [1:nringsall]
    ! sines of the polar angle for each iso-ring (redundant)
    real(s2hat_dp), dimension(:), pointer :: sth       ! ranges [1:nringsall]
  end type s2hat_pixeltype
  
  type s2hat_scandef     ! defines the scan parameters needed for the s2hat transforms
     ! a total number of observed pixels ! all 8byte int to simplify the C-interfacing
     integer(s2hat_i8b) :: npixsobs
     ! a total number of observed rings
     integer(s2hat_i8b) :: nringsobs
     ! observed flags:
     ! northern hemisphere (includes equator)
     integer(s2hat_i8b), dimension(:), pointer :: nfl ! ranges [1:nringsall]
     ! southern hemisphere 
     integer(s2hat_i8b), dimension(:), pointer :: sfl ! ranges [1:nringsall]
     ! either northern or southern (redundant)
     integer(s2hat_i8b), dimension(:), pointer :: fl  ! ranges [1:nringsall]
  end type s2hat_scandef

end module s2hat_defs
