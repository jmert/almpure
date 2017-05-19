
program test

  use s2hat_defs
  use s2hat_toolbox_mod
  use s2hat_pixelization
  use s2hat_map2alm_mod
  use s2hat_alm2map_mod

  implicit none

  ! include "mpif.h"

  integer(s2hat_i4b) :: myid, ierr, numprocs, root, comm
  integer(s2hat_i4b) :: nside, npix, lmax, nmaps, i, j, l, m, nstokes
  integer(s2hat_i4b) :: nmmax, nmvals, first_ring, last_ring, map_size, lda, nrings
  integer(s2hat_i4b) :: nringsall, precompute_plms, mapnum, nspec
  integer(s2hat_i8b) :: nplm
  character(len=128) :: param
  real(s2hat_dp)     :: t1, t2, t3, t4
  real(s2hat_dp),     pointer,     dimension(:,:)   :: local_w8ring, local_plms, map, mapout, w8ring
  real(s2hat_dp),     pointer,     dimension(:,:)   :: cls  
  real(s2hat_dp),     pointer,     dimension(:,:,:) :: local_map
  complex(s2hat_dpc), pointer,     dimension(:,:,:) :: alms
  complex(s2hat_dpc), pointer,     dimension(:,:,:,:) :: local_alms
  integer(s2hat_i4b), dimension(MPI_STATUS_SIZE)    :: status
  integer(s2hat_i4b), pointer,     dimension(:)     :: mask, mvals
  real(s2hat_dp),                  dimension(2)     :: zbounds
  type(s2hat_pixparameters) :: pixpars
  type(s2hat_pixeltype)     :: pixelization
  type(s2hat_scandef)       :: scan

  call cpu_time(t1)

  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  root = 0
  comm = MPI_COMM_WORLD

  nside = 128
  
  lmax            = 2*nside
  npix            = 12*nside**2
  nmaps           = 1
  nstokes         = 3
  lda             = lmax         ! s2hat convention for alm matrix arrangement
  ! lda             = nstokes    ! that would be HEALPix like
  precompute_plms = 0            ! Plms are not to be precomputed
  mapnum          = 1            ! # of maps to be processed simultanously
  nmmax           = lmax

  zbounds(:) = 0.d0              ! full map to be processed - no equatorial cut

  pixpars%par1 = nside           ! for HEALPIX the second field of pixpars is irrelevant

  call set_pixelization(PIXCHOICE_HEALPIX, pixpars, pixelization)  ! precomputes the pixelization structure assuming the HEALPIX pixelization
  
  nringsall = pixelization%nringsall  ! # of rings in the Northern hemisphere plus equator

  call zbounds2scan(zbounds, pixelization, scan)  ! computes the scan structure given 'zbounds'

  ! computes the object sizes needed for memory allocation etc
  call get_local_data_sizes(precompute_plms, pixelization, scan, lmax, nmmax, myid, numprocs, nmvals, &
          & first_ring, last_ring, map_size, nplm, root, comm)

  ! allocate storage for the values of the m numbers processed on this proc
  allocate(mvals(nmvals))

  call find_mvalues(myid, numprocs, nmmax, nmvals, mvals)  ! define the m-values subset

  nrings = last_ring-first_ring+1 ! # of equatorial rings o be processed by this proc

  allocate(local_w8ring(nrings,nstokes))    ! to store the ring quadrature weights - one per ring
  allocate(local_map(0:map_size-1,nstokes,nmaps))    ! to store the submap processed by this proc

  local_w8ring(1:nrings,1:nstokes) = 1.d0   ! i.e. no ringw weighting

  do i = 0, map_size-1
     local_map(i,1:nstokes,1:nmaps) = pixelization%fpix(first_ring)+i   ! fake map
  end do

  if( lda == nstokes) then
      allocate(local_alms(nstokes,0:lmax,0:nmvals-1,nmaps))    ! healpix like arrangement
  else
      allocate(local_alms(0:lmax,0:nmvals-1,nstokes,nmaps))    ! s2hat native arrangement
  end if     

  local_alms = 0.d0

  ! map-to-alm transform

  call s2hat_map2alm(0, pixelization, scan, lmax, nmmax, nmvals, mvals, nmaps, nstokes, &
          & first_ring, last_ring, local_w8ring, map_size, local_map, lda, local_alms, &
          & nplm, local_plms, numprocs, myid, comm)

  nspec = 4  ! get TT, EE, BB, TE spectra only

  if( myid == root) allocate( cls(0:lmax,1:nspec))                     ! only on the proc 'root'

  call collect_cls( nmaps, 0, nstokes, lmax, nmvals, mvals, lmax, local_alms,  &
                  & nspec, cls, myid, numprocs, root, comm)

  ! alm-to-map
     
  call s2hat_alm2map(precompute_plms, pixelization, scan, lmax, nmmax, nmvals, &
          & mvals, nmaps, nstokes, first_ring, last_ring, map_size, local_map, lda, local_alms, &
          & nplm, local_plms, numprocs, myid, comm)

  deallocate( local_map)
  deallocate( local_alms)
  deallocate( local_w8ring)
  deallocate( mvals)

  if( myid == root) deallocate( cls)

  ! And exit
  call mpi_finalize(ierr)

end program test
