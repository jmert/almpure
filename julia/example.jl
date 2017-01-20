using s2hat

MPI.Init()

myrank = Cint(MPI.Comm_rank(MPI.COMM_WORLD))
nprocs = Cint(MPI.Comm_size(MPI.COMM_WORLD))

nside = 128
pixelization = set_pixelization(s2hat.PIXCHOICE_HEALPIX, PixelParameters(nside,0))
scan = zbounds2scan([-1.0,1.0], pixelization)

nlmax = Cint(700)
nmmax = Cint(700)
nstokes = Cint(3)
delta = 80
seed = 0

sizes = get_local_data_sizes(Cint(0), pixelization, scan,
            nlmax, nmmax, myrank, nprocs, nstokes,
            MPI.COMM_WORLD)
nmvals = sizes[:nmvals]

nrings = sizes[:last_ring] - sizes[:first_ring] + 1
nalms  = (nlmax+1)*nmvals*nstokes

mvals = find_mvalues(myrank, nprocs, nmmax, nmvals);

# Make a delta map at ell == 80

if myrank == 0
    # In HEALPIX conventions:
    #   Dim 1 -> Stokes parameters (should be size 1 or 3)
    #   Dim 2 -> ell (size nlmax)
    #   Dim 3 -> m (size nmmax)
    alms = zeros(Complex128, nstokes, nlmax+1, nmmax+1);

    srand(seed)
    # m == 0 is special and must be real valued
    alms[1:nstokes,delta,1] = complex(randn(nstokes), 0)
    for mm in 1:delta
        alms[1:nstokes,delta,mm+1] = complex.(randn(nstokes)/sqrt(2), randn(nstokes)/sqrt(2))
    end
else
    alms = C_NULL
end

# Local alms adds on extra dimension for the number of maps. We're only using
# one here, so it'll be a singleton
local_alms = zeros(Complex128, nstokes, nlmax+1, nmmax+1, 1);
distribute_alms(nlmax, nmmax, Cint(1), Cint(0), nstokes, nmvals, mvals,
    nstokes, local_alms, alms, myrank, nprocs, Cint(0), MPI.COMM_WORLD)

local_map = zeros(Cdouble, sizes[:map_size], nstokes, 1);
alm2map(Cint(0), pixelization, scan, nlmax, nmmax, nmvals, mvals, Cint(1),
    nstokes, sizes[:first_ring], sizes[:last_ring], sizes[:map_size],
    local_map, nstokes, local_alms, Clonglong(0), C_NULL, nprocs, myrank,
    MPI.COMM_WORLD)

if myrank == 0
    map = zeros(Cdouble, pixelization.npixsall, nstokes);
else
    map = C_NULL
end
collect_map(pixelization, Cint(1), Cint(0), nstokes, map, sizes[:first_ring],
    sizes[:last_ring], sizes[:map_size], local_map, myrank, nprocs, Cint(0),
    MPI.COMM_WORLD)

nph = unsafe_wrap(Array, pixelization.nph, (pixelization.nringsall,), false);
fpix = unsafe_wrap(Array, pixelization.fpix, (pixelization.nringsall,), false);
kphi = unsafe_wrap(Array, pixelization.kphi, (pixelization.nringsall,), false);
cth = unsafe_wrap(Array, pixelization.cth, (pixelization.nringsall,), false);

destroy_scan(scan)
destroy_pixelization(pixelization)
MPI.finalize()

