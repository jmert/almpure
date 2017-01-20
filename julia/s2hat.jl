module s2hat
    using MPI
    export PixelParameters, PixelType, ScanDef
    export
       set_pixelization, destroy_pixelization,
       zbounds2scan, zbounds2mask, mask2scan, destroy_scan,
       get_local_data_sizes, find_mvalues, nummvalues, nummmodes,
       distribute_map, distribute_w8ring, distribute_mask, collect_map,
       distribute_alms, collect_alms,
       alm2map, map2alm, map2almpure

    include("./libpath.jl")
    const libs2hat = S2HAT_LIBPATH

    # From s2hat_defs.h

    const PIXCHOICE_HEALPIX = 0
    const PIXCHOICE_GLESP   = 1
    const PIXCHOICE_ECP     = 2
    const PIXCHOICE_GLCP    = 3

    const SPIN_CONV_SIGN = -1

    """
    Structure storing parameters defining any of the pixelization/sky gridding
    schemes pre-defined in S2HAT.
    """
    immutable PixelParameters
        par1::Int32
        par2::Int32
    end

    type CPixelType
        typeid::Int64
        npixsall::Int64
        nringsall::Int64
        nphmx::Int64
        fpix::Ptr{Int64}
        nph::Ptr{Int64}
        kphi::Ptr{Float64}
        qwght::Ptr{Float64}
        pixphi::Ptr{Float64}
        parea::Ptr{Float64}
        cth::Ptr{Float64}
        sth::Ptr{Float64}
    end
    type CScanDef
        npixsobs::Int64
        nringobs::Int64
        nfl::Ptr{Int64}
        sfl::Ptr{Int64}
        fl::Ptr{Int64}
    end


    """
    A structure storing all the pixelization characteristics needed by the
    transforms.
    """
    type PixelType
        """Pixelization type (e.g. `PIXCHOICE_HEALPIX`)"""
        typeid::Int64

        """Total number of pixels in the full sky map"""
        npixsall::Int64

        """Total number of iso-latitude rings covering the northern hemisphere
        plus equatorial ring (if present)"""
        nringsall::Int64

        """Maximum number of pixels per iso-ring"""
        nphmx::Int64

        """Vector containing, for each ring, the first pixel number in the
        global, ring-based pixel numbering scheme."""
        fpix::Vector{Int64}

        """Vector containing the number of pixels for each iso-ring."""
        nph::Vector{Int64}

        """Vector containing the phase shift info for the FFTs, i.e.  the
        azimuthal position (in radians) of the first pixel of each ring."""
        kphi::Vector{Float64}

        """Vector containing the quadrature weights for each iso-ring."""
        qwght::Vector{Float64}

        """Vector containing the pixel azimuthal separation for each iso-ring.
        **UNUSED:** The pixel separation is **always** assumed to be
        `2π/nph[ring number]`."""
        pixphi::Vector{Float64}

        """Vector containing the pixel areas for each iso-ring (and assumed to
        be constant for each iso-ring)."""
        parea::Vector{Float64}

        """Vector containing a value of the cosine for the polar angle measured
        from the North Pole for each iso-ring."""
        cth::Vector{Float64}

        """Vector containing a value of the sine for the polar angle measured
        from the North Pole for each iso-ring."""
        sth::Vector{Float64}
    end

    """
    A S2HAT structure storing all the info needed to define the sky patch to be
    analyzed.
    """
    type ScanDef
        """The number of observed pixels."""
        npixsobs::Int64

        """The number of observed rings, including all partially observed ones,
        observed in either the North or South. This should be the number of
        non-zero entries of the vector `scan.fl`."""
        nringobs::Int64

        """Vector storing info about each ring of the Northern hemisphere (plus
        the equatorial ring, if present). The stored values are either `1` or
        `0` if a ring is or is not observed. The rings are counted from the
        North Pole."""
        nfl::Vector{Int64}

        """Vector storing info about each ring of the Southern hemisphere. The
        stored values are either `1` or `0` if a ring is or is not observed.
        The rings are counted from the North Pole."""
        sfl::Vector{Int64}

        """Vector storing the info about each ring. The stored values are `1`
        if either a ring or its symmetric counterpart (or both) is observed,
        otherwise is `0`. This is the logical sum of `nfl` and `sfl`."""
        fl::Vector{Int64}
    end

    CPixelType() = CPixelType(-1,0,0,0, C_NULL,C_NULL,
            C_NULL,C_NULL,C_NULL,C_NULL,C_NULL,C_NULL)
    CScanDef() = CScanDef(0,0, C_NULL,C_NULL,C_NULL)

    # Provide conversions from C types to Julia types and back again.
    function PixelType(P::CPixelType)
        const N = P.nringsall
        return PixelType(
                P.typeid,
                P.npixsall,
                P.nringsall,
                P.nphmx,
                unsafe_wrap(Array, P.fpix,   (N,), false),
                unsafe_wrap(Array, P.nph,    (N,), false),
                unsafe_wrap(Array, P.kphi,   (N,), false),
                unsafe_wrap(Array, P.qwght,  (N,), false),
                unsafe_wrap(Array, P.pixphi, (N,), false),
                unsafe_wrap(Array, P.parea,  (N,), false),
                unsafe_wrap(Array, P.cth,    (N,), false),
                unsafe_wrap(Array, P.sth,    (N,), false)
            )
    end
    function CPixelType(P::PixelType)
        const N = P.nringsall
        return CPixelType(
                P.typeid,
                P.npixsall,
                P.nringsall,
                P.nphmx,
                pointer(P.fpix,  ),
                pointer(P.nph,   ),
                pointer(P.kphi,  ),
                pointer(P.qwght, ),
                pointer(P.pixphi,),
                pointer(P.parea, ),
                pointer(P.cth,   ),
                pointer(P.sth,   )
            )
    end

    function ScanDef(S::CScanDef, P::Union{PixelType,CPixelType})
        const N = P.nringsall
        return ScanDef(
                S.npixsobs,
                S.nringobs,
                unsafe_wrap(Array, S.nfl, (N,), false),
                unsafe_wrap(Array, S.sfl, (N,), false),
                unsafe_wrap(Array, S.fl,  (N,), false)
            )
    end

    function set_pixelization(typeid::Int, param::PixelParameters)
        pix = Ref(CPixelType())
        ccall((:set_pixelization, libs2hat), Void,
            (Cint, PixelParameters, Ref{CPixelType}),
            typeid, param, pix)
        return pix[]
    end

    function destroy_pixelization(pixtype::CPixelType)
        ccall((:destroy_pixelization, libs2hat), Void,
            (CPixelType,), pixtype)
        return nothing
    end

    function zbounds2scan(bounds::Vector{Float64}, pixelization::CPixelType)
        length(bounds) == 2 || throw(DimensionMismatch("bounds must be length 2"))
        scan = Ref{CScanDef}(CScanDef())
        ccall((:zbounds2scan, libs2hat), Void,
            (Ref{Float64}, CPixelType, Ref{CScanDef}),
            bounds, pixelization, scan)
        return scan[]
    end

    function zbounds2mask(bounds::Vector{Float64}, pixelization::CPixelType)
        length(bounds) == 2 || throw(DimensionMismatch("bounds must be length 2"))
        mask = Vector{Int32}(pixelization.npixsall)
        ccall((:zbounds2mask, libs2hat), Void,
            (Ref{Float64}, CPixelType, Ref{Int32}),
            bounds, pixelization, mask)
        return mask
    end

    function mask2scan(mask::Vector{Int32}, pixelization::CPixelType, scan::Ref{CScanDef})
        ccall((:mask2scan, libs2hat), Void,
            (Ref{Int32}, CPixelType, Ref{ScanDef}),
            mask, pixelization, scan)
    end

    function destroy_scan(scan::CScanDef)
        ccall((:destroy_scan, libs2hat), Void,
            (CScanDef,), scan)
        return nothing
    end

    function fft_setup(pixelization::CPixelType, nmmax, opt)
        ccall((:fft_setup, libs2hat), Void,
            (CPixelType, Int32, Int32), pixelization, nmmax, opt)
    end

    function fft_setup(pixelization::CPixelType, nmmax)
        ccall((:fft_mc_setup, libs2hat), Void,
            (CPixelType, Int32), pixelization, nmmax)
    end

    function fft_clean()
        ccall((:fft_mc_clean, libs2hat), Void, ())
    end

    function get_local_data_sizes(precompute_plms::Int32,
            pixelization::CPixelType, scan::CScanDef, nlmax::Int32,
            nmmax::Int32, myid::Int32, numprocs::Int32, root::Int32,
            comm::MPI.Comm)
        nmvals = Ref{Int32}(0)
        first_ring = Ref{Int32}(0)
        last_ring = Ref{Int32}(0)
        map_size = Ref{Int32}(0)
        nplm = Ref{Int64}(0)
        ccall((:get_local_data_sizes, libs2hat), Void,
            (Int32, CPixelType, CScanDef, Int32, Int32, Int32, Int32,
            Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Int64}, Int32,
            MPI.CComm), precompute_plms, pixelization, scan, nlmax, nmmax,
            myid, numprocs, nmvals, first_ring, last_ring, map_size, nplm,
            root, comm)
        return Dict(
            :nmvals => nmvals[],
            :first_ring => first_ring[],
            :last_ring => last_ring[],
            :map_size => map_size[],
            :nplm => nplm[])
    end

    function find_mvalues(myid::Int32, numprocs::Int32, nmmax::Int32,
            nmvals::Int32)
        mvals = Vector{Int32}(nmvals)
        ccall((:find_mvalues, libs2hat), Void,
            (Int32, Int32, Int32, Int32, Ref{Int32}),
            myid, numprocs, nmmax, nmvals, mvals)
        return mvals
    end

    function nummvalues(myid::Int32, numprocs::Int32, nmmax::Int32)
        return ccall((:nummvalues, libs2hat), Int32,
            (Int32, Int32, Int32),
            myid, numprocs, nmmax)
    end

    function nummmodes(lmax::Int32, nmvals::Int32, mvals::Vector{Int32})
        return ccall((:nummmodes, libs2hat), Int32,
            (Int32, Int32, Ref{Int32}),
            lmax, nmvals, mvals)
    end

    function find_scan_ring_range(pixelization::CPixelType, scan::ScanDef,
            mmax::Int32, myid::Int32, numprocs::Int32)
        first_ring = Ref{Int32}(0)
        last_ring = Ref{Int32}(0)
        outerror = Ref{Int32}(0)
        ccall((:find_scan_ring_range, libs2hat), Void,
            (CPixelType, CScanDef, Int32, Int32, Int32, Ref{Int32},
            Ref{Int32}, Ref{Int32}),
            pixelization, scan, mmax, myid, numprocs, first_ring, last_ring,
            outerror)
        return Dict(
                :first_ring => first_ring[],
                :last_ring  => last_ring[],
                :outerror   => outerror[])
    end

    function distribute_map(pixelization::CPixelType, nmaps::Int32,
            mapnum::Int32, nstokes::Int32, first_ring::Int32,
            last_ring::Int32, map_size::Int32, local_map::Array{Float64,3},
            map::Union{Ptr{Void},Matrix{Float64}}, myid::Int32,
            numprocs::Int32, root::Int32, comm::MPI.Comm)
        ccall((:distribute_map, libs2hat), Void,
            (CPixelType, Int32, Int32, Int32, Int32, Int32, Int32,
            Ref{Float64}, Ptr{Void}, Int32, Int32, Int32, MPI.CComm),
            pixelization, nmaps, mapnum, nstokes, first_ring, last_ring,
            map_size, local_map, map, myid, numprocs, root, comm)
    end

    function distribute_w8ring(npol::Int32, first_ring::Int32,
            last_ring::Int32, local_w8ring::Matrix{Float64}, nringsall::Int32,
            w8ring::Union{Ptr{Void},Matrix{Float64}}, myid::Int32,
            numprocs::Int32, root::Int32, comm::MPI.Comm)
        ccall((:distribute_w8ring, libs2hat), Void,
            (Int32, Int32, Int32, Ref{Float64}, Int32, Ptr{Void}, Int32,
            Int32, Int32, MPI.CComm),
            npol, first_ring, last_ring, local_w8ring, nringsall, w8ring,
            myid, numprocs, root, comm)
    end

    function distribute_mask(pixelization::CPixelType, nmasks::Int32,
            masknum::Int32, nstokes::Int32, first_ring::Int32,
            last_ring::Int32, mask_size::Int32, local_mask::Array{Float64,3},
            mask::Union{Ptr{Void},Matrix{Float64}}, myid::Int32,
            numprocs::Int32, root::Int32, comm::MPI.Comm)
        ccall((:distribute_mask, libs2hat), Void,
            (CPixelType, Int32, Int32, Int32, Int32, Int32, Int32,
            Ref{Float64}, Ptr{Void}, Int32, Int32, Int32, MPI.CComm),
            pixelization, nmasks, masknum, nstokes, first_ring, last_ring,
            mask_size, local_mask, mask, myid, numprocs, root, comm)
    end

    function collect_map(pixelization::CPixelType, nmaps::Int32, mapnum::Int32,
            nstokes::Int32, map::Union{Ptr{Void},Matrix{Float64}},
            first_ring::Int32, last_ring::Int32, map_size::Int32,
            local_map::Array{Float64,3}, myid::Int32, numprocs::Int32,
            root::Int32, comm::MPI.Comm)
        ccall((:collect_map, libs2hat), Void,
            (CPixelType, Int32, Int32, Int32, Ptr{Void}, Int32, Int32, Int32,
            Ref{Float64}, Int32, Int32, Int32, MPI.CComm),
            pixelization, nmaps, mapnum, nstokes, map, first_ring, last_ring,
            map_size, local_map, myid, numprocs, root, comm)
    end

    function collect_alms(lmax::Int32, mmax::Int32, nmaps::Int32,
            mapnum::Int32, nstokes::Int32, nmvals::Int32,
            mvals::Vector{Int32}, lda::Int32, local_alm::Array{Complex128,4},
            alms::Union{Ptr{Void},Array{Complex128,3}}, myid::Int32,
            numprocs::Int32, root::Int32, comm::MPI.Comm)
        ccall((:collect_alms, libs2hat), Void,
            (Int32, Int32, Int32, Int32, Int32, Int32, Ref{Int32}, Int32,
            Ref{Complex128}, Ptr{Void}, Int32, Int32, Int32, MPI.CComm),
            lmax, mmax, nmaps, mapnum, nstokes, nmvals, mvals, lda, local_alm,
            alms, myid, numprocs, root, comm)
    end

    function distribute_alms(nlmax::Int32, nmmax::Int32, nmaps::Int32,
            mapnum::Int32, nstokes::Int32, nmvals::Int32, mvals::Vector{Int32},
            lda::Int32, local_alm::Array{Complex128,4},
            alms::Union{Ptr{Void},Array{Complex128,3}}, myid::Int32,
            numprocs::Int32, root::Int32, comm::MPI.Comm)
        ccall((:distribute_alms, libs2hat), Void,
            (Int32, Int32, Int32, Int32, Int32, Int32, Ref{Int32}, Int32,
            Ref{Complex128}, Ptr{Void}, Int32, Int32, Int32, MPI.CComm),
            nlmax, nmmax, nmaps, mapnum, nstokes, nmvals, mvals, lda, local_alm,
            alms, myid, numprocs, root, comm)
    end

    function alm2map(precompute_plms::Int32, pixelization::CPixelType,
            scan::CScanDef, nlmax::Int32, nmmax::Int32, nmvals::Int32,
            mvals::Vector{Int32}, nmaps::Int32, nstokes::Int32,
            first_ring::Int32, last_ring::Int32, map_size::Int32,
            local_map::Array{Float64,3}, lda::Int32,
            local_alm::Array{Complex128,4}, nplm::Int64,
            local_plm::Union{Ptr{Void},Matrix{Float64}}, numprocs::Int32,
            myid::Int32, comm::MPI.Comm)
        ccall((:s2hat_alm2map, libs2hat), Void,
            (Int32, CPixelType, CScanDef, Int32, Int32, Int32, Ref{Int32},
            Int32, Int32, Int32, Int32, Int32, Ref{Float64}, Int32,
            Ref{Complex128}, Int64, Ptr{Void}, Int32, Int32, MPI.CComm),
            precompute_plms, pixelization, scan, nlmax, nmmax, nmvals, mvals,
            nmaps, nstokes, first_ring, last_ring, map_size, local_map, lda,
            local_alm, nplm, local_plm, numprocs, myid, comm)
    end

    function map2alm(precompute_plms::Int32, pixelization::CPixelType,
            scan::CScanDef, lmax::Int32, mmax::Int32, nmvals::Int32,
            mvals::Vector{Int32}, nmaps::Int32, nstokes::Int32,
            first_ring::Int32, last_ring::Int32,
            local_w8ring::Matrix{Float64}, map_size::Int32,
            local_map::Array{Float64,3}, lda::Int32,
            local_alm::Array{Complex128,4}, nplm::Int64,
            local_plm::Union{Ptr{Void},Matrix{Float64}}, numprocs::Int32,
            myid::Int32, comm::MPI.Comm)
        ccall((:s2hat_map2alm, libs2hat), Void,
            (Int32, CPixelType, CScanDef, Int32, Int32, Int32, Ref{Int32},
            Int32, Int32, Int32, Int32, Ref{Float64}, Int32, Ref{Float64},
            Int32, Ref{Complex128}, Int64, Ptr{Void}, Int32, Int32,
            MPI.CComm),
            precompute_plms, pixelization, scan, lmax, mmax, nmvals, mvals,
            nmaps, nstokes, first_ring, last_ring, local_w8ring, map_size,
            local_map, lda, local_alm, nplm, local_plm, numprocs, myid, comm)
    end

    function map2purealm(pixelization::CPixelType, scan::CScanDef,
            lmax::Int32, mmax::Int32, nmvals::Int32, mvals::Vector{Int32},
            nmaps::Int32, first_ring::Int32, last_ring::Int32,
            local_w8ring::Vector{Int32}, mapsize::Int32, maskstride::Int32,
            local_mask::Matrix{Int32}, winstride::Int32,
            local_swindow::Matrix{Float64}, local_cmbmap::Array{Float64,3},
            lda::Int32, local_apurelm::Array{Complex128,4}, my_rank::Int32,
            numprocs::Int32, my_comm::MPI.Comm)::Int32
        return ccall((:s2hat_map2purealm, libs2hat), Int32,
            (CPixelType, CScanDef, Int32, Int32, Int32, Ref{Int32}, Int32,
            Int32, Int32, Ref{Int32}, Int32, Int32, Ref{Int32}, Int32,
            Ref{Int32}, Ref{Int32}, Int32, Ref{Complex128}, Int32, Int32,
            MPI.CComm),
            pixelization, scan, lmax, mmax, nmvals, mvals, nmaps, first_ring,
            last_ring, local_w8ring, mapsize, maskstride, local_mask,
            winstride, local_swindow, local_cmbmap, lda, local_apurelm,
            my_rank, numprocs, my_comm)
    end
end

