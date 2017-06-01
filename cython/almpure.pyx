import numpy as np
cimport numpy as np
import cython
cimport cython

cimport almpure

@cython.boundscheck(False)
@cython.wraparound(False)
def map2alm(np.ndarray[double, ndim=2] hmap, np.ndarray[double, ndim=1] apmask,
            int lmax, int mmax, np.ndarray[double, ndim=1] qwghts):

    cdef int npix = hmap.shape[0]
    cdef int nstokes = hmap.shape[1]
    cdef int nmaps = 1

    cdef int nside = np.sqrt(npix/12)

    cdef np.ndarray aphmap = np.copy(hmap, order='F')
    aphmap = hmap * apmask[:,np.newaxis]
    cdef double[::1,:] aphmap_view = aphmap

    cdef np.ndarray alms = np.zeros([nstokes, lmax+1, mmax+1], order='F',
                                    dtype=[('re',np.double),('im',np.double)])
    cdef s2hat_dcomplex[::1,:,:] alms_view = alms

    cdef double[:] qwghts_view = qwghts

    map2alm_c(&alms_view[0,0,0], nstokes, lmax, mmax, nmaps,
              &aphmap_view[0,0], nside, &qwghts_view[0])

    return alms

@cython.boundscheck(False)
@cython.wraparound(False)
def map2almpure(np.ndarray[double, ndim=2] hmap, np.ndarray[double, ndim=1] apmask,
            int lmax, int mmax, np.ndarray[double, ndim=1] qwghts):

    cdef int npix = hmap.shape[0]
    cdef int nstokes = hmap.shape[1]
    cdef int nmaps = 1

    cdef int nside = np.sqrt(npix/12)

    cdef double[::1,:] hmap_view = hmap
    cdef double[::1] apmask_view = apmask

    cdef np.ndarray alms = np.zeros([nstokes, lmax+1, mmax+1], order='F',
                                    dtype=[('re',np.double),('im',np.double)])
    cdef s2hat_dcomplex[::1,:,:] alms_view = alms

    cdef double[:] qwghts_view = qwghts

    map2almpure_c(&alms_view[0,0,0], nstokes, lmax, mmax, nmaps,
              &hmap_view[0,0], &apmask_view[0], nside, &qwghts_view[0])

    return alms

@cython.boundscheck(False)
@cython.wraparound(False)
def alm2map(int nside, np.ndarray[s2hat_dcomplex, ndim=3] alms):

    cdef int npix = 12 * nside * nside
    cdef int nstokes = alms.shape[0]
    cdef int nlmax = alms.shape[1] - 1
    cdef int nmmax = alms.shape[2] - 1
    cdef int nmaps = 1

    cdef s2hat_dcomplex[::1,:,:] alms_view = alms

    cdef np.ndarray hmap = np.zeros([npix,nstokes], dtype=np.double, order='F')
    cdef double[::1,:] hmap_view = hmap

    alm2map_c(&hmap_view[0,0], nside, &alms_view[0,0,0], nstokes, nlmax, nmmax, nmaps)

    return hmap

