#line 1 "s2hat_map2almpure_c.c" /* workaround for mex compiler using full path */
/* These are needed before s2hat.h can be included since it doesn't take care
 * of loading it's own requirements, apparently... */
#include <sys/types.h>
#include <mpi.h>
/* */
#include <s2hat.h>
#include <s2hat_pure.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "debug.h"

#include <mex.h>
#include <matrix.h>

int map2almpure(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
                int nmaps, double* map, double* apmask, int nside,
                double* qwghts);

/*
 * alms = map2almpure(map, apmask, lmax, mmax);
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    /* Inputs from Matlab */
    const mxArray* ml_map = NULL;
    const mxArray* ml_apmask = NULL;
    const mxArray* ml_lmax = NULL;
    const mxArray* ml_mmax = NULL;
    const mxArray* ml_qwghts = NULL;

    /* Outputs back to Matlab */
    mxArray* ml_alms = NULL;

    /* Intermediaries between Matlab and s2hat */
    mwSize tmp_ndim = 0;
    const mwSize* tmp_dims = NULL;
    double* ml_alms_r = NULL;
    double* ml_alms_i = NULL;
    mwSize ml_alms_dims[4] = { 0 };

    /* Internal s2hat interfaces */
    const int nstokes = 2; /* always Q/U -> E/B for map2almpure */
    int lmax = 0;
    int mmax = 0;
    int nmaps = 0;
    int nside = 0;
    int npix = 0;
    s2hat_dcomplex* alms = NULL;
    double* map = NULL;
    double* apmask = NULL;
    double* qwghts = NULL;

    /* Initialize MPI if necessary */
    mexCallMATLAB(0, NULL, 0, NULL, "mpihelper");

    /* Validate MATLAB inputs */
    dbglog("Validating MATLAB inputs...\n");
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("map2almpure:args:nrhs",
                "Five input arguments are required");
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("map2almpure:args:nlhs",
                "One output argument is required");
    }

    ml_map = prhs[0];
    if (!mxIsDouble(ml_map) || mxIsComplex(ml_map)) {
        mexErrMsgIdAndTxt("map2almpure:args:notDouble",
                "Input map must be real and of type double");
    }
    tmp_ndim = mxGetNumberOfDimensions(ml_map);
    if (tmp_ndim != 2 && tmp_ndim != 3) {
        mexErrMsgIdAndTxt("map2almpure:args:dims",
                "Input map must be have 2 or 3 dimensions");
    }

    ml_apmask = prhs[1];
    if (!mxIsDouble(ml_apmask) || mxIsComplex(ml_apmask)) {
        mexErrMsgIdAndTxt("map2almpure:args:notDouble",
                "Input apmask must be real and of type double");
    }
    tmp_ndim = mxGetNumberOfDimensions(ml_apmask);
    if (tmp_ndim != 2) {
        mexErrMsgIdAndTxt("map2almpure:args:dims",
                "Input apmask must have 2 dimensions");
    }

    ml_lmax = prhs[2];
    if (!mxIsInt32(ml_lmax) || mxGetNumberOfElements(ml_lmax)!=1) {
        mexErrMsgIdAndTxt("map2almpure:args:notInt32",
                "Input lmax must be a scalar of type int32");
    }

    ml_mmax = prhs[3];
    if (!mxIsInt32(ml_mmax) || mxGetNumberOfElements(ml_mmax)!=1) {
        mexErrMsgIdAndTxt("map2almpure:args:notInt32",
                "Input mmax must be a scalar of type int32");
    }

    ml_qwghts = prhs[4];
    if (!mxIsDouble(ml_qwghts) || mxIsComplex(ml_qwghts)) {
        mexErrMsgIdAndTxt("map2almpure:args:notDouble",
                "Input qwghts must be real and of type double");
    }
    tmp_ndim = mxGetNumberOfDimensions(ml_qwghts);
    if (tmp_ndim != 2) {
        mexErrMsgIdAndTxt("map2almpure:args:dims",
                "Input qwghts must have 2 dimensions");
    }

    /* Retrieve information from the given inputs */

    /* Determine the number of components, nside, and number of maps from the
     * size of the input map array */
    tmp_ndim = mxGetNumberOfDimensions(ml_map);
    tmp_dims = mxGetDimensions(ml_map);
    npix = tmp_dims[0];
    if (tmp_dims[1] != nstokes) {
        mexErrMsgIdAndTxt("map2almpure:args:sizeMismatch",
                "Expected 2 Stokes parameters, got %d", tmp_dims[1]);
    }
    if (tmp_ndim > 2) {
        nmaps = tmp_dims[2];
    } else {
        nmaps = 1;
    }
    nside = (int)(sqrt((double)npix/12.0));
    dbglog("Using parameters given shape of input map:\n"
           "    npix    = %d\n"
           "    nstokes = %d\n"
           "    nmaps   = %d\n"
           "    nside   = %d\n",
           npix, nstokes, nmaps, nside);

    if (12*nside*nside != npix) {
        mexErrMsgIdAndTxt("map2almpure:args:dims",
                "map is not a valid nside. Got nside = %0.3lf",
                sqrt((double)npix/12.0));
    }

    /* Make sure the apodization mask has correct dimensions */
    tmp_ndim = mxGetNumberOfDimensions(ml_apmask);
    tmp_dims = mxGetDimensions(ml_apmask);
    if (tmp_dims[0] != npix) {
        mexErrMsgIdAndTxt("map2almpure:args:sizeMismatch",
                "Inputs map and apmask must have same number of pixels");
    }
    if (tmp_dims[1] != nmaps) {
        mexErrMsgIdAndTxt("map2almpure:args:sizeMismatch",
                "Inputs map and apmask must have same number of maps");
    }

    /* Make sure the quadrature weights vector has correct dimensions */
    tmp_ndim = mxGetNumberOfDimensions(ml_qwghts);
    tmp_dims = mxGetDimensions(ml_qwghts);
    if (tmp_dims[0] != 2*nside || tmp_dims[1] != 1) {
        mexErrMsgIdAndTxt("map2almpure:args:sizeMismatch",
                "Input qwghts must be an array of size (%d, 1). Got"
                "size (%d, %d).", 2*nside, tmp_dims[0], tmp_dims[1]);
    }

    /* Get the nlmax and nmmax for the output alms */
    lmax = *((int32_t*)mxGetData(ml_lmax));
    mmax = *((int32_t*)mxGetData(ml_mmax));

    /* Allocate the buffer required to output the alms */
    ml_alms_dims[0] = nstokes;
    ml_alms_dims[1] = lmax + 1;
    ml_alms_dims[2] = mmax + 1;
    ml_alms_dims[3] = nmaps;
    ml_alms = mxCreateNumericArray(4, ml_alms_dims, mxDOUBLE_CLASS, mxCOMPLEX);
    plhs[0] = ml_alms;

    /* For a plain double array, we can re-use the same buffer that Matlab
     * has allocated */
    map    = (double*)mxGetData(ml_map);
    apmask = (double*)mxGetData(ml_apmask);
    qwghts = (double*)mxGetData(ml_qwghts);

    /* Allocate the s2hat compatible complex array */
    tmp_ndim = nstokes * (lmax+1) * (mmax+1);
    alms = (s2hat_dcomplex*)calloc(tmp_ndim, sizeof(s2hat_dcomplex));
    if (alms == NULL) {
        mexErrMsgIdAndTxt("map2almpure:mem:error",
                "error allocating alms: %s", strerror(errno));
    }

    dbglog("Running map2almpure...\n");
    int ret = map2almpure(alms, nstokes, lmax, mmax, nmaps, map, apmask,
            nside, qwghts);

    /* Convert from s2hat complex numbers to Matlab format */
    ml_alms_r = mxGetPr(ml_alms);
    ml_alms_i = mxGetPi(ml_alms);
    dbglog("Converting from s2hat to Matlab complex format...\n");
    tmp_ndim = nstokes * (lmax+1) * (mmax+1);
    for (size_t ii=0; ii<tmp_ndim; ++ii) {
        /* Add 0*tmp_ndim to handle nmaps > 1 in the future */
        ml_alms_r[ii + 0*tmp_ndim] = alms[ii].re;
        ml_alms_i[ii + 0*tmp_ndim] = alms[ii].im;
    }

    /* Free the s2hat formatted alms array */
    if (alms) free(alms); alms = NULL;
}

int map2almpure(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
                int nmaps, double* map, double* apmask, int nside,
                double* qwghts)
{
    int ret = 0;
    /* Process communications */
    int myrank;
    int nprocs;
    /* Map description */
    int first_ring = -1;
    int last_ring = -1;
    int map_size = -1;
    s2hat_pixeltype pixel = { 0 };
    s2hat_scandef scan = { 0 };
    s2hat_pixparameters param;
    double bounds[2];
    int* mask = NULL;
    double* local_map = NULL;
    double* local_apmask = NULL;
    int32_t* local_mask = NULL;
    double* local_w8ring = NULL;
    /* Power spectrum/alms */
    int nmvals = -1;
    int* mvals = NULL;
    s2hat_dcomplex* local_alms = NULL;
    /* other */
    long nplm = -1;
    size_t tmp;


    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Setup to describe the type of map we want to analyze. */
    param.par1 = nside;
    param.par2 = 0;
    set_pixelization(PIXCHOICE_HEALPIX, param, &pixel);
    /* Automatically determine the observation region by building a mask from
     * the non-zero/inf/nan entries in the apodization mask */
    if (myrank == 0) {
        mask = (int*)calloc(pixel.npixsall, sizeof(int));
        if (mask == NULL) {
            perror("alloc mask");
            ret = -1;
            goto cleanup;
        }
        for (size_t ii=0; ii<pixel.npixsall; ++ii) {
            double ap = apmask[ii];
            if (ap == 0.0 || isnan(ap) || isinf(ap)) {
                mask[ii] = 0;
            } else {
                mask[ii] = 1;
            }
        }
        mask2scan(mask, pixel, &scan);
        /* Free the mask now that we know the boundaries */
        free(mask); mask = NULL;
        dbglog("Operating on %d of %d rings (%0.2lf%%)\n",
                scan.nringsobs, pixel.nringsall,
                100.0*scan.nringsobs/pixel.nringsall);
    }
    /* All processes need the information, so synchronize */
    MPI_scanBcast(pixel, &scan, 0, myrank, MPI_COMM_WORLD);

    /* To work in a distributed fashion, we need to know the sizes of the
     * various data sets' local size. */
    get_local_data_sizes(0, pixel, scan, lmax, mmax, myrank, nprocs,
            &nmvals, &first_ring, &last_ring, &map_size, &nplm,
            0, MPI_COMM_WORLD);

    mvals = (int*)calloc(nmvals, sizeof(int));
    if (mvals == NULL) {
        perror("alloc mvals");
        ret = -1;
        goto cleanup;
    }
    find_mvalues(myrank, nprocs, mmax, nmvals, mvals);

    /* Make enough space for us hold the local portion of the map,
     * apodization mask, binary mask, and ring weights */
    local_map = (double*)calloc(map_size*nstokes*nmaps, sizeof(double));
    if (local_map == NULL) {
        perror("alloc local_map");
        ret = -1;
        goto cleanup;
    }
    local_apmask = (double*)calloc(map_size*1*nmaps, sizeof(double));
    if (local_apmask == NULL) {
        perror("alloc local_apmask");
        ret = -1;
        goto cleanup;
    }
    local_mask = (int32_t*)calloc(map_size*1*nmaps, sizeof(int32_t));
    if (local_mask == NULL) {
        perror("alloc local_mask");
        ret = -1;
        goto cleanup;
    }
    tmp = (last_ring-first_ring+1)*nstokes;
    local_w8ring = (double*)calloc(tmp, sizeof(double));
    if (local_w8ring == NULL) {
        perror("alloc local_w8ring");
        ret = -1;
        goto cleanup;
    }

    /* Spread the map and apmask across the workers */
    distribute_map(pixel, nmaps, 0, nstokes, first_ring, last_ring, map_size,
            local_map, map, myrank, nprocs, 0, MPI_COMM_WORLD);
    distribute_map(pixel, nmaps, 0, 1, first_ring, last_ring, map_size,
            local_apmask, apmask, myrank, nprocs, 0, MPI_COMM_WORLD);
    /* We only accept a single ring weight vector that applies to all Stokes
     * parameters, but s2hat_map2purealm wants a buffer for each parameter.
     * Therefore, duplicate as necessary. */
    for (int ii=0; ii<nstokes; ++ii) {
        /* local_w8ring stores as (rings,stokes), compute a pointer to the
         * first element in each stokes parameters (assuming column-major
         * storage). */
        double* stokes_w8 = local_w8ring + ii*(last_ring-first_ring+1);
        /* Then broadcast the single weight ring into the local weights */
        distribute_w8ring(1, first_ring, last_ring, stokes_w8,
                pixel.nringsall, qwghts, myrank, nprocs, 0, MPI_COMM_WORLD);
    }

    /* Then allocate enough space for the local part of the alm array to work
     * with */
    local_alms = (s2hat_dcomplex*) calloc(nstokes*(lmax+1)*(mmax+1)*(nmaps),
            sizeof(s2hat_dcomplex));
    if (local_alms == NULL) {
        perror("alloc local_alms");
        ret = -1;
        goto cleanup;
    }

    /* Determine the binary mask by scanning the local apodization mask. Any
     * entries which are identically zero, NaN, or Inf are considered
     * unobserved. */
    tmp = map_size * nmaps;
    for (size_t ii=0; ii<tmp; ++ii) {
        double ap = local_apmask[ii];
        if (ap == 0.0 || isnan(ap) || isinf(ap)) {
            local_mask[ii] = 0;
        } else {
            local_mask[ii] = 1;
        }
    }

    s2hat_map2purealm(pixel, scan, lmax, mmax, nmvals, mvals, nmaps, first_ring,
            last_ring, local_w8ring, map_size, 0, local_mask, 0, local_apmask,
            local_map, nstokes, local_alms, myrank, nprocs, MPI_COMM_WORLD);

    /* All can now unconditionally free local map and weights */
    free(local_map); local_map = NULL;
    free(local_w8ring); local_w8ring = NULL;

    /* Collect alms onto root processor */
    collect_alms(lmax, mmax, nmaps, 0, nstokes, nmvals, mvals, nstokes,
            local_alms, alms, myrank, nprocs, 0, MPI_COMM_WORLD);

    /* All can now unconditionally free local alms */
    free(local_alms); local_alms = NULL;

cleanup:
    /* Cleanup */
    if (local_w8ring != NULL) free(local_w8ring);
    if (local_mask != NULL) free(local_mask);
    if (local_apmask != NULL) free(local_apmask);
    if (local_map != NULL) free(local_map);
    if (mask != NULL) free(mask);
    if (mvals != NULL) free(mvals);
    if (local_alms != NULL) free(local_alms);
    destroy_scan(scan);
    destroy_pixelization(pixel);
    return ret;
}
