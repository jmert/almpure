#line 1 "s2hat_alm2map_c.c" /* workaround for mex compiler using full path */
/* These are needed before s2hat.h can be included since it doesn't take care
 * of loading it's own requirements, apparently... */
#include <sys/types.h>
#include <mpi.h>
/* */
#include <s2hat.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "debug.h"

#include <mex.h>
#include <matrix.h>

int alm2map(double* map, int nside, s2hat_dcomplex* alms, int nstokes,
            int nlmax, int nmmax, int nmaps);

void alm2map_atexit() {
    dbglog("Finalizing MPI...\n");
    MPI_Finalize();
}

/*
 * map = alm2map(alms, nside);
 */
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    /* Inputs from Matlab */
    const mxArray* ml_alms = NULL;
    const mxArray* ml_nside = NULL;

    /* Outputs back to Matlab */
    mxArray* ml_map = NULL;

    /* Intermediaries between Matlab and s2hat */
    mwSize tmp_ndim = 0;
    const mwSize* tmp_dims = NULL;
    double* ml_alms_r = NULL;
    double* ml_alms_i = NULL;
    mwSize ml_map_dims[3] = { 0 };

    /* Internal s2hat interfaces */
    int nstokes = 0;
    int nlmax = 0;
    int nmmax = 0;
    int nmaps = 0;
    int nside = 0;
    int npix = 0;
    s2hat_dcomplex* alms = NULL;
    double* map = NULL;

    /* Initialize MPI if necessary */
    mexCallMATLAB(0, NULL, 0, NULL, "mpihelper");

    /* Validate MATLAB inputs */
    dbglog("Validating MATLAB inputs...\n");
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("alm2map:args:nrhs",
                "Two input arguments are required");
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("alm2map:args:nlhs",
                "One output argument is required");
    }

    ml_alms = prhs[0];
    if (!mxIsDouble(ml_alms) || !mxIsComplex(ml_alms)) {
        mexErrMsgIdAndTxt("alm2map:args:notComplex",
                "Input alms must be a complex");
    }
    tmp_ndim = mxGetNumberOfDimensions(ml_alms);
    if (tmp_ndim != 3 && tmp_ndim != 4) {
        mexErrMsgIdAndTxt("alm2map:args:dims",
                "Input alms must be have 3 or 4 dimensions");
    }

    ml_nside = prhs[1];
    if (!mxIsInt32(ml_nside) || mxGetNumberOfElements(ml_nside)!=1) {
        mexErrMsgIdAndTxt("alm2map:args:notInt32",
                "Input nside must be a scalar of type int32");
    }

    /* Retrieve information from the given inputs */

    /* Determine the number of components, lmax, mmax, and number of maps
     * from the size of the input alms array */
    tmp_ndim = mxGetNumberOfDimensions(ml_alms);
    tmp_dims = mxGetDimensions(ml_alms);
    nstokes = tmp_dims[0];
    nlmax = tmp_dims[1] - 1;
    nmmax = tmp_dims[2] - 1;
    if (tmp_ndim > 3) {
        nmaps = tmp_dims[3];
    } else {
        nmaps = 1;
    }

    /* Get the nside for the output map */
    nside = *((int32_t*)mxGetData(ml_nside));

    dbglog("Using parameters given shape of input alms:\n"
           "    nstokes = %d\n"
           "    nlmax   = %d\n"
           "    nmmax   = %d\n"
           "    nmaps   = %d\n",
           nstokes, nlmax, nmmax, nmaps);

    dbglog("Converting from Matlab to s2hat complex format...\n");
    /* Convert from Matlab complex numbers to s2hat format */
    ml_alms_r = mxGetPr(ml_alms);
    ml_alms_i = mxGetPi(ml_alms);
    tmp_ndim = mxGetNumberOfElements(ml_alms);
    alms = (s2hat_dcomplex*)calloc(tmp_ndim, sizeof(s2hat_dcomplex));
    /* Since both MATLAB and S2HAT use column-major ordering, we don't have
     * to convert indices and can just copy elements in order */
    for (int ii=0; ii<tmp_ndim; ++ii) {
        alms[ii].re = ml_alms_r[ii];
        alms[ii].im = ml_alms_i[ii];
    }

    /* Then allocate the buffer required to output the maps */
    npix = 12 * nside*nside;
    ml_map_dims[0] = npix;
    ml_map_dims[1] = nstokes;
    ml_map_dims[2] = nmaps;
    ml_map = mxCreateNumericArray(3, ml_map_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[0] = ml_map;

    /* For a plain double array, we can re-use the same buffer that Matlab
     * has allocated */
    map = (double*)mxGetData(ml_map);

    dbglog("Running alm2map...\n");
    /* Now execute the s2hat algorithm */
    int ret = alm2map(map, nside, alms, nstokes, nlmax, nmmax, nmaps);

    /* Free the s2hat formated alms array */
    if (alms) free(alms);
}

int alm2map(double* map, int nside, s2hat_dcomplex* alms, int nstokes,
            int nlmax, int nmmax, int nmaps)
{
    int ret = 0;
    /* Process communications */
    int myrank;
    int nprocs;
    /* Power spectrum/alms */
    int nmvals = -1;
    int* mvals = NULL;
    s2hat_dcomplex* local_alms = NULL;
    /* Map description and synthesis */
    int first_ring = -1;
    int last_ring = -1;
    int map_size = -1;
    s2hat_pixeltype pixel = { 0 };
    s2hat_scandef scan = { 0 };
    s2hat_pixparameters param;
    double bounds[2];
    double* local_map = NULL;
    /* other */
    long nplm = -1;


    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Setup to describe the type of map we want to synthesize. We'll do
     * a full-sky simulation at the selected nside. */
    param.par1 = nside;
    param.par2 = 0;
    set_pixelization(PIXCHOICE_HEALPIX, param, &pixel);
    bounds[0] = -1.0; /* full sky */
    bounds[1] =  1.0;
    zbounds2scan(bounds, pixel, &scan);

    /* To work in a distrubted fashion, we need to know the sizes of the
     * various data sets' local size. */
    get_local_data_sizes(0, pixel, scan, nlmax, nmmax, myrank, nprocs,
            &nmvals, &first_ring, &last_ring, &map_size, &nplm,
            0, MPI_COMM_WORLD);

    mvals = (int*)calloc(nmvals, sizeof(int));
    if (mvals == NULL) {
        perror("alloc mvals");
        ret = -1;
        goto cleanup;
    }
    find_mvalues(myrank, nprocs, nmmax, nmvals, mvals);

    /* Then allocate enough space for the local part of the alm array to work
     * with */
    local_alms = (s2hat_dcomplex*) calloc(nstokes*(nlmax+1)*(nmmax+1)*(nmaps),
            sizeof(s2hat_dcomplex));
    if (local_alms == NULL) {
        perror("alloc local_alms");
        ret = -1;
        goto cleanup;
    }

    /* Spread the alms across the workers */
    distribute_alms(nlmax, nmmax, nmaps, 0, nstokes, nmvals, mvals,
            nstokes, local_alms, alms, myrank, nprocs, 0, MPI_COMM_WORLD);

    /* Make enough space for us synthesize the local portion of the map */
    local_map = (double*)calloc(map_size*nstokes*nmaps, sizeof(double));
    if (local_map == NULL) {
        perror("alloc local_map");
        ret = -1;
        goto cleanup;
    }

    s2hat_alm2map(0, pixel, scan, nlmax, nmmax, nmvals, mvals, nmaps,
            nstokes, first_ring, last_ring, map_size, local_map, nstokes,
            local_alms, 0, NULL, nprocs, myrank, MPI_COMM_WORLD);

    /* All can now unconditionally free local alms */
    free(local_alms); local_alms = NULL;

    /* Collect map onto root processor */
    collect_map(pixel, nmaps, 0, nstokes, map, first_ring, last_ring,
            map_size, local_map, myrank, nprocs, 0, MPI_COMM_WORLD);

    /* All can now unconditionally free local maps */
    free(local_map); local_map = NULL;

cleanup:
    /* Cleanup */
    if (local_map != NULL) free(local_map);
    if (mvals != NULL) free(mvals);
    if (local_alms != NULL) free(local_alms);
    destroy_scan(scan);
    destroy_pixelization(pixel);
    return ret;
}
