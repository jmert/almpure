#line 1 "map2alm.c" /* workaround for mex compiler using full path */
#include <sys/types.h>
#include <mpi.h>
#include <s2hat.h>
#include <s2hat_pure.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "debug.h"

int map2alm_c(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
            int nmaps, double* map, int nside, double* qwghts)
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
    double* local_map = NULL;
    double* local_w8ring = NULL;
    /* Power spectrum/alms */
    int nmvals = -1;
    int* mvals = NULL;
    s2hat_dcomplex* local_alms = NULL;
    /* other */
    long nplm = -1;
    size_t tmp;


    /* Initialize MPI if necessary */
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        dbglog("Initializing MPI...\n");
        MPI_Init(0, NULL);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Setup to describe the type of map we want to analyze. We'll assume
     * full-sky at the selected nside. */
    param.par1 = nside;
    param.par2 = 0;
    set_pixelization(PIXCHOICE_HEALPIX, param, &pixel);
    bounds[0] = -1.0; /* full sky */
    bounds[1] =  1.0;
    zbounds2scan(bounds, pixel, &scan);

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

    /* Make enough space for us hold the local portion of the map and ring
     * weights */
    local_map = (double*)calloc(map_size*nstokes*nmaps, sizeof(double));
    if (local_map == NULL) {
        perror("alloc local_map");
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

    /* Spread the map across the workers */
    distribute_map(pixel, nmaps, 0, nstokes, first_ring, last_ring, map_size,
            local_map, map, myrank, nprocs, 0, MPI_COMM_WORLD);
    /* We only accept a single ring weight vector that applies to all Stokes
     * parameters, but s2hat_map2alm wants a buffer for each parameter.
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

    s2hat_map2alm(0, pixel, scan, lmax, mmax, nmvals, mvals, nmaps,
            nstokes, first_ring, last_ring, local_w8ring, map_size, local_map,
            nstokes, local_alms, 0, NULL, nprocs, myrank, MPI_COMM_WORLD);

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
    if (local_map != NULL) free(local_map);
    if (mvals != NULL) free(mvals);
    if (local_alms != NULL) free(local_alms);
    destroy_scan(scan);
    destroy_pixelization(pixel);
    return ret;
}

int map2almpure_c(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
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


    /* Initialize MPI if necessary */
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        dbglog("Initializing MPI...\n");
        MPI_Init(0, NULL);
    }

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

int alm2map_c(double* map, int nside, s2hat_dcomplex* alms, int nstokes,
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


    /* Initialize MPI if necessary */
    int mpi_init;
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        dbglog("Initializing MPI...\n");
        MPI_Init(0, NULL);
    }

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
