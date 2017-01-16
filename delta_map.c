/* These are needed before s2hat.h can be included since it doesn't take care
 * of loading it's own requirements, apparently... */
#include "sys/types.h"
#include "mpi.h"
/* */
#include "s2hat/s2hat.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* Hard code several parameters for ease of use right now. */
const int nlmax = 700;
const int nmmax = 700;
const int nstokes = 3;
const int nmaps = 1;
/* Choose the HEALPix nside to synthesize */
const int nside = 128;
/* Choose the delta ell to generate */
const int delta_ell = 80;
/* Choose the random seed to use */
const int seed = 0;

double randn();
void gen_delta_alms(s2hat_dcomplex* alms, int nstokes, int nlmax, int nmmax,
                    int delta_ell);

int main(int argc, char** argv) {
    int ret = 0;
    /* Process communications */
    int myrank;
    int nprocs;
    /* Power spectrum/alms */
    int nmvals = -1;
    int* mvals = NULL;
    s2hat_dcomplex* alms = NULL;
    s2hat_dcomplex* local_alms = NULL;
    /* Map description and synthesis */
    int first_ring = -1;
    int last_ring = -1;
    int map_size = -1;
    s2hat_pixeltype pixel = { 0 };
    s2hat_scandef scan = { 0 };
    s2hat_pixparameters param = { nside, 0 };
    double bounds[2] = {-1.0, 1.0}; /* full sky */
    double* local_map = NULL;
    double* map = NULL;
    /* other */
    long nplm = -1;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (myrank == 0) {
        /* Allocate enough space for alms for nstokes-many spectra */
        alms = (s2hat_dcomplex*) calloc(nstokes*(nlmax+1)*(nmmax+1),
                sizeof(s2hat_dcomplex));
        if (alms == NULL) {
            perror("alloc alms");
            ret = -1;
            goto cleanup;
        }
    }

    /* Setup to describe the type of map we want to synthesize. We'll do
     * a full-sky simulation at the selected nside. */
    set_pixelization(PIXCHOICE_HEALPIX, param, &pixel);
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

    /* Fill in a random delta realization */
    srand(seed);
    gen_delta_alms(alms, nstokes, nlmax, nmmax, delta_ell);
    /* Then spread across the workers */
    distribute_alms(nlmax, nmmax, nmaps, 0, nstokes, nmvals, mvals,
            nstokes, local_alms, alms, myrank, nprocs, 0, MPI_COMM_WORLD);
    /* No longer need the entire array, so save some RAM */
    if (myrank == 0) {
        free(alms);
        alms = NULL;
    }

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

    if (myrank == 0) {
        map = (double*)calloc(pixel.npixsall*nstokes, sizeof(double));
        if (map == NULL) {
            perror("alloc map");
            ret = -1;
            goto cleanup;
        }
    }
    collect_map(pixel, nmaps, 0, nstokes, map, first_ring, last_ring,
            map_size, local_map, myrank, nprocs, 0, MPI_COMM_WORLD);

    /* All can now unconditionally free local maps */
    free(local_map); local_map = NULL;

    FILE* hfile = fopen("delta_map.dat", "w+");
    int npix = pixel.npixsall;
    for (int ii=0; ii<npix; ++ii) {
        fprintf(hfile, "%i\t%020.15lf\t%020.15lf\t%020.15lf\n", ii,
                map[ii + 0*npix],
                map[ii + 1*npix],
                map[ii + 2*npix]);
    }
    fclose(hfile);

cleanup:
    /* Cleanup */
    if (map != NULL) free(map);
    if (local_map != NULL) free(local_map);
    if (mvals != NULL) free(mvals);
    if (local_alms != NULL) free(local_alms);
    if (alms != NULL) free(alms);
    destroy_scan(scan);
    destroy_pixelization(pixel);
    MPI_Finalize();
    return ret;
}

void gen_delta_alms(s2hat_dcomplex* alms, int nstokes, int nlmax, int nmmax,
                    int delta_ell)
{
    /* Note that the array is assumed to be in column-major order by the
     * underlying library, so take note that this disagrees with the typical
     * C convention. */

    /* m == 0 must be real, so handle it before we start the loop. */
    for (int ii=0; ii<nstokes; ++ii) {
        s2hat_dcomplex v;
        v.re = randn();
        v.im = 0;
        alms[ii + delta_ell*(nstokes) + 0] = v;
    }
    /* All other m values are then complex */
    for (int mm=1; mm<delta_ell; ++mm) {
        for (int ii=0; ii<nstokes; ++ii) {
            s2hat_dcomplex v;
            v.re = randn() / sqrt(2);
            v.im = randn() / sqrt(2);
            alms[ii + delta_ell*(nstokes) + mm*(nstokes*nlmax)] = v;
        }
    }
}

/* This isn't as efficient as it could be since the return value U2*mult
 * is never used and just thrown away, but whatever... this was simpler. */
double randn() {
    double U1;
    double U2;
    double W;
    double mult;

    do {
        U1 = -1 + 2 * ((double)rand() / RAND_MAX);
        U2 = -1 + 2 * ((double)rand() / RAND_MAX);
        W = pow(U1,2) + pow(U2,2);
    } while (W >= 1 || W == 0);

    mult = sqrt( (-2*log(W)) / W);
    return U1 * mult;
}
