#line 1 "mpihelper.c" /* workaround for mex compiler using full path */
#include <mpi.h>
#include <mex.h>
#include "debug.h"

void mpihelper_atexit() {
    dbglog("Finalizing MPI...\n");
    MPI_Finalize();
}

void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
    int mpi_init;

    /* Initialize MPI if necessary */
    MPI_Initialized(&mpi_init);
    if (mpi_init == 0) {
        dbglog("Initializing MPI...\n");
        MPI_Init(0, NULL);

        mexAtExit(&mpihelper_atexit);
    }
}

