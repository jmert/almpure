
/* have to be 'included' earlier ...
  #include "sys/types.h"
  #include "mpi.h"
  #include "malloc.h"
*/

#define s2hat_int4   int
#define s2hat_int8   int64_t /* long long int */

#define s2hat_flt4   float
#define s2hat_flt8   double

#define s2hat_mpi_int4  MPI_INT
#define s2hat_mpi_int8  MPI_LONG_LONG_INT
#define s2hat_mpi_flt4  MPI_FLOAT
#define s2hat_mpi_flt8  MPI_DOUBLE
#define s2hat_mpi_2int4 MPI_2INT

/* define pixel choices here */
#define PIXCHOICE_HEALPIX 0    /* for HEALPIX pixelization */
#define PIXCHOICE_GLESP   1    /* for GLESP pixelization */
#define PIXCHOICE_ECP     2    /* for ECP gridding */
#define PIXCHOICE_GLCP    3    /* for GLCP gridding */

/* define the spin convention */

#define SPIN_CONV_SIGN   -1    /* i.e. HEALPIX choice */

typedef struct {
  s2hat_flt8   re;
  s2hat_flt8   im;
} s2hat_dcomplex;

typedef struct {
  s2hat_flt4   re;
  s2hat_flt4   im;
} s2hat_fcomplex;

typedef struct {
  s2hat_int4 par1;
  s2hat_int4 par2;
} s2hat_pixparameters;

typedef struct {    /* pixelization is assumed to be iso-latitudal with poixels evenly spaced for each latitude */
					/* and a symmetry present between northern and southern hemispheres.                        */
  /* a total number of pixels - 8 byte ints to avoid padding problems in the f90/C-interface */
  s2hat_int8 type;
  /* a total number of pixels */  
  s2hat_int8 npixsall;
  /* a total number of iso-latitude rings in the north hemisphere (including equator) */
  s2hat_int8 nringsall;
  /* a total maximum number of pixels per iso-ring */  
  s2hat_int8 nphmx;
  /* a number of the first pixel for each ring [0, nringsall-1] */
  s2hat_int8 *fpix;
  /* a number of pixel for each iso-ring [0, nringsall-1] */
  s2hat_int8 *nph;
  /* an offset of the 1st pixel of each ring wrt the meridian zero (in radians) [0, nringsall-1] */
  s2hat_flt8 *kphi;
  /* quadrature weights [0, nringsall-1] */
  s2hat_flt8 *qwght;
  /* pixel center separations for each iso-ring [0, nringsall-1] */
  s2hat_flt8 *pixphi;
  /* pixel area (assumed to be constant) for each iso-ring [0, nringsall-1] */
  s2hat_flt8 *parea;
  /* cosines of the polar angle for each iso-ring [0, nringsall-1] */
  s2hat_flt8 *cth;
  /* sines of the polar angle for each iso-ring (redundant) [0, nringsall-1] */
  s2hat_flt8 *sth; 
} s2hat_pixeltype;

typedef struct { /* defines the scan parameters needed for the s2hat transforms */
  /* a total number of observed pixels ! all 8byte int to simplify the f90/C-interfacing */
  s2hat_int8 npixsobs; 
  /* a total number of observed rings */
  s2hat_int8 nringsobs;
  /* observation flags: */
  /* northern hemisphere (includes equator) [0:nringsall-1] */
  s2hat_int8 *nfl;
  /* southern hemisphere [ 0:nringsall-1] */
  s2hat_int8 *sfl;
  /* either northern or southern (redundant) [0:nringsall-1] */
  s2hat_int8 *fl;
} s2hat_scandef;
