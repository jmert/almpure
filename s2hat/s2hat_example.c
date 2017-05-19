
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "sys/types.h"
#include "sys/stat.h"
#include "sys/param.h"
#include "mpi.h"

#include "s2hat.h"

/* this code computes first map2alm transform of some map, collects the alm coefficients in the memory       *
 * of the proc root and outputs them as a binary file 'test_alm.bin'. Then it performs alm2map and collects  *
 * the map in the memory of the proc root. It outputs the map as 'test_map.bin'. Then it computes the power  *
 * spectrum out of the precomputed alm coefficients. It is written as a file 'test_cl.bin'. Note that though *
 * in principle the code performs the 'full loop' map2alm followed by alm2map the outputted map will in not  *
 * agree perfectly with the input due to the bandwidth issues.                                               *
 *                                                                                                           *
 *                                                                                      - rs@apc, 2010/02/10 *
 * Works correctly for a multiple maps and any accepted number of Stokes parameters                          *
 *                                                                                                           *
 *                                                                                      - rs@apc, 2010/09/07 */

s2hat_int4 main(s2hat_int4 argc, char * argv[])
{

  s2hat_int4 pixchoice = PIXCHOICE_HEALPIX;   /* use HEALPIX pixelization */

  s2hat_int4 nside = 64;
  s2hat_int4 nmaps = 1;
  s2hat_int4 lmax = 4*nside-1;
  s2hat_int4 npix=12*nside*nside;
  s2hat_int4 root=0;

  s2hat_int4 myrank, nprocs;

  s2hat_int4 nlmax=lmax, nmmax=lmax;
  s2hat_int4 nstokes=3;  /* no of Stokes params 1, 2 or 3 */
  s2hat_int4 nspec = 1;  /* no of map spectra to be computed - check the rules in the manual */
  s2hat_int4 i, j, lout, mout;
  s2hat_int4 plms=0;     /* i.e., do not assume that plm are precomputed */
  s2hat_int4 nmvals, first_ring, last_ring, map_size;
  s2hat_int4 *mvals;
  s2hat_int8 nplm;
  s2hat_int4 one = 1, zero = 0;

  s2hat_flt8 *local_map, *map;
  s2hat_int4 nrings;
  s2hat_flt8 *local_w8ring;

  s2hat_flt8 *cls;

  FILE *fout;

  s2hat_int4 nalms;
  s2hat_dcomplex *local_alm, *alms;

  s2hat_pixeltype cpixelization;
  s2hat_pixparameters pixpar;
  s2hat_scandef cscan;

  s2hat_flt8 zbounds[2];

  /* MPI initialization */

  MPI_Init( &argc, &argv);

  MPI_Comm_rank( MPI_COMM_WORLD, &myrank);
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

  /* define the maps parameters and data distribution */

  /* - define the pixelization structure -> cpixelization */
  pixpar.par1 = nside; pixpar.par2 = 0;
  set_pixelization( pixchoice, pixpar, &cpixelization); 

  /* - define the scan structure -> cscan : corresponding to the observed sky given by zbounds[0] >= cos theta >= zbounds[1] */

  zbounds[0] = 0.0; zbounds[1] = 0.0;
  zbounds2scan( zbounds, cpixelization, &cscan);

  /* computes the data distribution over the processors */

  get_local_data_sizes( plms, cpixelization, cscan, nlmax, nmmax, myrank, nprocs,
			&nmvals, &first_ring, &last_ring, &map_size, &nplm, root, MPI_COMM_WORLD);  

  /* number of equilatitude rings on this processor */
  nrings = last_ring-first_ring+1;

  /* m-values to be stored on this proc */
  mvals = (int *) calloc( nmvals, sizeof( int));
  find_mvalues( myrank, nprocs, nmmax, nmvals, mvals);

  /* number of alm coefficients stored on this proc */
  nalms = (nlmax+1)*nmvals*nstokes;

  /* to store part of the map */
  local_map = (s2hat_flt8 *)calloc( nmaps*map_size*nstokes, sizeof( s2hat_flt8));

  /* to assign some arbitrary values to the map NB. these depend on a number of pixs used */
  for( j=0; j<nmaps; j++) for( i=0; i<map_size*nstokes; i++) local_map[i+j*nstokes*map_size] = 1.0*i;

  /* quadrature weights */
  local_w8ring = (s2hat_flt8 *)calloc( nrings*nstokes, sizeof( s2hat_flt8));

  /* i.e. use no weighting */
  for( i=0; i<nrings*nstokes; i++) local_w8ring[i] = 1.0;

  local_alm = (s2hat_dcomplex *) calloc( nalms*nstokes, sizeof( s2hat_dcomplex));

  /* calculate the alm coefficients */

  s2hat_map2alm( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals,
	         nmaps, nstokes, first_ring, last_ring,
	         local_w8ring, map_size, local_map, nlmax,
	         local_alm, nplm, NULL, nprocs, myrank, MPI_COMM_WORLD);

  free( local_w8ring);

  /* gather all computed alm in the memory of the proc 'root' */

  if( myrank == root)
  {
     alms = (s2hat_dcomplex *)calloc( nstokes*(nlmax+1)*(nmmax+1), sizeof( s2hat_dcomplex));
     fout = fopen( "test_alm.bin","w"); fclose( fout);
  }

  /* do a loop over maps */
  for( i=0; i<nmaps; i++)
  {
     collect_alms( nlmax, nmmax, nmaps, i, nstokes, nmvals, mvals, nlmax, local_alm, alms, myrank, nprocs, root, MPI_COMM_WORLD);

     if( myrank == root) 
     {
        fout = fopen( "test_alm.bin", "a"); fwrite( alms, sizeof( s2hat_dcomplex), nstokes*(nmmax+1)*(nlmax+1), fout); fclose( fout);
     }
  }

  if( myrank == root) free( alms);

  /* go back from alms to map */

  s2hat_alm2map( plms, cpixelization, cscan, nlmax, nmmax, nmvals, mvals, nmaps, nstokes, first_ring, last_ring,
		 map_size, local_map, nlmax, local_alm, nplm, NULL, nprocs, myrank, MPI_COMM_WORLD);

  /* gather complete map in the memory of the proc root */

  if( myrank == root)
  {
      map = (s2hat_flt8 *)calloc( nstokes*npix, sizeof( s2hat_flt8));
      fout = fopen( "test_map.bin","w"); fclose( fout);
  }

  /* do the loop over maps */
  for( i=0; i<nmaps; i++)
  {
     collect_map(cpixelization,nmaps,i,nstokes,map,first_ring,last_ring,map_size,local_map,myrank,nprocs,root,MPI_COMM_WORLD);

     if( myrank == root)
     {
         fout=fopen( "test_map.bin","a"); fwrite( map, sizeof( double), npix*nstokes, fout); fclose( fout);
     }
  }

  if( myrank == root) free( map);

  free( local_map);

  /* and compute (auto)spectra just for the first map */

  cls = (s2hat_flt8 *)calloc( (nlmax+1)*nspec, sizeof( s2hat_flt8));
  collect_cls( nmaps, zero, nstokes, nlmax, nmvals, mvals, nlmax, local_alm, 
               nspec, cls, myrank, nprocs, root, MPI_COMM_WORLD);

  free( local_alm);

  if( myrank == root) 
  {
     fout = fopen( "test_cl.bin", "w"); fwrite( cls, sizeof( double), (nlmax+1)*nspec, fout); fclose( fout);
  }

  free( cls);

  /* conclude */

  destroy_scan( cscan);
  destroy_pixelization( cpixelization);

  /* and finish up */

  MPI_Finalize();
  return( 0);

}
