
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "malloc.h"
#include "sys/types.h"
#include "mpi.h"

#include "s2hat.h"

/* MPI communication routines for S2HAT structures */

s2hat_int4 MPI_scanBcast( s2hat_pixeltype cpixelization, s2hat_scandef *cscan, s2hat_int4 root, s2hat_int4 my_rank, MPI_Comm current_comm)

{

    MPI_Bcast( &((*cscan).npixsobs), 1, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( &((*cscan).nringsobs), 1, s2hat_mpi_int8, root, current_comm);

    if( my_rank != root)
    {
       (*cscan).nfl = (s2hat_int8 *)calloc( cpixelization.nringsall, sizeof( s2hat_int8));
       (*cscan).sfl = (s2hat_int8 *)calloc( cpixelization.nringsall, sizeof( s2hat_int8));
       (*cscan).fl = (s2hat_int8 *)calloc( cpixelization.nringsall, sizeof( s2hat_int8));
    }

    MPI_Bcast( (*cscan).nfl, (s2hat_int4 )cpixelization.nringsall, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( (*cscan).sfl, (s2hat_int4 )cpixelization.nringsall, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( (*cscan).fl, (s2hat_int4 )cpixelization.nringsall, s2hat_mpi_int8, root, current_comm);

    return( 1);
}

s2hat_int4 MPI_pixelizationBcast( s2hat_pixeltype* cpixelization, s2hat_int4 root, s2hat_int4 my_rank, MPI_Comm current_comm)

{

    MPI_Bcast( &((*cpixelization).type), 1, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( &((*cpixelization).npixsall), 1, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( &((*cpixelization).nringsall), 1, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( &((*cpixelization).nphmx), 1, s2hat_mpi_int8, root, current_comm);

    if( my_rank != root)
    {
       (*cpixelization).fpix = (s2hat_int8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_int8));
       (*cpixelization).nph = (s2hat_int8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).kphi = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).qwght = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).pixphi = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).parea = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).cth = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));
       (*cpixelization).sth = (s2hat_flt8 *)calloc( (*cpixelization).nringsall, sizeof( s2hat_flt8));

    }

    MPI_Bcast( (*cpixelization).fpix, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_int8, root, current_comm);
    MPI_Bcast( (*cpixelization).nph, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_int8, root, current_comm);

    MPI_Bcast( (*cpixelization).kphi, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);
    MPI_Bcast( (*cpixelization).qwght, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);
    MPI_Bcast( (*cpixelization).pixphi, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);
    MPI_Bcast( (*cpixelization).parea, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);
    MPI_Bcast( (*cpixelization).cth, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);
    MPI_Bcast( (*cpixelization).sth, (s2hat_int4 )(*cpixelization).nringsall, s2hat_mpi_flt8, root, current_comm);

    return( 1);

}

/* c-wrappers on the main routines */

void  set_pixelization( s2hat_int4 pixchoice, s2hat_pixparameters pixpars, s2hat_pixeltype *cpixelization)

{
    s2hat_pixeltype f90pixelization;
    s2hat_int4 *parint, ntwo = 2;

    parint = (s2hat_int4 *)calloc( ntwo, sizeof( s2hat_int4));

    parint[0] = pixpars.par1;
    parint[1] = pixpars.par2;

    c_set_pixelization( &pixchoice, &ntwo, parint, &f90pixelization);

    free( parint);

    f2c_pixelization( &f90pixelization, cpixelization);

    c_destroy_pixelization( &f90pixelization);
}

void destroy_pixelization( s2hat_pixeltype cpixelization)

{

    free( cpixelization.fpix);
    free( cpixelization.nph);
    free( cpixelization.kphi);
    free( cpixelization.pixphi);
    free( cpixelization.parea);
    free( cpixelization.cth);
    free( cpixelization.sth);

}

void destroy_scan( s2hat_scandef cscan)

{

    free( cscan.nfl);
    free( cscan.sfl);
    free( cscan.fl);

}

void zbounds2scan( s2hat_flt8 *zbounds, s2hat_pixeltype cpixelization, s2hat_scandef *cscan)

{

    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_zbounds2scan( zbounds, &f90pixelization, &f90scan);

    f2c_scan( &f90pixelization, &f90scan, cscan);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void zbounds2mask( s2hat_flt8 *zbounds, s2hat_pixeltype cpixelization, s2hat_int4 *mask)

{
    s2hat_pixeltype f90pixelization;

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_zbounds2mask( zbounds, &f90pixelization, mask);

    c_destroy_pixelization( &f90pixelization);

}

void mask2scan( s2hat_int4 *mask, s2hat_pixeltype cpixelization, s2hat_scandef *cscan)

{
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_mask2scan( mask, &f90pixelization, &f90scan);

    f2c_scan( &f90pixelization, &f90scan, cscan);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void fft_setup( s2hat_pixeltype cpixelization, s2hat_int4 nmmax, s2hat_int4 opt)

{
    s2hat_pixeltype f90pixelization;

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_fft_setup( &f90pixelization, &nmmax, &opt);
}

void fft_mc_setup( s2hat_pixeltype cpixelization, s2hat_int4 nmmax)

{
    s2hat_pixeltype f90pixelization;

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_fft_mc_setup( &f90pixelization, &nmmax);
}

void fft_mc_clean()

{
    c_fft_mc_clean();
}


void distribute_local_data_objects_map2alm( s2hat_int4 precompute_plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, 
                                            s2hat_int4 nstokes, s2hat_flt8 *map, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int8 nplm, s2hat_int4 lda,
                                            s2hat_flt8 *local_plm, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_flt8 *w8ring, s2hat_int4 myid, 
                                            s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;

    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    /* no plms precomputed - default */
    ldim1_dwn = 1; ldim1_up = 1;
    ldim2_dwn = 1; ldim2_up = 1;

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
    {

      if( precompute_plms == 1)
      {
          ldim1_dwn = 1; ldim1_up = 1;
          ldim2_dwn = 0; ldim2_up = nplm-1;
      }
      else
      {
	if( precompute_plms == 2)
        {
            ldim1_dwn = 1; ldim1_up = nstokes;
            ldim2_dwn = 0; ldim2_up = nplm-1;
	}
      }

    } 
    else
    {
        if( precompute_plms == 1)
	{
            ldim1_dwn = 0; ldim1_up = nplm-1;
            ldim2_dwn = 1; ldim2_up = 1;
        }
        else
	{
	    if( precompute_plms == 2)
	    {
                ldim1_dwn = 0; ldim1_up = nplm-1;
                ldim2_dwn = 1; ldim2_up = nstokes;
	    } 
	}
    }

    c2f_scan( &cpixelization, &cscan, &f90scan);
    c2f_pixelization( &cpixelization, &f90pixelization);

    c_distribute_local_data_objects_map2alm( &precompute_plms, &f90pixelization, &f90scan, &nlmax, &nmmax, &nmaps, &nstokes, map, &map_size, local_map,
					     &nmvals, mvals, &nplm, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, local_plm, &first_ring, &last_ring, 
                                              local_w8ring, w8ring, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void distribute_local_data_objects_alm2map( s2hat_int4 precompute_plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, 
                                            s2hat_int4 nstokes, s2hat_int4 lda, s2hat_dcomplex *alms, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_dcomplex *local_alm,
		       			    s2hat_int8 nplm, s2hat_flt8 *local_plm, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;

    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up, ldim7_dwn, ldim7_up, ldim8_dwn, ldim8_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    /* no plms precomputed - default */
    ldim7_dwn = 1; ldim7_up = 1;
    ldim8_dwn = 1; ldim8_up = 1;

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmmax;

        ldim4_dwn = 1; ldim4_up = nstokes;
        ldim5_dwn = 0; ldim5_up = nlmax;
        ldim6_dwn = 0; ldim6_up = nmvals-1;


        if( precompute_plms == 1)
        {
            ldim7_dwn = 1; ldim7_up = 1;
            ldim8_dwn = 0; ldim8_up = nplm-1;
        }
        else
        {
	   if( precompute_plms == 2)
           {
               ldim7_dwn = 1; ldim7_up = nstokes;
               ldim8_dwn = 0; ldim8_up = nplm-1;
	   }
        }

    } 
    else
    {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmmax;
        ldim3_dwn = 1; ldim3_up = nstokes;

        ldim4_dwn = 0; ldim4_up = nlmax;
        ldim5_dwn = 0; ldim5_up = nmvals-1;
        ldim6_dwn = 1; ldim6_up = nstokes;

        if( precompute_plms == 1)
	{
            ldim7_dwn = 0; ldim7_up = nplm-1;
            ldim8_dwn = 1; ldim8_up = 1;
        }
        else
	{
	    if( precompute_plms == 2)
	    {
                ldim7_dwn = 0; ldim7_up = nplm-1;
                ldim8_dwn = 1; ldim8_up = nstokes;
	    } 
	}

    }

    c2f_scan( &cpixelization, &cscan, &f90scan);
    c2f_pixelization( &cpixelization, &f90pixelization);

    c_distribute_local_data_objects_alm2map( &precompute_plms, &f90pixelization, &f90scan, &nlmax, &nmmax, &nmaps, &nstokes, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, 
                                             &ldim3_dwn, &ldim3_up, alms,  &nmvals, mvals, &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, &ldim6_dwn, &ldim6_up, local_alm, 
                                             &nplm, &ldim7_dwn, &ldim7_up, &ldim8_dwn, &ldim8_up, local_plm, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void get_local_data_sizes( s2hat_int4 precompute_plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 myid, 
                           s2hat_int4 numprocs, s2hat_int4 *nmvals, s2hat_int4 *first_ring, s2hat_int4 *last_ring, s2hat_int4 *map_size, s2hat_int8 *nplm, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_scan( &cpixelization, &cscan, &f90scan);
    c2f_pixelization( &cpixelization, &f90pixelization);

    c_get_local_data_sizes( &precompute_plms, &f90pixelization, &f90scan, &nlmax, &nmmax, &myid, &numprocs, nmvals, first_ring, 
                            last_ring, map_size, nplm, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void find_mvalues( s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals)

{

    c_find_mvalues( &myid, &numprocs, &nmmax, &nmvals, mvals);

}

s2hat_int4 nummvalues( s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 nmmax)

{
  s2hat_int4 nummvals;

  nummvals = c_nummvalues( &myid, &numprocs, &nmmax);

  return( nummvals);
}

s2hat_int4 nummmodes( s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals)

{
    s2hat_int4 nmodes;

    nmodes = c_nummmodes( &nlmax, &nmvals, mvals);

    return( nmodes);

}

void find_scan_ring_range( s2hat_pixeltype pixelization, s2hat_scandef scan, s2hat_int4 nmmax, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 *first_ring, s2hat_int4 *last_ring, s2hat_int4 *outerror)

{

  c_find_scan_ring_range( &pixelization, &scan, &nmmax, &myid, &numprocs, first_ring, last_ring, outerror);

}

void distribute_map( s2hat_pixeltype cpixelization, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 map_size, s2hat_flt8 *local_map, 
                     s2hat_flt8 *map, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_distribute_map( &f90pixelization, &nmaps, &mapnum, &nstokes, &first_ring, &last_ring, &map_size, local_map, map, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);

}

void distribute_w8ring(s2hat_int4 npol, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring, s2hat_int4 nringsall, s2hat_flt8 *w8ring, s2hat_int4 myid, s2hat_int4 numprocs, 
                       s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c_distribute_w8ring( &npol, &first_ring, &last_ring, local_w8ring, &nringsall, w8ring, &myid, &numprocs, &root, &fcomm);

}

void distribute_mask( s2hat_pixeltype cpixelization, s2hat_int4 nmasks, s2hat_int4 masknum, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 mask_size, s2hat_int4 *local_mask, 
                      s2hat_int4 *mask, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_distribute_mask( &f90pixelization, &nmasks, &masknum, &nstokes, &first_ring, &last_ring, &mask_size, local_mask, mask, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);

}

void distribute_partialmap( s2hat_pixeltype cpixelization, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 localmap_size, s2hat_flt8 *local_map, 
                            s2hat_int8 firstPix, s2hat_int4 map_size, s2hat_flt8 *map, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_distribute_partialmap( &f90pixelization, &nmaps, &mapnum, &nstokes, &first_ring, &last_ring, &localmap_size, local_map, &firstPix, &map_size, map, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);

}

void distribute_alms( s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, 
                      s2hat_dcomplex *alms, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;

    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
      {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

        ldim4_dwn = 1; ldim4_up = nstokes;
        ldim5_dwn = 0; ldim5_up = nlmax;
        ldim6_dwn = 0; ldim6_up = nmmax;
      } 
    else
      {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

        ldim4_dwn = 0; ldim4_up = nlmax;
        ldim5_dwn = 0; ldim5_up = nmmax;
        ldim6_dwn = 1; ldim6_up = nstokes;
      }


    c_distribute_alms( &nlmax, &nmmax, &nmaps, &mapnum, &nstokes, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, 
                       local_alm, &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, &ldim6_dwn, &ldim6_up, alms, &myid, &numprocs, &root, &fcomm);

}

void distribute_partialalms( s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, 
                      s2hat_int4 mmin, s2hat_int4 mmax, s2hat_dcomplex *alms, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;

    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
      {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

        ldim4_dwn = 1; ldim4_up = nstokes;
        ldim5_dwn = 0; ldim5_up = nlmax;
        ldim6_dwn = 0; ldim6_up = nmmax;
      } 
    else
      {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

        ldim4_dwn = 0; ldim4_up = nlmax;
        ldim5_dwn = 0; ldim5_up = nmmax;
        ldim6_dwn = 1; ldim6_up = nstokes;
      }

    c_distribute_partialalms( &nlmax, &nmmax, &nmaps, &mapnum, &nstokes, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, 
                              local_alm, &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, &ldim6_dwn, &ldim6_up, &mmin, &mmax, alms, &myid, &numprocs, &root, &fcomm);

}

void collect_map( s2hat_pixeltype cpixelization, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_flt8 *map, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 map_size, 
                  s2hat_flt8 *local_map, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_collect_map( &f90pixelization, &nmaps, &mapnum, &nstokes, map, &first_ring, &last_ring, &map_size, local_map, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);

}

void collect_alms( s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, 
                   s2hat_dcomplex *alms, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
      {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

        ldim4_dwn = 1; ldim4_up = nstokes;
        ldim5_dwn = 0; ldim5_up = nlmax;
        ldim6_dwn = 0; ldim6_up = nmmax;
      } 
    else
      {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

        ldim4_dwn = 0; ldim4_up = nlmax;
        ldim5_dwn = 0; ldim5_up = nmmax;
        ldim6_dwn = 1; ldim6_up = nstokes;
      }

     c_collect_alms( &nlmax, &nmmax, &nmaps, &mapnum, &nstokes, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, 
                     &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, &ldim6_dwn, &ldim6_up, alms, &myid, &numprocs, &root, &fcomm);

}

void collect_partialmap( s2hat_pixeltype cpixelization, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int8 firstPix, s2hat_int4 map_size, s2hat_flt8 *map, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 localmap_size,
                         s2hat_flt8 *local_map, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);

    c_collect_partialmap( &f90pixelization, &nmaps, &mapnum, &nstokes, &firstPix, &map_size, map, &first_ring, &last_ring, &localmap_size, local_map, &myid, &numprocs, &root, &fcomm);

    c_destroy_pixelization( &f90pixelization);

}

void collect_partialalms( s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, 
						  s2hat_int4 mmin, s2hat_int4 mmax, s2hat_dcomplex *alms, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up, ldim6_dwn, ldim6_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
      {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

        ldim4_dwn = 1; ldim4_up = nstokes;
        ldim5_dwn = 0; ldim5_up = nlmax;
        ldim6_dwn = 0; ldim6_up = nmmax;
      } 
    else
      {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

        ldim4_dwn = 0; ldim4_up = nlmax;
        ldim5_dwn = 0; ldim5_up = nmmax;
        ldim6_dwn = 1; ldim6_up = nstokes;
      }

     c_collect_partialalms( &nlmax, &nmmax, &nmaps, &mapnum, &nstokes, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, 
                            &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, &ldim6_dwn, &ldim6_up, &mmin, &mmax, alms, &myid, &numprocs, &root, &fcomm);

}

void collect_cls( s2hat_int4 nmaps, s2hat_int4 mapnum, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_int4 nspec, s2hat_flt8 *cl, 
                  s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;
    } 
    else
    {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;
    }

    c_collect_cls( &nmaps, &mapnum, &nstokes, &nlmax, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, &nspec, cl, &myid, &numprocs, &root, &fcomm);

}

void collect_xls( s2hat_int4 nmaps1, s2hat_int4 mapnum1, s2hat_int4 nmaps2, s2hat_int4 mapnum2, s2hat_int4 nstokes, s2hat_int4 nlmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_dcomplex *local_alm1, 
                  s2hat_dcomplex *local_alm2, s2hat_int4 nspec, s2hat_flt8 *xcl, s2hat_int4 myid, s2hat_int4 numprocs, s2hat_int4 root, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    if( (lda == 1) || (lda == 2) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;
    } 
    else
      {
        ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;
      }


    c_collect_xls( &nmaps1, &mapnum1, &nmaps2, &mapnum2, &nstokes, &nlmax, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up,
                    local_alm1, local_alm2, &nspec, xcl, &myid, &numprocs, &root, &fcomm);

}

void s2hat_alm2map( s2hat_int4 precompute_plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps, s2hat_int4 nstokes, 
                    s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_int8 nplm, s2hat_flt8 *local_plm, s2hat_int4 numprocs, 
                    s2hat_int4 myid, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    /* no plms precomputed - default */
    ldim4_dwn = 1; ldim4_up = 1;
    ldim5_dwn = 1; ldim5_up = 1;

    if( (lda == 1) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

        if( precompute_plms == 1)
	{
            ldim4_dwn = 1; ldim4_up = 1;
            ldim5_dwn = 0; ldim5_up = nplm-1;
        }
        else
        {
	    if( precompute_plms == 2)
            {
		ldim4_dwn = 1; ldim4_up = nstokes;
		ldim5_dwn = 0; ldim5_up = nplm-1;
            }
        }

    } 
    else
    {
	ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

        if( precompute_plms == 1)
	{
            ldim4_dwn = 0; ldim4_up = nplm-1;
            ldim5_dwn = 1; ldim5_up = 1;
        }
        else
	{
            if( precompute_plms == 2)
	    {
                ldim4_dwn = 0; ldim4_up = nplm-1;
                ldim5_dwn = 1; ldim5_up = nstokes;
	    }
        }

    }

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    c_s2hat_alm2map( &precompute_plms, &f90pixelization, &f90scan, &nlmax, &nmmax, &nmvals, mvals, &nmaps, &nstokes, &first_ring, &last_ring,
         	     &map_size, local_map, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, 
                     &nplm, &ldim4_dwn, &ldim4_up, &ldim5_dwn, &ldim5_up, local_plm, &numprocs, &myid, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void s2hat_alm2map_spin( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 spin, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps, 
                         s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_int4 numprocs, 
                         s2hat_int4 myid, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    /* define the shape of the local_alm matrix */

    if( lda != 2)   /* i.e., s2hat convention */
    {
	ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = 2;
    }
    else            /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = 2;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;
    }  

    c_s2hat_alm2map_spin( &f90pixelization, &f90scan, &spin, &nlmax, &nmmax, &nmvals, mvals, &nmaps, &first_ring, &last_ring,
                          &map_size, local_map, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, &numprocs, &myid, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void s2hat_alm2mapderv_spin( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 spin, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps, 
                         s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_flt8 *local_derv, s2hat_int4 lda, s2hat_dcomplex *local_alm, 
                         s2hat_int4 numprocs, s2hat_int4 myid, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    /* define the shape of the local_alm matrix */

    if( lda != 2)   /* i.e., s2hat convention */
    {
	ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = 2;
    }
    else            /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = 2;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;
    }  

    c_s2hat_alm2mapderv_spin( &f90pixelization, &f90scan, &spin, &nlmax, &nmmax, &nmvals, mvals, &nmaps, &first_ring, &last_ring,
                              &map_size, local_map, local_derv, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, &numprocs, &myid, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void s2hat_map2alm( s2hat_int4 precompute_plms, s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps, 
                    s2hat_int4 nstokes, s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring,
		    s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_int8 nplm, s2hat_flt8 *local_plm, s2hat_int4 numprocs, s2hat_int4 myid, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up, ldim4_dwn, ldim4_up, ldim5_dwn, ldim5_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    /* define the shape of the local_alm matrix */

    /* no plms precomputed - default */
    ldim4_dwn = 1; ldim4_up = 1;
    ldim5_dwn = 1; ldim5_up = 1;

    if( (lda == 1) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = nstokes;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;

	if( precompute_plms == 1)
	{
            ldim4_dwn = 1; ldim4_up = 1;
            ldim5_dwn = 0; ldim5_up = nplm-1;
	}
        else
	{
	    if( precompute_plms == 2)
	    {
		ldim4_dwn = 1; ldim4_up = nstokes;
		ldim5_dwn = 0; ldim5_up = nplm-1;
	    }
	}

    } 
    else
    {
	ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = nstokes;

	if( precompute_plms == 1)
	{
            ldim4_dwn = 0; ldim4_up = nplm-1;
            ldim5_dwn = 1; ldim5_up = 1;
	}
        else
	{
            if( precompute_plms == 2)
            {
                ldim4_dwn = 0; ldim4_up = nplm-1;
                ldim5_dwn = 1; ldim5_up = nstokes;
            }
        }

    }

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    c_s2hat_map2alm( &precompute_plms, &cpixelization, &cscan, &nlmax, &nmmax, &nmvals, mvals, &nmaps, &nstokes, &first_ring, &last_ring, local_w8ring,
		     &map_size, local_map, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, &nplm, &ldim4_dwn, &ldim4_up, 
                     &ldim5_dwn, &ldim5_up, local_plm, &numprocs, &myid, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void s2hat_map2alm_spin( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 spin, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 nmaps, 
                         s2hat_int4 first_ring, s2hat_int4 last_ring, s2hat_flt8 *local_w8ring,
         		 s2hat_int4 map_size, s2hat_flt8 *local_map, s2hat_int4 lda, s2hat_dcomplex *local_alm, s2hat_int4 numprocs, s2hat_int4 myid, MPI_Comm comm)

{

    MPI_Fint fcomm;
    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up, ldim3_dwn, ldim3_up;

#if (defined(LAM_MPI) || defined(OPEN_MPI))
    fcomm = MPI_Comm_c2f( comm);
#else
    fcomm = comm;
#endif

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    /* define the shape of the local_alm matrix */

    if( lda != 2)   /* i.e., s2hat convention */
    {
	ldim1_dwn = 0; ldim1_up = nlmax;
        ldim2_dwn = 0; ldim2_up = nmvals-1;
        ldim3_dwn = 1; ldim3_up = 2;
    }
    else            /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = 2;
        ldim2_dwn = 0; ldim2_up = nlmax;
        ldim3_dwn = 0; ldim3_up = nmvals-1;
    }  

    c_s2hat_map2alm_spin( &cpixelization, &cscan, &spin, &nlmax, &nmmax, &nmvals, mvals, &nmaps, &first_ring, &last_ring, local_w8ring,
		          &map_size, local_map, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &ldim3_dwn, &ldim3_up, local_alm, &numprocs, &myid, &fcomm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

void  plm_mvalues_gen( s2hat_pixeltype cpixelization, s2hat_scandef cscan, s2hat_int4 npols, s2hat_int4 nlmax, s2hat_int4 nmmax, s2hat_int4 nmvals, s2hat_int4 *mvals, s2hat_int4 lda, s2hat_int8 nplm, s2hat_flt8 *local_plm)

{

    s2hat_pixeltype f90pixelization;
    s2hat_scandef f90scan;
    s2hat_int4 ldim1_dwn, ldim1_up, ldim2_dwn, ldim2_up;

    /* define the shape of the local_plm matrix */

    if( (lda == 1) || (lda == 3))           /* i.e., healpix convention */
    {
        ldim1_dwn = 1; ldim1_up = npols;
        ldim2_dwn = 0; ldim2_up = nplm-1;
    } 
    else
    {
        ldim1_dwn = 0; ldim1_up = nplm-1;
        ldim2_dwn = 1; ldim2_up = npols;
    }

    c2f_pixelization( &cpixelization, &f90pixelization);
    c2f_scan( &cpixelization, &cscan, &f90scan);

    c_plm_mvalues_gen( &f90pixelization, &f90scan, &npols, &nlmax, &nmmax, &nmvals, mvals, &ldim1_dwn, &ldim1_up, &ldim2_dwn, &ldim2_up, &nplm, local_plm);

    c_destroy_pixelization( &f90pixelization);
    c_destroy_scan( &f90scan);

}

/* F90 - to - C conversion routines */

s2hat_int4 f2c_pixelization( s2hat_pixeltype *f90pixelization, s2hat_pixeltype *cpixelization)

{
  s2hat_int8 ivec[4], nringsall, i;
  s2hat_int8 *ringInts;
  s2hat_flt8 *ringDoubs;

  nringsall = totalringsnumber( f90pixelization);

  ringInts = (s2hat_int8 *)calloc( 2*(s2hat_int4 )nringsall, sizeof( s2hat_int8));
  ringDoubs = (s2hat_flt8 *)calloc( 6*(s2hat_int4 )nringsall, sizeof( s2hat_flt8));

  pixelization2vectors( f90pixelization, &ivec[0], &ringInts[0], &ringDoubs[0]);
  
  (*cpixelization).type = ivec[0];
  (*cpixelization).npixsall = ivec[1];
  (*cpixelization).nringsall = ivec[2];
  (*cpixelization).nphmx = ivec[3];

  (*cpixelization).fpix = (s2hat_int8 *)calloc( nringsall, sizeof( s2hat_int8));
  (*cpixelization).nph = (s2hat_int8 *)calloc( nringsall, sizeof( s2hat_int8));
  (*cpixelization).kphi = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));
  (*cpixelization).qwght = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));

  (*cpixelization).pixphi = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));
  (*cpixelization).parea = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));
  (*cpixelization).cth = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));
  (*cpixelization).sth = (s2hat_flt8 *)calloc( nringsall, sizeof( s2hat_flt8));

  for( i = 0; i < nringsall; i++) 
  {
     (*cpixelization).fpix[i] = ringInts[i];
     (*cpixelization).nph[i] = ringInts[i+nringsall];

     (*cpixelization).kphi[i] = ringDoubs[i+4*nringsall];
     (*cpixelization).qwght[i] = ringDoubs[i+5*nringsall];

     (*cpixelization).pixphi[i] = ringDoubs[i];    
     (*cpixelization).parea[i] = ringDoubs[i+nringsall];
     (*cpixelization).cth[i] = ringDoubs[i+2*nringsall];
     (*cpixelization).sth[i] = ringDoubs[i+3*nringsall];
  }

  free( ringInts); free( ringDoubs);

  return( 1);
}

s2hat_int4 f2c_scan( s2hat_pixeltype *f90pixelization, s2hat_scandef *f90scan, s2hat_scandef *cscan)

{
  s2hat_int8 ivec[4], nringsall, i;
  s2hat_int8 *ringInts;

  nringsall = totalringsnumber( f90pixelization);

  ringInts = (s2hat_int8 *)calloc( 3*(s2hat_int4 )nringsall, sizeof( s2hat_int8));

  scan2vectors( f90pixelization, f90scan, &ivec[0], &ringInts[0]);
  
  (*cscan).npixsobs = ivec[0];
  (*cscan).nringsobs = ivec[1];

  (*cscan).nfl = (s2hat_int8 *)calloc( nringsall, sizeof( s2hat_int8));
  (*cscan).sfl = (s2hat_int8 *)calloc( nringsall, sizeof( s2hat_int8));
  (*cscan).fl = (s2hat_int8 *)calloc( nringsall, sizeof( s2hat_int8));

  for( i = 0; i < nringsall; i++) 
  {
     (*cscan).nfl[i] = ringInts[i];
     (*cscan).sfl[i] = ringInts[i+nringsall];
     (*cscan).fl[i] = ringInts[i+2*nringsall];
  }

  free( ringInts);

  return( 1);
}

/* C - to - F90 conversion routines */

s2hat_int4 c2f_pixelization( s2hat_pixeltype *cpixelization, s2hat_pixeltype *f90pixelization)

{
  s2hat_int8 ivec[4], nringsall, i;
  s2hat_int8 *ringInts;
  s2hat_flt8 *ringDoubs;

  nringsall = (*cpixelization).nringsall;

  ringInts = (s2hat_int8 *)calloc( 2*(s2hat_int4 )nringsall, sizeof( s2hat_int8));
  ringDoubs = (s2hat_flt8 *)calloc( 6*(s2hat_int4 )nringsall, sizeof( s2hat_flt8));
  
  ivec[0] = (*cpixelization).type;
  ivec[1] = (*cpixelization).npixsall;
  ivec[2] = (*cpixelization).nringsall;
  ivec[3] = (*cpixelization).nphmx;

  for( i = 0; i < nringsall; i++) 
  {
     ringInts[i] = (*cpixelization).fpix[i];
     ringInts[i+nringsall] = (*cpixelization).nph[i];

     ringDoubs[i+4*nringsall] = (*cpixelization).kphi[i];
     ringDoubs[i+5*nringsall] = (*cpixelization).qwght[i];

     ringDoubs[i] = (*cpixelization).pixphi[i];
     ringDoubs[i+nringsall] = (*cpixelization).parea[i];
     ringDoubs[i+2*nringsall] = (*cpixelization).cth[i];
     ringDoubs[i+3*nringsall] = (*cpixelization).sth[i];
  }

  vectors2pixelization( &nringsall, &ivec[0], &ringInts[0], &ringDoubs[0], f90pixelization);

  free( ringInts); free( ringDoubs);

  return( 1);
}

s2hat_int4 c2f_scan( s2hat_pixeltype *cpixelization, s2hat_scandef *cscan, s2hat_scandef *f90scan)

{
  s2hat_int8 ivec[4], nringsall, i;
  s2hat_int8 *ringInts;

  nringsall = (*cpixelization).nringsall;

  ringInts = (s2hat_int8 *)calloc( 3*(s2hat_int4 )nringsall, sizeof( s2hat_int8));

  ivec[0] = (*cscan).npixsobs;
  ivec[1] = (*cscan).nringsobs;

  for( i = 0; i < nringsall; i++) 
  {
     ringInts[i] = (*cscan).nfl[i];
     ringInts[i+nringsall] = (*cscan).sfl[i];
     ringInts[i+2*nringsall] = (*cscan).fl[i];
  }

  vectors2scan( &nringsall, &ivec[0], &ringInts[0], f90scan);
  
  free( ringInts);

  return( 1);
}

