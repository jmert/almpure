
/* definitions for the C-interface         */
/* Radek Stompor (APC) February, 2007      */
/* spin transforms added rs@apc, Aug, 2007 */

#include "s2hat_f2c.h"

/* s2hat C routines */

/* MPI communication routines */

s2hat_int4 MPI_scanBcast( s2hat_pixeltype, s2hat_scandef*, s2hat_int4, s2hat_int4, MPI_Comm);
s2hat_int4 MPI_pixelizationBcast( s2hat_pixeltype*, s2hat_int4, s2hat_int4, MPI_Comm);

/* data sizes */

s2hat_int4 nummmodes( s2hat_int4, s2hat_int4, s2hat_int4*);
s2hat_int4 nummvalues( s2hat_int4, s2hat_int4, s2hat_int4);
void find_mvalues( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*);
void find_scan_ring_range( s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4*, s2hat_int4*);
void get_local_data_sizes( s2hat_int4, s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                           s2hat_int8*, s2hat_int4, MPI_Comm);

/* pixelization/scan routines */

void set_pixelization( s2hat_int4, s2hat_pixparameters, s2hat_pixeltype*);
void zbounds2scan( s2hat_flt8*, s2hat_pixeltype, s2hat_scandef*);
void zbounds2mask( s2hat_flt8*, s2hat_pixeltype, s2hat_int4*);
void mask2scan(s2hat_int4*, s2hat_pixeltype, s2hat_scandef*);
void destroy_pixelization( s2hat_pixeltype);
void destroy_scan( s2hat_scandef);

/* fft precomputation routines */

void fft_setup( s2hat_pixeltype, s2hat_int4, s2hat_int4);
void fft_mc_setup( s2hat_pixeltype, s2hat_int4);
void fft_mc_clean();

/* data distribution/collection */

void distribute_local_data_objects_map2alm( s2hat_int4,  s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*,  s2hat_int4, s2hat_flt8*,
                                            s2hat_int4,  s2hat_int4*,  s2hat_int8, s2hat_int4, s2hat_flt8*,  s2hat_int4,  s2hat_int4, s2hat_flt8*, s2hat_flt8*,  s2hat_int4, 
                                            s2hat_int4,  s2hat_int4, MPI_Comm); 
void distribute_local_data_objects_alm2map( s2hat_int4, s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_dcomplex*, s2hat_int4, 
                                            s2hat_int4*, s2hat_dcomplex*, s2hat_int8, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void distribute_map( s2hat_pixeltype, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void distribute_partialmap( s2hat_pixeltype, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int8, s2hat_int4, s2hat_flt8*, s2hat_int4, 
                            s2hat_int4, s2hat_int4, MPI_Comm);
void distribute_w8ring( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void distribute_alms( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_dcomplex*, s2hat_int4, s2hat_int4, 
                      s2hat_int4, MPI_Comm);
void distribute_partialalms( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_int4, s2hat_int4, s2hat_dcomplex*, 
                             s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void collect_map( s2hat_pixeltype, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void collect_partialmap( s2hat_pixeltype, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int8, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, 
                         s2hat_int4, s2hat_int4, MPI_Comm);
void collect_alms( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_dcomplex*, s2hat_int4, s2hat_int4, 
                   s2hat_int4, MPI_Comm);
void collect_partialalms( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_int4, s2hat_int4, s2hat_dcomplex*, 
                          s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);
void collect_cls( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_int4, 
                  s2hat_int4, MPI_Comm);
void collect_xls( s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_dcomplex*, s2hat_dcomplex*, s2hat_int4, s2hat_flt8*, 
                  s2hat_int4, s2hat_int4, s2hat_int4, MPI_Comm);

/* transform routines */

void s2hat_alm2map( s2hat_int4, s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, 
                    s2hat_int4, s2hat_dcomplex*, s2hat_int8, s2hat_flt8*, s2hat_int4, s2hat_int4, MPI_Comm);
void s2hat_alm2map_spin( s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, 
                    s2hat_dcomplex*, s2hat_int4, s2hat_int4, MPI_Comm);
void s2hat_alm2mapderv_spin( s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, 
		    s2hat_flt8*, s2hat_int4, s2hat_dcomplex*, s2hat_int4, s2hat_int4, MPI_Comm);
void s2hat_map2alm( s2hat_int4, s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, 
                    s2hat_flt8*, s2hat_int4, s2hat_dcomplex*, s2hat_int8, s2hat_flt8*, s2hat_int4, s2hat_int4, MPI_Comm);
void s2hat_map2alm_spin( s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_flt8*, s2hat_int4, s2hat_flt8*, 
                         s2hat_int4, s2hat_dcomplex*, s2hat_int4, s2hat_int4, MPI_Comm);

/* plm calculation */

void plm_mvalues_gen( s2hat_pixeltype, s2hat_scandef, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4, s2hat_int4*, s2hat_int4, s2hat_int8, s2hat_flt8*);
