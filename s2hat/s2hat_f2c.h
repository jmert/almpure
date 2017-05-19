
/* definitions for the F2C-interface       *
 * Radek Stompor (APC) February, 2007      *
 * spin transforms added rs@apc, Sep, 2007 *
 * some headers added rs@apc, Sep, 2008    *
 * extra flags defs added, rs@apc, Aug, 09 *
* vectors2scan fixed, rs@pc, Jun, 10       */

#include "s2hat_defs.h"

/* define the optmization */

#if defined(SP3OPT)

   #define OPTIMIZATION_SP3

#else

   #if defined( OPTERONOPT)

      #define OPTIMIZATION_OPTERON

   #else

      #if defined(SP3)

         #define OPTIMIZATION_SP3

      #else

         #define OPTIMIZATION_OPTERON       /* i.e., a default for all but SP3 like platforms */

      #endif

   #endif

#endif

/* define the mpi version to be used */

#if !defined(MPICH)
   #define OPEN_MPI   /* default */
#endif

/* the keywords definitions below and the function f77func taken from Chris Cantalupo's MADCAP3dataFortran.h */

#if defined(SP3)

#elif defined(T3E)

  #define F77_UPPERCASE_NAMES 

#elif defined(SX5)

  #define F77_APPEND_UNDERSCORE

#elif defined(O2K)

  #define F77_APPEND_UNDERSCORE

#elif defined(powermac)

#elif defined(opteron)

  #define F77_APPEND_UNDERSCORE

#elif defined(x86_32)

  #define F77_APPEND_UNDERSCORE

#elif defined(x86_64)

  #define F77_APPEND_UNDERSCORE

#endif

#ifndef f77func
#  if defined (F77_APPEND_UNDERSCORE)
#    if defined (F77_UPPERCASE_NAMES)
#      define f77func(f, F) F##_
#    else
#      define f77func(f, F) f##_
#    endif
#  else
#    if defined (F77_UPPERCASE_NAMES)
#      define f77func(f, F) F
#    else
#      define f77func(f, F) f
#    endif
#  endif
#endif

/* functional interface goes below here */

/* - internal conversion rooutines */

s2hat_int4 f2c_pixelization( s2hat_pixeltype*, s2hat_pixeltype*);
s2hat_int4 f2c_scan( s2hat_pixeltype*, s2hat_scandef*, s2hat_scandef*);
s2hat_int4 f2c_pixelization( s2hat_pixeltype*, s2hat_pixeltype*);
s2hat_int4 f2c_scan( s2hat_pixeltype*, s2hat_scandef*, s2hat_scandef*);
s2hat_int4 c2f_pixelization( s2hat_pixeltype *cpixelization, s2hat_pixeltype *f90pixelization);
s2hat_int4 c2f_scan( s2hat_pixeltype *cpixelization, s2hat_scandef *cscan, s2hat_scandef *f90scan);

#define pixelization2vectors f77func( pixelization2vectors, PIXELIZATION2VECTORS)
#define vectors2pixelization f77func( vectors2pixelization, VECTORS2PIXELIZATION)
#define scan2vectors f77func( scan2vectors, SCAN2VECTORS)
#define vectors2scan f77func( vectors2scan, VECTORS2SCAN)
#define totalringsnumber f77func( totalringsnumber, TOTALRINGSNUMBER)

extern void pixelization2vectors( void*, s2hat_int8*, s2hat_int8*, s2hat_flt8*);
extern void vectors2pixelization( s2hat_int8*, s2hat_int8*, s2hat_int8*, s2hat_flt8*, void*);
extern void scan2vectors( void*, void*, s2hat_int8*, s2hat_int8*);
extern void vectors2scan( s2hat_int8*, s2hat_int8*, s2hat_int8*, void*);
extern s2hat_int4 totalringsnumber( void*);

/* - Fortran interface */

/* -- pixelization/scan definitions and conversions */

#define c_set_pixelization f77func( f2c_set_pixelization, F2C_SET_PIXELIZATION)
#define c_zbounds2scan f77func( f2c_zbounds2scan, F2C_ZBOUNDS2SCAN)
#define c_zbounds2mask f77func( f2c_zbounds2mask, F2C_ZBOUNDS2MASK)
#define c_mask2scan f77func( f2c_mask2scan, F2C_MASK2SCAN)
#define c_destroy_pixelization f77func( f2c_destroy_pixelization, F2C_DESTROY_PIXELIZATION)
#define c_destroy_scan f77func( f2c_destroy_scan, F2C_DESTROY_SCAN)

extern void c_set_pixelization( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_pixeltype*);
extern void c_zbounds2scan( s2hat_flt8*, s2hat_pixeltype*, s2hat_scandef*);
extern void c_zbounds2mask( s2hat_flt8*, s2hat_pixeltype*, s2hat_int4*);
extern void c_mask2scan( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*);
extern void c_destroy_scan( s2hat_scandef*);
extern void c_destroy_pixelization( s2hat_pixeltype*);

/* -- FFT precomputations */

#define c_fft_setup f77func( f2c_fft_setup, F2C_FFT_SETUP)
#define c_fft_mc_setup f77func( f2c_fft_mc_setup, F2C_FFT_MC_SETUP)
#define c_fft_mc_clean f77func( f2c_fft_mc_clean, F2C_FFT_MC_CLEAN)

extern void c_fft_setup( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*);
extern void c_fft_mc_setup( s2hat_pixeltype*, s2hat_int4*);
extern void c_fft_mc_clean();

/* -- data size estimation routines */

#define c_nummmodes f77func( f2c_nummmodes, F2C_NUMMMODES)
#define c_nummvalues f77func( f2c_nummvalues, F2C_NUMMVALUES)
#define c_find_mvalues f77func( f2c_find_mvalues, F2C_FIND_MVALUES)
#define c_find_scan_ring_range f77func( f2c_find_scan_ring_range, F2C_FIND_SCAN_RING_RANGE)
#define c_get_local_data_sizes f77func( f2c_get_local_data_sizes, F2C_GET_LOCAL_DATA_SIZES)

extern s2hat_int4 c_nummmodes( s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern s2hat_int4 c_nummvalues( s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_find_mvalues( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_find_scan_ring_range( s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_get_local_data_sizes( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                    s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int8*, s2hat_int4*, s2hat_int4*);

/* -- data distribution/collection routine */

#define c_distribute_local_data_objects_map2alm f77func( f2c_distribute_local_data_objects_map2alm, F2C_DISTRIBUTE_LOCAL_DATA_OBJECTS_MAP2ALM)
#define c_distribute_local_data_objects_alm2map f77func( f2c_distribute_local_data_objects_alm2map, F2C_DISTRIBUTE_LOCAL_DATA_OBJECTS_ALM2ALM)

#define c_distribute_map f77func( f2c_distribute_map, F2C_DISTRIBUTE_MAP)
#define c_distribute_mask f77func( f2c_distribute_mask, F2C_DISTRIBUTE_MASK)
#define c_distribute_partialmap f77func( f2c_distribute_partialmap, F2C_DISTRIBUTE_PARTIALMAP)
#define c_distribute_alms f77func( f2c_distribute_alms, F2C_DISTRIBUTE_ALMS)
#define c_distribute_partialalms f77func( f2c_distribute_partialalms, F2C_DISTRIBUTE_PARTIALALMS)
#define c_distribute_w8ring f77func( f2c_distribute_w8ring, F2C_DISTRIBUTE_W8RING)
#define c_collect_map f77func( f2c_collect_map, F2C_COLLECT_MAP)
#define c_collect_partialmap f77func( f2c_collect_partialmap, F2C_COLLECT_PARTIALMAP)
#define c_collect_alms f77func( f2c_collect_alms, F2C_COLLECT_ALMS)
#define c_collect_partialalms f77func( f2c_collect_partialalms, F2C_COLLECT_PARTIALALMS)
#define c_collect_cls f77func( f2c_collect_cls, F2C_COLLECT_CLS)
#define c_collect_xls f77func( f2c_collect_xls, F2C_COLLECT_XLS)

/* Fortran interfaces called by the c-wrappers */

extern void c_distribute_local_data_objects_map2alm( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*,  s2hat_int4*, 
                                                     s2hat_flt8*, s2hat_int4*,  s2hat_int4*,  s2hat_int8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*,  
                                                     s2hat_int4*,  s2hat_int4*, s2hat_flt8*, s2hat_flt8*,  s2hat_int4*, s2hat_int4*,  s2hat_int4*, s2hat_int4*); 
extern void c_distribute_local_data_objects_alm2map( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                                     s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                                     s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int8*, s2hat_int4*, s2hat_int4*, 
                                                     s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_map( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_flt8*, s2hat_int4*, 
                              s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_mask( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                               s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_w8ring( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_alms( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                               s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, 
                               s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_partialmap( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int8*, s2hat_int4*, 
                                     s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_distribute_partialalms( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                      s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                      s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);						   
extern void c_collect_map( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, 
                           s2hat_int4*, s2hat_int4*);
extern void c_collect_alms( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                            s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, 
                            s2hat_int4*, s2hat_int4*);
extern void c_collect_partialmap( s2hat_pixeltype*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int8*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                  s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_collect_partialalms( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                   s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                   s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_collect_cls( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                           s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_collect_xls( s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                           s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_dcomplex*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*);

/* plm precomputation */

#define c_plm_mvalues_gen f77func( f2c_plm_mvalues_gen, F2C_PLM_MVALUES_GEN)
extern void c_plm_mvalues_gen( s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                               s2hat_int4*, s2hat_int8*, s2hat_flt8*);

/* transforms */

#define c_s2hat_alm2map f77func( f2c_s2hat_alm2map, F2C_S2HAT_ALM2MAP)
#define c_s2hat_alm2map_spin f77func( f2c_s2hat_alm2map_spin, F2C_S2HAT_ALM2MAP_SPIN)
#define c_s2hat_alm2mapderv_spin f77func( f2c_s2hat_alm2mapderv_spin, F2C_S2HAT_ALM2MAPDERV_SPIN)
#define c_s2hat_map2alm f77func( f2c_s2hat_map2alm, F2C_S2HAT_MAP2ALM)
#define c_s2hat_map2alm_spin f77func( f2c_s2hat_map2alm_spin, F2C_S2HAT_MAP2ALM_SPIN)

extern void c_s2hat_alm2map( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                             s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                             s2hat_dcomplex*, s2hat_int8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_s2hat_alm2map_spin( s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                  s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, 
                                  s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_s2hat_alm2mapderv_spin( s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                                  s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, 
                                  s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_s2hat_map2alm( s2hat_int4*, s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, 
                             s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, 
                             s2hat_int8*,s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*);
extern void c_s2hat_map2alm_spin( s2hat_pixeltype*, s2hat_scandef*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_flt8*, 
                                  s2hat_int4*, s2hat_flt8*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_int4*, s2hat_dcomplex*, s2hat_int4*, s2hat_int4*, s2hat_int4*);

