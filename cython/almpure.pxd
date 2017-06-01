cdef extern from "<s2hat.h>":
    ctypedef struct s2hat_dcomplex:
        double re
        double im

cdef extern from "almpure_c.c":
    int alm2map_c(double* map, int nside, s2hat_dcomplex* alms, int nstokes,
                int nlmax, int nmmax, int nmaps)

    int map2alm_c(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
                int nmaps, double* map, int nside, double* qwghts)

    int map2almpure_c(s2hat_dcomplex* alms, int nstokes, int lmax, int mmax,
                    int nmaps, double* map, double* apmask, int nside,
                    double* qwghts)
