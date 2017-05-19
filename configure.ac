AC_INIT(almpure, version-0.1)
AC_CONFIG_HEADERS([config.h])

m4_include([m4/ax_mpi.m4])

AC_PROG_CPP
AC_PROG_RANLIB

dnl Take some tricks from
dnl   https://github.com/edrosten/autoconf_tutorial

    dnl Force the compiler to run once and do all the basic checks
    dnl if you don't, it will do the test on the first invocation
    dnl below and so your pretty section titles won't work as well
    AC_CHECK_HEADERS(string.h)

    dnl Small macro to print out a nice, pretty section title.
    define(SECTION_TITLE,
    [
        echo >&AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
        echo '   $1' >&AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
    ])


SECTION_TITLE([Configuring Fortran MPI])
AC_LANG(Fortran)
AX_MPI([], [AC_ERROR([A Fortran MPI compiler is required])])



SECTION_TITLE([Configuring C MPI])
AC_LANG(C)
AX_MPI([], [AC_ERROR([A C MPI compiler is required])])



SECTION_TITLE([Configuring FFT libraries])
AC_ARG_WITH([healpix],
    [AS_HELP_STRING([--with-healpix],[use the HEALPix FFT library @<:@default=yes@:>@])],
    [],
    [with_healpix=yes])
AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3],[use the FFTW3 FFT library @<:@default=yes@:>@])],
    [],
    [with_fftw3=yes])

# Check for appropriate libraries given the choice of FFT package requested.
# By default we check for all available options, but prefer to use libfftw3
# with full complex support by overriding in order each time a check
# succeeds.

# The S2HAT preprocessor define to use
WITH_FFT=
# Appropriate flags for each choice of FFT library
FFT_CFLAGS=
FFT_LIBS=

have_healpix=no
AS_IF([test x"$with_healpix" != xno],
    [AC_CHECK_LIB([healpix],[healpix_fft],
        [have_healpix=yes
         FFT_LIBS=-lhealpix
         WITH_FFT=HEALPIX_fft
         AC_DEFINE([HAVE_HEALPIX], [1], [Define if you have libhealpix])
        ],
        [])],
    [])
have_fftw3=no
have_fftw3_c2r=no
have_fftw3_hc2r=no
AS_IF([test x"$with_fftw3" != xno],
    [
     AC_CHECK_LIB([fftw3],[fftw_plan_dft_hc2r],
        [have_fftw3=yes
         have_fftw3_hc2r=yes
         FFT_LIBS=-lfftw3
         WITH_FFT=FFTW3_HC2R
         AC_DEFINE([HAVE_FFTW3], [1], [Define if you have libfftw3])
        ],
        [])
     AC_CHECK_LIB([fftw3],[fftw_plan_dft_c2r],
        [have_fftw3=yes
         have_fftw3_c2r=yes
         FFT_LIBS=-lfftw3
         WITH_FFT=FFTW3_C2R
         AC_DEFINE([HAVE_FFTW3], [1], [Define if you have libfftw3])
        ],
        [])
    ],
    [])

AS_IF([test x"$have_healpix" == xno -a x"$have_fftw3" == xno],
    [AC_MSG_ERROR([no FFT library enabled and-or found])],
    [])
AC_SUBST([FFT_CFLAGS])
AC_SUBST([FFT_LIBS])
AC_SUBST([WITH_FFT])


SECTION_TITLE([Configuration results])
AC_CONFIG_FILES([s2hat/Makefile])
AC_OUTPUT

echo "CC=$CC" >&AS_MESSAGE_FD
echo "MPICC=$MPICC" >&AS_MESSAGE_FD
echo "CFLAGS=$CFLAGS $FFT_CFLAGS" >&AS_MESSAGE_FD
echo "FC=$FC" >&AS_MESSAGE_FD
echo "MPIFC=$MPIFC" >&AS_MESSAGE_FD
echo "FCFLAGS=$FCFLAGS" >&AS_MESSAGE_FD
echo "MPILIBS=$MPILIBS" >&AS_MESSAGE_FD
echo "LDFLAGS=$LDFLAGS" >&AS_MESSAGE_FD
echo "LIBS=$FFT_LIBS" >&AS_MESSAGE_FD

