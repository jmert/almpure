AC_INIT(almpure, version-0.1)

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



SECTION_TITLE([Configuration results])
AC_CONFIG_FILES([s2hat.makefile])
AC_OUTPUT

echo "CC=$CC" >&AS_MESSAGE_FD
echo "MPICC=$MPICC" >&AS_MESSAGE_FD
echo "CFLAGS=$CFLAGS" >&AS_MESSAGE_FD
echo "FC=$FC" >&AS_MESSAGE_FD
echo "MPIFC=$MPIFC" >&AS_MESSAGE_FD
echo "FCFLAGS=$FCFLAGS" >&AS_MESSAGE_FD
echo "MPILIBS=$MPILIBS" >&AS_MESSAGE_FD
echo "LDFLAGS=$LDFLAGS" >&AS_MESSAGE_FD
echo "LIBS=$LIBS" >&AS_MESSAGE_FD

