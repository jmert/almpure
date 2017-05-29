AC_INIT([almpure], [0.1])
AC_CONFIG_AUX_DIR([config])

dnl Take some tricks from
dnl   https://github.com/edrosten/autoconf_tutorial

    dnl Small macro to print out a nice, pretty section title.
    define(SECTION_TITLE,
    [
        echo >&AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
        echo '   $1' >&AS_MESSAGE_FD
        echo '   $1   ' | sed -e's/./-/g' >&AS_MESSAGE_FD
    ])



SECTION_TITLE([Configuring build environment])
AC_CONFIG_SRCDIR([s2hat/s2hat.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_SILENT_RULES([yes])
AM_PROG_AR
PKG_PROG_PKG_CONFIG


SECTION_TITLE([Configuring C])
AC_LANG_PUSH(C)
AC_PROG_CC
AC_LANG_POP(C)



SECTION_TITLE([Configuring Fortran])
AC_LANG_PUSH(Fortran)
AC_PROG_FC
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])
AC_FC_PP_DEFINE
AC_FC_MODULE_FLAG
AC_FC_MODULE_EXTENSION
AC_FC_LIBRARY_LDFLAGS
AM_CONDITIONAL([KNOWN_MODEXT], [test x"$FC_MODEXT" != x""])
AC_LANG_POP(Fortran)



AC_ARG_WITH([mpi-stub],
    [AS_HELP_STRING([--with-mpi-stub],[use a stub MPI implementation @<:@default=no@:>@])],
    [],
    [with_mpi_stub=no])
AM_CONDITIONAL([USE_MPISTUB], [test x"$with_mpi_stub" == xyes])
AS_IF([test x"$with_mpi_stub" == xyes],
    [
        WITH_MPISTUB=true
    ],
    [
        SECTION_TITLE([Configuring MPI])
        WITH_MPISTUB=false

        AC_LANG_PUSH(C)
        AX_MPI([], [AC_ERROR([A C MPI compiler is required])])
        CC="$MPICC"
        AC_LANG_POP(C)

        AC_LANG_PUSH(Fortran)
        AX_MPI([], [AC_ERROR([A Fortran MPI compiler is required])])
        FC="$MPIFC"
        AC_LANG_POP(Fortran)
    ])


dnl Configure libtool after the MPI compilers have been setup. Because of how
dnl libtool works, the linking compiler is set during configuration rather
dnl than as an environment variable during make. Therefore, the only way to
dnl get the correct MPI libraries linked into the shared libs2hat is to set
dnl the linking compilers *after* MPI setup.
SECTION_TITLE([Configuring libtool linking])
LT_INIT([shared static pic-only])



dnl For the rest of the configuring, go back to using a C compiler
AC_LANG(C)

SECTION_TITLE([Configuring FFT libraries])
PKG_CHECK_MODULES([FFTW3], [fftw3])
AC_SUBST([FFTW3_CFLAGS])
AC_SUBST([FFTW3_LIBS])
dnl For at least gfortran, a default include search path isn't used, so even
dnl if fftw3.f is in the system include directory, the following checks can
dnl still fail.
dnl
dnl FFTW3_CFLAGS will be empty if it's installed to the standard location, so
dnl in this case manually create a FFTW3_FCFLAGS variable that just gets the
dnl prefix from pkg-config. Otherwise, copy from FFTW3_CFLAGS.
AC_ARG_VAR([FFTW3_FCFLAGS], [Fortran compiler flags for FFTW3])
AS_IF([test -z "$FFTW3_FCFLAGS"], [
    dnl Try to construct FFTW3_FCFLAGS manually
    AS_IF([test -z "$FFTW3_CFLAGS"],
        [
            dnl Get the prefix from pkg-config and append "/include" to the
            dnl path as an initial guess.
            FFTW3_FCFLAGS=`$PKG_CONFIG --variable=prefix fftw3 2>/dev/null`
            FFTW3_FCFLAGS="$ac_cv_fc_module_flag$FFTW3_FCFLAGS/include"
        ],
        [
            dnl This probably isn't quite right since the Fortran compiler
            dnl doesn't necessarily share the flag format with the C compiler.
            FFTW3_FCFLAGS="$FFTW3_CFLAGS"
        ])
    ])
AC_SUBST([FFTW3_FCFLAGS])
AC_LANG_PUSH(Fortran)

dnl Make sure fftw3.f can be included
AC_MSG_CHECKING([for fftw3.f])
save_FCFLAGS="$FCFLAGS"
FCFLAGS="$save_FCFLAGS $FFTW3_FCFLAGS"
have_fftw3_f=no
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([], [include "fftw3.f"])],
    [have_fftw3_f=yes],
    [have_fftw3_f=no])
AC_MSG_RESULT([$have_fftw3_f])
AS_IF([test x"$have_fftw3_f" != xyes],
    [AC_MSG_ERROR([Could not find fftw3.f])])
FCFLAGS="$save_FCFLAGS"

dnl Now figure out which set of API bindings s2hat has available to work with.
WITH_FFT=unknown
with_fft_name="full complex"
AC_CHECK_LIB([fftw3],[dfftw_plan_dft_c2r],
    [WITH_FFT=FFTW3_C2R],
    [
        AC_CHECK_LIB([fftw3],[dfftw_plan_dft_hc2r],
            [WITH_FFT=FFTW3_HC2R])
        with_fft_name="half complex"
    ])
AS_IF([test x"$WITH_FFT" == x"unknown"],
    AC_MSG_ERROR([Unable to determine half- or full-complex support in libfftw3])])
AC_LANG_POP(Fortran)
AC_SUBST([WITH_FFT])



AC_ARG_ENABLE([matlab],
    [AS_HELP_STRING([--enable-matlab],[build Matlab MEX @<:@default=no@:>@])],
    [enable_matlab=$enableval], [enable_matlab=])
AS_IF([test -n "$enable_matlab"],
    [
        SECTION_TITLE([Configure Matlab])
        WITH_MATLAB=true
        AC_ARG_VAR([MEX], [Matlab MEX compiler command])
        AC_PATH_PROG([MEX], [mex], [not found])
        AS_IF([test x"$MEX" == x"not found"], [WITH_MATLAB=false])
    ],
    [
        WITH_MATLAB=false
    ])
AM_CONDITIONAL([HAVE_MATLAB], [$WITH_MATLAB])
AS_IF([$WITH_MATLAB], [mex_summary=yes], [mex_summary=no])

SECTION_TITLE([Configuration results])
AC_CONFIG_FILES([
    Makefile
    s2hat/Makefile
    ])
AS_IF([$WITH_MPISTUB],[AC_CONFIG_FILES([mpistub/Makefile])])
AS_IF([$WITH_MATLAB], [AC_CONFIG_FILES([matlab/Makefile])])
AC_OUTPUT
AS_ECHO("") >&AS_MESSAGE_FD

AS_ECHO(" ---- Compiler summary ----") >&AS_MESSAGE_FD
AS_ECHO("C compiler             = $CC") >&AS_MESSAGE_FD
AS_ECHO("Fortran compiler       = $FC") >&AS_MESSAGE_FD
AS_ECHO("C compiler flags       = $CFLAGS") >&AS_MESSAGE_FD
AS_ECHO("Fortran compiler flags = $FCFLAGS") >&AS_MESSAGE_FD
AS_ECHO("Linker flags           = $LDFLAGS") >&AS_MESSAGE_FD
AS_ECHO("Fortran linking flags  = $FCLIBS") >&AS_MESSAGE_FD
AS_ECHO("") >&AS_MESSAGE_FD

AS_ECHO(" ---- FFTW3 summary ----") >&AS_MESSAGE_FD
AS_ECHO("FFTW3 compiler flags = $FFTW3_CFLAGS") >&AS_MESSAGE_FD
AS_ECHO("FFTW3 Fortran flags  = $FFTW3_FCFLAGS") >&AS_MESSAGE_FD
AS_ECHO("FFTW3 linker flags   = $FFTW3_LIBS") >&AS_MESSAGE_FD
AS_ECHO("FFTW3 bindings style = $with_fft_name") >&AS_MESSAGE_FD
AS_ECHO("") >&AS_MESSAGE_FD

AS_ECHO(" ---- Bindings summary ----") >&AS_MESSAGE_FD
AS_ECHO("Building Matlab bindings = $mex_summary")
AS_ECHO("") >&AS_MESSAGE_FD
