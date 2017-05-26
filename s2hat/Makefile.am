# Because of Fortran's mod files not being dependency trackable, we cannot
# support parallel builds in this directory. Disable parallel building.
.NOTPARALLEL:

# Generic build options
AM_CPPFLAGS =
AM_CFLAGS = @FFTW3_CFLAGS@
AM_FCFLAGS = @FFTW3_FCFLAGS@ -fno-second-underscore
AM_LDFLAGS = @FFTW3_LIBS@
# S2HAT options
AM_CFLAGS += -Dx86_64 -D@WITH_FFT@ -DOPEN_MPI
AM_FCFLAGS += $(DD)x86_64 $(DD)@WITH_FFT@ $(DD)OPEN_MPI

if USE_MPISTUB
AM_CPPFLAGS += -I$(top_builddir)/mpistub
AM_FCFLAGS += -I$(top_builddir)/mpistub
AM_LDFLAGS += -L$(top_builddir)/mpistub -lmpi_stubs
endif #USE_MPISTUB

# This is always OK since we set MPICC and MPIFC to CC and FC, respectively,
# for the MPI stub case.
CC = @MPICC@
FC = @MPIFC@
# Shortcut for how to pass define flags to the Fortran compiler
DD = @FC_DEFINE@

lib_LTLIBRARIES = libs2hat.la

# Install s2hat C headers and Fortran module files using custom target
# and path. (Fortran only if the module extension is known.)
cincdir = $(includedir)/s2hat
fmoddir = $(includedir)/s2hat


cinc_HEADERS = \
	s2hat.h \
	s2hat_defs.h \
	s2hat_f2c.h \
	s2hat_types.h

if KNOWN_MODEXT
MOD = @FC_MODEXT@
nodist_fmod_HEADERS = \
	s2hat_alm2map_mod.$(MOD) \
	s2hat_defs.$(MOD) \
	s2hat_map2alm_mod.$(MOD) \
	s2hat_pixelization.$(MOD) \
	s2hat_toolbox_mod.$(MOD) \
	s2hat_types_internal.$(MOD)

mostlyclean-local:: mostlyclean-mods
mostlyclean-mods:
	rm -f *.$(MOD)

endif #KNOWN_MODEXT

libs2hat_la_DEPENDENCIES = $(cinc_HEADERS)

libs2hat_la_SOURCES = \
	s2hat_defs.f90 \
	s2hat_types_internal.F90 \
	s2hat_pixelization.f90 \
	s2hat_toolbox.F90 \
	s2hat_alm2map.f90 \
	s2hat_map2alm.F90 \
	s2hat_c_interface.f90 \
	s2hat_c_wrappers.c \
	s2hat_map2purealm.c

# Generate links for the built libraries automatically.
include $(top_srcdir)/convenience-link.mk
