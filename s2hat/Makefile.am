CC = @MPICC@
FC = @MPIFC@
DD = @FC_DEFINE@

# Generic build options
AM_CFLAGS =
AM_FCFLAGS = -fno-second-underscore
AM_LDFLAGS = @FCLIBS@
# S2HAT options
AM_CFLAGS += -Dx86_64 -D@WITH_FFT@ -DOPEN_MPI
AM_FCFLAGS += $(DD)x86_64 $(DD)@WITH_FFT@ $(DD)OPEN_MPI

lib_LTLIBRARIES = libs2hat.la

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
endif

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

