# Because of Fortran's mod files not being dependency trackable, we cannot
# support parallel builds in this directory. Disable parallel building.
.NOTPARALLEL:

# The libraries to build
noinst_LTLIBRARIES = libmpi_stubs.la
noinst_HEADERS = \
        mpi.h

# The files to add to the library and to the source distribution
libmpi_stubs_la_SOURCES = \
        $(libmpi_stubs_a_HEADERS) \
        mpi_stubs_f90.f90 \
        mpi_stubs_c.c

if KNOWN_MODEXT
clean-local::
	rm -f *.@FC_MODEXT@
endif #KNOWN_MODEXT

