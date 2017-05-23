# The libraries to build
noinst_LIBRARIES = libmpi_stubs.a
noinst_HEADERS = \
        mpi.h

# The files to add to the library and to the source distribution
libmpi_stubs_a_SOURCES = \
        $(libmpi_stubs_a_HEADERS) \
        mpi_stubs_f90.f90 \
        mpi_stubs_c.c

clean:
	rm -f *.o *.mod *.a
