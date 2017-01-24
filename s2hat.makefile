# Choose CPU platform:
#   SP3, T3E, SX5, O2K, opteron, powermac, x86_32, x86_64
PLATFORM := x86_64
# Choose FFT package:
#   ACML_fft, HEALPIX_fft, FFTW3_HC2R, FFTW3_C2R
FFTPKG := FFTW3_C2R
# Choose the MPI library: (Optional)
#   MPICH, OPEN_MPI, LAM_MPI
MPILIB := OPEN_MPI

AR = ar
FC = mpif90
CC = mpicc

# Set compiler defines for the choices set above.
DFLAGS = -g -D$(PLATFORM) -D$(FFTPKG)
ifdef MPILIB
	DFLAGS += -D$(MPILIB)
endif

# Set the build flags used in compile rules below
FFLAGS += -O3 -std=gnu -fPIC $(DFLAGS) \
	-fno-second-underscore -ffixed-line-length-none -ffree-line-length-none
CFLAGS += -O3 -fPIC $(DFLAGS)
LDFLAGS += -lgfortran

# If one of the FFTW packages was chosen, automatically insert the correct
# compilation flags retrieved from pkgconfig
ifneq (,$(filter $(FFTPKG),FFTW3_HC2R FFTW3_C2R))
	FFLAGS += $(shell pkg-config --cflags fftw3)
	CFLAGS += $(shell pkg-config --cflags fftw3)
	LDFLAGS += $(shell pkg-config --libs fftw3)
	FFTW_STATIC = $(shell $(CC) $(LDFLAGS) --print-file-name=libfftw3.a)
endif

# Ordered list of objects to compile. The order is important since the Fortran
# module files from a prior compile may be required by a later compilation
# unit.
OBJECTS = \
	s2hat_defs.o s2hat_types_internal.o s2hat_pixelization.o \
	s2hat_toolbox.o s2hat_alm2map.o s2hat_map2alm.o \
	s2hat_c_interface.o s2hat_c_wrappers.o \
	s2hat_map2purealm.o

.PHONY: clean

all: libs2hat.a libs2hat.so

libs2hat.a: $(OBJECTS)
	$(AR) rcs $@ $^

# Target required for Matlab's mex functions which includes the FFTW functions
# statically so that the Matlab-internal functions aren't used instead.
libs2hat_fftw.a: libs2hat.a
	echo -e "CREATE $@\nADDLIB $^\nADDLIB $(FFTW_STATIC)\nSAVE\nEND" | ar -M
	ranlib $@

libs2hat.so: $(OBJECTS)
	$(CC) -o $@ -shared $(CFLAGS) $^ $(LDFLAGS)

clean:
	rm -f *.so *.a *.o *.mod

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

