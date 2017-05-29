SUFFIXES = .mexa64

MEXCC = @CC@
MEXFC = @FC@
MEXLD = @FC@
if USE_MPISTUB
MPIFLAGS = -I$(top_builddir)/mpistub -L$(top_builddir)/mpistub -lmpi_stubs
else
MPIFLAGS =
endif #USE_MPISTUB

.c.mexa64:
	$(MEX) -largeArrayDims -output $@ $^ \
		CC='$(MEXCC)' \
		FC='$(MEXFC)' \
		LD='$(MEXLD)' \
		CFLAGS='$$CFLAGS -std=gnu99' \
		CLIBS='$$CLIBS -Wl,-Bstatic -ls2hat @FFTW3_LIBS@ -Wl,-Bdynamic' \
		$(MPIFLAGS) \
		-I$(top_builddir)/s2hat -L$(top_builddir)/s2hat \
		@FFTW3_CFLAGS@

mexsrcs = \
	debug.h \
	mpihelper.c \
	s2hat_alm2map_c.c \
	s2hat_map2alm_c.c \
	s2hat_map2almpure_c.c
mexfuncs = \
	mpihelper.mexa64 \
	s2hat_alm2map_c.mexa64 \
	s2hat_map2alm_c.mexa64 \
	s2hat_map2almpure_c.mexa64
mfiles = \
	s2hat_alm2map.m \
	s2hat_map2alm.m \
	s2hat_map2almpure.m

# To keep the executable bits set during install, the MEX files are installed
# using the SCRIPTS type. We make ${mexdir} == ${pgkdatadir} so that the MEX
# and M-files will live side-by-side.
mexdir = $(pkgdatadir)
mex_SCRIPTS = $(mexfuncs)
# The M-files go to the standard package data directory.
pkgdata_DATA = $(mfiles)

# None of the targets yet get the sources included in the distributable, so
# manually tell Automake to include the sources.
EXTRA_DIST = $(mfiles) $(mexsrcs)

all-local: $(mexfuncs)

clean-local::
	rm -f *.mexa64

