ORIGIN := $(shell pwd)
TARGET := release

CFLAGS += -std=gnu99
ifeq ($(TARGET),release)
	MFLAGS = -O
else
	MFLAGS = -g
endif
MFLAGS += -largeArrayDims

S2HATPATH := $(shell readlink -m $(ORIGIN)/../s2hat)
CFLAGS += -I$(S2HATPATH)
LINKLIBS += -Wl,-Bstatic -L$(S2HATPATH) -ls2hat_fftw -Wl,-Bdynamic
LINKLIBS += -lgfortran

# Add in MPI flags
CFLAGS += $(shell mpifort --showme:compile)
LINKLIBS += $(shell mpifort --showme:link)

.PHONY: clean

all: mpihelper.mexa64 s2hat_alm2map_c.mexa64 s2hat_map2alm_c.mexa64 \
     s2hat_map2almpure_c.mexa64

%.mexa64: %.c
	mex $(MFLAGS) -output $@ CFLAGS='$$CFLAGS $(CFLAGS)' LINKLIBS='$$LINKLIBS $(LINKLIBS)' $^

clean:
	rm -f *.mexa64

