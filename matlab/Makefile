ORIGIN := $(shell pwd)
TARGET := release

CFLAGS += -std=gnu99
ifeq ($(TARGET),release)
	MFLAGS = -O
else
	MFLAGS = -g
endif

S2HATPATH := $(shell readlink -m $(ORIGIN)/../s2hat)
CFLAGS += -I$(S2HATPATH)
LDFLAGS += -L$(S2HATPATH) -ls2hat -Wl,-rpath=$(S2HATPATH)

# Add in MPI flags
CFLAGS += $(shell mpifort --showme:compile)
LDFLAGS += $(shell mpifort --showme:link)

.PHONY: clean

all: alm2map_c.mexa64 map2alm_c.mexa64 map2almpure_c.mexa64

%.mexa64: %.c
	mex $(MFLAGS) -output $@ CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' $^

clean:
	rm -f *.mexa64

