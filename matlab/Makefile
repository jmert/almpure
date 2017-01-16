ORIGIN := $(shell pwd)

CFLAGS += -g -O3 -std=gnu99

S2HATPATH := $(shell readlink -m $(ORIGIN)/../s2hat)
CFLAGS += -I$(S2HATPATH)
LDFLAGS += -L$(S2HATPATH) -ls2hat -Wl,-rpath=$(S2HATPATH)

# Add in MPI flags
CFLAGS += $(shell mpifort --showme:compile)
LDFLAGS += $(shell mpifort --showme:link)

.PHONY: clean

all: alm2map_c.mexa64

alm2map_c.mexa64: alm2map_c.c
	mex -output $@ CFLAGS='$$CFLAGS $(CFLAGS)' LDFLAGS='$$LDFLAGS $(LDFLAGS)' $^

clean:
	rm -f *.mexa64

