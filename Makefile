ORIGIN := $(shell pwd)

CFLAGS += -g -O0 -I./s2hat
LDFLAGS += -lm -L./s2hat -ls2hat -Wl,-rpath=$(ORIGIN)/s2hat

# Since libs2hat uses fortran, we need to pull in extra fortran requirements
# when compiling the C code for the link stage to work correctly.
CFLAGS += $(shell mpifort --showme:compile)
LDFLAGS += $(shell mpifort --showme:link)

all: delta_map

.PHONY: s2hat clean cleanall

s2hat:
	./setup.sh

delta_map: delta_map.c | s2hat
	mpicc $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f delta_map

cleanall: clean
	make -C s2hat clean

