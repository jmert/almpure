ORIGIN := $(shell pwd)

CFLAGS += -g -O0 -I./s2hat
LDFLAGS += -lm -L./s2hat -ls2hat -Wl,-rpath=$(ORIGIN)/s2hat

# Since libs2hat uses fortran, we need to pull in extra fortran requirements
# when compiling the C code for the link stage to work correctly.
CFLAGS += $(shell mpifort --showme:compile)
LDFLAGS += $(shell mpifort --showme:link)

all: s2hat delta_map

.PHONY: s2hat matlab clean cleanall

s2hat:
	@./setup.sh
	make -C s2hat

delta_map: delta_map.c | s2hat
	mpicc $(CFLAGS) -o $@ $^ $(LDFLAGS)

matlab: | s2hat
	make -C s2hat libs2hat_fftw.a
	make -C matlab

julia: | s2hat
	make -C julia

clean:
	rm -f delta_map
	make -C matlab clean
	make -C julia clean

cleanall: clean
	make -C s2hat clean

