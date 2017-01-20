ORIGIN = $(shell cd .. && pwd)
S2HAT_LIBPATH = $(ORIGIN)/s2hat/libs2hat.so

.PHONY: clean

all: libpath.jl

libpath.jl:
	echo "const S2HAT_LIBPATH = \"$(S2HAT_LIBPATH)\"" > libpath.jl

clean:
	rm -f libpath.jl

