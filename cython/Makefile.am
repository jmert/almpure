EXTRA_DIST = \
	almpure_c.c \
	debug.h \
	almpure.pxd \
	almpure.pyx

.PHONY: almpure-cython

all-local: almpure-cython

almpure-cython:
	$(PYTHON) setup.py build_ext --inplace
