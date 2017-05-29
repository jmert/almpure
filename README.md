# almpure

`almpure` is a set of bindings to the [S2HAT][1] and [PS2HAT][2] libraries.

## Getting Started

### Downloading

Either download a release tarball from
https://github.com/jmert/almpure/releases or clone the git repository:

```bash
$ git clone https://github.com/jmert/almpure
```

### Prerequisites

If starting from git sources, the following will be required to bootstrap
the build system:

- The GNU Autotools
- `pkg-config`

Generic build dependencies then include:

- A `make` build tool.
- An MPI library implementation.
- The FFTW3 library.
- C and Fortran compilers.

> At present, only OpenMPI has been tested, though other compliant MPI
> implementations should work as well.

Optional components include:

- Matlab to build MEX “all-in-one” function wrappers.
  * A static FFTW3 library is required for MEX compilation.

### Compiling

To bootstrap the build system from git sources,

```bash
$ autoreconf -vi
```

Next (or for a release tarball), configure and build the library. For example,
for home directory install without the Matlab MEX bindings:

```bash
$ ./configure --without-matlab --prefix=$HOME/local/almpure
$ make && make install
```

More options can be discovered via

```bash
$ ./configure --help
```

which will show other available configuration options and variables.

[1]: http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat.html
[2]: http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/pures2hat/pureS2HAT.html

