# almpure

`almpure` is a set of bindings to the [S2HAT][1] and [PS2HAT][2] libraries.

## Getting Started

### Downloading
Clone this repository to a local directory:

```bash
$ git clone https://github.com/jmert/almpure
```

### Compiling

An MPI compiler is required to compile the (P)S2HAT library. At present, only
OpenMPI has been tested, though other compliant MPI implementations should
work as well.

On some systems (such as Harvard's Odyssey cluster), this software is provided
in loadable modules. The support directory contains small shell scripts which
can be sourced in the shell to load up the correct environment. For example,
on Odyssey,

```bash
$ cd almpure
$ source support/odyssey_env.sh
```

Then configure and compile the library:

```bash
$ autoreconf -vi
$ ./configure
$ make
```

### Compiling Matlab bindings

To build the Matlab MEX bindings, the following conditions must be met:

1. Matlab available with a MEX compiler.
2. A static FFTW library. (The dynamic FFTW library used when compiling
   the PS2HAT library conflicts with Matlab's internal library, so use
   of the static library instead works around this.)

The bindings can be built with

```bash
$ make matlab
```

> At this time, only use with Matlab 2016b has been tested.

[1]: http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/s2hat/s2hat.html
[2]: http://www.apc.univ-paris7.fr/APC_CS/Recherche/Adamis/MIDAS09/software/pures2hat/pureS2HAT.html

