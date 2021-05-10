## g_lomepro ##

Local membrane property analysis

## Installation ##

1. Try if the statically linked version `g_lomepro_static` works for you. If so, you don't need to compile the software yourself; simply use `g_lomepro_static`

2. If you nevertheless need to build the binary yourself, firstly, adjust the paths in the Makefile. There are 5 paths to be entered:
 - SRCTOOLS: gromacs tools
 - INCLUDE: gromacs headers
 - GMXLIB: compiled gromacs libraries
 - FFTW: compiled fft library path
 - FFTWINCL: fft headers

Afterwards, run the command
```
make
```

Note, to compile `g_lomepro` it will be necessary to use gmx4.6 or older version, yet this is only a requirement for compilation: the software itself is compatible with the latest gromacs version (see below).

## Is `g_lomepro` compatible with the latest gromacs version? ##
Yes.

Simply, when supplying the structure file (option -s) use a .pdb instead of .tpr.
