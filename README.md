# libis
iSys Library.

## Description

This contains the functions using MPFR library, multi precision floating
point arithmetic. For examples, linear equation solver, eigenvalues
solver, machine interval operations, etc...

## Installation
1. Install the latest version of GMP Library from https://gmplib.org.
     # tar xvfj gmp-version.tar.bz2
     # cd gmp-version
     # ./configure -prefix=/usr/local
     # make
     # sudo make install
2. Install the latest version of The GNU MPFR Library from http://www.mpfr.org.
     # tar xvfj mpfr-version.tar.bz2
     # cd mpr-version
     # ./configure -prefix=/usr/local
     # make
     # sudo make install
3. Compile libis.
     # cd libis/src
     # make
     # make install
4. Compile libs/prog.
     # cd libis/prog
     # make
     # make install
5. Compile libs/matlab
     # cd matlab
     # mex -I/usr/local/include -I../include -L/usr/local/lib -L../MACOSX -lis -lmpfr -lgmp multi_mex.c

#EOF