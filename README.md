# libis
iSys Library.
This contains the functions using MPFR library, multi precision floating
point arithmetic. For examples, linear equation solver, eigenvalues
solver, machine interval operations, etc...

1. Install the latest version of GMP Library from https://gmplib.org.
     # tar xvfj gmp-ver.tar.bz2
     # cd gmp-ver
     # ./configure -prefix=/usr/local
     # make
     # sudo make install
2. Install the latest version of The GNU MPFR Library from http://www.mpfr.org.
     # tar xvfj mpfr-ver.tar.bz2
     # cd mpr-ver
     # ./configure -prefix=/usr/local
     # make
     # sudo make install
3. Edit src/Makefile.
4. cd src ; make ; make install
5. Edit prog/Makefile.
6. cd prog ; make 
7. Edit matlab/@multi/Makefile
8. cd matlab/@multi ; make
