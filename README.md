# libis
iSys Library.
This contains the functions using MPFR library, multi precision floating
point arithmetic. For examples, linear equation solver, eigenvalues
solver, machine interval operations, etc...

1. Install the latest version of GMP Library from https://gmplib.org.
2. Install the latest version of The GNU MPFR Library from http://www.mpfr.org.
3. Edit src/Makefile.
4. cd src ; make ; make install
5. Edit prog/Makefile.
6. cd prog ; make 
7. Edit matlab/@rmulti/private/Makefile
8. cd matlab/@rmulti/private ; make
9. Edit matlab/@cmulti/private/Makefile
10. cd matlab/@cmulti/private ; make
