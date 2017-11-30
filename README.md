# iSystem Library (libis)
iSys Library provides some functions, which runs on multi precision floating point arithmetic.

## Description
This contains the functions using MPFR library, multi precision floating point arithmetic. For examples, linear equation solver, eigenvalues solver, machine interval operations, etc...

## Requirement
* The latest version of GMP Library from <https://gmplib.org>.
* The latest version of The GNU MPFR Library from <http://www.mpfr.org>.

## Installation
1. Install the latest version of GMP Library from <https://gmplib.org>.  
   `# tar xvfj gmp-version.tar.bz2`  
   `# cd gmp-version`  
   `# ./configure --prefix=/usr/local`  
   `# make`  
   `# sudo make install`  
2. Install the latest version of The GNU MPFR Library from <http://www.mpfr.org>.  
   `# tar xvfj mpfr-version.tar.bz2`  
   `# cd mpr-version`  
   `# ./configure --prefix=/usr/local`  
   `# make`  
   `# sudo make install`  
3. Compile libis.  
   `# cd libis/src`  
   `# make`  
   `# make install`  
4. Compile libs/prog if you need them.  
   `# cd libis/prog`  
   `# make`  
   `# make install`  
5. Compile libs/matlab if you need them.  
   Exec MATLAB.  
   Change directory to libis/matlab.  
   Run the following commands on MATLAB.  
   `# mex -setup C`  
   `# make`  
## Easy Installation
   Run iSystem Toolbox.mltbx
Operating Environment
1 macOS Sierra Ver.10.12.6
  Mathematica R2017b
