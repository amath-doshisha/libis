#!/bin/bash

if [ -d /System ] ; then
  ./check_mpfr.sh || exit
fi

echo "====================================================== cleanig..."
make clean

echo "====================================================== making..."
make

echo "====================================================== installing..."
make install

echo "====================================================== done."
exit

#EOF
