#!/bin/bash

OS=""
if [ -d /System ] ; then
  OS=MACOSX
else
  OS=LINUX
fi

echo "====================================================== creating Makefile for ${OS}..."
rm -v Makefile.inc
ln -sv Makefile.inc.$OS++ Makefile.inc

echo "====================================================== cleanig..."
make clean

echo "====================================================== making..."
make -j 4

echo "====================================================== installing..."
make install

echo "====================================================== done."
exit

#EOF
