#!/bin/sh

message()
{
  echo "Re-install gmp and mpfr"
  exit 1
}

[ ! -f /usr/local/lib/libgmp.dylib ] && message
[ ! -f /usr/local/lib/libmpfr.dylib ] && message

ver_gmp_ok="14.2.0"
ver_mpfr_ok="6.6.0"

ver_gmp=`otool -L /usr/local/lib/libmpfr.dylib | grep libgmp | sed 's/[()]//g' | awk '{print$7}'`
ver_mpfr=`otool -L /usr/local/lib/libmpfr.dylib | grep libmpfr | tail -n 1 | sed 's/[()]//g' | awk '{print$7}'`

if [ "$ver_gmp" = "$ver_gmp_ok" ] ; then
  echo "[OK] GMP Version = $ver_gmp"
else
  echo "[NG] GMP Version = $ver_gmp < $ver_gmp_ok"
  message
fi

if [ "$ver_mpfr" = "$ver_mpfr_ok" ] ; then
  echo "[OK] MPFR Version = $ver_mpfr"
else
  echo "[NG] MPFR Version = $ver_mpfr < $ver_mpfr_ok"
  message
fi

