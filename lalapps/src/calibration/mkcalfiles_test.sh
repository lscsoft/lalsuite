#!/bin/sh

srcdir=${srcdir:-.}

if test -x lalapps_mkcalfac; then
  ./lalapps_mkcalfac -run S1 -ifo L1 $srcdir/S1-L1-CAL-FACTORS.txt || exit 1
fi
if test -x lalapps_mkcalref; then
  ./lalapps_mkcalref -run S1 -time 715388533 -ifo L1 $srcdir/S1-L1-CAL-RESPONSE.txt $srcdir/S1-L1-CAL-CAV_GAIN.txt || exit 1
fi
