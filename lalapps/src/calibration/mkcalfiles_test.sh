#!/bin/sh

srcdir=${srcdir:-.}

if test -x lalapps_mkcalfac; then
  ./lalapps_mkcalfac -run S1 -ifo L1 $srcdir/S1-L1-CAL-FACTORS.txt || exit 1
  ./lalapps_mkcalfac -run S1 -ifo H1 $srcdir/S1-H1-CAL-FACTORS.txt -ifo H2 $srcdir/S1-H2-CAL-FACTORS.txt || exit 1
fi
if test -x lalapps_mkcalref; then
  ./lalapps_mkcalref -run S1 -time 715388533 -ifo L1 $srcdir/S1-L1-CAL-RESPONSE.txt $srcdir/S1-L1-CAL-CAV_GAIN.txt || exit 1
  ./lalapps_mkcalref -run S1 -time 715156759 -ifo H1 $srcdir/S1-H1-CAL-RESPONSE.txt $srcdir/S1-H1-CAL-CAV_GAIN.txt || exit 1
  ./lalapps_mkcalref -run S1 -time 714644464 -ifo H2 $srcdir/S1-H2-CAL-RESPONSE.txt $srcdir/S1-H2-CAL-CAV_GAIN.txt || exit 1
fi
