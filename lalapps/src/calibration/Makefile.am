## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_mkcalfac lalapps_mkcalref
else
FRAMEPROGRAMS = 
endif
bin_PROGRAMS = lalapps_ZPG2Transfer lalapps_BuildTransfer $(FRAMEPROGRAMS)
lalapps_ZPG2Transfer_SOURCES = ZPG2Transfer.c
lalapps_BuildTransfer_SOURCES = BuildTransfer.c
if FRAME
lalapps_mkcalfac_SOURCES = mkcalfac.c series.c series.h
lalapps_mkcalref_SOURCES = mkcalref.c series.c series.h
endif
INCLUDES = -I$(top_srcdir)/src
pkgdata_DATA = \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H2-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt
EXTRA_DIST = mkcalfac.c mkcalref.c series.c series.h mkcalfiles_test.sh \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt
MOSTLYCLEANFILES = *.ilwd *.gwf
TESTS = mkcalfiles_test.sh
