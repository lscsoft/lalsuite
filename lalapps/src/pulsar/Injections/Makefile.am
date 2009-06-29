## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/misc
#LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID


bin_PROGRAMS = lalapps_makefakedata_v2 lalapps_Makefakedata_v4 lalapps_compareSFTs lalapps_compareTS \
		lalapps_dumpSFT lalapps_makefakedatabinary lalapps_makefakedata_test

EXTRA_PROGRAMS = psinject

lalapps_Makefakedata_v4_SOURCES = makefakedata_v4.c sft_extra.c sft_extra.h
lalapps_Makefakedata_v4_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_makefakedata_v2_SOURCES =  makefakedata_v2.c sft_extra.c sft_extra.h
lalapps_makefakedata_v2_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_compareSFTs_SOURCES = compareSFTs.c
lalapps_compareSFTs_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_compareTS_SOURCES = compareTS.c
lalapps_compareTS_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_dumpSFT_SOURCES = dumpSFT.c sft_extra.c sft_extra.h
lalapps_dumpSFT_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_makefakedatabinary_SOURCES = makefakedatabinary.c
lalapps_makefakedatabinary_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

lalapps_makefakedata_test_SOURCES =  makefakedata_test.c sft_extra.c sft_extra.h
lalapps_makefakedata_test_LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## extra-programs
psinject_SOURCES = psinject.c

EXTRA_DIST = testMakefakedata.sh testMFDv4.sh testMFDv4-HW.sh testT8_1800
TESTS = testMakefakedata.sh testMFDv4.sh testMFDv4-HW.sh

CLEANFILES =
