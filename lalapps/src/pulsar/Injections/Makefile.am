## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
frameprogs = lalapps_sw_inj_frames
else
frameprogs =
endif

bin_PROGRAMS = lalapps_makefakedata_v2 lalapps_Makefakedata_v4 lalapps_compareSFTs lalapps_compareTS \
	lalapps_dumpSFT lalapps_makefakedatabinary lalapps_makefakedata_test

EXTRA_PROGRAMS = eah_Makefakedata_v4 $(frameprogs)

lalapps_Makefakedata_v4_SOURCES = makefakedata_v4.c sft_extra.c sft_extra.h

# force linking in C++ mode for E@H
eah_Makefakedata_v4_SOURCES = $(lalapps_Makefakedata_v4_SOURCES) ../hough/src2/EinsteinAtHome/win_lib.cpp

lalapps_makefakedata_v2_SOURCES =  makefakedata_v2.c sft_extra.c sft_extra.h

lalapps_compareSFTs_SOURCES = compareSFTs.c

lalapps_compareTS_SOURCES = compareTS.c

lalapps_dumpSFT_SOURCES = dumpSFT.c sft_extra.c sft_extra.h

lalapps_makefakedatabinary_SOURCES = makefakedatabinary.c

lalapps_makefakedata_test_SOURCES =  makefakedata_test.c sft_extra.c sft_extra.h

## ----- extra-programs
lalapps_sw_inj_frames_SOURCES = sw_inj_frames.c

EXTRA_DIST = testMakefakedata.sh testMFDv4.sh testMFDv4-HW.sh testT8_1800
TESTS = testMakefakedata.sh testMFDv4.sh testMFDv4-HW.sh

TESTS_ENVIRONMENT = LALPULSAR_DATADIR="$(LALPULSAR_DATADIR)"


## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
