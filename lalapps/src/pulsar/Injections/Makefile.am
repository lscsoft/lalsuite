## Process this file with automake to produce Makefile.in

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_makefakedata_v2 lalapps_Makefakedata_v4 lalapps_Makefakedata_v5 lalapps_compareSFTs lalapps_compareTS \
	lalapps_dumpSFT lalapps_makefakedatabinary lalapps_makefakedata_test

EXTRA_PROGRAMS = eah_Makefakedata_v4 eah_Makefakedata_v5

if LALFRAME
EXTRA_PROGRAMS += lalapps_sw_inj_frames
endif

lalapps_Makefakedata_v4_SOURCES = makefakedata_v4.c sft_extra.c sft_extra.h

lalapps_Makefakedata_v5_SOURCES = makefakedata_v5.c sft_extra.c sft_extra.h

# force linking in C++ mode for E@H
eah_Makefakedata_v4_SOURCES = $(lalapps_Makefakedata_v4_SOURCES) ../hough/src2/EinsteinAtHome/win_lib.cpp

eah_Makefakedata_v5_SOURCES = $(lalapps_Makefakedata_v5_SOURCES) ../hough/src2/EinsteinAtHome/win_lib.cpp

lalapps_makefakedata_v2_SOURCES =  makefakedata_v2.c sft_extra.c sft_extra.h

lalapps_compareSFTs_SOURCES = compareSFTs.c

lalapps_compareTS_SOURCES = compareTS.c

lalapps_dumpSFT_SOURCES = dumpSFT.c sft_extra.c sft_extra.h

lalapps_makefakedatabinary_SOURCES = makefakedatabinary.c

lalapps_makefakedata_test_SOURCES =  makefakedata_test.c sft_extra.c sft_extra.h

## ----- extra-programs
lalapps_sw_inj_frames_SOURCES = sw_inj_frames.c

EXTRA_DIST = testMakefakedata.sh testMFDv4.sh testMFDv5.sh testMFDv4-HW.sh testT8_1800
TESTS = testMakefakedata.sh testMFDv4.sh testMFDv4-HW.sh testMFDv5.sh

TESTS_ENVIRONMENT = \
	LAL_DEBUG_LEVEL=msglvl1,memdbg; export LAL_DEBUG_LEVEL; \
	LAL_DATA_PATH="$(LAL_DATA_PATH)"; export LAL_DATA_PATH;


## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
