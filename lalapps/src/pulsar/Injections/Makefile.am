MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_Makefakedata_v4 \
	lalapps_Makefakedata_v5 \
	lalapps_compareTS \
	$(END_OF_LIST)

if LALFRAME
EXTRA_PROGRAMS = lalapps_sw_inj_frames
endif

lalapps_Makefakedata_v4_SOURCES = makefakedata_v4.c

lalapps_Makefakedata_v5_SOURCES = makefakedata_v5.c

lalapps_compareTS_SOURCES = compareTS.c

lalapps_sw_inj_frames_SOURCES = sw_inj_frames.c

EXTRA_DIST = \
	testMFDv4.sh \
	testMFDv5.sh \
	testT8_1800 \
	$(END_OF_LIST)

TESTS = \
	testMFDv4.sh \
	testMFDv5.sh \
	$(END_OF_LIST)
