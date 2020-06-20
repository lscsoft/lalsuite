BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = \
	lalapps_Makefakedata_v4 \
	lalapps_Makefakedata_v5 \
	$(END_OF_LIST)

if LALFRAME
EXTRA_PROGRAMS = lalapps_sw_inj_frames
endif

lalapps_Makefakedata_v4_SOURCES = makefakedata_v4.c

lalapps_Makefakedata_v5_SOURCES = makefakedata_v5.c

lalapps_sw_inj_frames_SOURCES = sw_inj_frames.c

# Add shell test scripts to this variable
test_scripts += testMakefakedata_v4.sh
test_scripts += testMakefakedata_v4_hwinjection.sh
test_scripts += testMakefakedata_v5.sh

# Add any helper programs required by tests to this variable
test_helpers += compareTS
