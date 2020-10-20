BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

# Makefile for fscan

bin_PROGRAMS = lalapps_spec_avg lalapps_spec_avg_long

lalapps_spec_avg_SOURCES = spec_avg.c
lalapps_spec_avg_long_SOURCES = spec_avg_long.c

# Add shell test scripts to this variable
test_scripts += testFscan.sh
