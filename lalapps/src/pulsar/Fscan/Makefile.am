BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_spec_avg \
	lalapps_spec_avg_long \
	$(END_OF_LIST)

lalapps_spec_avg_SOURCES = spec_avg.c
lalapps_spec_avg_long_SOURCES = spec_avg_long.c

if HAVE_PYTHON
pybin_scripts = \
	lalapps_Fscan_checkCoherencePath \
	lalapps_Fscan_coherenceFromSFTs \
	lalapps_Fscan_findCombs \
	lalapps_Fscan_fscanDriver \
	lalapps_Fscan_plotSpecAvgOutput \
	$(END_OF_LIST)
endif

# Add shell test scripts to this variable
test_scripts += testFscan.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
