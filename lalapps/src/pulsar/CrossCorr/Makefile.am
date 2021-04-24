BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic

bin_PROGRAMS = \
	lalapps_pulsar_crosscorr \
	lalapps_pulsar_crosscorr_v2 \
	$(END_OF_LINE)

lalapps_pulsar_crosscorr_SOURCES = \
	pulsar_crosscorr.c \
	pulsar_crosscorr.h \
	$(END_OF_LINE)

lalapps_pulsar_crosscorr_v2_SOURCES = \
	pulsar_crosscorr_v2.c \
	CrossCorrToplist.c \
	CrossCorrToplist.h \
	HeapToplist.c \
	$(END_OF_LINE)

if HAVE_PYTHON
pybin_scripts = \
	lalapps_combine_crosscorr_toplists \
	lalapps_run_pulsar_crosscorr_v2 \
	$(END_OF_LINE)
endif

# Add shell test scripts to this variable
test_scripts += test_pulsar_crosscorr_v2.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
