BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
man1_MANS =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_crosscorr_v2 \
	$(END_OF_LINE)

noinst_PROGRAMS = \
	lalpulsar_crosscorr \
	$(END_OF_LINE)

lalpulsar_crosscorr_SOURCES = \
	pulsar_crosscorr.c \
	pulsar_crosscorr.h \
	$(END_OF_LINE)

lalpulsar_crosscorr_v2_SOURCES = \
	CrossCorrToplist.c \
	CrossCorrToplist.h \
	pulsar_crosscorr_v2.c \
	$(END_OF_LINE)

if HAVE_PYTHON

pybin_scripts = $(END_OF_LINE)

endif

# Add shell test scripts to this variable
test_scripts += test_crosscorr_v2.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
