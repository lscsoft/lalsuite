BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_spec_avg \
	lalpulsar_spec_avg_long \
	lalpulsar_spec_coherence \
	$(END_OF_LIST)

lalpulsar_spec_avg_SOURCES = \
	spec_avg.c \
	fscanutils.c \
	fscanutils.h \
	$(END_OF_LIST)

lalpulsar_spec_avg_long_SOURCES = \
	spec_avg_long.c \
	fscanutils.c \
	fscanutils.h \
	$(END_OF_LIST)

lalpulsar_spec_coherence_SOURCES = \
	spec_coherence.c \
	fscanutils.c \
	fscanutils.h \
	$(END_OF_LIST)

# Add shell test scripts to this variable
test_scripts += testFscan.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
