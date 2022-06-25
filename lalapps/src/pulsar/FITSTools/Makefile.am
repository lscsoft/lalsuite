BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps_pulsar.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

AM_CPPFLAGS += $(PAGER_CPPFLAGS)

if CFITSIO

bin_PROGRAMS = \
	lalpulsar_fits_header_getval \
	lalpulsar_fits_header_list \
	lalpulsar_fits_overview \
	lalpulsar_fits_table_list \
	$(END_OF_LIST)

lalpulsar_fits_header_getval_SOURCES = \
	fits_header_getval.c \
	$(END_OF_LIST)

lalpulsar_fits_header_list_SOURCES = \
	fits_header_list.c \
	$(END_OF_LIST)

lalpulsar_fits_overview_SOURCES = \
	fits_overview.c \
	$(END_OF_LIST)

lalpulsar_fits_table_list_SOURCES = \
	fits_table_list.c \
	$(END_OF_LIST)

endif

# Add shell test scripts to this variable
test_scripts += testfits.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=

# All tests require CFITSIO
if !CFITSIO

skip_tests += $(test_scripts)

endif
