BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += $(PAGER_CPPFLAGS)

if CFITSIO

bin_PROGRAMS = \
	lalapps_fits_header_getval \
	lalapps_fits_header_list \
	lalapps_fits_overview \
	lalapps_fits_table_list \
	$(END_OF_LIST)

lalapps_fits_header_getval_SOURCES = fits_header_getval.c
lalapps_fits_header_list_SOURCES = fits_header_list.c
lalapps_fits_overview_SOURCES = fits_overview.c
lalapps_fits_table_list_SOURCES = fits_table_list.c

endif

# Add shell test scripts to this variable
test_scripts += testfits.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# All tests require CFITSIO
if !CFITSIO
skip_tests += $(test_scripts)
endif
