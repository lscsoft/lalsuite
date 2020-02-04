MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
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
