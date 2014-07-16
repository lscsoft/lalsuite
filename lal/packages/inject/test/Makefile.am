include $(top_srcdir)/gnuscripts/lalsuite_header_links.am

TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

check_PROGRAMS = \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

TESTS = \
	$(check_PROGRAMS)

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h
