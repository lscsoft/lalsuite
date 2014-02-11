AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include @AM_CPPFLAGS@
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	GeocentricGeodeticTest \
	SkyCoordinatesTest

TESTS = \
	$(check_PROGRAMS)

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h
