AM_CFLAGS = @AM_CFLAGS@ $(LALSUPPORT_CFLAGS) -DTEST_DATA_DIR='"$(abs_srcdir)/"'
AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalburst.la

check_PROGRAMS = \
	CLRTest \
	TfrPswvTest \
	TfrRspTest \
	TfrSpTest \
	TfrWvTest \
	ThresholdsTest \
	TrackSearchTest
TESTS = \
	$(check_PROGRAMS)

noinst_DATA = \
	a.pgm \
	CLRindata.asc

EXTRA_DIST = \
	$(DATA)

MOSTLYCLEANFILES = CLRoutdata.asc
