AM_CFLAGS = $(LALSUPPORT_CFLAGS) -DDATADIR='"$(abs_srcdir)/"'
LDADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalburst.la

if LALSUPPORT
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
else
check_PROGRAMS =
TESTS =
endif

noinst_DATA = \
	a.pgm \
	CLRindata.asc

EXTRA_DIST = \
	$(DATA)

MOSTLYCLEANFILES = CLRoutdata.asc
