AM_CFLAGS = $(LALSUPPORT_CFLAGS)
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

TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(srcdir)

noinst_DATA = \
	a.pgm \
	CLRindata.asc

EXTRA_DIST = \
	$(DATA)

MOSTLYCLEANFILES = CLRoutdata.asc
