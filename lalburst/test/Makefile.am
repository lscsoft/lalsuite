AMCFLAGS = $(LALSUPPORT_CFLAGS)

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

CLRTest_SOURCES = CLRTest.c
CLRTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
TfrPswvTest_SOURCES = TfrPswvTest.c
TfrPswvTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
TfrRspTest_SOURCES = TfrRspTest.c
TfrRspTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
TfrSpTest_SOURCES = TfrSpTest.c
TfrSpTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
TfrWvTest_SOURCES = TfrWvTest.c
TfrWvTest_LDADD = \
$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
ThresholdsTest_SOURCES = ThresholdsTest.c
ThresholdsTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la
TrackSearchTest_SOURCES = TrackSearchTest.c
TrackSearchTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalburst.la

MOSTLYCLEANFILES = CLRoutdata.asc
