TEST_CPPFLAGS = -DTEST_DATA_DIR='"$(abs_srcdir)/"'

FRAMEFILES = \
	F-TEST-600000000-60.gwf \
	F-TEST-600000060-60.gwf \
	F-TEST-600000120-60.gwf \
	H-CAL_FAC_V03-729273600-5094000.gwf \
	H-CAL_REF_V03-734073939-64.gwf \
	L-CAL_FAC_V03-729273600-5094000.gwf \
	L-CAL_REF_V03-731488397-64.gwf

#if PTHREAD
#PTHREADPRGS = \
#	TestLowLatencyData1 \
#	TestLowLatencyData2 \
#	TestLowLatencyData3
#
# FIXME: technically codes that link against libpthread need to build
# using the PTHREAD_CC compiler and not CC, this is only a different
# compiler on AIX
#
# set appropriate pthread flags
#TestLowLatencyData1_LDADD = $(PTHREAD_LIBS) $(LDADD)
#TestLowLatencyData1_CFLAGS = $(PTHREAD_CFLAGS) $(AM_CFLAGS)
#TestLowLatencyData2_LDADD = $(PTHREAD_LIBS) $(LDADD)
#TestLowLatencyData2_CFLAGS = $(PTHREAD_CFLAGS) $(AM_CFLAGS)
#TestLowLatencyData3_LDADD = $(PTHREAD_LIBS) $(LDADD)
#TestLowLatencyData3_CFLAGS = $(PTHREAD_CFLAGS) $(AM_CFLAGS)
#endif

check_PROGRAMS = \
	AggregationTest \
	LALFrSeriesTest \
	MakeFrames
TESTS = \
	LALFrSeriesTest

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

noinst_DATA = catalog catalog.test $(FRAMEFILES)

EXTRA_DIST = \
	$(FRAMEFILES) \
	catalog \
	mkcatalog.sh \
	H1-CAL-V03-729273600-734367600.cache \
	L1-CAL-V03-729273600-734367600.cache \
	TestLowLatencyData1.c \
	TestLowLatencyData2.c \
	TestLowLatencyData3.c

MOSTLYCLEANFILES = *.out *.[0-9][0-9][0-9] catalog.out Response*.txt
CLEANFILES = catalog catalog.test H-H1_LSC_AS_Q-600000120-60.gwf

catalog.test:
	@echo "creating" $@
	@$(srcdir)/mkcatalog.sh catalog.test $(srcdir) $(srcdir)
	@echo "X EMPTY 0 0 file:///dev/null" >> catalog.test

catalog:
	@echo "creating" $@
	@$(srcdir)/mkcatalog.sh catalog $(srcdir) $(pkgdatadir)
