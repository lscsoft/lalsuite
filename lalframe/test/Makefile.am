EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += LALFrSeriesTest
#test_programs += AggregationTest

# Add shell, Python, etc. test scripts to this variable
test_scripts +=

# Add any helper programs required by tests to this variable
test_helpers += MakeFrames

FRAMEFILES = \
	F-TEST-600000000-60.gwf \
	F-TEST-600000060-60.gwf \
	F-TEST-600000120-60.gwf \
	H-CAL_FAC_V03-729273600-5094000.gwf \
	H-CAL_REF_V03-734073939-64.gwf \
	L-CAL_FAC_V03-729273600-5094000.gwf \
	L-CAL_REF_V03-731488397-64.gwf \
	$(END_OF_LIST)

BUILT_SOURCES = \
	catalog \
	catalog.test \
	$(END_OF_LIST)

catalog:
	$(AM_V_GEN)$(srcdir)/mkcatalog.sh $@ $(srcdir) $(pkgdatadir)

catalog.test:
	$(AM_V_GEN)$(srcdir)/mkcatalog.sh $@ $(srcdir) $(srcdir); \
	echo "X EMPTY 0 0 file:///dev/null" >> $@

MOSTLYCLEANFILES = \
	*.[0-9][0-9][0-9] \
	*.out \
	H-H1_LSC_AS_Q-600000120-60.gwf \
	Response*.txt \
	catalog \
	catalog.out \
	catalog.test \
	$(END_OF_LIST)

EXTRA_DIST += \
	$(FRAMEFILES) \
	H1-CAL-V03-729273600-734367600.cache \
	L1-CAL-V03-729273600-734367600.cache \
	TestLowLatencyData1.c \
	TestLowLatencyData2.c \
	TestLowLatencyData3.c \
	mkcatalog.sh \
	$(END_OF_LIST)

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
