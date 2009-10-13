AMCFLAGS = $(LALSUPPORT_CFLAGS)

FRAMEFILES = \
	F-TEST-600000000-60.gwf \
	F-TEST-600000060-60.gwf \
	F-TEST-600000120-60.gwf \
	H-CAL_FAC_V03-729273600-5094000.gwf \
	H-CAL_REF_V03-734073939-64.gwf \
	L-CAL_FAC_V03-729273600-5094000.gwf \
	L-CAL_REF_V03-731488397-64.gwf

if LALSUPPORT
check_PROGRAMS = \
	AggregationTest \
	FrameCacheTest \
	FrameDataTest \
	FrameStreamTest \
	MakeFrames \
	Read40mData \
	ReadGEOData \
	ReadLHOData
TESTS = \
	FrameCacheTest.sh \
	FrameStreamTest
else
check_PROGRAMS =
TESTS =
endif

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

pkgdata_DATA = catalog $(FRAMEFILES)
noinst_DATA = catalog.test

EXTRA_DIST = \
	$(FRAMEFILES) \
	mkcatalog.sh \
	H1-CAL-V03-729273600-734367600.cache \
	L1-CAL-V03-729273600-734367600.cache \
	FrameCacheTest.sh

AggregationTest_SOURCES = AggregationTest.c
AggregationTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
FrameCacheTest_SOURCES = FrameCacheTest.c
FrameCacheTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
FrameDataTest_SOURCES = FrameDataTest.c
FrameDataTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
FrameStreamTest_SOURCES = FrameStreamTest.c
FrameStreamTest_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
MakeFrames_SOURCES = MakeFrames.c
MakeFrames_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
Read40mData_SOURCES = Read40mData.c
Read40mData_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
ReadGEOData_SOURCES = ReadGEOData.c
ReadGEOData_LDADD = \
$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la
ReadLHOData_SOURCES = ReadLHOData.c
ReadLHOData_LDADD = \
	$(LALSUPPORT_LIBS) \
	$(top_builddir)/src/liblalframe.la

MOSTLYCLEANFILES = *.out *.[0-9][0-9][0-9] catalog.out Response*.txt
CLEANFILES = catalog catalog.test H-H1_LSC_AS_Q-600000120-60.gwf

catalog.test:
	$(srcdir)/mkcatalog.sh catalog.test $(srcdir) $(srcdir)
	echo "X EMPTY 0 0 file:///dev/null" >> catalog.test

catalog:
	$(srcdir)/mkcatalog.sh catalog $(srcdir) $(pkgdatadir)
