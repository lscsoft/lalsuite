## Process this file with automake to produce Makefile.in
FRAMEFILES = F-TEST-600000000-60.gwf F-TEST-600000060-60.gwf F-TEST-600000120-60.gwf \
H-CAL_FAC_V03-729273600-5094000.gwf H-CAL_REF_V03-734073939-64.gwf \
L-CAL_FAC_V03-729273600-5094000.gwf L-CAL_REF_V03-731488397-64.gwf

if FRAME
#check_PROGRAMS   = FrameStreamTest Read40mData ReadGEOData ReadLHOData MakeFrames FrameCacheTest FrameCalibrationTest
check_PROGRAMS   = FrameStreamTest Read40mData ReadGEOData ReadLHOData MakeFrames FrameCacheTest AggregationTest
#TESTS            = FrameCacheTest.sh FrameStreamTest FrameCalibrationTest
TESTS            = FrameCacheTest.sh FrameStreamTest
else
check_PROGRAMS   =
TESTS            = NoTest
endif

TESTS_ENVIRONMENT = LAL_FRAME_PATH=$(srcdir)

pkgdata_DATA = catalog $(FRAMEFILES)
noinst_DATA = catalog.test

EXTRA_DIST = NoTest FrameDataTest.c SpecBufferTest.c DataBufferTest.c \
        mkcatalog.sh FrameCacheTest.sh $(FRAMEFILES) \
        H1-CAL-V03-729273600-734367600.cache L1-CAL-V03-729273600-734367600.cache \
	FrameCalibrationTest.c

EXTRA_TEXSOURCES = FrameCalibrationTest.c

LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la

FrameStreamTest_SOURCES = FrameStreamTest.c
Read40mData_SOURCES = Read40mData.c
ReadGEOData_SOURCES = ReadGEOData.c
ReadLHOData_SOURCES = ReadLHOData.c
MakeFrames_SOURCES = MakeFrames.c
#FrameCalibrationTest_SOURCES = FrameCalibrationTest.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep *.out *.[0-9][0-9][0-9] catalog.out Response*.txt
CLEANFILES = catalog catalog.test H-H1_LSC_AS_Q-600000120-60.gwf
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep

catalog.test:
	$(srcdir)/mkcatalog.sh catalog.test $(srcdir) $(srcdir)
	echo "X EMPTY 0 0 file:///dev/null" >> catalog.test

catalog:
	$(srcdir)/mkcatalog.sh catalog $(srcdir) $(pkgdatadir)
@DVIMKFILE@
