## Process this file with automake to produce Makefile.in
if FRAME
NEED_FRAME = ResampleTimeSeriesTest ComputeCalibrationFactorsTest
else
NEED_FRAME = 
endif
check_PROGRAMS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse ComputeTransferTest CoherentEstimationTest $(NEED_FRAME)
UnitsTest_SOURCES = UnitsTest.c
UnitsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
DetectorSiteTest_SOURCES = DetectorSiteTest.c
DetectorSiteTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALTestDetResponse0_SOURCES = LALTestDetResponse0.c
LALTestDetResponse0_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
LALIndependentTestDetResponse_SOURCES = LALIndependentTestDetResponse.c 
LALIndependentTestDetResponse_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ComputeTransferTest_SOURCES = ComputeTransferTest.c
ComputeTransferTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
CoherentEstimationTest_SOURCES = CoherentEstimationTest.c
CoherentEstimationTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
if FRAME
ResampleTimeSeriesTest_SOURCES = ResampleTimeSeriesTest.c
ResampleTimeSeriesTest_LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
ComputeCalibrationFactorsTest_SOURCES = ComputeCalibrationFactorsTest.c
ComputeCalibrationFactorsTest_LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
endif
EXTRA_DIST = ResampleTimeSeriesTest.c ComputeCalibrationFactorsTest.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
if FRAME
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse ComputeTransferTest ComputeCalibrationFactorsTest
else
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse ComputeTransferTest
endif
MOSTLYCLEANFILES = .dvi-dep PrintVector.* transfer.out
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
