## Process this file with automake to produce Makefile.in
if FRAME
#NEED_FRAME = ResampleTimeSeriesTest ComputeCalibrationFactorsTest
NEED_FRAME = ResampleTimeSeriesTest
#EXTRA_TEXSOURCES =
EXTRA_TEXSOURCES = ComputeCalibrationFactorsTest.c
else
NEED_FRAME = 
EXTRA_TEXSOURCES = ResampleTimeSeriesTest.c ComputeCalibrationFactorsTest.c
endif
check_PROGRAMS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse ComputeTransferTest CoherentEstimationTest SegmentsTest $(NEED_FRAME) SequenceTest TimeSeriesTest
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
SegmentsTest_SOURCES = SegmentsTest.c
SegmentsTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
SequenceTest_SOURCES = SequenceTest.c
SequenceTest_LDADD = $(top_builddir)/lib/liblal.la
TimeSeriesTest_SOURCES = TimeSeriesTest.c
TimeSeriesTest_LDADD = $(top_builddir)/lib/liblal.la
if FRAME
ResampleTimeSeriesTest_SOURCES = ResampleTimeSeriesTest.c
ResampleTimeSeriesTest_LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
#ComputeCalibrationFactorsTest_SOURCES = ComputeCalibrationFactorsTest.c
#ComputeCalibrationFactorsTest_LDADD = $(top_builddir)/packages/framedata/src/liblalframe.la $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
endif
EXTRA_DIST = ResampleTimeSeriesTest.c ComputeCalibrationFactorsTest.c LALIndTestDetRes.sh indTestDetResLHO.cfg indTestDetResLLO.cfg indTestDetResGEO.cfg
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
if FRAME
#TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse ComputeTransferTest ComputeCalibrationFactorsTest SequenceTest TimeSeriesTest
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse LALIndTestDetRes.sh ComputeTransferTest SequenceTest TimeSeriesTest
else
TESTS = UnitsTest DetectorSiteTest LALTestDetResponse0 LALIndependentTestDetResponse LALIndTestDetRes.sh ComputeTransferTest SequenceTest TimeSeriesTest
endif
TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null
MOSTLYCLEANFILES = .dvi-dep PrintVector.* transfer.out \
	ff_diff_plus.txt ff_rms_diff_plus_vs_fudge.txt \
	ff_rms_diff_cros_vs_fudge.txt ff_local_cros.txt ff_lal_cros.txt \
	ff_local_plus.txt ff_lal_plus.txt ff_diff_cros.txt

DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
