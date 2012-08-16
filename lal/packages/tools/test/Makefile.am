INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include

check_PROGRAMS = \
	CoherentEstimationTest \
	ComputeTransferTest \
	DetectorSiteTest \
	FrequencySeriesTest \
	LALBarycenterTest \
	LALIndependentTestDetResponse \
	LALTestDetResponse0 \
	SegmentsTest \
	SequenceTest \
	SkymapTest \
	CubicSplineTriggerInterpolantTest \
	LanczosTriggerInterpolantTest \
	NearestNeighborTriggerInterpolantTest \
	QuadraticFitTriggerInterpolantTest \
	TimeSeriesTest \
	UnitsTest

CoherentEstimationTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

ComputeTransferTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

DetectorSiteTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

FrequencySeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALBarycenterTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALIndependentTestDetResponse_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LALTestDetResponse0_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SegmentsTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SequenceTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

SkymapTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

TimeSeriesTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

UnitsTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

CubicSplineTriggerInterpolantTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

LanczosTriggerInterpolantTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

NearestNeighborTriggerInterpolantTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

QuadraticFitTriggerInterpolantTest_LDADD = \
	$(top_builddir)/lib/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

EXTRA_DIST = \
	ComputeCalibrationFactorsTest.c \
	LALIndTestDetRes.sh \
	ResampleTimeSeriesTest.c \
	earth98.dat \
	indTestDetResGEO.cfg \
	indTestDetResLHO.cfg \
	indTestDetResLLO.cfg \
	sun98.dat

TESTS = \
	ComputeTransferTest \
	DetectorSiteTest \
	FrequencySeriesTest \
	LALBarycenterTest \
	LALIndTestDetRes.sh \
	LALTestDetResponse0 \
	SequenceTest \
	SkymapTest \
	TimeSeriesTest \
	CubicSplineTriggerInterpolantTest \
	LanczosTriggerInterpolantTest \
	NearestNeighborTriggerInterpolantTest \
	QuadraticFitTriggerInterpolantTest \
	UnitsTest

TESTS_ENVIRONMENT = 1>/dev/null 2>/dev/null LAL_DATA_PATH=$(srcdir)

MOSTLYCLEANFILES = \
	PrintVector.* \
	ff_diff_cros.txt \
	ff_diff_plus.txt \
	ff_lal_cros.txt \
	ff_lal_plus.txt \
	ff_local_cros.txt \
	ff_local_plus.txt \
	ff_rms_diff_cros_vs_fudge.txt \
	ff_rms_diff_plus_vs_fudge.txt \
	transfer.out
