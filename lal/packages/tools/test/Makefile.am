AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -DTEST_DATA_DIR='"$(abs_srcdir)/"'
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	CoherentEstimationTest \
	ComputeTransferTest \
	DetectorSiteTest \
	FrequencySeriesTest \
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

EXTRA_DIST = \
	ComputeCalibrationFactorsTest.c \
	LALIndTestDetRes.sh \
	ResampleTimeSeriesTest.c \
	indTestDetResGEO.cfg \
	indTestDetResLHO.cfg \
	indTestDetResLLO.cfg

TESTS = \
	ComputeTransferTest \
	DetectorSiteTest \
	FrequencySeriesTest \
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
