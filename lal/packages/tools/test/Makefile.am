AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include -DTEST_DATA_DIR='"$(abs_srcdir)/"' @AM_CPPFLAGS@
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
	circ_series.txt \
	cross.txt \
	cross_at_0_0.txt \
	cross_at_2_10.txt \
	cross_at_4_15.txt \
	cross_at_m4_15.txt \
	cross_series.txt \
	cross_sq_avg.txt \
	cross_sq_time_avg.txt \
	ff_diff_cros.txt \
	ff_diff_plus.txt \
	ff_lal_cros.txt \
	ff_lal_plus.txt \
	ff_local_cros.txt \
	ff_local_plus.txt \
	ff_rms_diff_cros_vs_fudge.txt \
	ff_rms_diff_plus_vs_fudge.txt \
	phi.txt \
	plus.txt \
	plus_at_0_0.txt \
	plus_at_2_10.txt \
	plus_at_4_15.txt \
	plus_at_m4_15.txt \
	plus_series.txt \
	plus_sq_avg.txt \
	plus_sq_time_avg.txt \
	scalar_series.txt \
	sum_of_sq_time_avg.txt \
	sum_series.txt \
	sum_sq.txt \
	sum_sq_avg.txt \
	theta.txt \
	transfer.out
