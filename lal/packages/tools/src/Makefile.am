## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libtools.la
libtools_la_SOURCES = \
  UnitDefs.c UnitNormalize.c UnitRaise.c UnitMultiply.c UnitCompare.c \
  CreateDetector.c DetResponse.c ComputeTransfer.c CoherentEstimation.c \
  ResampleTimeSeries.c SnglInspiralUtils.c SnglBurstUtils.c \
  ComputeCalibrationFactors.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
