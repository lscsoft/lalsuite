## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libpulsar.la
libpulsar_la_SOURCES = \
  GetEarthTimes.c TBaryPtolemaic.c TSpin.c TComp.c CoherentMetric.c      \
  StackMetric.c ProjectMetric.c FlatMesh.c CreateResampleRules.c         \
  DestroyResampleRules.c ApplyResampleRules.c PolycoToTimingDifference.c \
  RulesToTimingDifference.c LALDemod.c ComputeSky.c FoldAmplitudes.c     \
  TwoDMesh.c TwoDMeshInternal.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
