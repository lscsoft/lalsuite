## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libutilities.la
libutilities_la_SOURCES = \
  Random.c MersenneRandom.c FindRoot.c Integrate.c Interpolate.c HeapSort.c ODE.c \
  CoarseGrainFrequencySeries.c Dirichlet.c LALRunningMedian.c
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
