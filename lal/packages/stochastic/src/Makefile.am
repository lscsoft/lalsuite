## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libstochastic.la
libstochastic_la_SOURCES = StochasticCrossCorrelation.c \
  ZeroPadAndFFT.c StochasticOptimalFilter.c OverlapReductionFunction.c \
  StochasticOmegaGW.c StochasticInverseNoise.c Dirichlet.c

#  FalseAlarm.c FalseDismissal.c FDismissalVsFAlarm.c

INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
