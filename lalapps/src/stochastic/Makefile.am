## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_olapredfcn lalapps_stochasticPipeline
lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c
lalapps_stochasticPipeline_SOURCES = stochasticPipeline.c \
                                     stochasticPipeline.h \
                                     stochasticPipelineAux.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
