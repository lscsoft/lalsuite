## Process this file with automake to produce Makefile.in
# $Id$
#if LALFRAME
#LALFRAMEPROGRAMS = lalapps_stochasticPipeline
#else
LALFRAMEPROGRAMS =
#endif
bin_PROGRAMS = lalapps_olapredfcn $(LALFRAMEPROGRAMS)
lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c
#if LALFRAME
#lalapps_stochasticPipeline_SOURCES = stochasticPipeline.c \
                                     stochasticPipeline.h \
                                     stochasticPipelineAux.c
#endif
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
