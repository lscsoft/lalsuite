## Process this file with automake to produce Makefile.in
# $Id$
if LALFRAME
LALFRAMEPROGRAMS = lalapps_stochastic
else
LALFRAMEPROGRAMS =
endif
bin_PROGRAMS = lalapps_olapredfcn $(LALFRAMEPROGRAMS)
lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c
if LALFRAME
lalapps_stochastic_SOURCES = stochastic.c stochastic.h
endif
bin_SCRIPTS = lalapps_stochastic_pipe
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = stochastic.py

lalapps_stochastic_pipe: stochastic_pipe.in
	sed s+@PYTHONLIBDIR@+${pythonlibdir}+g stochastic_pipe.in > lalapps_stochastic_pipe

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
