## Process this file with automake to produce Makefile.in
# $Id$

if LALFRAME
LALFRAMEPROGRAMS = lalapps_stochastic lalapps_stochastic_dev
if LALMETAIO
LALMETAIOFRAMEPROGRAMS = lalapps_stopp
else
LALMETAIOFRAMEPROGRAMS =
endif
else
LALFRAMEPROGRAMS =
endif

bin_PROGRAMS = lalapps_olapredfcn $(LALFRAMEPROGRAMS) \
${LALMETAIOFRAMEPROGRAMS}

lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c

if LALFRAME
lalapps_stochastic_SOURCES = stochastic.c stochastic.h
lalapps_stochastic_dev_SOURCES = stochastic_dev.c stochastic_dev.h
endif

if LALMETAIO
lalapps_stopp_SOURCES = stopp.c
endif

bin_SCRIPTS = lalapps_stochastic_pipe
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = stochastic.py

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_stochastic_pipe: stochastic_pipe.in python_config.sed
	sed -f python_config.sed stochastic_pipe.in > lalapps_stochastic_pipe
	chmod +x lalapps_stochastic_pipe

CLEANFILES = python_config.sed lalapps_stochastic_pipe
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
