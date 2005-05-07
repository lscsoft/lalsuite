## Process this file with automake to produce Makefile.in
# $Id$

if FRAME
if LALFRAME
FRAMEPROGRAMS = lalapps_stochastic
lalapps_stochastic_SOURCES = stochastic.c
endif
endif
if LALMETAIO
METAIOPROGRAMS = lalapps_stopp
lalapps_stopp_SOURCES = stopp.c
endif
bin_SCRIPTS = lalapps_stochastic_pipe lalapps_stochastic_bayes
bin_PROGRAMS = lalapps_olapredfcn ${FRAMEPROGRAMS} ${METAIOPROGRAMS}
lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = stochastic.py
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
EXTRA_DIST = stochastic.py stochastic_pipe.in stochastic_bayes.in
MOSTLYCLEANFILES = *STOCHASTIC*.xml *CCSPECTRA*.gwf
CLEANFILES = python_config.sed lalapps_stochastic_pipe \
             lalapps_stochastic_bayes

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed

lalapps_stochastic_pipe: stochastic_pipe.in python_config.sed
	sed -f python_config.sed stochastic_pipe.in > lalapps_stochastic_pipe
	chmod +x lalapps_stochastic_pipe

lalapps_stochastic_bayes: stochastic_bayes.in python_config.sed
	sed -f python_config.sed stochastic_bayes.in > lalapps_stochastic_bayes
	chmod +x lalapps_stochastic_bayes
