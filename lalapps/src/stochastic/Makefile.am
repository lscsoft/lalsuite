## Process this file with automake to produce Makefile.in
# $Id$

if FRAME
if LALFRAME
STOCHPROGRAM = lalapps_stochastic
endif
endif
if LALMETAIO
STOPPPROGRAM = lalapps_stopp
endif
bin_SCRIPTS = lalapps_stochastic_pipe lalapps_stochastic_bayes
bin_PROGRAMS = lalapps_olapredfcn $(STOCHPROGRAM) ${STOPPPROGRAM}
if FRAME
if LALFRAME
lalapps_stochastic_SOURCES = stochastic.c helper.c stochastic.h
endif
endif
if LALMETAIO
lalapps_stopp_SOURCES = stopp.c
endif
lalapps_olapredfcn_SOURCES = olapredfcn.c olapredfcn.h olapredfcn_args.c \
                             olapredfcn_bars.c
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = stochastic.py
pkgdata_DATA = stochastic_S3_H1L1.ini stochastic_S3_H1L1_bayes.ini
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
MOSTLYCLEANFILES = *stochastic*.xml *ccSpectra*.gwf
CLEANFILES = python_config.sed lalapps_stochastic_pipe \
						 lalapps_stochastic_bayes stochastic_S3_H1L1.ini \
						 stochastic_S3_H1L1_bayes.ini

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@EXECPATH@+${bindir}+g" >> python_config.sed

lalapps_stochastic_pipe: stochastic_pipe.in python_config.sed
	sed -f python_config.sed stochastic_pipe.in > lalapps_stochastic_pipe
	chmod +x lalapps_stochastic_pipe

lalapps_stochastic_bayes: stochastic_bayes.in python_config.sed
	sed -f python_config.sed stochastic_bayes.in > lalapps_stochastic_bayes
	chmod +x lalapps_stochastic_bayes

stochastic_S3_H1L1.ini: stochastic_S3_H1L1.in python_config.sed
	sed "s+@EXECPATH@+${bindir}+" stochastic_S3_H1L1.in > stochastic_S3_H1L1.ini
	
stochastic_S3_H1L1_bayes.ini: stochastic_S3_H1L1_bayes.in python_config.sed
	sed "s+@EXECPATH@+${bindir}+" stochastic_S3_H1L1_bayes.in > stochastic_S3_H1L1_bayes.ini
