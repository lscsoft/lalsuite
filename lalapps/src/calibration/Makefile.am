## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_mkcalfac lalapps_mkcalref
else
FRAMEPROGRAMS = 
endif
bin_PROGRAMS = $(FRAMEPROGRAMS) lalapps_ComputeStrain lalapps_ComputeFactors lalapps_ComputeStrainDriver
if FRAME
lalapps_mkcalfac_SOURCES = mkcalfac.c
lalapps_mkcalref_SOURCES = mkcalref.c
endif
bin_SCRIPTS = lalapps_strain_pipe
lalapps_ComputeStrain_SOURCES = ComputeStrain.c
lalapps_ComputeStrainDriver_SOURCES = ComputeStrainDriver.c
lalapps_ComputeFactors_SOURCES = ComputeFactors.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
pythonlibdir = ${prefix}/lib/python
pythonlib_SCRIPTS = mkcal.py strain.py
pkgdata_DATA = \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H2-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt
EXTRA_DIST = mkcalfac.c mkcalref.c mkcalfiles_test.sh mkcal.py \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H2-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt \
  filters-H1-S2.h filters-H1-S3.h filters-H2-S2.h filters-H2-S3.h \
  filters-L1-S2.h filters-L1-S2.h strain_pipe.in
MOSTLYCLEANFILES = *.ilwd *.gwf
TESTS = mkcalfiles_test.sh

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+${pythonlibdir}+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_strain_pipe: strain_pipe.in python_config.sed
	sed -f python_config.sed strain_pipe.in > lalapps_strain_pipe 
	chmod +x lalapps_strain_pipe

