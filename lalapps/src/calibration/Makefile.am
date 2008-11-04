## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_mkcalfac lalapps_mkcalref lalapps_dumpcalframe
lalapps_mkcalfac_SOURCES = mkcalfac.c
lalapps_mkcalref_SOURCES = mkcalref.c
lalapps_dumpcalframe_SOURCES = dumpcalframe.c
if LALFRAME
LALFRAMEPROGRAMS = lalapps_makecalframe
lalapps_makecalframe_SOURCES = makecalframe.c LALASCIIFileRead.c LALASCIIFileRead.h
else
LALFRAMEPROGRAMS =
endif
else
LALFRAMEPROGRAMS =
FRAMEPROGRAMS =
endif
if LALFRAME
LALFRAMEONLYPROGRAMS = lalapps_readcalframe lalapps_getresp
lalapps_readcalframe_SOURCES = readcalframe.c
lalapps_getresp_SOURCES = getresp.c
else
LALFRAMEONLYPROGRAMS =
endif
lalapps_ComputeStrain_SOURCES = ComputeStrain.c
lalapps_ComputeStrainDriver_SOURCES = ComputeStrainDriver.c LALASCIIFileRead.c LALASCIIFileRead.h
lalapps_ComputeStrainDriverDMT_SOURCES = ComputeStrainDriverDMT.c
lalapps_ComputeFactors_SOURCES = ComputeFactors.c
lalapps_NoiseComparison_SOURCES = NoiseComparison.c
bin_PROGRAMS = $(FRAMEPROGRAMS) $(LALFRAMEPROGRAMS) $(LALFRAMEONLYPROGRAMS) lalapps_NoiseComparison lalapps_ComputeStrain lalapps_ComputeFactors lalapps_ComputeStrainDriver lalapps_ComputeStrainDriverDMT
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if HAVE_PYTHON
bin_SCRIPTS = lalapps_strain_pipe lalapps_noise_comp_pipe
pkgpython_DATA = mkcal.py strain.py
endif

pkgdata_DATA = \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H2-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt
EXTRA_DIST = mkcalfac.c mkcalref.c mkcalfiles_test.sh mkcal.py \
	     simplereadcalframe.c \
  S1-L1-CAL-FACTORS.txt S1-L1-CAL-RESPONSE.txt S1-L1-CAL-CAV_GAIN.txt \
  S1-H1-CAL-FACTORS.txt S1-H1-CAL-RESPONSE.txt S1-H1-CAL-CAV_GAIN.txt \
  S1-H2-CAL-FACTORS.txt S1-H2-CAL-RESPONSE.txt S1-H2-CAL-CAV_GAIN.txt \
  filters-H1-S2.h filters-H1-S3.h filters-H2-S2.h filters-H2-S3.h \
  filters-L1-S2.h filters-L1-S2.h strain_pipe.in noise_comp_pipe.in \
  mkcal.py strain.py
MOSTLYCLEANFILES = *.ilwd *.gwf
CLEANFILES = python_config.sed lalapps_strain_pipe lalapps_noise_comp_pipe

python_config.sed:
	-rm -f python_config.sed
	echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_strain_pipe: $(srcdir)/strain_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/strain_pipe.in > lalapps_strain_pipe 
	chmod +x lalapps_strain_pipe

lalapps_noise_comp_pipe: $(srcdir)/noise_comp_pipe.in python_config.sed
	sed -f python_config.sed $(srcdir)/noise_comp_pipe.in > lalapps_noise_comp_pipe 
	chmod +x lalapps_noise_comp_pipe

