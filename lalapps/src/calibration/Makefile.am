BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalapps.am

lalapps_ComputeStrainDriver_SOURCES = ComputeStrainDriver.c
lalapps_ComputeFactors_SOURCES = ComputeFactors.c
lalapps_NoiseComparison_SOURCES = NoiseComparison.c

bin_PROGRAMS = \
	lalapps_ComputeFactors \
	lalapps_ComputeStrainDriver \
	lalapps_NoiseComparison

if HAVE_PYTHON
pybin_SCRIPTS = lalapps_strain_pipe lalapps_noise_comp_pipe
pkgpython_PYTHON = strain.py
endif

MOSTLYCLEANFILES += \
	*.ilwd \
	*.gwf
