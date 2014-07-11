BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/lalapps.am

lalapps_ComputeStrainDriver_SOURCES = ComputeStrainDriver.c
lalapps_ComputeFactors_SOURCES = ComputeFactors.c
lalapps_NoiseComparison_SOURCES = NoiseComparison.c

bin_PROGRAMS = \
	lalapps_ComputeFactors \
	lalapps_ComputeStrainDriver \
	lalapps_NoiseComparison

if HAVE_PYTHON
bin_SCRIPTS = lalapps_strain_pipe lalapps_noise_comp_pipe
pkgpython_PYTHON = strain.py
endif

EXTRA_DIST = \
	strain_pipe.in \
	noise_comp_pipe.in \
	strain.py

MOSTLYCLEANFILES = \
	*.ilwd \
	*.gwf

CLEANFILES = \
	python_config.sed \
	lalapps_strain_pipe \
	lalapps_noise_comp_pipe

python_config.sed:
	@-rm -f python_config.sed
	@echo "s+@PYTHONLIBDIR@+@pkgpythondir@+g" > python_config.sed
	@echo "s+@PYTHONPROG@+@PYTHON@+g" >> python_config.sed
	@echo "s+@RINGEXECPATH@+${bindir}+g" >> python_config.sed

lalapps_strain_pipe: $(srcdir)/strain_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/strain_pipe.in > lalapps_strain_pipe
	@chmod +x lalapps_strain_pipe

lalapps_noise_comp_pipe: $(srcdir)/noise_comp_pipe.in python_config.sed
	$(AM_V_GEN)sed -f python_config.sed $(srcdir)/noise_comp_pipe.in > lalapps_noise_comp_pipe
	@chmod +x lalapps_noise_comp_pipe
