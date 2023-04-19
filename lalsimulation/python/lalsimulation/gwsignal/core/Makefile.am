BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/gwsignal/core

pymodule_PYTHON = \
	__init__.py \
	waveform.py \
	waveform_conditioning.py \
	utils.py \
	parameter_conventions.py \
	gw.py \
	errors.py \
	conditioning_subroutines.py \
	$(END_OF_LIST)

endif
