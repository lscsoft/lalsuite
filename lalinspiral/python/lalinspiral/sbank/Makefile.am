BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/sbank

pymodule_PYTHON = \
	__init__.py \
	bank.py \
	overlap_cuda.py \
	psds.py \
	tau0tau3.py \
	waveforms.py \
	$(END_OF_LIST)

endif
