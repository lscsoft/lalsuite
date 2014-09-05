BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

pkgdata_DATA = \
	laldetchar-idq.ini

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/idq

pymodule_PYTHON = \
	__init__.py \
	auxmvc.py \
	auxmvc_utils.py \
	event.py \
	idq_gdb_utils.py \
	idq.py \
	idq_summary_plots.py \
	idq_tables.py \
	ovl.py \
	pdf_estimation.py \
	svmkit.py \
	$(END_OF_LIST)

endif
