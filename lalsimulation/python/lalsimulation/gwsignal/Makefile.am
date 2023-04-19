BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = \
	core \
	models \
	$(END_OF_LIST)

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/gwsignal

pymodule_PYTHON = \
	__init__.py \
	$(END_OF_LIST)

endif
