BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = events

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/io

pymodule_PYTHON = \
	__init__.py \
	fits.py \
	hdf5.py \
	$(END_OF_LIST)

endif
