BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	bayestar \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	fits.py \
	plot.py \
	git_version.py \
	$(END_OF_LIST)

dist_pymodule_DATA = \
	ne_simplified_coastline.json \
	$(END_OF_LIST)
