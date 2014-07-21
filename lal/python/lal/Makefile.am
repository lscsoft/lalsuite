BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	spectrum \
	utils \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	git_version.py \
	gpstime.py \
	series.py \
	$(END_OF_LIST)
