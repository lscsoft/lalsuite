BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	utils \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	antenna.py \
	git_version.py \
	gpstime.py \
	iterutils.py \
	pipeline.py \
	rate.py \
	series.py \
	$(END_OF_LIST)

endif
