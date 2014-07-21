BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	hveto \
	triggers \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	dqsegs.py \
	git_version.py \
	$(END_OF_LIST)
