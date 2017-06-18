BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	bayestar \
	imrtgr \
	io \
	plot \
	popprior \
	rapid_pe \
	tiger \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	healpix_tree.py \
	lalinference_pipe_utils.py \
	nest2pos.py \
	git_version.py \
	$(END_OF_LIST)

endif
