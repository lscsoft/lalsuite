BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

SUBDIRS = lalstochastic

# workaround autotools not cleaning up __pycache__
clean-local:
	-find . -name "__pycache__" -type d -exec rm -r "{}" \;
