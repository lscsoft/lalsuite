BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# The Numpy 1.7 headers contain some static function definitions that
# themselves are not '-Werror'-safe. Turn off the '-Werror' (treat warnings
# as errors) flag, for just targets declared in this file.
CFLAGS += -Wno-error

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/tiger

pymodule_PYTHON = \
	__init__.py \
	make_injtimes.py \
	postproc.py \
	omegascans_dag.py

endif
