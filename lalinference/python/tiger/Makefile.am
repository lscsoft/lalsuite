# The Numpy 1.7 headers contain some static function definitions that
# themselves are not '-Werror'-safe. Turn off the '-Werror' (treat warnings
# as errors) flag, for just targets declared in this file.
CFLAGS += -Wno-error

tigerdir = $(pkgpythondir)/tiger

tiger_PYTHON = \
	__init__.py \
	make_injtimes.py \
	postproc.py \
	omegascans_dag.py
