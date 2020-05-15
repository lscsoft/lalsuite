BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- C programs -------------

bin_PROGRAMS = \
	lalmetaio_version \
	$(END_OF_LIST)

lalmetaio_version_SOURCES = version.c

TESTS += \
	lalmetaio_version \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif
