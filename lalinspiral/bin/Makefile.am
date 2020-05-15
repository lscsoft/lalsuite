BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- C programs -------------

bin_PROGRAMS = \
	lalinspiral_version \
	$(END_OF_LIST)

lalinspiral_version_SOURCES = version.c

TESTS += \
	lalinspiral_version \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif
