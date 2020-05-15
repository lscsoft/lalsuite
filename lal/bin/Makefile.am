BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# -- C programs -------------

bin_PROGRAMS = \
	lal_simd_detect \
	lal_version \
	$(END_OF_LIST)

lal_simd_detect_SOURCES = simd_detect.c
lal_version_SOURCES = version.c

TESTS += \
	lal_simd_detect \
	lal_version \
	$(END_OF_LIST)

# -- Python programs --------

if HAVE_PYTHON

pybin_scripts = \
	$(END_OF_LIST)

TESTS += \
	$(END_OF_LIST)

endif
