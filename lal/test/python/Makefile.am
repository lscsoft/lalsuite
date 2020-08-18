include $(top_srcdir)/gnuscripts/lalsuite_test_python.am

EXTRA_DIST += \
	874000000-20000.cache \
	$(END_OF_LIST)

CLEANFILES = \
	iterutils_randindex_histogram.png \
	$(END_OF_LIST)

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts += \
	test_antenna.py \
	test_cache.py \
	test_gpstime.py \
	test_iterutils.py \
	test_rate.py \
	utils_cache_verify.py \
	$(END_OF_LIST)

# Add any helper programs required by tests to this variable
test_helpers +=
