BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_test_python.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

EXTRA_DIST += \
	test_thinca_correct_min1.xml.gz \
	test_thinca_correct_min2.xml.gz \
	test_thinca_input.xml.gz \
	$(END_OF_LIST)

noinst_pybin_scripts = \
	thinca_brute_force_coinc \
	thinca_compare_coinc_lists \
	$(END_OF_LIST)

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts += \
	test_thinca.sh \
	$(END_OF_LIST)

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES += \
	brute_min1.xml.gz \
	brute_min2.xml.gz \
	thinca_min1.xml.gz \
	thinca_min2.xml.gz \
	$(END_OF_LIST)

endif
