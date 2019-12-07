EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts += \
	test_python_imports.py \
	test_python_scripts.py \
	test_frread.py \
	$(END_OF_LIST)

# Add any helper programs required by tests to this variable
test_helpers +=

# Clean up byte-compiled files from tests
DISTCLEANFILES = \
    __pycache__/ \
    .pytest_cache/ \
    $(END_OF_LIST)
