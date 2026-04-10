include $(top_srcdir)/gnuscripts/lalsuite_test_python.am

# Run tests against each configured frame library
test_script_compiler = $(top_srcdir)/gnuscripts/lalframe_test_compiler $(FRAME_LIBRARIES_AVAILABLE) --

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
test_scripts += \
	test_frread.py \
	test_utils_frtools.py \
	$(END_OF_LIST)

# Add any helper programs required by tests to this variable
test_helpers +=
