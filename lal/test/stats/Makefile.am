EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += LALMomentTest
test_programs += XLALChisqTest

# Add shell, Python, etc. test scripts to this variable
test_scripts +=
if SWIG_BUILD_PYTHON
test_scripts += LALMarcumQTest.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

EXTRA_DIST += \
	CheckStatus.h \
	$(END_OF_LIST)
