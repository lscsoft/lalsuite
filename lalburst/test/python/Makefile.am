XFAIL_TESTS =
include $(top_srcdir)/gnuscripts/lalsuite_test_python.am

# Add shell, Python, etc. test scripts to this variables
test_scripts += \
	date_verify.py \
	snglcoinc_verify.py \
	$(END_OF_LIST)

# with python-ligo-lw installed, and scipy<0.19.0
# snglcoinc_verify.py will fail, see lalburst's configure.ac and
# https://git.ligo.org/lscsoft/lalsuite/issues/208
# for more details
if !SNGLCOINC_XPASS
XFAIL_TESTS += snglcoinc_verify.py
endif
