EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs += BurstFileReadTest

# Add shell, Python, etc. test scripts to this variable
if SWIG_BUILD_PYTHON
test_scripts += \
	date_verify.py \
	snglcoinc_verify.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	$(END_OF_LIST)

EXTRA_DIST += \
    ad_hoc_test.xml.gz \
    test_inj.txt \
	a.pgm \
	CLRindata.asc \
	$(END_OF_LIST)
