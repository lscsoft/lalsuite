EXTRA_DIST =
SUBDIRS =
include $(top_srcdir)/gnuscripts/lalsuite_test.am

# Add compiled test programs to this variable
test_programs +=

# Add shell, Python, etc. test scripts to this variable
if HAVE_PYTHON
SUBDIRS += python

test_scripts += \
	test_bucluster.sh \
	$(END_OF_LIST)
test_extra_files = \
	H1-POWER_S5-830046410-50274.xml.gz \
	$(END_OF_LIST)
endif

# Add any helper programs required by tests to this variable
test_helpers +=

MOSTLYCLEANFILES = \
	bucluster_test.xml.gz \
	$(END_OF_LIST)

EXTRA_DIST += \
	H1-POWER_S5-830046410-50274.xml.gz \
	$(END_OF_LIST)
