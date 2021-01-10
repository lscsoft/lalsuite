BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pybin_scripts = \
	lalpulsar_PiecewiseSearch \
	$(END_OF_LINE)

noinst_pybin_scripts = \
	lalpulsar_PiecewiseSearchTemplateBank \
	$(END_OF_LINE)

endif

# Add shell test scripts to this variable
if HAVE_PYTHON
test_scripts +=
endif

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files +=
