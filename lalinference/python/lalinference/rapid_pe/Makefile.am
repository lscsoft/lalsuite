BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/rapid_pe

pymodule_PYTHON = \
	__init__.py \
    factored_likelihood.py \
    lalsimutils.py \
    mcsampler.py \
    statutils.py \
    xmlutils.py \
    common_cl.py \
	$(END_OF_LIST)

endif
