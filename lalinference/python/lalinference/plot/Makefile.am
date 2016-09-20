BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/plot

pymodule_PYTHON = \
	__init__.py \
	backend_pdf_vcs_friendly.py \
	cmap.py \
	marker.py \
	$(END_OF_LIST)

endif

dist_pymodule_DATA = \
	cylon.csv \
	ne_simplified_coastline.json \
	$(END_OF_LIST)

EXTRA_DIST += \
	cylon.py \
	$(END_OF_LIST)

#
# Note: the following script generates a color map RGB data file.
# The output file is distributed to avoid a build dependency on the Python
# `colormath` package (https://pypi.python.org/pypi/colormath).
#
# $(srcdir)/cylon.csv: $(srcdir)/cylon.py
# 	$(AM_V_GEN)$(PYTHON) $< > $@
