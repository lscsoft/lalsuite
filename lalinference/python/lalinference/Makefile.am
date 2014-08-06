BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST = cylon.py
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	bayestar \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	cmap.py \
	fits.py \
	plot.py \
	git_version.py \
	$(END_OF_LIST)

endif

dist_pymodule_DATA = \
	cylon.csv \
	ne_simplified_coastline.json \
	$(END_OF_LIST)

#
# Note: the following script generates a color map RGB data file.
# The output file is distributed to avoid a build dependency on the Python
# `colormath` package (https://pypi.python.org/pypi/colormath).
#
$(srcdir)/cylon.csv: $(srcdir)/cylon.py
	$(AM_V_GEN)$(PYTHON) $< > $@
