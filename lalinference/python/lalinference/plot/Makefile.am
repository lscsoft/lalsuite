BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/plot

pymodule_PYTHON = \
	__init__.py \
	allsky.py \
	backend_pdf_vcs_friendly.py \
	cmap.py \
	healpix.py \
	marker.py \
	poly.py \
	pp.py \
	$(END_OF_LIST)

dist_pymodule_DATA = \
	cylon.csv \
	ne_simplified_coastline.json \
	$(END_OF_LIST)

all-local: data-all-local
data-all-local: data-mostlyclean-local
	$(AM_V_at)if test "x$(builddir)" != "x$(srcdir)"; then \
		for file in $(dist_pymodule_DATA); do \
			if test -r "$(srcdir)/$${file}"; then \
				$(LN_S) "$(srcdir)/$${file}" "$(builddir)/$${file}" || exit 1; \
			fi; \
		done; \
	fi

mostlyclean-local: data-mostlyclean-local
data-mostlyclean-local:
	$(AM_V_at)if test "x$(builddir)" != "x$(srcdir)"; then \
		for file in $(dist_pymodule_DATA); do \
			if test -r "$(srcdir)/$${file}"; then \
				rm -f "$(builddir)/$${file}" || exit 1; \
			fi; \
		done; \
	fi

endif

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
