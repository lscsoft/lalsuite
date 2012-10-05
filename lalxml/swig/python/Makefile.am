# this Makefile.am should NOT be modified; all SWIG configuration should
# be performed in lalsuite_swig.m4, swiglal_common.am, or swiglal_python.am

CPPFLAGS = $(SWIG_CPPFLAGS)
CFLAGS = $(SWIG_CFLAGS)
CXXFLAGS = $(SWIG_CXXFLAGS)
LDFLAGS = $(SWIG_LDFLAGS)
LIBS = $(SWIG_LIBS)

TEST_EXTENSIONS = .py

include $(top_srcdir)/swig/python/swiglal_python.am

if SWIG_BUILD_PYTHON

pkgpyexec_LTLIBRARIES = _lalxml.la
_lalxml_la_SOURCES = $(swig_lt_sources)
nodist__lalxml_la_SOURCES = $(swig_lt_nodist_sources)
_lalxml_la_CPPFLAGS = $(swig_lt_cppflags)
_lalxml_la_CFLAGS = $(swig_lt_cflags)
_lalxml_la_CXXFLAGS = $(swig_lt_cxxflags)
_lalxml_la_LDFLAGS = $(swig_lt_ldflags)
_lalxml_la_LIBADD = $(swig_lt_libadd)
_lalxml_la_LIBTOOLFLAGS = $(swig_lt_libtoolflags)

endif # SWIG_BUILD_PYTHON
