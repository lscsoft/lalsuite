# this Makefile.am should NOT be modified; all SWIG configuration should
# be performed in lalsuite_swig.m4, swiglal_common.am, or swiglal_octave.am

CPPFLAGS = $(SWIG_CPPFLAGS)
CFLAGS = $(SWIG_CFLAGS)
CXXFLAGS = $(SWIG_CXXFLAGS)
LDFLAGS = $(SWIG_LDFLAGS)
LIBS = $(SWIG_LIBS)

TEST_EXTENSIONS = .m

include $(top_srcdir)/swig/octave/swiglal_octave.am

if SWIG_BUILD_OCTAVE

octexec_LTLIBRARIES = laldetchar.la
laldetchar_la_SOURCES = $(swig_lt_sources)
nodist_laldetchar_la_SOURCES = $(swig_lt_nodist_sources)
laldetchar_la_CPPFLAGS = $(swig_lt_cppflags)
laldetchar_la_CFLAGS = $(swig_lt_cflags)
laldetchar_la_CXXFLAGS = $(swig_lt_cxxflags)
laldetchar_la_LDFLAGS = $(swig_lt_ldflags)
laldetchar_la_LIBADD = $(swig_lt_libadd)
laldetchar_la_LIBTOOLFLAGS = $(swig_lt_libtoolflags)

endif # SWIG_BUILD_OCTAVE
