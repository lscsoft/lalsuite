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

octexec_LTLIBRARIES = lal.la
lal_la_SOURCES = $(swig_lt_sources)
nodist_lal_la_SOURCES = $(swig_lt_nodist_sources)
lal_la_CPPFLAGS = $(swig_lt_cppflags)
lal_la_CFLAGS = $(swig_lt_cflags)
lal_la_CXXFLAGS = $(swig_lt_cxxflags)
lal_la_LDFLAGS = $(swig_lt_ldflags)
lal_la_LIBADD = $(swig_lt_libadd)
lal_la_LIBTOOLFLAGS = $(swig_lt_libtoolflags)

endif # SWIG_BUILD_OCTAVE
