# this Makefile.am should NOT be modified; all SWIG configuration should
# be performed in lalsuite_swig.m4, swiglal_common.am, or swiglal_octave.am

CPPFLAGS = $(SWIG_CPPFLAGS)
CFLAGS = $(SWIG_CFLAGS)
CXXFLAGS = $(SWIG_CXXFLAGS)
LDFLAGS = $(SWIG_LDFLAGS)
LIBS = $(SWIG_LIBS)

TEST_EXTENSIONS = .m

if SWIG_BUILD
if SWIG_BUILD_OCTAVE

include $(top_srcdir)/swig/swiglal_common.am
include $(top_srcdir)/swig/octave/swiglal_octave.am

octexec_LTLIBRARIES = lalframe.la
lalframe_la_SOURCES = $(swig_lt_sources)
nodist_lalframe_la_SOURCES = $(swig_lt_nodist_sources)
lalframe_la_CPPFLAGS = $(swig_lt_cppflags)
lalframe_la_CFLAGS = $(swig_lt_cflags)
lalframe_la_CXXFLAGS = $(swig_lt_cxxflags)
lalframe_la_LDFLAGS = $(swig_lt_ldflags)
lalframe_la_LIBADD = $(swig_lt_libadd)
lalframe_la_LIBTOOLFLAGS = $(swig_lt_libtoolflags)

endif # SWIG_BUILD_OCTAVE
endif # SWIG_BUILD
