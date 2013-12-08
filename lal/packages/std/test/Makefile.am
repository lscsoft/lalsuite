AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	LALConstantsTest \
	LALGSLTest \
	LALMallocTest \
	LALVersionTest \
	StringConvertTest

TESTS = \
	$(check_PROGRAMS)
