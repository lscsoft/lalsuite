AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/lib/lalsupport/src/liblalsupport.la

check_PROGRAMS = \
	LALMomentTest \
	LALCorrelationTest

LALMomentTest_SOURCES = \
	LALMomentTest.c \
	CheckStatus.c

TESTS = \
	$(check_PROGRAMS)

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h
