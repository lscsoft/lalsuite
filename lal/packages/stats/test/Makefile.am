AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

check_PROGRAMS = \
	LALMomentTest \
	LALCorrelationTest

LALMomentTest_SOURCES = \
	LALMomentTest.c \
	CheckStatus.c

TESTS = \
	$(check_PROGRAMS)

TESTS_ENVIRONMENT = 1>/dev/null

EXTRA_DIST = \
	CheckStatus.c \
	CheckStatus.h
