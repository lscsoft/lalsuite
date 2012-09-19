AM_CPPFLAGS = -I$(top_builddir)/include -I$(srcdir)/../include
LDADD = \
	$(top_builddir)/lib/lal/liblal.la \
	$(top_builddir)/packages/support/src/liblalsupport.la

check_PROGRAMS = \
	VectorOpsTest \
	VectorIndexRangeTest

TESTS = \
	$(check_PROGRAMS)
