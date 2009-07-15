## Process this file with automake to produce Makefile.in
check_PROGRAMS = TfrPswvTest TfrRspTest TfrSpTest TfrWvTest
TfrPswvTest_SOURCES = TfrPswvTest.c
TfrPswvTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TfrRspTest_SOURCES = TfrRspTest.c
TfrRspTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TfrSpTest_SOURCES = TfrSpTest.c
TfrSpTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
TfrWvTest_SOURCES = TfrWvTest.c
TfrWvTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
TESTS = TfrPswvTest TfrRspTest TfrSpTest TfrWvTest
TESTS_ENVIRONMENT = 1>/dev/null
EXTRA_DIST = TfrPswvTest.c TfrRspTest.c TfrSpTest.c TfrWvTest.c
@DVIMKFILE@
