# Process this file with automake to produce Makefile.in
check_PROGRAMS = setMaskTest clusterTest LALGetLayerWaveletTest 
LALGetLayerWaveletTest_SOURCES = LALGetLayerWaveletTest.c 
LALGetLayerWaveletTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
setMaskTest_SOURCES = setMaskTest.c 
setMaskTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
clusterTest_SOURCES = clusterTest.c 
clusterTest_LDADD = $(top_builddir)/packages/support/src/liblalsupport.la $(top_builddir)/lib/liblal.la
INCLUDES = -I$(top_builddir)/include -I$(srcdir)/../include
#TESTS =  setMaskTest clusterTest LALGetLayerWaveletTest
TESTS =  clusterTest LALGetLayerWaveletTest
TESTS_ENVIRONMENT = 1>/dev/null LAL_DATA_PATH=$(srcdir)
EXTRA_DIST = test_getLayer.dat wavelet1.dat params_getLayer0.dat layer0.dat  \
        test_setMask.dat setMask_1_0_2.in params_1_0_2.dat setMask_1_0_2.out \
        test_cluster.dat params cluster_1_0_2.in cluster_1_0_2.out
MOSTLYCLEANFILES = .dvi-dep
DISTCLEANFILES = *.tex
dvi-local: .dvi-dep
@DVIMKFILE@
