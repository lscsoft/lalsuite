pkginclude_HEADERS = \
	AVFactories.h \
	SeqFactories.h \
	StringVector.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =

MAINTAINERCLEANFILES = \
	SeqFactories.h

EXTRA_DIST = \
	SeqFactoriesH.m4 \
	SeqFactoriesBaseH.m4 \
	ArraySeqFactoriesBaseH.m4

all-local: include-link

SeqFactories.h: SeqFactoriesH.m4 SeqFactoriesBaseH.m4 ArraySeqFactoriesBaseH.m4
	cd $(srcdir) && @M4@ SeqFactoriesH.m4 > SeqFactories.h

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swig-header.am
