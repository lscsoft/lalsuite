pkginclude_HEADERS = \
	AVFactories.h \
	SeqFactories.h \
	StringVector.h

CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =

all-local: include-link

@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
