## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALMalloc.h LALError.h PrintVector.h

all-local: $(top_srcdir)/include/LALMalloc.h \
           $(top_srcdir)/include/LALError.h \
           $(top_srcdir)/include/PrintVector.h

$(top_srcdir)/include/LALMalloc.h: LALMalloc.h
	cp -f LALMalloc.h $(top_srcdir)/include

$(top_srcdir)/include/LALError.h: LALError.h
	cp -f LALError.h $(top_srcdir)/include

$(top_srcdir)/include/PrintVector.h: PrintVector.h
	cp -f PrintVector.h $(top_srcdir)/include

