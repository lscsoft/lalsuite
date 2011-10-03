pkginclude_HEADERS = \
	LALNoiseModels.h
CLEANFILES =
MOSTLYCLEANFILES =
DISTCLEANFILES =
all-local: include-link
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
CLEANFILES += $(swig_clean)
@SWIG_HEADER_MK@
