## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALConfig.h
EXTRA_DIST = LALConfig.h.in
DISTCLEANFILES = *.h *.i
all-local:
	@for file in *.h *.i; do \
		test -r $$file || rm -f $$file; \
	done; \
	set fnord $(MAKEFLAGS); amf=$$2; \
		d=$(top_builddir)/packages; \
		(cd $$d && $(MAKE) $(AM_MAKEFLAGS) include-links) \
			|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
		test -z "$$fail"
		@d=$(top_builddir)/lib/lal; \
		(cd $$d && $(MAKE) $(AM_MAKEFLAGS) include-link) \
			|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
		test -z "$$fail"

include $(top_srcdir)/swig/swiglal_header.am
