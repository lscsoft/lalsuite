## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALConfig.h
EXTRA_DIST = LALConfig.h.in
DISTCLEANFILES = *.h
all-local:
	@set fnord $(MAKEFLAGS); amf=$$2; \
		d=$(top_builddir)/packages; \
		echo "Making include-links in $$d"; \
		(cd $$d && $(MAKE) $(AM_MAKEFLAGS) include-links) \
			|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
		test -z "$$fail"
#		d=$(top_builddir)/lib; \
#		(cd $$d && $(MAKE) $(AM_MAKEFLAGS) include-link) \
#			|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
#		test -z "$$fail"
