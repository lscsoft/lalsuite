## Process this file with automake to produce Makefile.in
pkginclude_HEADERS = LALConfig.h
EXTRA_DIST = LALConfig.h.in
DISTCLEANFILES = *.h
PKGLIST = @PKGLIST@
all-local:
	@set fnord $(MAKEFLAGS); amf=$$2; \
        for pkg in $(PKGLIST); do \
          subdir=$(top_builddir)/packages/$$pkg/include; \
          echo "Making all in $$subdir"; \
          (cd $$subdir && $(MAKE) $(AM_MAKEFLAGS) all) \
           || case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
        done; test -z "$$fail"
