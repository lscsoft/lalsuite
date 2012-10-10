DISTCLEANFILES = *.h
all-local:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	d=$(top_builddir)/src; \
	echo "Making header-links in $$d"; \
	(cd $$d && $(MAKE) $(AM_MAKEFLAGS) header-links) \
		|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
        test -z "$$fail"
