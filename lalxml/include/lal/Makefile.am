DISTCLEANFILES = *.h *.i
all-local:
	@for file in *.h *.i; do \
		test -r $$file || rm -f $$file; \
	done; \
	set fnord $(MAKEFLAGS); amf=$$2; \
	d=$(top_builddir)/src; \
	echo "Making header-links in $$d"; \
	(cd $$d && $(MAKE) $(AM_MAKEFLAGS) header-links) \
		|| case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
        test -z "$$fail"
