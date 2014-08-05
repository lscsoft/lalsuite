.PHONY: force-header-link-generation

EXTRA_DIST = header-link-stamp
BUILT_SOURCES = header-link-stamp
DISTCLEANFILES = header-link-stamp *.h *.i

header-link-stamp: force-header-link-generation
	$(AM_V_at)for file in $(DISTCLEANFILES); do \
		test -r "$$file" || rm -f "$$file"; \
	done; \
	if test $(PACKAGE) = lal; then \
		srcdirs="packages lib"; \
	else \
		srcdirs="src"; \
	fi; \
	for dir in $$srcdirs swig; do \
		( cd "$(top_builddir)/$$dir" && $(MAKE) $(AM_MAKEFLAGS) header-links ) || exit 1; \
	done; \
	echo "$@" > $@
