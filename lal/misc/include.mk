## This file gets automatically appended to the end of the Makefile
include-link:
	@for file in $(HEADERS) ; do \
		targetdir="$(abs_top_builddir)/include/lal"; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			$(LN_S) $(abs_srcdir)/$$file $$targetdir; \
		fi; \
	done
