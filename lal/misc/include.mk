## This file gets automatically appended to the end of the Makefile
include-link:
	@for file in $(pkginclude_HEADERS) ; do \
		sourcedir=$(abs_srcdir); \
		builddir=$(abs_builddir); \
		targetdir=$(abs_top_builddir)/include/lal; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			if test -f $$sourcedir/$$file; then \
				$(LN_S) $$sourcedir/$$file $$targetdir; \
			elif test -f $$builddir/$$file; then \
				$(LN_S) $$builddir/$$file $$targetdir; \
			fi; \
		fi; \
	done
