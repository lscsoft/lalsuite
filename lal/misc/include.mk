## This file gets automatically appended to the end of the Makefile
include-link:
	@for file in $(pkginclude_HEADERS) ; do \
		sourcedir=$(abs_srcdir); \
		targetdir=$(abs_top_builddir)/include/lal; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			$(LN_S) $$sourcedir/$$file $$targetdir; \
		fi; \
	done
	@for file in LALVCSInfo.h ; do \
		if [ -f $$file ]; then \
			d=$(abs_builddir) ; \
			targetdir=$(abs_top_builddir)/include/lal ; \
			if test ! -r $$targetdir/$$file ; then \
				rm -f $$targetdir/$$file; \
				$(LN_S) $$d/$$file $$targetdir; \
			fi; \
		fi; \
	done
