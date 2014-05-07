all-local: header-links
header-links:
	@for file in $(pkginclude_HEADERS) ; do \
		sourcedir=$(abs_srcdir) ; \
		targetdir=$(abs_top_builddir)/include/lal ; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			$(LN_S) $$sourcedir/$$file $$targetdir; \
		fi; \
	done
	@for file in LALFrameConfig.h LALFrameVCSInfo.h ; do \
		d=$(abs_builddir) ; \
			targetdir=$(abs_top_builddir)/include/lal ; \
			if test ! -r $$targetdir/$$file ; then \
				rm -f $$targetdir/$$file; \
				$(LN_S) $$d/$$file $$targetdir; \
			fi; \
	done
