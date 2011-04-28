all-local: header-links
header-links:
	@for file in $(lalinspiralinclude_HEADERS) ; do \
		sourcedir=`cd $(srcdir) && pwd`; \
		targetdir=`cd $(top_builddir)/include/lal && pwd`; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			$(LN_S) $$sourcedir/$$file $$targetdir; \
		fi; \
	done
	@for file in LALInspiralConfig.h LALInspiralVCSInfo.h ; do \
		d=`pwd`; \
			targetdir=`cd $(top_builddir)/include/lal && pwd`; \
			if test ! -r $$targetdir/$$file ; then \
				rm -f $$targetdir/$$file; \
				$(LN_S) $$d/$$file $$targetdir; \
			fi; \
	done
