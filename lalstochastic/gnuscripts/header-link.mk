all-local: header-links
header-links:
	@for file in $(HEADERS) ; do \
		sourcedir=`cd $(srcdir) && pwd`; \
		targetdir=`cd $(top_builddir)/include/lal && pwd`; \
		if test ! -r $$targetdir/$$file ; then \
			rm -f $$targetdir/$$file; \
			$(LN_S) $$sourcedir/$$file $$targetdir; \
		fi; \
	done
