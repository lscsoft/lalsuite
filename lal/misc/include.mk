## This file gets automatically appended to the end of the Makefile
include-link:
	@for file in $(HEADERS) ; do \
          sourcedir=`cd $(srcdir) && pwd`; \
          targetdir=`cd $(top_builddir)/include/lal && pwd`; \
          if test ! -f $$targetdir/$$file ; then \
            @LN_S@ $$sourcedir/$$file $$targetdir; \
            echo "@LN_S@ $$sourcedir/$$file $$targetdir"; \
          fi; \
        done
