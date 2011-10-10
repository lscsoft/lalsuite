@SWIG_COMMON_MK@

SUBDIRS = octave python

if SWIG_BUILD

# install SWIG interface headers
EXTRA_HEADERS =
swigincludedir = $(includedir)/lal
swiginclude_HEADERS = $(swig_ifacefile) \
                      $(SWIG_HEADERS)

# create links to SWIG interface headers
.PHONY : iface-links
iface-links :
	@for file in $(swiginclude_HEADERS); do \
		basename=`echo "/$$file" | $(SED) 's|^.*/\([^/]*\)$$|\1|'`; \
		srcdir=`cd $(top_srcdir)/$(subdir) && pwd`; \
		builddir=`cd $(top_builddir)/$(subdir) && pwd`; \
		incldir=`cd $(top_builddir)/include/lal && pwd`; \
		if test ! -r $$incldir/$$basename; then \
			rm -f $$incldir/$$basename; \
			if test -f $$srcdir/$$file; then \
				$(LN_S) $$srcdir/$$file $$incldir; \
			elif test -f $$builddir/$$file; then \
				$(LN_S) $$builddir/$$file $$incldir; \
			fi; \
		fi; \
	done

# extra files for distribution
EXTRA_DIST = swig-common.mk \
             swig-header.mk

endif # SWIG_BUILD
