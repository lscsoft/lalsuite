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

# generic environment setup for SWIG wrapping module
$(swig_envfile).shell ::
	@echo '# source to access $(PACKAGE_NAME) SWIG wrapping modules for:' > $@; \
	for subdir in $(SUBDIRS); do \
		( cd $$subdir && $(MAKE) $@ ); \
		test ! -f $$subdir/$@ || cat $$subdir/$@ >> $@; \
	done

# generate environmental setups for specific shells
$(swig_envfile).sh : $(swig_envfile).shell
	@${SED} 's|^\([A-Z_][A-Z_]*\)=|export \1=|' $< > $@
$(swig_envfile).csh : $(swig_envfile).shell
	@${SED} 's|^\([A-Z_][A-Z_]*\)=|[ $${?\1} ] \&\& setenv \1; setenv \1 |' $< > $@
install-data-local : $(swig_envfile).sh $(swig_envfile).csh
	$(MKDIR_P) $(sysconfdir)
	$(INSTALL_DATA) $(swig_envfile).sh $(sysconfdir)
	$(INSTALL_DATA) $(swig_envfile).csh $(sysconfdir)

# cleanup
CLEANFILES = $(swig_envfile).shell \
             $(swig_envfile).sh \
             $(swig_envfile).csh

# extra files for distribution
EXTRA_DIST = swig-common.mk \
             swig-header.mk

else # !SWIG_BUILD

# do not generate environmental setups for specific shells
$(swig_envfile).sh :
$(swig_envfile).csh :

endif # SWIG_BUILD
