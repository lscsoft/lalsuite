include $(top_srcdir)/swig/swig-common.am

if SWIG_BUILD
if SWIG_BUILD_OCTAVE

# set language
swig_language = octave

# name of SWIG wrapping module object
swig_wrap = $(swig_wrapname).oct

# start from scratch if build script changes
swig_wrap_depends += $(srcdir)/mkoct.m

# build SWIG wrapping module
$(swig_wrap) : $(swig_wrap_depends)
	$(OCTAVE) $(srcdir)/mkoct.m

# test SWIG wrapping module
swig-check : $(swig_wrap)
	test -f $(swig_wrapcheck).m && env \
	$(SWIG_LD_LIBPATH_NAME)="$(swig_ldlibpath):$${$(SWIG_LD_LIBPATH_NAME)}" \
	OCTAVE_PATH="$${PWD}:$${OCTAVE_PATH}" \
	$(OCTAVE) $(swig_wrapcheck).m

# install SWIG wrapping module
swig-install : $(swig_wrap)
	$(MKDIR_P) $(DESTDIR)$(octfiledir)
	$(INSTALL) $(swig_wrap) $(DESTDIR)$(octfiledir)

# uninstall SWIG wrapping module
swig-uninstall :
	rm -f $(DESTDIR)$(octfiledir)/$(swig_wrap)

# extra distribution files
EXTRA_DIST += \
	$(swig_wrapcheck).m \
	$(srcdir)/mkoct.m

# cleanup
CLEANFILES += \
	$(swig_wrapname).o \
	$(swig_wrap)

endif # SWIG_BUILD_OCTAVE
endif # SWIG_BUILD
