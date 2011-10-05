@SWIG_COMMON_MK@

if SWIG_BUILD_OCTAVE

swig_language = octave
swig_wrapobj = $(swig_wrapname).oct

# start from scratch if build script changes
$(swig_wrapobj) : $(srcdir)/mkoct.m

# build SWIG wrapping code
$(swig_wrapobj) :
	$(OCTAVE) $(srcdir)/mkoct.m

# test SWIG wrapping module
check .PHONY : check-swig
check-swig : $(swig_wrapobj)
	test -f $(swig_wrapcheck).m && env \
	$(swig_ldlibpathname)="$(swig_ldlibpath):$${$(swig_ldlibpathname)}" \
	OCTAVE_PATH="$${PWD}:$${OCTAVE_PATH}" \
	$(OCTAVE) $(swig_wrapcheck).m

# install SWIG wrapping module
install-exec-local :
	$(MKDIR_P) $(octfiledir)
	$(INSTALL) $(swig_wrapobj) $(octfiledir)

# uninstall SWIG wrapping module
uninstall-local :
	rm -f $(octfiledir)/$(swig_wrapobj)

# cleanup
CLEANFILES = $(swig_cleanfiles) \
             $(swig_wrapname).o \
             $(swig_wrapobj)

EXTRA_DIST = mkoct.m

endif # SWIG_BUILD_OCTAVE
