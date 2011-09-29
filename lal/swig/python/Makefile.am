@SWIG_COMMON_MK@

if SWIG_BUILD_PYTHON

swig_language = python
swig_language_opts = -modern
swig_wrapobj = _$(swig_wrapname).so

# start from scratch if build script changes
$(swig_wrapobj) : $(srcdir)/setup.py

# build SWIG wrapping code
$(swig_wrapobj) :
	$(PYTHON) $(srcdir)/setup.py build
	rm -f $(swig_wrapobj) && $(LN_S) build/lib/$(swig_wrapobj) $(swig_wrapobj)

# test SWIG wrapping module
check .PHONY : check-swig
check-swig : $(swig_wrapobj)
	test -f $(swig_wrapcheck).py && \
	$(swig_ldlibpath) \
	PYTHONPATH="$${PWD}:$${PYTHONPATH}" \
	$(PYTHON) $(swig_wrapcheck).py

# install SWIG wrapping module
install-exec-local :
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(prefix) --record .install.txt

# uninstall SWIG wrapping module
uninstall-local :
	-test -f .install.txt && rm -f `cat .install.txt`

# cleanup
clean-local :
	$(PYTHON) $(srcdir)/setup.py clean --all

CLEANFILES = $(swig_cleanfiles) \
             $(swig_wrapname).py \
             $(swig_wrapname).pyc \
             $(swig_wrapobj)

DISTCLEANFILES = .install.txt

EXTRA_DIST = setup.py

endif # SWIG_BUILD_PYTHON
