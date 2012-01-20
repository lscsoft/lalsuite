include $(top_srcdir)/swig/swig-common.am

if SWIG_BUILD
if SWIG_BUILD_PYTHON

# set language
swig_language = python
swig_language_opts = -modern

# name of SWIG wrapping module object
swig_wrap = _$(swig_wrapname).so

# start from scratch if build script changes
swig_wrap_depends += $(srcdir)/setup.py

# build SWIG wrapping module
$(swig_wrap) : $(swig_wrap_depends)
	$(PYTHON) $(srcdir)/setup.py build
	rm -f $(swig_wrap) && $(LN_S) build/lib/$(swig_wrap) $(swig_wrap)

# test SWIG wrapping module
swig-check : $(swig_wrap)
	test -f $(swig_wrapcheck).py && env \
	$(SWIG_LD_LIBPATH_NAME)="$(swig_ldlibpath):$${$(SWIG_LD_LIBPATH_NAME)}" \
	PYTHONPATH="$${PWD}:$${PYTHONPATH}" \
	$(PYTHON) $(swig_wrapcheck).py

# install SWIG wrapping module
swig-install : $(swig_wrap)
	$(PYTHON) $(srcdir)/setup.py install --prefix=$(DESTDIR)$(prefix) --record .install.txt

# uninstall SWIG wrapping module
swig-uninstall :
	-test -f .install.txt && rm -f `cat .install.txt`

# language-specific headers (with LAL only)
swiginclude_HEADERS += $(SWIG_PYTHON_HEADERS)

# extra distribution files
EXTRA_DIST += \
	$(swig_wrapcheck).py \
	$(srcdir)/setup.py

# cleanup
clean-local :
	$(PYTHON) $(srcdir)/setup.py clean --all

CLEANFILES += \
	$(swig_wrapname).py \
	$(swig_wrapname).pyc \
	$(swig_wrap)

DISTCLEANFILES = .install.txt

endif # SWIG_BUILD_PYTHON
endif # SWIG_BUILD
