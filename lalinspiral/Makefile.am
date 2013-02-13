ACLOCAL_AMFLAGS = -I gnuscripts

# Override automake's default rule for rebuilding ./configure, to ensure that
# it is successfully recreated after certain build system changes which also
# require other autotools (e.g. autoheader) to be run. See Redmine issue #728.
$(top_srcdir)/configure: $(am__configure_deps)
	$(am__cd) $(srcdir) && $(SHELL) ./00boot

EXTRA_DIST = \
	COPYING \
	lalinspiral.pc.in \
	lalinspiral.spec.in \
	lalinspiral.spec \
	aclocal.m4 \
	doc

SUBDIRS = \
	gnuscripts \
	debian \
	include \
	src \
	python \
	swig \
	test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalinspiral.pc

DISTCLEANFILES = config.lt

if HAVE_GIT_REPO
vcs-clean: maintainer-clean
	-cd $(top_srcdir) && $(GIT) clean -Xdf
endif

all-local:
	@echo "";\
	echo "==================================================";\
	echo "LALInspiral has now been successfully built.";\
	echo "Now run 'make install' to install LALInspiral.";\
	echo "==================================================";\
	echo ""

install-exec-local:
	@echo "";\
	echo "==================================================";\
	echo "LALInspiral has now been successfully installed.";\
	echo "Now you may want to setup your environment:";\
	echo "";\
	echo "* Bourne shell [bash] users: please add the line";\
	echo "  . $(sysconfdir)/lalinspiral-user-env.sh";\
	echo "  to your .profile file.";\
	echo "";\
	echo "* C-shell [tcsh] users: please add the line";\
	echo "  source $(sysconfdir)/lalinspiral-user-env.csh";\
	echo "  to your .login file.";\
	echo "==================================================";\
	echo ""
