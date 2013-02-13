ACLOCAL_AMFLAGS = -I gnuscripts

# Override automake's default rule for rebuilding ./configure, to ensure that
# it is successfully recreated after certain build system changes which also
# require other autotools (e.g. autoheader) to be run. See Redmine issue #728.
$(top_srcdir)/configure: $(am__configure_deps)
	$(am__cd) $(srcdir) && $(SHELL) ./00boot

EXTRA_DIST = \
	COPYING \
	doc \
	lalstochastic.pc.in \
	lalstochastic.spec.in \
	lalstochastic.spec \
	aclocal.m4

SUBDIRS = \
	gnuscripts \
	debian \
	include \
	src \
	test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalstochastic.pc

DISTCLEANFILES = config.lt

if HAVE_GIT_REPO
vcs-clean: maintainer-clean
	-cd $(top_srcdir) && $(GIT) clean -Xdf
endif

all-local:
	@echo "";\
	echo "==================================================";\
	echo "LALStochastic has now been successfully built.";\
	echo "Now run 'make install' to install LALStochastic.";\
	echo "==================================================";\
	echo ""

install-exec-local:
	@echo "";\
	echo "==================================================";\
	echo "LALStochastic has now been successfully installed.";\
	echo "Now you may want to setup your environment:";\
	echo "";\
	echo "* Bourne shell [bash] users: please add the line";\
	echo "  . $(sysconfdir)/lalstochastic-user-env.sh";\
	echo "  to your .profile file.";\
	echo "";\
	echo "* C-shell [tcsh] users: please add the line";\
	echo "  source $(sysconfdir)/lalstochastic-user-env.csh";\
	echo "  to your .login file.";\
	echo "==================================================";\
	echo ""
