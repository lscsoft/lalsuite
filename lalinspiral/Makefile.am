ACLOCAL_AMFLAGS = -I gnuscripts

if HAVE_GIT_REPO
# Override automake's default rule for rebuilding aclocal.m4, to ensure that
# ./configure is successfully recreated after certain build system changes which
# require other autotools (e.g. autoheader) to be re-run. See Redmine issue #728.
$(ACLOCAL_M4): $(am__aclocal_m4_deps)
	$(am__cd) $(srcdir) && $(SHELL) ./00boot
endif # HAVE_GIT_REPO

EXTRA_DIST = \
	doc \
	lalinspiral.spec \
	$(DATA)

SUBDIRS = \
	debian \
	include \
	src \
	python \
	test \
	swig

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalinspiral.pc
noinst_DATA = lalinspiral-uninstalled.pc
BUILT_SOURCES = lalinspiral-uninstalled.pc

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
