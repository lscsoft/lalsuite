ACLOCAL_AMFLAGS = -I gnuscripts

EXTRA_DIST = \
	COPYING \
	lalinference.pc.in \
	lalinference.spec.in \
	lalinference.spec \
	aclocal.m4

SUBDIRS = \
	gnuscripts \
	debian \
	include \
	src \
	test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalinference.pc

DISTCLEANFILES = config.lt

if HAVE_GIT_REPO
vcs-clean: maintainer-clean
	-cd $(top_srcdir) && $(GIT) clean -Xdf
endif

all-local:
	@echo "";\
	echo "==================================================";\
	echo "LALInference has now been successfully built.";\
	echo "Now run 'make install' to install LALInference.";\
	echo "==================================================";\
	echo ""

install-exec-local:
	@echo "";\
	echo "==================================================";\
	echo "LALInference has now been successfully installed.";\
	echo "Now you may want to setup your environment:";\
	echo "";\
	echo "* Bourne shell [bash] users: please add the line";\
	echo "  . $(sysconfdir)/lalinference-user-env.sh";\
	echo "  to your .profile file.";\
	echo "";\
	echo "* C-shell [tcsh] users: please add the line";\
	echo "  source $(sysconfdir)/lalinference-user-env.csh";\
	echo "  to your .login file.";\
	echo "==================================================";\
	echo ""
