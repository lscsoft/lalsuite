ACLOCAL_AMFLAGS = -I gnuscripts

EXTRA_DIST = \
	COPYING \
	lalmetaio.pc.in \
	lalmetaio.spec.in \
	lalmetaio.spec \
	aclocal.m4

SUBDIRS = \
	gnuscripts \
	debian \
	include \
	src

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalmetaio.pc

DISTCLEANFILES = config.lt

if HAVE_GIT_REPO
vcs-clean: maintainer-clean
	-cd $(top_srcdir) && $(GIT) clean -Xdf
endif

all-local:
	@echo "";\
	echo "==================================================";\
	echo "LALMetaIO has now been successfully built.";\
	echo "Now run 'make install' to install LALMetaIO.";\
	echo "==================================================";\
	echo ""

install-exec-local:
	@echo "";\
	echo "==================================================";\
	echo "LALMetaIO has now been successfully installed.";\
	echo "Now you may want to setup your environment:";\
	echo "";\
	echo "* Bourne shell [bash] users: please add the line";\
	echo "  . $(sysconfdir)/lalmetaio-user-env.sh";\
	echo "  to your .profile file.";\
	echo "";\
	echo "* C-shell [tcsh] users: please add the line";\
	echo "  source $(sysconfdir)/lalmetaio-user-env.csh";\
	echo "  to your .login file.";\
	echo "==================================================";\
	echo ""
