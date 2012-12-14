ACLOCAL_AMFLAGS = -I misc

EXTRA_DIST = \
	COPYING \
	README \
	README.bugs \
	README.contrib \
	lalapps.spec.in \
	lalapps.spec

SUBDIRS = \
	debian \
	misc \
	src \
	doc \
	man

DISTCLEANFILES = config.lt

if HAVE_GIT_REPO
vcs-clean: maintainer-clean
	-cd $(top_srcdir) && $(GIT) clean -Xdf
endif

all-local:
	@echo "";\
	echo "==================================================";\
	echo "LALApps has now been successfully built.";\
	echo "Now run 'make install' to install LALApps.";\
	echo "==================================================";\
	echo ""

install-exec-local:
	@echo "";\
	echo "==================================================";\
	echo "LALApps has now been successfully installed.";\
	echo "Now you may want to setup your environment:";\
	echo "";\
	echo "* Bourne shell [bash] users: please add the line";\
	echo "  . $(sysconfdir)/lalapps-user-env.sh";\
	echo "  to your .profile file.";\
	echo "";\
	echo "* C-shell [tcsh] users: please add the line";\
	echo "  source $(sysconfdir)/lalapps-user-env.csh";\
	echo "  to your .login file.";\
	echo "==================================================";\
	echo ""
