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
distcleancheck_listfiles = find . -type f \! -name "LALInferenceVCSInfo.h" -print

# vcs-clean: more than maintainer-clean, this cleans everything not in
# the vcs archive
vcs-clean: maintainer-clean
	-rm -f aclocal.m4
	-rm -f configure
	-rm -f Makefile.in
	-rm -f include/config.h.in
	-rm -f include/Makefile.in
	-rm -f include/lal/Makefile.in
	-rm -r src/Makefile.in
	-rm -f test/Makefile.in

all-local:
	@echo "";\
		echo "=================================================================";\
		echo "";\
		echo "        LAL Inference has now been successfully built.";\
		echo "        Now run 'make install' to install LAL Inference.";\
		echo "";\
		echo "================================================================="

install-exec-local:
	@echo "";\
		echo "=================================================================";\
		echo "";\
		echo "        LAL Inference has now been successfully installed.";\
		echo "";\
		echo "        Now you may want to setup your environment:";\
		echo "";\
		echo "        Bourne shell [bash] users: please add the";\
		echo "        following lines to your .profile file";\
		echo "";\
		echo '            . $(sysconfdir)/lalinference-user-env.sh';\
		echo "";\
		echo "        C-shell [tcsh] users: please add the following";\
		echo "        lines to your .login file";\
		echo "";\
		echo '            source $(sysconfdir)/lalinference-user-env.csh';\
		echo "";\
		echo "================================================================="
