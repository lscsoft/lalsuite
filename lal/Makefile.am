ACLOCAL_AMFLAGS = -I misc

SUBDIRS = \
	debian \
	misc \
	include \
	packages \
	lib \
	swig

MAINTAINERCLEANFILES = lal.spec
DISTCLEANFILES = config.lt
distcleancheck_listfiles = find . -type f \! -name "LALVCSInfo.h" -print

EXTRA_DIST = \
	COPYING \
	README.bugs \
	README.contrib \
	README.install \
	lal.pc.in \
	lalsupport.pc.in \
	aclocal.m4 \
	lal.spec.in \
	lal.spec

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = \
	lal.pc \
	lalsupport.pc

# vcs-clean: more than maintainer-clean, this cleans
# everything not in the vcs archive.
vcs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f */Makefile.in
	rm -f */*/Makefile.in
	rm -f */*/*/Makefile.in
	rm -f include/stamp-h*.in
	rm -f include/config.h.in
	rm -f include/lal/stamp-h*.in
	rm -f configure
	rm -f aclocal.m4
	rm -f acinclude.m4
	rm -f acconfig.h
	rm -f testscript.log
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f config.h.in

# map cvs-clean to vcs-clean for backwards compatibility
cvs-clean: vcs-clean
	@echo "";\
		echo "=================================================================";\
		echo "";\
		echo "        cvs-clean is deprecated, use vcs-clean instead";\
		echo "";\
		echo "================================================================="

all-local:
	@echo "";\
		echo "=================================================================";\
		echo "";\
		echo "        LAL has now been successfully built.";\
		echo "        Now run 'make install' to install LAL.";\
		echo "";\
		echo "================================================================="

install-exec-local:
	@echo "";\
		echo "=================================================================";\
		echo "";\
		echo "        LAL has now been successfully installed.";\
		echo "";\
		echo "        Now you may want to setup your environment:";\
		echo "";\
		echo "        Bourne shell [bash] users: please add the";\
		echo "        following lines to your .profile file";\
		echo "";\
		echo '            . $(sysconfdir)/lal-user-env.sh';\
		echo "";\
		echo "        C-shell [tcsh] users: please add the following";\
		echo "        lines to your .login file";\
		echo "";\
		echo '            source $(sysconfdir)/lal-user-env.csh';\
		echo "";\
		echo "================================================================="
