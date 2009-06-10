## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS = -I misc
SUBDIRS = src doc man misc debian
EXTRA_DIST = README README.bugs README.contrib \
	lalapps.spec.in lalapps.spec

# cvs-clean: more than maintainer-clean, this cleans
# everything not in the cvs archive.
cvs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f */Makefile.in
	rm -f src/*/Makefile.in
	rm -f src/pulsar/*/Makefile.in
	rm -f src/pulsar/*/*/Makefile.in
	rm -f src/stamp-h.in
	rm -f src/config.h.in
	rm -f configure
	rm -f aclocal.m4
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f lalapps.spec

all-local:
	@echo "";\
echo "=================================================================";\
echo "";\
echo "	LALApps has now been successfully built.";\
echo "	Now run 'make install' to install LALApps.";\
echo "";\
echo "================================================================="

install-exec-local:
	@echo "";\
echo "=================================================================";\
echo "";\
echo "	LALApps has now been successfully installed.";\
echo "";\
echo "	Now you may want to setup your environment:";\
echo "";\
echo "	Bourne shell [bash] users: please add the";\
echo "	following lines to your .profile file";\
echo "";\
echo '		. $(sysconfdir)/lalapps-user-env.sh';\
echo "";\
echo "	C-shell [tcsh] users: please add the following";\
echo "	lines to your .login file";\
echo "";\
echo '	source $(sysconfdir)/lalapps-user-env.csh';\
echo "";\
echo "================================================================="
