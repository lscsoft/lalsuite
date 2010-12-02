ACLOCAL_AMFLAGS = -I gnuscripts

SUBDIRS = \
	doxygen \
	gnuscripts \
	@subdirs@

distcleancheck_listfiles = find . -type f \! -name "LAL*VCSInfo.h" -print

vcs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f configure
	rm -f aclocal.m4
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f gnuscripts/Makefile.in

cvs-clean: vcs-clean
