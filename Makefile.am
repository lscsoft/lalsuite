ACLOCAL_AMFLAGS = -I lal/misc
SUBDIRS = @subdirs@
distcleancheck_listfiles = find . -type f \! -name "LAL*VCSInfo.h" -print
vcs-clean: cvs-clean
cvs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f configure
	rm -f aclocal.m4
	rm -f config.cache
	rm -f config.status
	rm -f config.log
