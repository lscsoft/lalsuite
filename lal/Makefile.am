## Process this file with automake to produce Makefile.in
SUBDIRS = include packages doc lib
CLEANFILES = $(srcdir)/include/*.h
EXTRA_DIST = lalspec.pdf mklibmake

#disable autoheader
AUTOHEADER=echo

# cvs-clean: more than maintainer-clean, this cleans
# everything not in the cvs archive.
cvs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f */Makefile.in
	rm -f */*/Makefile.in
	rm -f */*/*/Makefile.in
	rm -f lib/Makefile.am
	rm -f include/stamp-h.in
	rm -f configure
	rm -f aclocal.m4
