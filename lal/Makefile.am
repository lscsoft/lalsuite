## Process this file with automake to produce Makefile.in
SUBDIRS = include packages doc lib
CLEANFILES = $(top_builddir)/include/*.h
EXTRA_DIST = README.contrib lalspec.pdf acconfig.h

# cvs-clean: more than maintainer-clean, this cleans
# everything not in the cvs archive.
cvs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f */Makefile.in
	rm -f */*/Makefile.in
	rm -f */*/*/Makefile.in
	rm -f include/stamp-h*.in
	rm -f include/config.h.in
	rm -f include/lal/stamp-h*.in
	rm -f configure
	rm -f aclocal.m4
