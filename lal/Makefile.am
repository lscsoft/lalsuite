## Process this file with automake to produce Makefile.in
SUBDIRS = include packages doc lib
#CLEANFILES = $(top_builddir)/include/*.h
EXTRA_DIST = README.bugs README.contrib acconfig.h \
  misc/doc.mk misc/dvi.mk misc/include.mk misc/main.tex.in


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
	rm -f configure.in
	rm -f aclocal.m4
	rm -f testscript.log
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f config.h.in
