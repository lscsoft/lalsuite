## Process this file with automake to produce Makefile.in
SUBDIRS = misc include packages doc lib
MAINTAINERCLEANFILES = README.install
EXTRA_DIST = README.bugs README.contrib README.install \
  lal-config.in lal.pc.in acinclude.m4 aclocal.m4 \
  misc/doc.mk misc/dvi.mk misc/include.mk misc/main.tex.in \
  misc/lal.m4 misc/lallibtool.m4 misc/lalam15patch.m4 \
  misc/acconfig.h
BUILT_SOURCES = README.install
bin_SCRIPTS = lal-config
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lal.pc

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
	rm -f acinclude.m4
	rm -f acconfig.h
	rm -f testscript.log
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f config.h.in
	rm -f configure.ac

# make README.install
README.install:
	@set fnord $(MAKEFLAGS); amf=$$2; \
	echo "Making README.install in $(top_builddir)/doc"; \
	(cd $(top_builddir)/doc && $(MAKE) $(AM_MAKEFLAGS) README.install) \
	  || case "$$amf" in *=*) exit1;; *k*) fail=yes;; *) exit 1;; esac; \
	test -z "$$fail"
