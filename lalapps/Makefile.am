## Process this file with automake to produce Makefile.in
SUBDIRS = src doc man misc
EXTRA_DIST = acconfig.h lalapps.spec.in lalapps.spec

# cvs-clean: more than maintainer-clean, this cleans
# everything not in the cvs archive.
cvs-clean: maintainer-clean
	rm -f Makefile.in
	rm -f */Makefile.in
	rm -f src/*/Makefile.in
	rm -f src/pulsar/*/Makefile.in
	rm -f src/pulsar/hough/src/Makefile.in
	rm -f src/stamp-h.in
	rm -f src/config.h.in
	rm -f configure
	rm -f aclocal.m4
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f lalapps.spec
