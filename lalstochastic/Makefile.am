ACLOCAL_AMFLAGS = -I gnuscripts
EXTRA_DIST = doc lalstochastic.pc.in lalstochastic.spec.in \
	lalstochastic.spec aclocal.m4
SUBDIRS = gnuscripts debian include src test

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalstochastic.pc

maintainer-clean-local:
	-rm -f aclocal.m4
	-rm -f configure
	-rm -f Makefile.in
	-rm -f include/config.h.in
	-rm -f include/Makefile.in
	-rm -f include/lal/Makefile.in
	-rm -r src/Makefile.in
	-rm -f test/Makefile.in
