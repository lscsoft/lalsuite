ACLOCAL_AMFLAGS = -I gnuscripts

EXTRA_DIST = \
	lalmetaio.pc.in \
	lalmetaio.spec.in \
	lalmetaio.spec \
	aclocal.m4

SUBDIRS = \
	gnuscripts \
	debian \
	include \
	src

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lalmetaio.pc

maintainer-clean-local:
	-rm -f aclocal.m4
	-rm -f configure
	-rm -f Makefile.in
	-rm -f include/config.h.in
	-rm -f include/Makefile.in
	-rm -f include/lal/Makefile.in
	-rm -r src/Makefile.in
	-rm -f test/Makefile.in
