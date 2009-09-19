ACLOCAL_AMFLAGS = -I gnuscripts
EXTRA_DIST = doc lalstochastic.pc.in lalstochastic.spec.in \
	lalstochastic.spec aclocal.m4 \
	gnuscripts/config.guess \
	gnuscripts/config.sub \
	gnuscripts/depcomp \
	gnuscripts/header-link.mk \
	gnuscripts/install-sh \
	gnuscripts/lal_c99.m4 \
	gnuscripts/lalstochastic.m4 \
	gnuscripts/libtool.m4 \
	gnuscripts/ltmain.sh \
	gnuscripts/lt~obsolete.m4 \
	gnuscripts/ltoptions.m4 \
	gnuscripts/ltsugar.m4 \
	gnuscripts/ltversion.m4 \
	gnuscripts/missing \
	gnuscripts/pkg.m4
SUBDIRS = debian include src test

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
