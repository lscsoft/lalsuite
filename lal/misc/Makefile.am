EXTRA_DIST = \
	ax_pthread.m4 \
	config.guess \
	config.sub \
	depcomp \
	doc.mk \
	dvi.mk \
	include.mk \
	install-sh \
	lal.m4 \
	lalsuite_c99.m4 \
	libtool.m4 \
	ltmain.sh \
	ltoptions.m4 \
	ltsugar.m4 \
	ltversion.m4 \
	lt~obsolete.m4 \
	main.tex.in \
	missing \
	mkinstalldirs \
	pkg.m4

CLEANFILES = lal-user-env.sh lal-user-env.csh

sysconf_DATA = lal-user-env.sh lal-user-env.csh
BUILT_SOURCES = lal-user-env.sh lal-user-env.csh

all-local: include-link

lal-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL" >$@ ;\
	echo "LAL_PREFIX=\"$(prefix)\"" >>$@ ;\
	echo "export LAL_PREFIX" >>$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "LD_LIBRARY_PATH=\"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "DYLD_LIBRARY_PATH=\"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "PKG_CONFIG_PATH=\"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "export PATH MANPATH LD_LIBRARY_PATH DYLD_LIBRARY_PATH PKG_CONFIG_PATH" >>$@

lal-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LAL" >$@ ;\
	echo "setenv LAL_PREFIX \"$(prefix)\"" >>$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?LD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir):\$${LD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv LD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?DYLD_LIBRARY_PATH ) then" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir):\$${DYLD_LIBRARY_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv DYLD_LIBRARY_PATH \"$(libdir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?PKG_CONFIG_PATH ) then" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig:\$${PKG_CONFIG_PATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "	setenv PKG_CONFIG_PATH \"$(libdir)/pkgconfig\"" >>$@ ;\
	echo "endif" >>$@

@INCLUDEMKFILE@
