EXTRA_DIST = \
	ax_pthread.m4 \
	config.guess \
	config.sub \
	cuda.common \
	depcomp \
	generate_user_env.awk \
	include.mk \
	install-sh \
	lal.m4 \
	lalsuite_build.m4 \
	lalsuite_gccflags.m4 \
	libtool.m4 \
	ltmain.sh \
	ltoptions.m4 \
	ltsugar.m4 \
	ltversion.m4 \
	lt~obsolete.m4 \
	missing \
	mkinstalldirs \
	test-driver

CLEANFILES = \
	lal-user-env.shell \
	lal-user-env.sh \
	lal-user-env.csh

sysconf_DATA = \
	lal-user-env.sh \
	lal-user-env.csh

BUILT_SOURCES = \
	lal-user-env.sh \
	lal-user-env.csh

all-local: include-link
@INCLUDEMKFILE@

lal-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LAL_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir) $(pythondir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@printf "$(SWIG_USER_ENV)" >>$@

lal-user-env.csh lal-user-env.sh: lal-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }
