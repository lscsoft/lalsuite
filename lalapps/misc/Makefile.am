EXTRA_DIST = \
	ax_mpi.m4 \
	compile \
	config.guess \
	config.sub \
	cuda.common \
	depcomp \
	generate_user_env.awk \
	install-sh \
	lalapps.m4 \
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
	test-driver \
	py-compile \
	vcsID.common

CLEANFILES = \
	lalapps-user-env.shell \
	lalapps-user-env.sh \
	lalapps-user-env.csh

sysconf_DATA = \
	lalapps-user-env.sh \
	lalapps-user-env.csh

BUILT_SOURCES = \
	lalapps-user-env.sh \
	lalapps-user-env.csh

if HAVE_PYTHON
lalapps_pythonpath = $(pyexecdir) $(pythondir)
else
lalapps_pythonpath =
endif

lalapps-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PYTHONPATH $(lalapps_pythonpath)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

lalapps-user-env.csh lalapps-user-env.sh: lalapps-user-env.shell $(srcdir)/generate_user_env.awk
	@echo "creating $@"
	@$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }
