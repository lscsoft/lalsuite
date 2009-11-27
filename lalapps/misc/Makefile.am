EXTRA_DIST = \
	config.guess \
	config.sub \
	depcomp \
	install-sh \
	lalapps.m4 \
	lalapps_c99.m4 \
	lalappslibtool.m4 \
	lalsuite_build.m4 \
	ltmain.sh \
	missing \
	mkinstalldirs \
	pkg.m4 \
	py-compile

sysconf_DATA = lalapps-user-env.sh lalapps-user-env.csh
BUILT_SOURCES = lalapps-user-env.sh lalapps-user-env.csh
CLEANFILES = lalapps-user-env.sh lalapps-user-env.csh
MAINTAINERCLEANFILES = __lalappsGitID.h lalappsGitID.h

noinst_HEADERS = lalappsGitID.h
lalappsGitID.h: gitID

.PHONY: gitID
gitID:
	if test -f $(top_srcdir)/misc/generateGitID.sh; then \
	$(SHELL) $(top_srcdir)/misc/generateGitID.sh $(top_srcdir) lalapps; \
	fi


if HAVE_PYTHON
lalapps-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	if test "$(pyexecdir)" != "$(pythondir)" ; then \
		echo "PYTHONPATH=\"$(pyexecdir):$(pythondir):\$${PYTHONPATH}\"" >>$@ ;\
	else \
		echo "PYTHONPATH=\"$(pythondir):\$${PYTHONPATH}\"" >>$@ ;\
	fi ;\
	echo "export PATH MANPATH PYTHONPATH" >>$@
lalapps-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@ ;\
	echo "if ( \$$?PYTHONPATH ) then" >>$@ ;\
	if test "$(pyexecdir)" != "$(pythondir)" ; then \
		echo "  setenv PYTHONPATH \"$(pyexecdir):$(pythondir):\$${PYTHONPATH}\"" >>$@ ;\
	else \
		echo "  setenv PYTHONPATH \"$(pythondir):\$${PYTHONPATH}\"" >>$@ ;\
	fi ;\
	echo "else" >>$@ ;\
	if test "$(pyexecdir)" != "$(pythondir)" ; then \
		echo "  setenv PYTHONPATH \"$(pyexecdir):$(pythondir)\"" >>$@ ;\
	else \
		echo "  setenv PYTHONPATH \"$(pythondir)\"" >>$@ ;\
	fi ;\
	echo "endif" >>$@
else
lalapps-user-env.sh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "export PATH MANPATH" >>$@
lalapps-user-env.csh: Makefile
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@
endif
