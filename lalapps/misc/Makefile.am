## Process this file with automake to produce Makefile.in
sysconf_DATA = lalapps-user-env.sh lalapps-user-env.csh
BUILT_SOURCES = lal-user-env.sh lal-user-env.csh
CLEANFILES = lal-user-env.sh lal-user-env.csh

lalapps-user-env.sh:
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "PATH=\"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "MANPATH=\"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "export PATH MANPATH" >>$@

lalapps-user-env.csh:
	@echo "creating $@" ;\
	echo "# Source this file to access LALApps" >$@ ;\
	echo "setenv PATH \"$(bindir):\$${PATH}\"" >>$@ ;\
	echo "if ( \$$?MANPATH ) then" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir):\$${MANPATH}\"" >>$@ ;\
	echo "else" >>$@ ;\
	echo "  setenv MANPATH \"$(mandir)\"" >>$@ ;\
	echo "endif" >>$@
