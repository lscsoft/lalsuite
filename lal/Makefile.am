## Process this file with automake to produce Makefile.in
SUBDIRS = include packages doc lib misc @subdirs@
MAINTAINERCLEANFILES = README.install lal.spec
EXTRA_DIST = README.bugs README.contrib README.install \
  lal.pc.in lalsupport.pc.in lalframe.pc.in lalmetaio.pc.in lalmpi.pc.in \
  lal-config.in acinclude.m4 aclocal.m4 \
  misc/doc.mk misc/dvi.mk misc/include.mk misc/main.tex.in \
  misc/lal.m4 misc/lallibtool.m4 misc/lalam15patch.m4 \
  misc/acconfig.h lal.spec.in lal.spec
BUILT_SOURCES = README.install
bin_SCRIPTS = lal-config

if FRAME
LALFRAMEPC = lalframe.pc
else
LALFRAMEPC =
endif
if METAIO
LALMETAIOPC = lalmetaio.pc
else
LALMETAIOPC =
endif
if MPI
LALMPIPC = lalmpi.pc
else
LALMPIPC =
endif
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lal.pc lalsupport.pc $(LALFRAMEPC) $(LALMETAIOPC) $(LALMPIPC)

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

all-local:
	@echo "";\
echo "=================================================================";\
echo "";\
echo "        LAL has now been successfully built.";\
echo "        Now run 'make install' to install LAL.";\
echo "";\
echo "================================================================="

install-exec-local:
	@echo "";\
echo "=================================================================";\
echo "";\
echo "        LAL has now been successfully installed.";\
echo "";\
echo "        Now you may want to setup your environment:";\
echo "";\
echo "        Bourne shell [bash] users: please add the";\
echo "        following lines to your .profile file";\
echo "";\
echo "            LAL_LOCATION=$(prefix)";\
echo "            export LAL_LOCATION";\
echo '            . $${LAL_LOCATION}/etc/lal-user-env.sh';\
echo "";\
echo "        C-shell [tcsh] users: please add the following";\
echo "        lines to your .login file";\
echo "";\
echo "            setenv LAL_LOCATION $(prefix)";\
echo '            source $${LAL_LOCATION}/etc/lal-user-env.csh';\
echo "";\
echo "================================================================="

