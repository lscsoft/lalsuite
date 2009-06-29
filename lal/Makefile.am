## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS = -I misc
SUBDIRS = include packages doc lib misc debian @subdirs@
MAINTAINERCLEANFILES = lal.spec
distcleancheck_listfiles = find . -type f \! -name "LALVCSInfo.h" -print
EXTRA_DIST = README.bugs README.contrib README.install \
  lal.pc.in lalsupport.pc.in lalframe.pc.in lalmetaio.pc.in lalxml.pc.in \
  lal-config.in aclocal.m4 \
  lal.spec.in lal.spec
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
if XML
LALXMLPC = lalxml.pc
else
LALXMLPC =
endif
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = lal.pc lalsupport.pc $(LALFRAMEPC) $(LALMETAIOPC) $(LALXMLPC)

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
	rm -f aclocal.m4
	rm -f acinclude.m4
	rm -f acconfig.h
	rm -f testscript.log
	rm -f config.cache
	rm -f config.status
	rm -f config.log
	rm -f config.h.in

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
echo '            . $(sysconfdir)/lal-user-env.sh';\
echo "";\
echo "        C-shell [tcsh] users: please add the following";\
echo "        lines to your .login file";\
echo "";\
echo '            source $(sysconfdir)/lal-user-env.csh';\
echo "";\
echo "================================================================="

