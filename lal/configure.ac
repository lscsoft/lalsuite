dnl Process this file with autoconf to produce a configure script.
AC_INIT([lal],[6.2],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_HEADERS([include/config.h] [include/lal/LALConfig.h])
AC_CONFIG_SRCDIR([include/config.h.in])
AC_CONFIG_AUX_DIR([misc])
AC_CONFIG_FILES([ \
                 lal.spec \
                 lal.pc \
                 lalframe.pc \
                 lalmetaio.pc \
                 lalxml.pc \
                 lalstochastic.pc \
                 lalsupport.pc \
                 lal-config \
                 Makefile \
                 debian/Makefile \
                 doc/Makefile \
                 doc/laldoc/Makefile \
                 doc/autodoc/Makefile \
                 include/Makefile \
                 include/lal/Makefile \
                 lib/git_version
                 lib/Makefile \
                 misc/Makefile \
                 packages/Makefile \
                 packages/std/Makefile \
                 packages/std/doc/Makefile \
                 packages/std/include/Makefile
                 packages/std/src/Makefile \
                 packages/std/test/Makefile \
                 packages/tools/Makefile \
                 packages/tools/doc/Makefile \
                 packages/tools/include/Makefile \
                 packages/tools/src/Makefile \
                 packages/tools/test/Makefile \
                 packages/factories/Makefile \
                 packages/factories/doc/Makefile \
                 packages/factories/include/Makefile \
                 packages/factories/src/Makefile \
                 packages/factories/test/Makefile \
                 packages/vectorops/Makefile \
                 packages/vectorops/doc/Makefile \
                 packages/vectorops/include/Makefile \
                 packages/vectorops/src/Makefile \
                 packages/vectorops/test/Makefile \
                 packages/utilities/Makefile \
                 packages/utilities/doc/Makefile \
                 packages/utilities/include/Makefile \
                 packages/utilities/src/Makefile \
                 packages/utilities/test/Makefile \
                 packages/stats/Makefile \
                 packages/stats/doc/Makefile \
                 packages/stats/include/Makefile \
                 packages/stats/src/Makefile \
                 packages/stats/test/Makefile \
                 packages/inject/Makefile \
                 packages/inject/doc/Makefile \
                 packages/inject/include/Makefile \
                 packages/inject/src/Makefile \
                 packages/inject/test/Makefile \
                 packages/date/Makefile \
                 packages/date/doc/Makefile \
                 packages/date/include/Makefile \
                 packages/date/src/Makefile \
                 packages/date/test/Makefile \
                 packages/tdfilter/Makefile \
                 packages/tdfilter/doc/Makefile \
                 packages/tdfilter/include/Makefile \
                 packages/tdfilter/src/Makefile \
                 packages/tdfilter/test/Makefile \
                 packages/window/Makefile \
                 packages/window/doc/Makefile \
                 packages/window/include/Makefile \
                 packages/window/src/Makefile \
                 packages/window/test/Makefile \
                 packages/fft/Makefile \
                 packages/fft/doc/Makefile \
                 packages/fft/include/Makefile \
                 packages/fft/src/Makefile \
                 packages/fft/test/Makefile \
                 packages/clremoval/Makefile \
                 packages/clremoval/doc/Makefile \
                 packages/clremoval/include/Makefile \
                 packages/clremoval/src/Makefile \
                 packages/clremoval/test/Makefile \
                 packages/inspiral/Makefile \
                 packages/inspiral/doc/Makefile \
                 packages/inspiral/include/Makefile \
                 packages/inspiral/src/Makefile \
                 packages/inspiral/test/Makefile \
                 packages/noisemodels/Makefile \
                 packages/noisemodels/doc/Makefile \
                 packages/noisemodels/include/Makefile \
                 packages/noisemodels/src/Makefile \
                 packages/noisemodels/test/Makefile \
                 packages/bank/Makefile \
                 packages/bank/doc/Makefile \
                 packages/bank/include/Makefile \
                 packages/bank/src/Makefile \
                 packages/bank/test/Makefile \
                 packages/ring/Makefile \
                 packages/ring/doc/Makefile \
                 packages/ring/include/Makefile \
                 packages/ring/src/Makefile \
                 packages/ring/test/Makefile \
                 packages/findchirp/Makefile \
                 packages/findchirp/doc/Makefile \
                 packages/findchirp/include/Makefile \
                 packages/findchirp/src/Makefile \
                 packages/findchirp/test/Makefile \
                 packages/pulsar/Makefile \
                 packages/pulsar/doc/Makefile \
                 packages/pulsar/include/Makefile \
                 packages/pulsar/src/Makefile \
                 packages/pulsar/test/Makefile \
                 packages/houghpulsar/Makefile \
                 packages/houghpulsar/doc/Makefile \
                 packages/houghpulsar/include/Makefile \
                 packages/houghpulsar/src/Makefile \
                 packages/houghpulsar/test/Makefile \
                 packages/burstsearch/Makefile \
                 packages/burstsearch/doc/Makefile \
                 packages/burstsearch/include/Makefile \
                 packages/burstsearch/src/Makefile \
                 packages/burstsearch/test/Makefile \
                 packages/tracksearch/Makefile \
                 packages/tracksearch/doc/Makefile \
                 packages/tracksearch/include/Makefile \
                 packages/tracksearch/src/Makefile \
                 packages/tracksearch/test/Makefile \
                 packages/timefreq/Makefile \
                 packages/timefreq/doc/Makefile \
                 packages/timefreq/include/Makefile \
                 packages/timefreq/src/Makefile \
                 packages/timefreq/test/Makefile \
                 packages/stochastic/Makefile \
                 packages/stochastic/doc/Makefile \
                 packages/stochastic/include/Makefile \
                 packages/stochastic/src/Makefile \
                 packages/stochastic/test/Makefile \
                 packages/support/Makefile \
                 packages/support/doc/Makefile \
                 packages/support/include/Makefile \
                 packages/support/src/Makefile \
                 packages/support/test/Makefile \
                 packages/framedata/Makefile \
                 packages/framedata/doc/Makefile \
                 packages/framedata/include/Makefile \
                 packages/framedata/src/Makefile \
                 packages/framedata/test/Makefile
])
AM_INIT_AUTOMAKE([])
AH_TOP([
#ifndef CONFIG_H__LAL
#define CONFIG_H__LAL
])
AH_TEMPLATE([LAL_SIZEOF_DOUBLE], [The number of bytes in a double.])
AH_TEMPLATE([LAL_SIZEOF_FLOAT], [The number of bytes in a float.])
AH_BOTTOM([
#endif /* CONFIG_H__LAL */
])

AC_CANONICAL_HOST

LAL_WITH_CC
LAL_WITH_CFLAGS
LAL_WITH_EXTRA_CFLAGS
LAL_WITH_EXTRA_CPPFLAGS
LAL_WITH_EXTRA_LDFLAGS
LAL_WITH_EXTRA_LIBS
LAL_WITH_GCC_FLAGS

LAL_ENABLE_INTELFFT
LAL_ENABLE_FRAME
LAL_ENABLE_METAIO
LAL_ENABLE_XML

LAL_ENABLE_DEBUG
LAL_ENABLE_MACROS
LAL_ENABLE_PTHREAD_LOCK
LAL_ENABLE_NIGHTLY

AC_PREFIX_DEFAULT([/opt/lscsoft/lal])

lal_version_major=`echo "$VERSION" | cut -d. -f1`
lal_version_minor=`echo "$VERSION" | cut -d. -f2`
lal_version_micro=`echo "$VERSION" | cut -d. -f3-`
test -z "$lal_version_micro" && lal_version_micro=0
lal_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
# This sed script is more complex than it needs to be in order to prevent CVS
# from corrupting it.
lal_cvs_tag=`echo '$Name$' | sed -e 's? \\$\$??' -e 's?^\\$Name: ??'`
AC_DEFINE_UNQUOTED(LAL_VERSION, "$VERSION", [LAL Version])
AC_DEFINE_UNQUOTED(LAL_VERSION_MAJOR, $lal_version_major, [LAL Version Major Number] )
AC_DEFINE_UNQUOTED(LAL_VERSION_MINOR, $lal_version_minor, [LAL Version Minor Number] )
AC_DEFINE_UNQUOTED(LAL_VERSION_MICRO, $lal_version_micro, [LAL Version Micro Number] )
AC_DEFINE_UNQUOTED(LAL_CONFIGURE_ARGS, "$ac_configure_args", [LAL Configure Arguments])
AC_DEFINE_UNQUOTED(LAL_CONFIGURE_DATE, "$lal_configure_date", [LAL Configure Date])
AC_DEFINE_UNQUOTED(LAL_CVS_TAG, "$lal_cvs_tag", [LAL CVS Tag])
test "$prefix"="NONE" && lal_prefix="$ac_default_prefix" || lal_prefix="$prefix"
AC_DEFINE_UNQUOTED(LAL_PREFIX, "$lal_prefix", [LAL Install Prefix])

lal_pthread_available="yes"
ACX_PTHREAD( , [lal_pthread_available="no"])
if test "$lal_pthread_lock" = "true" ; then
  if test "$lal_pthread_available" = "yes" ; then
    LIBS="$PTHREAD_LIBS $LIBS"
    CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
    CC="$PTHREAD_CC"
  else
    AC_MSG_ERROR([do not know how to compile posix threads])
  fi
fi

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99], [AC_PROG_CC_C99], [LAL_AC_PROG_CC_C99])

# use silent build rules if appropriate
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

dnl Checks for programs.
AC_PROG_AWK
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_CC_C_O
AM_PROG_LIBTOOL
AM_PATH_PYTHON([2.4],,[AC_MSG_ERROR([Python-2.4, or higher, is required])])
AC_CHECK_PROGS(M4, gm4 m4, m4)
AC_CHECK_PROGS(LATEX, pdflatex latex, echo)
AC_CHECK_PROGS(MKIND, makeindex, echo)
AC_CHECK_PROGS(DVIPS, dvips, echo)
if test -z "$TEXFLAGS"
then
AC_MSG_CHECKING(if $LATEX accepts --interaction=batchmode)
cat > conflatextest.tex <<\CEOF
\documentclass{article}
\begin{document}
\end{document}
CEOF
if $LATEX conflatextest.tex --interaction=batchmode >/dev/null 2>&1
then
TEXFLAGS="--interaction=batchmode"
AC_MSG_RESULT(yes)
else
TEXFLAGS=""
AC_MSG_RESULT(no)
fi
rm -f conflatextest.*
fi
AC_SUBST(TEXFLAGS)

dnl Set needed compiler flags
AC_MSG_CHECKING([for additional compiler flags])
case "${host_os}" in
  solaris*) CFLAGS="$CFLAGS -D_POSIX_C_SOURCE=200112L"
  AC_MSG_RESULT([-D_POSIX_C_SOURCE=200112L]);;
  *) AC_MSG_RESULT(none);;
esac

dnl Checks for libraries.
LALSUPPORT_LIBADD=""
lal_save_LIBS="$LIBS"
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(z, compress, [ LALSUPPORT_LIBADD="$LALSUPPORT_LIBADD -lz" ])
dnl if solaris, check for sunmath library
case "${host_os}" in
  solaris*) AC_CHECK_LIB(sunmath, sincosp);;
esac
PKG_CHECK_MODULES(GSL, gsl, [ true ], [ false ])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$LIBS $GSL_LIBS"
AC_CHECK_LIB(gslcblas, main)
AC_CHECK_LIB(gsl, gsl_strerror, , [AC_MSG_ERROR(could not find GSL library)])
LAL_CHECK_GSL_VERSION(1.9)

dnl Checks for fft libraries.
if test "${intelfft}" = "false" ; then
  fftw3="true"
  PKG_CHECK_MODULES(FFTW3, fftw3 fftw3f, [ true ], [ false ])
  CPPFLAGS="$CPPFLAGS $FFTW3_CFLAGS"
  LIBS="$LIBS $FFTW3_LIBS"
  AC_CHECK_LIB(fftw3f, fftwf_execute_dft, , [AC_MSG_ERROR(could not find fftw3f library)], -lm)
  AC_CHECK_LIB(fftw3, fftw_execute_dft, , [AC_MSG_ERROR(could not find fftw3 library)], -lm)
else
  AC_MSG_WARN(Using Intel FFT routines)
  if test "x${qthread}" = "xtrue" ; then
    LAL_INTEL_MKL_QTHREAD_WARNING
  else
    AC_CHECK_LIB(pthread, pthread_create, , [AC_MSG_ERROR(could not find pthread library)])
  fi
  AC_CHECK_LIB(guide, __kmp_abort, , [AC_MSG_ERROR(could not find guide library included with Intel MKL)])
  AC_CHECK_HEADERS(pthread.h, , [AC_MSG_ERROR(could not find pthread header)])
  if test "${enable_static}" = "yes"; then
    if test "${enable_shared}" = "yes" ; then
      LAL_INTEL_FFT_LIBS_MSG_ERROR
    else
      AC_CHECK_LIB(mkl_ia32, DftiCreateDescriptor, , [AC_MSG_ERROR(could not find the Intel FFT library)])
    fi
  else
    AC_CHECK_LIB(mkl, DftiCreateDescriptor, , [AC_MSG_ERROR(could not find the Intel FFT library)])
  fi
  AC_CHECK_HEADERS(mkl_fft.h, , [AC_MSG_ERROR(could not find Intel FFT header)])
fi

dnl Checks for optional FrameL, metaio, and xml2 libraries.
if test "${metaio}" = "true"; then
  PKG_CHECK_MODULES(METAIO, libmetaio, [ true ], [ false ])
  lal_pre_metaio_LIBS="$LIBS"
  LIBS="$LIBS $METAIO_LIBS"
  AC_SEARCH_LIBS(MetaioOpen, metaio dataflow, ,
                 [AC_MSG_WARN(disabling metaio support)]
                 [metaio="false"]
                 [LIBS="$lal_pre_metaio_LIBS"])
fi
if test "${frame}" = "true"; then
  PKG_CHECK_MODULES(FRAME, libframe, [ true ], [ false ])
  lal_pre_frame_LIBS="$LIBS"
  LIBS="$LIBS $FRAME_LIBS"
  AC_SEARCH_LIBS(FrLibIni, Frame, ,
                 [AC_MSG_WARN(disabling frame library)]
                 [frame="false"]
                 [LIBS="$lal_pre_frame_LIBS"])
fi
if test "${xml}" = "true"; then
  PKG_CHECK_MODULES(LIBXML2, libxml-2.0 >= 2.6, [ true ], [ false ])
  lal_pre_xml2_LIBS="$LIBS"
  LIBS="$LIBS $LIBXML2_LIBS"
  AC_SEARCH_LIBS(xmlInitParser, xml2, ,
                 [AC_MSG_WARN(disabling xml support)]
                 [xml="false"]
                 [LIBS="$lal_pre_xml2_LIBS"])
fi


dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(sys/time.h unistd.h getopt.h malloc.h zlib.h regex.h glob.h)
case "${host_os}" in
  solaris*) AC_CHECK_HEADERS(sunmath.h);;
esac
AC_CHECK_HEADERS(gsl/gsl_errno.h, , [AC_MSG_ERROR(could not find gsl/gsl_errno.h)])
if test "${intelfft}" = "false" ; then
  AC_CHECK_HEADERS(fftw3.h, , [AC_MSG_ERROR(could not find fftw3.h)])
fi
AC_CHECK_HEADER(stdint.h, , [AC_MSG_ERROR(could not find stdint.h)])
AC_CHECK_HEADER(inttypes.h, , [AC_MSG_ERROR(could not find inttypes.h)])

AC_MSG_CHECKING(for gethostname prototype in unistd.h)
AC_EGREP_HEADER(gethostname, unistd.h,AC_MSG_RESULT(yes)
AC_DEFINE(HAVE_GETHOSTNAME_PROTOTYPE, 1, [Define if gethostname prototype is in unistd.h]), AC_MSG_RESULT(no))

if test "${metaio}" = "true"; then
  lal_pre_metaio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
  AC_CHECK_HEADERS(metaio.h, ,
                   [AC_MSG_WARN(disabling metaio support)]
                   [metaio="false"]
                   [CPPFLAGS="$lal_pre_metaio_CPPFLAGS"])
fi
if test "${frame}" = "true"; then
  lal_pre_frame_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $FRAME_CFLAGS"
  AC_CHECK_HEADERS(FrameL.h, ,
                   [AC_MSG_WARN(disabling frame library support)]
                   [frame="false"]
                   [CPPFLAGS="$lal_pre_frame_CPPFLAGS"])
fi
if test "${xml}" = "true"; then
  lal_pre_xml2_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
  AC_CHECK_HEADERS(libxml/xmlversion.h, ,
                   [AC_MSG_WARN(disabling xml support)]
                   [xml="false"]
                   [CPPFLAGS="$lal_pre_xml2_CPPFLAGS"])
fi
if test "${lal_pthread_lock}" = "true"; then
  AC_CHECK_HEADERS(pthread.h, break)
fi


dnl Define variables describing what was enabled.
if test "${fftw3}" = "true"; then
  AC_DEFINE(LAL_FFTW3_ENABLED, 1, [Define if using fftw3 library])
fi
if test "${metaio}" = "true"; then
  AC_DEFINE(LAL_METAIO_ENABLED, 1, [Define if using MetaIo library])
fi
if test "${frame}" = "true"; then
  AC_DEFINE(LAL_FRAME_ENABLED, 1, [Define if using frame library])
fi
if test "${xml}" = "true"; then
  AC_DEFINE(LAL_XML_ENABLED, 1, [Define if using xml2 library])
fi


dnl Record the libraries needed for linking LAL.
LALLIBADD="$LIBS"
LIBS="$lal_save_LIBS"
AC_SUBST(LALLIBADD)
LAL_LIBADD="$LIBS $LALLIBADD"
AC_SUBST(LAL_LIBADD)
AC_SUBST(LALSUPPORT_LIBADD)


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_STRINGIZE
AC_TYPE_SIZE_T
AC_HEADER_TIME
LAL_CHECK_SIZEOF(float, 4)
LAL_CHECK_SIZEOF(double, 8)

dnl Checks for library functions.
AC_CHECK_FUNCS(getopt)
AC_CHECK_FUNCS(getopt_long)
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(vsnprintf)
AC_CHECK_FUNCS(putenv)

dnl Set package documentation install path
pkgdocdir="\${prefix}/share/doc/${PACKAGE}-${VERSION}"
AC_SUBST(pkgdocdir)

dnl Setup automake conditionals
AM_CONDITIONAL(FRAME, test x$frame = xtrue)
AM_CONDITIONAL(METAIO, test x$metaio = xtrue)
AM_CONDITIONAL(XML, test x$xml = xtrue)
AM_CONDITIONAL(INTELFFT, test x$intelfft = xtrue)
AM_CONDITIONAL(QTHREAD, test x$qthread = xtrue)


dnl Configure subdirectory lalapps if present
if test -d $srcdir/lalapps ; then
  lal_build_dir='$(top_builddir)/..'
  LAL_BUILD="true"
  LAL_LIBS=$lal_build_dir/lib/liblal.la
  LAL_CFLAGS=-I$lal_build_dir/include
  LALSUPPORT_LIBS=$lal_build_dir/packages/support/src/liblalsupport.la
  LALSUPPORT_CFLAGS=-I$lal_build_dir/include
  if test x$frame = xtrue ; then
    LALFRAME_LIBS=$lal_build_dir/packages/framedata/src/liblalframe.la
    LALFRAME_CFLAGS=-I$lal_build_dir/include
  else
    LALFRAME_LIBS=""
    LALFRAME_CFLAGS=""
  fi
  if test x$metaio = xtrue ; then
    LALMETAIO_LIBS=$lal_build_dir/packages/support/src/liblalmetaio.la
    LALMETAIO_CFLAGS=-I$lal_build_dir/include
  else
    LALMETAIO_LIBS=""
    LALMETAIO_CFLAGS=""
  fi
  if test x$xml = xtrue ; then
    LALXML_LIBS=$lal_build_dir/packages/support/src/liblalxml.la
    LALXML_CFLAGS=-I$lal_build_dir/include
  else
    LALXML_LIBS=""
    LALXML_CFLAGS=""
  fi
  export LAL_BUILD LAL_LIBS LAL_CFLAGS LALSUPPORT_LIBS LALSUPPORT_CFLAGS LALFRAME_LIBS LALFRAME_CFLAGS LALMETAIO_LIBS LALMETAIO_CFLAGS LALXML_LIBS LALXML_CFLAGS
  AC_CONFIG_SUBDIRS(lalapps)
fi


dnl Add gcc specific flags --- must be done at end or tests will break!
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS $lal_gcc_flags"
fi
INCLUDEMKFILE=$srcdir/misc/include.mk
DOCMKFILE=$srcdir/misc/doc.mk
DVIMKFILE=$srcdir/misc/dvi.mk
AC_SUBST_FILE(INCLUDEMKFILE)
AC_SUBST_FILE(DOCMKFILE)
AC_SUBST_FILE(DVIMKFILE)
AC_CONFIG_COMMANDS([default-1],[[chmod +x lal-config]],[[]])
AC_CONFIG_COMMANDS([default-2],[[chmod +x lib/git_version]],[[]])
AC_OUTPUT
frameenabled="`test x${frame} = xtrue && echo "ENABLED" || echo "DISABLED"`"
metaioenabled="`test x${metaio} = xtrue && echo "ENABLED" || echo "DISABLED"`"
xmlenabled="`test x${xml} = xtrue && echo "ENABLED" || echo "DISABLED"`"
echo "
================================================================

        LAL has now been successfully configured:

                LAL Frame library support is $frameenabled
                LAL MetaIo library support is $metaioenabled
                LAL XML library support is $xmlenabled

        Now run 'make' to build LAL
        and run 'make install' to install LAL
        LAL will be installed under the directory:

                $prefix

================================================================="
