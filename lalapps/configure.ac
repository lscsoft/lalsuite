dnl Process this file with autoconf to produce a configure script.
AC_INIT([lalapps],[6.2],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/config.h.in])
AC_CONFIG_AUX_DIR([misc])
AC_CONFIG_FILES([\
  Makefile \
  lalapps.spec \
  debian/Makefile \
  doc/Makefile \
  man/Makefile \
  misc/Makefile \
  src/Makefile \
  src/lalapps/Makefile \
  src/calibration/Makefile \
  src/frametools/Makefile \
  src/findchirp/Makefile \
  src/stochastic/Makefile \
  src/power/Makefile \
  src/ring/Makefile \
  src/string/Makefile \
  src/pulsar/Makefile \
  src/pulsar/FDS_isolated/Makefile \
  src/pulsar/FDS_isolated/ResamplingFstat/Makefile \
  src/pulsar/Injections/Makefile \
  src/pulsar/FDS_griphyn/Makefile \
  src/pulsar/hough/Makefile \
  src/pulsar/hough/src/Makefile \
  src/pulsar/FDS_binary/Makefile \
  src/pulsar/SideBand/Makefile \
  src/pulsar/templateBanks/Makefile \
  src/pulsar/TimingTests/Makefile \
  src/pulsar/TDS_isolated/Makefile \
  src/inspiral/Makefile \
  src/inspiral/posterior/Makefile \
  src/zmsearch/Makefile \
  src/detresponse/Makefile \
  src/pulsar/hough/src2/Makefile \
  src/pulsar/MakeSFTs/Makefile \
  src/pulsar/crosscorr/Makefile \
  src/pulsar/fscan/Makefile \
  src/tracksearch/Makefile \
  src/online/Makefile \
])
AM_INIT_AUTOMAKE([])
AC_ARG_VAR(LALSUITE_BUILD,[Set if part of lalsuite build])
AH_TOP([
#ifndef CONFIG_H
#define CONFIG_H])
AH_BOTTOM([
#endif /* CONFIG_H */
])

AC_CANONICAL_HOST

LALAPPS_WITH_CC
LALAPPS_WITH_LAL_PREFIX
LALAPPS_WITH_EXTRA_CPPFLAGS
LALAPPS_WITH_EXTRA_CFLAGS
LALAPPS_WITH_EXTRA_LDFLAGS
LALAPPS_WITH_EXTRA_LIBS
LALAPPS_WITH_GCC_FLAGS
LALAPPS_ENABLE_CONDOR
LALAPPS_ENABLE_FRAME
LALAPPS_ENABLE_METAIO
LALAPPS_ENABLE_XML
LALAPPS_ENABLE_CFITSIO
LALAPPS_ENABLE_NIGHTLY

AC_PREFIX_DEFAULT([/opt/lscsoft/lalapps])

dnl Put version/configure info in config.h header
lalapps_version_major=`echo "$VERSION" | cut -d. -f1`
lalapps_version_minor=`echo "$VERSION" | cut -d. -f2`
lalapps_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
lalapps_cvs_tag=`echo '$Name$' | sed 's/\\$//g' | sed 's/Name: //'`
AC_DEFINE_UNQUOTED(LALAPPS_VERSION, "$VERSION", [LALApps Version])
AC_DEFINE_UNQUOTED(LALAPPS_VERSION_MAJOR, $lalapps_version_major,
                   [LALApps Version Major Number])
AC_DEFINE_UNQUOTED(LALAPPS_VERSION_MINOR, $lalapps_version_minor,
                   [LALApps Version Minor Number])
AC_DEFINE_UNQUOTED(LALAPPS_CONFIGURE_ARGS, "$ac_configure_args",
                   [LALApps Configure Arguments])
AC_DEFINE_UNQUOTED(LALAPPS_CONFIGURE_DATE, "$lalapps_configure_date",
                   [LALApps Configure Date])
AC_DEFINE_UNQUOTED(LALAPPS_CVS_TAG, "$lalapps_cvs_tag", [LALApps CVS Tag])


dnl Remove -all-static from ldflags
allstatic=`echo $LDFLAGS | sed '/-all-static/!d'`
if test "x$allstatic" = "x" ; then
  allstatic="false"
else
  LDFLAGS=`echo $LDFLAGS | sed s/-all-static//g`
  allstatic="true"
fi

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99], [AC_PROG_CC_C99], [LAL_AC_PROG_CC_C99])

# use silent build rules if appropriate
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

dnl Checks for programs.
AC_PROG_AWK
AC_PROG_CPP

AC_ENABLE_SHARED
if test "x$condor" = "xtrue"; then
  AC_CHECK_PROGS(CONDOR_COMPILE, condor_compile)
  CC="$CONDOR_COMPILE $CC"
  case "$host" in
    *-*-linux* | *-*-osf* | *-*-hpux* ) CFLAGS="$CFLAGS -static" AC_DISABLE_SHARED;;
  esac
  AC_DEFINE(LALAPPS_CONDOR, 1, [LALApps is condor compiled])
fi
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS -W -Wall $lalapps_gcc_flags"
fi
AC_PROG_INSTALL
AC_PROG_LN_S
AM_PROG_LIBTOOL
AC_CHECK_PROGS(LATEX, pdflatex latex, echo)
AC_CHECK_PROGS(MKIND, makeindex, echo)
AC_CHECK_PROGS(DVIPS, dvips, echo)
AC_CHECK_PROGS(BIBTEX, bibtex, echo)
AM_PATH_PYTHON(2.4,,:)
AM_CONDITIONAL([HAVE_PYTHON],[test "$PYTHON" != :])

dnl Checks for libraries.
AC_CHECK_LIB(m, sin)
AC_CHECK_LIB(z, compress)

PKG_CHECK_MODULES(LAL, lal, [ true ], [ false ])
PKG_CHECK_MODULES(LALSUPPORT, lalsupport , [ true ], [ false ])
if test "$LALSUITE_BUILD" != "true" ; then
	CPPFLAGS="$CPPFLAGS $LAL_CFLAGS $LALSUPPORT_CFLAGS"
	LIBS="$LALSUPPORT_LIBS $LAL_LIBS $LIBS"
	AC_CHECK_LIB(lal, LALVersion, , AC_MSG_ERROR(could not find the LAL library))
	AC_CHECK_LIB(lalsupport, LALOpenDataFile, , AC_MSG_ERROR(could not find the LAL support library))
fi

dnl LALFrame, LALMetaIo, and LALXML (optional).
lalframe="true"
lalmetaio="true"
lalxml="true"
PKG_CHECK_MODULES(LALFRAME, lalframe, [ true ], [ false ])
PKG_CHECK_MODULES(LALMETAIO, lalmetaio, [ true ], [ false ])
PKG_CHECK_MODULES(LALXML, lalxml, [ true ], [ false ])
if test "$LALSUITE_BUILD" != "true" ; then
	CPPFLAGS="$CPPFLAGS $LALFRAME_CFLAGS $LALMETAIO_CFLAGS $LALXML_CFLAGS"
	LIBS="$LALFRAME_LIBS $LALMETAIO_LIBS $LALXML_LIBS $LIBS"
fi

dnl cfitsio
PKG_CHECK_MODULES(CFITSIO, cfitsio, [ true ], [ false ])
CPPFLAGS="$CPPFLAGS $CFITSIO_CFLAGS"
LIBS="$CFITSIO_LIBS $LIBS"

dnl GSL
PKG_CHECK_MODULES(GSL, gsl, [ true ], [ false ])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$GSL_LIBS $LIBS"
extra_CPPFLAGS="$extra_CPPFLAGS $GSL_CFLAGS"
extra_LIBS="$GSL_LIBS $extra_LIBS"
AC_CHECK_LIB(gslcblas, main)
AC_CHECK_LIB(gsl, gsl_version)

dnl FrameL (libframe)
if test "${frame}" = "true"; then
  PKG_CHECK_MODULES(FRAME, libframe, [ true ], [ false ])
  lal_pre_frame_LIBS="$LIBS"
  LIBS="$LIBS $FRAME_LIBS"
  AC_SEARCH_LIBS(FrLibIni, Frame, ,
                 [LALAPPS_DISABLE_FRAME]
                 [frame="false"]
                 [LIBS="$lal_pre_frame_LIBS"])
fi

dnl MetaIo (libmetaio)
if test "${metaio}" = "true"; then
  PKG_CHECK_MODULES(METAIO, libmetaio, [ true ], [ false ])
  lal_pre_metaio_LIBS="$LIBS"
  LIBS="$LIBS $METAIO_LIBS"
  AC_SEARCH_LIBS(MetaioOpen, metaio dataflow, ,
                 [AC_MSG_WARN(disabling metaio support)]
                 [metaio="false"]
                 [LIBS="$lal_pre_metaio_LIBS"])
fi

dnl XML2 (libxml2)
if test "${xml}" = "true"; then
  PKG_CHECK_MODULES(LIBXML2, libxml-2.0 >= 2.6, [ true ], [ false ])
  lal_pre_xml2_LIBS="$LIBS"
  LIBS="$LIBS $LIBXML2_LIBS"
  AC_SEARCH_LIBS(xmlInitParser, xml2, ,
                 [AC_MSG_WARN(disabling xml support)]
                 [xml="false"]
                 [LIBS="$lal_pre_xml2_LIBS"])
fi

dnl Check for libcfitsio
if test "${cfitsio}" = "true"; then
  AC_CHECK_LIB(cfitsio, ffopen, , cfitsio="false")
fi

dnl Now check for LALFrame and LALMetaIo
if test "$LALSUITE_BUILD" != "true" ; then
  AC_CHECK_LIB(lalframe, LALFrOpen, , lalframe="false")
  AC_CHECK_LIB(lalmetaio, LALSnglInspiralTableFromLIGOLw, , lalmetaio="false")
  AC_CHECK_LIB(lalxml, XLALXMLFilePrintElements, , lalxml="false")
fi

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(getopt.h)
AC_CHECK_HEADERS([gsl/gsl_fft_real.h])
if test "$LALSUITE_BUILD" != "true" ; then
  AC_CHECK_HEADERS([lal/LALStdlib.h])
fi
AC_CHECK_HEADERS(glob.h)
AC_CHECK_HEADERS(unistd.h)

if test "${frame}" = "true"; then
  lal_pre_frame_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $FRAME_CFLAGS"
  AC_CHECK_HEADERS(FrameL.h, ,
                   [LALAPPS_DISABLE_FRAME]
                   [frame="false"]
                   [CPPFLAGS="$lal_pre_frame_CPPFLAGS"])
fi

if test "${metaio}" = "true"; then
  lal_pre_metaio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
  AC_CHECK_HEADERS(metaio.h, ,
                   [AC_MSG_WARN(disabling metaio support)]
                   [metaio="false"]
                   [CPPFLAGS="$lal_pre_metaio_CPPFLAGS"])
fi

if test "${xml}" = "true"; then
  lal_pre_xml2_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
  AC_CHECK_HEADERS(libxml/xmlversion.h, ,
                   [AC_MSG_WARN(disabling xml support)]
                   [xml="false"]
                   [CPPFLAGS="$lal_pre_xml2_CPPFLAGS"])
fi

if test "${cfitsio}" = "true"; then
  lal_pre_cfitsio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $CFITSIO_CFLAGS"
  AC_CHECK_HEADER(fitsio.h, ,
                  [cfitsio="false"]
                  [CPPFLAGS="$lal_pre_cfitsio_CPPFLAGS"])
fi

LALAPPS_CHECK_QTHREAD

AC_MSG_CHECKING(for gethostname prototype in unistd.h)
AC_EGREP_HEADER(gethostname,unistd.h,AC_MSG_RESULT(yes)
AC_DEFINE(HAVE_GETHOSTNAME_PROTOTYPE,1,
[Define if gethostname prototype is in unistd.h]),
AC_MSG_RESULT(no))

AC_CHECK_FUNCS(setenv)

if test "$LALSUITE_BUILD" = "true"; then
  LIBS="$LALXML_LIBS $LALFRAME_LIBS $LALMETAIO_LIBS $LALSUPPORT_LIBS $LAL_LIBS $LIBS"
  CPPFLAGS="$CPPFLAGS $LAL_CFLAGS $LALSUPPORT_CFLAGS $LALFRAME_CFLAGS $LALMETAIO_CFLAGS $LALXML_LIBS"
  lalframe="false"
  lalmetaio="false"
  lalxml="false"
  test -n "$LALFRAME_LIBS" && lalframe="true"
  test -n "$LALMETAIO_LIBS" && lalmetaio="true"
  test -n "$LALXML_LIBS" && lalxml="true"
fi

AM_CONDITIONAL(FRAME, test x$frame = xtrue)
AM_CONDITIONAL(LALFRAME, test x$lalframe = xtrue)
AM_CONDITIONAL(LALMETAIO, test x$lalmetaio = xtrue)
AM_CONDITIONAL(LALXML, test x$lalxml = xtrue)
AM_CONDITIONAL(DATAFLOW, test x$metaio = xtrue)
AM_CONDITIONAL(CFITSIO, test x$cfitsio = xtrue)
AC_DEFINE_UNQUOTED(PREFIX, "$prefix", [Install prefix])

dnl restore all static to ldflags
if test "${allstatic}" = "true" ; then
  LDFLAGS="-all-static $LDFLAGS"
fi

AC_OUTPUT

lalframeenabled="`test x${lalframe} = xtrue && echo "ENABLED" || echo "DISABLED"`"
lalmetaioenabled="`test x${lalmetaio} = xtrue && echo "ENABLED" || echo "DISABLED"`"
lalxmlenabled="`test x${lalxml} = xtrue && echo "ENABLED" || echo "DISABLED"`"
cfitsioenabled="`test x${cfitsio} = xtrue && echo "ENABLED" || echo "DISABLED"`"
echo "
================================================================

        LALApps has now been successfully configured:

                LAL Frame library support is $lalframeenabled
                LAL MetaIo library support is $lalmetaioenabled
                LAL XML library support is $lalxmlenabled
                FITS library support is $cfitsioenabled

        Now run 'make' to build LALApps
        and run 'make install' to install LALApps
        LALApps will be installed under the directory:

                $prefix

================================================================="
