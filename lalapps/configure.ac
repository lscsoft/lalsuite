AC_INIT([lalapps],[6.3.2.1],[lal-discuss@gravity.phys.uwm.edu])
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
  src/lalapps/defaults.py \
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
  src/inspiral/posterior/SPINspiral/Makefile \
  src/zmsearch/Makefile \
  src/detresponse/Makefile \
  src/pulsar/hough/src2/Makefile \
  src/pulsar/GCT/Makefile \
  src/pulsar/MakeSFTs/Makefile \
  src/pulsar/crosscorr/Makefile \
  src/pulsar/fscan/Makefile \
  src/tracksearch/Makefile \
  src/online/Makefile \
])
AC_CONFIG_FILES([src/lalapps/git_version],[chmod +x src/lalapps/git_version])
AM_INIT_AUTOMAKE([foreign])
AC_ARG_VAR(LALSUITE_BUILD,[Set if part of lalsuite build])
AC_ARG_VAR(LAL_TOP_SRCDIR,[Set to top source directory of lal])
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
LALAPPS_ENABLE_GCC_FLAGS
LALAPPS_ENABLE_CONDOR
LALAPPS_ENABLE_FRAME
LALAPPS_ENABLE_METAIO
LALAPPS_ENABLE_CFITSIO
LALAPPS_ENABLE_LALFRAME
LALAPPS_ENABLE_LALMETAIO
LALAPPS_ENABLE_LALXML
LALAPPS_ENABLE_LALBURST
LALAPPS_ENABLE_LALSTOCHASTIC
LALAPPS_ENABLE_NIGHTLY

AC_PREFIX_DEFAULT([/opt/lscsoft/lalapps])

# put version/configure info in config.h header
lalapps_version_major=`echo "$VERSION" | cut -d. -f1`
lalapps_version_minor=`echo "$VERSION" | cut -d. -f2`
lalapps_version_micro=`echo "$VERSION" | cut -d. -f3`
lalapps_version_devel=`echo "$VERSION" | cut -d. -f4-`
lalapps_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALAPPS_VERSION],["$VERSION"],[LALApps Version])
AC_DEFINE_UNQUOTED([LALAPPS_VERSION_MAJOR],[$lalapps_version_major],
                   [LALApps Version Major Number])
AC_DEFINE_UNQUOTED([LALAPPS_VERSION_MINOR],[$lalapps_version_minor],
                   [LALApps Version Minor Number])
AC_DEFINE_UNQUOTED([LALAPPS_VERSION_MICRO],[$lalapps_version_micro],
                   [LALApps Version Micro Number])
AC_DEFINE_UNQUOTED([LALAPPS_VERSION_DEVEL],[$lalapps_version_devel],
                   [LALApps Version Devel Number])
AC_DEFINE_UNQUOTED([LALAPPS_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALApps Configure Arguments])
AC_DEFINE_UNQUOTED([LALAPPS_CONFIGURE_DATE],["$lalapps_configure_date"],
                   [LALApps Configure Date])


# remove -all-static from ldflags, breaks test
allstatic=`echo $LDFLAGS | sed '/-all-static/!d'`
if test "x$allstatic" = "x" ; then
  allstatic="false"
else
  LDFLAGS=`echo $LDFLAGS | sed s/-all-static//g`
  allstatic="true"
fi

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99],[AC_PROG_CC_C99],[LALSUITE_AC_PROG_CC_C99])

# use silent build rules if appropriate
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

# checks for programs
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S

# check for condor_compile
AC_ENABLE_SHARED
if test "x$condor" = "xtrue"; then
  AC_CHECK_PROGS([CONDOR_COMPILE],[condor_compile])
  CC="$CONDOR_COMPILE $CC"
  case "$host" in
    *-*-linux* | *-*-osf* | *-*-hpux* ) CFLAGS="$CFLAGS -static" AC_DISABLE_SHARED;;
  esac
  AC_DEFINE([LALAPPS_CONDOR],[1],[LALApps is condor compiled])
fi
CONDOR_ENABLE_VAL="`eval test x$condor = xtrue && echo "ENABLED" || echo "DISABLED"`"

# further program checks
AC_PROG_LIBTOOL
AC_CHECK_PROGS([LATEX],[pdflatex latex],[echo])
AC_CHECK_PROGS([MKIND],[makeindex],[echo])
AC_CHECK_PROGS([DVIPS],[dvips],[echo])
AC_CHECK_PROGS([BIBTEX],[bibtex],[echo])
AM_PATH_PYTHON([2.4],,[:])
AM_CONDITIONAL([HAVE_PYTHON],[test "$PYTHON" != :])

# check for supported mac os x version
if test "x$build_vendor" = "xapple"; then
  AC_CHECK_PROGS([SW_VERS],[sw_vers])
  if test "x$SW_VERS" != "x"; then
    AC_MSG_CHECKING([Mac OS X version])
    MACOSX_VERSION=`$SW_VERS -productVersion`
    AC_MSG_RESULT([$MACOSX_VERSION])
  fi
  case "$MACOSX_VERSION" in
    10.0*|10.1*|10.2*|10.3*)
      AC_MSG_ERROR([This version of Mac OS X is not supported])
      ;;
    10.4*|10.5*|10.6*)
      # supported version
      ;;
    *)
      AC_MSG_WARN([Unknown Mac OS X version])
      ;;
  esac
fi

# check for system libraries
AC_CHECK_LIB([m],[sin])
AC_CHECK_LIB([z],[compress])

# check for system headers
AC_HEADER_STDC
AC_CHECK_HEADERS([glob.h])
AC_CHECK_HEADERS([punistd.h])
AC_CHECK_HEADERS([getopt.h])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$GSL_LIBS $LIBS"
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for libFrame
if test "${frame}" = "true"; then
  PKG_CHECK_MODULES([FRAME],[libframe],[true],[false])
  lal_pre_frame_LIBS="$LIBS"
  LIBS="$LIBS $FRAME_LIBS"
  AC_SEARCH_LIBS([FrLibIni],[Frame],,
                 [AC_MSG_WARN([could not find the frame library])]
                 [frame="false"]
                 [lalframe="false"]
                 [LIBS="$lal_pre_frame_LIBS"])
fi

# check for libFrame headers
if test "${frame}" = "true"; then
  lal_pre_frame_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $FRAME_CFLAGS"
  AC_CHECK_HEADERS([FrameL.h],,
                   [AC_MSG_WARN([could not find the FrameL.h header])]
                   [frame="false"]
                   [lalframe="false"]
                   [CPPFLAGS="$lal_pre_frame_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([FRAME],[frame])

# check for libmetaio
if test "${metaio}" = "true"; then
  PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
  lal_pre_metaio_LIBS="$LIBS"
  LIBS="$LIBS $METAIO_LIBS"
  AC_SEARCH_LIBS([MetaioOpen],[metaio],,
                 [AC_MSG_WARN([could not find the metaio library])]
                 [metaio="false"]
                 [lalmetaio="false"]
                 [LIBS="$lal_pre_metaio_LIBS"])
fi

# check for libmetaio headers
if test "${metaio}" = "true"; then
  lal_pre_metaio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
  AC_CHECK_HEADERS([metaio.h],,
                   [AC_MSG_WARN([could not find the metaio.h header])]
                   [metaio="false"]
                   [lalmetaio="false"]
                   [CPPFLAGS="$lal_pre_metaio_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([METAIO],[metaio])

# check for libcfitsio
if test "${cfitsio}" = "true"; then
  PKG_CHECK_MODULES([CFITSIO],[cfitsio],[true],[false])
  lal_pre_cfitsio_libs="$LIBS"
  LIBS="$CFITSIO_LIBS $LIBS"
  AC_SEARCH_LIBS([ffopen],[cfitsio],,
                 [cfitsio="false"]
                 [LIBS="$lal_pre_cfitsio_LIBS"])
fi

# check for libcfitsio headers
if test "${cfitsio}" = "true"; then
  lal_pre_cfitsio_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $CFITSIO_CFLAGS"
  AC_CHECK_HEADERS([fitsio.h],,
                   [cfitsio="false"]
                   [CPPFLAGS="$lal_pre_cfitsio_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([CFITSIO],[cfitsio])

# check for lal libararies and headers
LALSUITE_CHECK_LIB([LAL],[lal],[6.3.2.1],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSUPPORT],[lalsupport],[6.3.2.1],[LALOpenDataFile],[lal/FileIO.h])
LALSUITE_CHECK_OPT_LIB([LALFRAME],[lalframe],[1.0],[LALFrOpen],[lal/LALFrameL.h])
LALSUITE_CHECK_OPT_LIB([LALMETAIO],[lalmetaio],[1.0],[LALSnglInspiralTableFromLIGOLw],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_OPT_LIB([LALXML],[lalxml],[1.0],[XLALXMLFilePrintElements],[lal/LALXML.h])
LALSUITE_CHECK_OPT_LIB([LALBURST],[lalburst],[1.0],[XLALEPSearch],[lal/EPSearch.h])
LALSUITE_CHECK_OPT_LIB([LALSTOCHASTIC],[lalstochastic],[1.0],[LALStochasticOptimalFilter],[lal/StochasticCrossCorrelation.h])

# check for qthread
LALAPPS_CHECK_QTHREAD

# check for gethostname prototype
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
  [AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
             [Define if gethostname prototype is in unistd.h])],
  [AC_MSG_RESULT([no])])

AC_CHECK_FUNCS([setenv])

# set prefix
AC_DEFINE_UNQUOTED([PREFIX],["$prefix"],[Install prefix])

# add gcc specific flags
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS $lalapps_gcc_flags"

  # add mac os x specific flags
  if test "x$MACOSX_VERSION" != "x"; then
    CFLAGS="$CFLAGS -mmacosx-version-min=10.4"
  fi
fi

# restore all static to ldflags
if test "${allstatic}" = "true" ; then
  LDFLAGS="-all-static $LDFLAGS"
fi

AC_OUTPUT

echo "
================================================================

        LALApps has now been successfully configured:

                Frame library support is $FRAME_ENABLE_VAL
                MetaIO library support is $METAIO_ENABLE_VAL
                LAL Frame library support is $LALFRAME_ENABLE_VAL
                LAL MetaIo library support is $LALMETAIO_ENABLE_VAL
                LAL XML library support is $LALXML_ENABLE_VAL
                LAL Burst library support is $LALBURST_ENABLE_VAL
                LAL Stochastic library support is $LALSTOCHASTIC_ENABLE_VAL
                FITS library support is $CFITSIO_ENABLE_VAL
                Condor support is $CONDOR_ENABLE_VAL

        Now run 'make' to build LALApps
        and run 'make install' to install LALApps
        LALApps will be installed under the directory:

                $prefix

================================================================="
