AC_PREREQ([2.63])
AC_INIT([lalapps],[6.13.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/lalapps/lalapps.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([\
  Makefile \
  lalapps.spec \
  debian/Makefile \
  doc/Makefile \
  gnuscripts/Makefile \
  man/Makefile \
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
  src/pulsar/FDS_isolated/OptimizedCFS/Makefile \
  src/pulsar/FDS_isolated/ResamplingFstat/Makefile \
  src/pulsar/Injections/Makefile \
  src/pulsar/hough/Makefile \
  src/pulsar/hough/src/Makefile \
  src/pulsar/FDS_binary/Makefile \
  src/pulsar/SideBand/Makefile \
  src/pulsar/CombSearch/Makefile \
  src/pulsar/templateBanks/Makefile \
  src/pulsar/TimingTests/Makefile \
  src/pulsar/TDS_isolated/Makefile \
  src/inspiral/Makefile \
  src/inspiral/bayestar/Makefile \
  src/inspiral/posterior/Makefile \
  src/inspiral/posterior/SPINspiral/Makefile \
  src/inspiral/posterior/mpi/Makefile \
  src/pulsar/hough/src2/Makefile \
  src/pulsar/GCT/Makefile \
  src/pulsar/MakeSFTs/Makefile \
  src/pulsar/crosscorr/Makefile \
  src/pulsar/fscan/Makefile \
  src/pulsar/TwoSpect/Makefile \
  src/pulsar/SFTReferenceLibrary/Makefile \
  src/tracksearch/Makefile \
  src/online/Makefile \
  src/hwinjection/Makefile \
])
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests])
AH_TOP([
#ifndef CONFIG_H
#define CONFIG_H])
AH_BOTTOM([
#endif /* CONFIG_H */
])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# check for git repo
LALSUITE_CHECK_GIT_REPO

LALAPPS_WITH_CC
LALAPPS_WITH_LAL_PREFIX
LALAPPS_WITH_EXTRA_CPPFLAGS
LALAPPS_WITH_EXTRA_CFLAGS
LALAPPS_WITH_EXTRA_LDFLAGS
LALAPPS_WITH_EXTRA_LIBS

LALSUITE_ENABLE_DEBUG

LALSUITE_WITH_CUDA

LALSUITE_ENABLE_GCC_FLAGS
LALSUITE_ENABLE_NIGHTLY
LALAPPS_ENABLE_CONDOR
LALSUITE_ENABLE_BOINC
LALAPPS_ENABLE_STATIC_BINARIES

LALAPPS_ENABLE_MPI

LALAPPS_ENABLE_FFTW
LALAPPS_ENABLE_FRAME
LALAPPS_ENABLE_METAIO
LALAPPS_ENABLE_CFITSIO
LALAPPS_ENABLE_PSS
LALAPPS_ENABLE_GDS

LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALXML
LALSUITE_ENABLE_LALSIMULATION
LALSUITE_ENABLE_LALBURST
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALINFERENCE
LALSUITE_ENABLE_LALSTOCHASTIC

# put version/configure info in config header
lalapps_version_major=`echo "$VERSION" | cut -d. -f1`
lalapps_version_minor=`echo "$VERSION" | cut -d. -f2`
lalapps_version_micro=`echo "$VERSION" | cut -d. -f3`
lalapps_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalapps_version_micro" && lalapps_version_micro=0
test -z "$lalapps_version_devel" && lalapps_version_devel=0
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

# remove -all-static from ldflags, breaks tests
allstatic=`echo $LDFLAGS | sed '/-all-static/!d'`
if test "x$allstatic" = "x" ; then
  allstatic="false"
else
  LDFLAGS=`echo $LDFLAGS | sed s/-all-static//g`
  allstatic="true"
fi

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for C/C++ compilers
LALSUITE_PROG_CC_CXX

# check for MPI compiler
if test "x$mpi" = "xtrue"; then
  AX_MPI([mpi=true],[mpi=false])
fi
LALSUITE_ENABLE_MODULE([MPI],[mpi])

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED

# check for condor_compile
AC_ENABLE_SHARED
if test "x$condor" = "xtrue"; then
  AC_CHECK_PROGS([CONDOR_COMPILE],[condor_compile],[false])
  if test "x$CONDOR_COMPILE" = "xfalse"; then
    AC_MSG_ERROR([condor_compile not found])
  fi

  # does condor_compile suffer from condorbug 19343
  AC_MSG_CHECKING([if condor_compile suffers from bug 19343])
  if $CONDOR_COMPILE $CC -print-prog-name=ld 1>/dev/null 2>/dev/null; then
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    LD=/usr/bin/ld
  fi
fi
CONDOR_ENABLE_VAL="`eval test x$condor = xtrue && echo "ENABLED" || echo "DISABLED"`"
AM_CONDITIONAL([CONDOR_ENABLED],[test "$CONDOR_ENABLE_VAL" = "ENABLED"])

# further program checks
LT_INIT
AC_CHECK_PROGS([LATEX],[pdflatex latex],[echo])
AC_CHECK_PROGS([MKIND],[makeindex],[echo])
AC_CHECK_PROGS([DVIPS],[dvips],[echo])
AC_CHECK_PROGS([BIBTEX],[bibtex],[echo])
AM_PATH_PYTHON([2.6],,[:])
AM_CONDITIONAL([HAVE_PYTHON],[test "$PYTHON" != :])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for system libraries
AC_CHECK_LIB([m],[sin])
AC_CHECK_LIB([z],[compress],,[AC_MSG_ERROR([could not find the zlib library])])

# check for system headers
AC_HEADER_STDC
AC_CHECK_HEADERS([zlib.h],,[AC_MSG_ERROR([could not find zlib.h])])
AC_CHECK_HEADERS([glob.h])
AC_CHECK_HEADERS([punistd.h])
AC_CHECK_HEADERS([getopt.h])

# check for boinc
if test "${boinc}" = "true" ; then
  if test -n "$BOINC_PREFIX}" ; then
    # do compilation checks with c++
    AC_LANG_PUSH([C++])

    # has lal build built with boinc support
    LALSUITE_CHECK_BOINC

    # point to boinc headers and libraries
    CPPFLAGS="$CPPFLAGS -I$BOINC_PREFIX/include"
    LDFLAGS="$LDFLAGS -L$BOINC_PREFIX/lib"

    # check for boinc libraries
    AC_CHECK_LIB([boinc],[boinc_fopen],[true],[boinc=false])
    AC_CHECK_LIB([boinc_api],[boinc_finish],[true],[boinc=false],[-lboinc])

    # check for boinc headers
    AC_CHECK_HEADERS([boinc/boinc_api.h],,[boinc=false])

    # go back to c
    AC_LANG_POP([C++])
  fi
fi
BOINC_ENABLE_VAL="`eval test x$boinc = xtrue && echo "ENABLED" || echo "DISABLED"`"

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$GSL_LIBS $LIBS"
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fftw libraries
if test "${fftw}" = "true"; then
  PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
  CPPFLAGS="$CPPFLAGS $FFTW3_CFLAGS"
  LIBS="$LIBS $FFTW3_LIBS"
  AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_WARN([could not find the fftw3f library])],[-lm])
  AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_WARN([could not find the fftw3 library])],[-lm])
fi

# check for fftw headers
if test "${fftw}" = "true"; then
  AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_WARN([could not find the fftw3.h header])])
fi
LALSUITE_ENABLE_MODULE([FFTW],[fftw])

# check for GDS system: we check if we can link against libawg which is in GDS
# GDS is a requirements for building 'online_psinject', which is used for hardware injections
if test "${gds}" = "true"; then
  AC_CHECK_LIB([awg],[awgGetChannelNames],,[gds="false"])
  AC_CHECK_HEADERS([gdstype.h],,[gds="false"])
fi
LALSUITE_ENABLE_MODULE([GDS],[gds])

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

# check for libxml2
if test "$lalxml" = "true"; then
  PKG_CHECK_MODULES([LIBXML2],[libxml-2.0 >= 2.6],[true],[false])
  lal_pre_libxml2_LIBS="$LIBS"
  LIBS="$LIBS $LIBXML2_LIBS"
  AC_SEARCH_LIBS([xmlInitParser],[xml2],,
                 [AC_MSG_WARN([cannot find the xml2 library])]
                 [lalxml="false"]
                 [LIBS="$lal_pre_libxml2_LIBS"])
fi

# check for libxml2 headers
if test "$lalxml" = "true"; then
  lal_pre_libxml2_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
  AC_CHECK_HEADERS([libxml/tree.h],,
                   [AC_MSG_WARN([could not the libxml/tree.h header])]
                   [lalxml="false"]
                   [CPPFLAGS="$lal_pre_libxml2_CPPFLAGS"])
fi

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

# check for pss
if test "${pss}" = "true"; then
  PKG_CHECK_MODULES([PSS],[pss],[true],[false])
  lal_pre_pss_LIBS="$LIBS"
  LIBS="$PSS_LIBS $LIBS"
  AC_SEARCH_LIBS([crea_sfdbheader],[psssfdb],,
                 [pss="false"]
                 [LIBS="$lal_pre_pss_LIBS"])
  AC_SEARCH_LIBS([pss_set_gw],[pss],,
                 [pss="false"]
                 [LIBS="$lal_pre_pss_LIBS"])
  AC_SEARCH_LIBS([PSS_utc2tdt],[pssastro],,
                 [pss="false"]
                 [LIBS="$lal_pre_pss_LIBS"])
  AC_SEARCH_LIBS([app_star],[pssnovas],,
                 [pss="false"]
                 [LIBS="$lal_pre_pss_LIBS"])
fi


# check for pss headers
if test "${pss}" = "true"; then
  lal_pre_pss_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$PSS_CFLAGS $CPPFLAGS"
  AC_CHECK_HEADERS([pss_sfdb.h],,
                   [pss="false"]
                   [CPPFLAGS="$lal_pre_pss_CPPFLAGS"])
  AC_CHECK_HEADERS([pss_gw.h],,
                   [pss="false"]
                   [CPPFLAGS="$lal_pre_pss_CPPFLAGS"])
  AC_CHECK_HEADERS([pss_astro.h],,
                   [pss="false"]
                   [CPPFLAGS="$lal_pre_pss_CPPFLAGS"])
  AC_CHECK_HEADERS([novas.h],,
                   [pss="false"]
                   [CPPFLAGS="$lal_pre_pss_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([PSS],[pss])

# use c++ mode for boinc build
if test "${boinc}" = "true" ; then
   AC_LANG_PUSH([C++])
fi

# check for lal libararies and headers
LALSUITE_CHECK_LIB([LAL],[6.10.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.10.0],[LALOpenDataFile],[lal/FileIO.h])
LALSUITE_CHECK_OPT_LIB([LALFrame],[1.0.10],[LALFrOpen],[lal/LALFrameL.h])
LALSUITE_CHECK_OPT_LIB([LALMetaIO],[1.1.0],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_OPT_LIB([LALXML],[1.1.3],[XLALXMLFilePrintElements],[lal/LALXML.h])
LALSUITE_CHECK_OPT_LIB([LALSimulation],[0.7.0],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_OPT_LIB([LALBurst],[1.2.0],[XLALEPSearch],[lal/EPSearch.h])
LALSUITE_CHECK_OPT_LIB([LALInspiral],[1.5.0],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_OPT_LIB([LALPulsar],[1.7.0],[XLALComputeFaFb],[lal/ComputeFstat.h])
LALSUITE_CHECK_OPT_LIB([LALInference],[1.1.0],[LALInferenceAddVariable],[lal/LALInference.h])
LALSUITE_CHECK_OPT_LIB([LALStochastic],[1.1.4],[LALStochasticOptimalFilter],[lal/StochasticCrossCorrelation.h])

# back to c mode
if test "${boinc}" = "true" ; then
   AC_LANG_POP([C++])
fi

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
LALSUITE_ADD_GCC_FLAGS

# define if CUDA is enabled
if test "${cuda}" = "true"; then
  AC_DEFINE([LALAPPS_CUDA_ENABLED],[1],[Define if using cuda library])
fi

# enable condor compile
if test "x${condor}" = "xtrue" ; then
  CC="$CONDOR_COMPILE $CC"
  AM_LDFLAGS="$AM_LDFLAGS -all-static"
  AC_DISABLE_SHARED
  AC_SUBST([AM_LDFLAGS])
  AC_DEFINE([LALAPPS_CONDOR],[1],[LALApps is condor compiled])
fi

# add -all-static to AM_LDFLAGS for --enable-static-binaries
if test "${static_binaries}" = "true"; then
  AM_LDFLAGS="$AM_LDFLAGS -all-static"
  AC_DISABLE_SHARED
  AC_SUBST([AM_LDFLAGS])
fi

# restore all static to ldflags
if test "${allstatic}" = "true" ; then
  LDFLAGS="-static -all-static $LDFLAGS"
fi

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

AC_OUTPUT
echo "
==================================================
LALApps has now been successfully configured:

* FFTW library support is $FFTW_ENABLE_VAL
* Frame library support is $FRAME_ENABLE_VAL
* MetaIO library support is $METAIO_ENABLE_VAL
* LALFrame library support is $LALFRAME_ENABLE_VAL
* LALMetaIO library support is $LALMETAIO_ENABLE_VAL
* LALXML library support is $LALXML_ENABLE_VAL
* LALSimulation library support is $LALSIMULATION_ENABLE_VAL
* LALBurst library support is $LALBURST_ENABLE_VAL
* LALInspiral library support is $LALINSPIRAL_ENABLE_VAL
* LALStochastic library support is $LALSTOCHASTIC_ENABLE_VAL
* LALPulsar library support is $LALPULSAR_ENABLE_VAL
* LALInference library support is $LALINFERENCE_ENABLE_VAL
* FITS library support is $CFITSIO_ENABLE_VAL
* PSS library support is $PSS_ENABLE_VAL
* Condor support is $CONDOR_ENABLE_VAL
* BOINC support is $BOINC_ENABLE_VAL
* MPI support is $MPI_ENABLE_VAL
* GDS support is $GDS_ENABLE_VAL
* CUDA support is $CUDA_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALApps,
and run 'make install' to install LALApps.
==================================================
"
