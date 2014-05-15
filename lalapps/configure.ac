AC_PREREQ([2.63])
AC_INIT([LALApps],[6.14.0.1],[lal-discuss@ligo.org])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_SRCDIR([src/lalapps/lalapps.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([\
  Makefile \
  lalapps.spec \
  debian/Makefile \
  doc/Makefile \
  man/Makefile \
  src/Makefile \
  src/lalapps/Makefile \
  src/lalapps/LALAppsVCSInfo.c \
  src/lalapps/LALAppsVCSInfo.h \
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
  src/pulsar/hough/Makefile \
  src/pulsar/hough/src/Makefile \
  src/pulsar/SideBand/Makefile \
  src/pulsar/CombSearch/Makefile \
  src/pulsar/templateBanks/Makefile \
  src/pulsar/TimingTests/Makefile \
  src/pulsar/TDS_isolated/Makefile \
  src/pulsar/CreateEphemeris/Makefile \
  src/pulsar/Tools/Makefile \
  src/inspiral/Makefile \
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
# FIXME: -Wno-unsupported is a hack which silences warnings regarding
# subdir-objects in automake-1.14+. this will be the default, and only,
# behaviour in automake-2.0+ so this is only a temporary measure
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests -Wno-unsupported])

AC_CANONICAL_HOST

LALSUITE_CHECK_GIT_REPO
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

LALSUITE_WITH_CUDA

LALSUITE_ENABLE_GCC_FLAGS
LALSUITE_ENABLE_NIGHTLY
LALAPPS_ENABLE_CONDOR
LALSUITE_ENABLE_BOINC
LALAPPS_ENABLE_STATIC_BINARIES

LALAPPS_ENABLE_MPI
LALAPPS_ENABLE_BAMBI

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

# define version/configure info
LALSUITE_VERSION_CONFIGURE_INFO

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# enable Fortran and C++ compilers for BAMBI
AS_IF([test "$bambi" = "true"], [
  LALSUITE_REQUIRE_F77
  LALSUITE_REQUIRE_CXX
])

# check for required compilers
LALSUITE_PROG_COMPILERS

# check for MPI compilers
bambimpi=false
if test "x$mpi" = "xtrue"; then
  AX_MPI([mpi=true],[mpi=false])
  AS_IF([test "$bambi" = "true"], [
    MPILIBS_SAVE="${MPILIBS}"
    AC_LANG([C++])
    AX_MPI([bambimpi=true],[])
    AC_LANG([Fortran])
    AX_MPI([],[bambimpi=false])
    AC_LANG([C])
    MPILIBS="${MPILIBS_SAVE}"
  ])
  AC_LANG([C])
fi
LALSUITE_ENABLE_MODULE([MPI])

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PATH_PYTHON([2.6],,[:])
AM_CONDITIONAL([HAVE_PYTHON],[test "$PYTHON" != :])

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
AC_CHECK_PROGS([LATEX],[pdflatex latex],[echo])
AC_CHECK_PROGS([MKIND],[makeindex],[echo])
AC_CHECK_PROGS([DVIPS],[dvips],[echo])
AC_CHECK_PROGS([BIBTEX],[bibtex],[echo])

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

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
LALSUITE_ADD_FLAGS([C],[${GSL_CFLAGS}],[${GSL_LIBS}])
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_version],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_CHECK_GSL_VERSION([1.13])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# check for fftw libraries
if test "${fftw}" = "true"; then
  PKG_CHECK_MODULES([FFTW3],[fftw3 fftw3f],[true],[false])
  LALSUITE_ADD_FLAGS([C],[${FFTW3_CFLAGS}],[${FFTW3_LIBS}])
  AC_CHECK_LIB([fftw3f],[fftwf_execute_dft],,[AC_MSG_WARN([could not find the fftw3f library])],[-lm])
  AC_CHECK_LIB([fftw3],[fftw_execute_dft],,[AC_MSG_WARN([could not find the fftw3 library])],[-lm])
fi

# check for fftw headers
if test "${fftw}" = "true"; then
  AC_CHECK_HEADERS([fftw3.h],,[AC_MSG_WARN([could not find the fftw3.h header])])
fi
LALSUITE_ENABLE_MODULE([FFTW])

# check for GDS system: we check if we can link against libawg which is in GDS
# GDS is a requirements for building 'online_psinject', which is used for hardware injections
if test "${gds}" = "true"; then
  AC_CHECK_LIB([awg],[awgGetChannelNames],,[gds="false"])
  AC_CHECK_HEADERS([gdstype.h],,[gds="false"])
fi
LALSUITE_ENABLE_MODULE([GDS])

# check for libFrame
if test "${frame}" = "true"; then
  PKG_CHECK_MODULES([FRAME],[libframe],[true],[false])
  lal_pre_frame_LIBS="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${FRAME_LIBS}])
  AC_SEARCH_LIBS([FrLibIni],[Frame],,
                 [AC_MSG_WARN([could not find the frame library])]
                 [frame="false"]
                 [lalframe="false"]
                 [LIBS="$lal_pre_frame_LIBS"])
fi

# check for libFrame headers
if test "${frame}" = "true"; then
  lal_pre_frame_CPPFLAGS="$CPPFLAGS"
  LALSUITE_ADD_FLAGS([C],[${FRAME_CFLAGS}],[])
  AC_CHECK_HEADERS([FrameL.h],,
                   [AC_MSG_WARN([could not find the FrameL.h header])]
                   [frame="false"]
                   [lalframe="false"]
                   [CPPFLAGS="$lal_pre_frame_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([FRAME])

# check for libmetaio
if test "${metaio}" = "true"; then
  PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
  lal_pre_metaio_LIBS="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${METAIO_LIBS}])
  AC_SEARCH_LIBS([MetaioOpen],[metaio],,
                 [AC_MSG_WARN([could not find the metaio library])]
                 [metaio="false"]
                 [lalmetaio="false"]
                 [LIBS="$lal_pre_metaio_LIBS"])
fi

# check for libmetaio headers
if test "${metaio}" = "true"; then
  lal_pre_metaio_CPPFLAGS="$CPPFLAGS"
  LALSUITE_ADD_FLAGS([C],[${METAIO_CFLAGS}],[])
  AC_CHECK_HEADERS([metaio.h],,
                   [AC_MSG_WARN([could not find the metaio.h header])]
                   [metaio="false"]
                   [lalmetaio="false"]
                   [CPPFLAGS="$lal_pre_metaio_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([METAIO])

# check for libxml2
if test "$lalxml" = "true"; then
  PKG_CHECK_MODULES([LIBXML2],[libxml-2.0 >= 2.6],[true],[false])
  lal_pre_libxml2_LIBS="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${LIBXML2_LIBS}])
  AC_SEARCH_LIBS([xmlInitParser],[xml2],,
                 [AC_MSG_WARN([cannot find the xml2 library])]
                 [lalxml="false"]
                 [LIBS="$lal_pre_libxml2_LIBS"])
fi

# check for libxml2 headers
if test "$lalxml" = "true"; then
  lal_pre_libxml2_CPPFLAGS="$CPPFLAGS"
  LALSUITE_ADD_FLAGS([C],[${LIBXML2_CFLAGS}],[])
  AC_CHECK_HEADERS([libxml/tree.h],,
                   [AC_MSG_WARN([could not the libxml/tree.h header])]
                   [lalxml="false"]
                   [CPPFLAGS="$lal_pre_libxml2_CPPFLAGS"])
fi

# check for libcfitsio
if test "${cfitsio}" = "true"; then
  PKG_CHECK_MODULES([CFITSIO],[cfitsio],[true],[false])
  lal_pre_cfitsio_libs="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${CFITSIO_LIBS}])
  AC_SEARCH_LIBS([ffopen],[cfitsio],,
                 [cfitsio="false"]
                 [LIBS="$lal_pre_cfitsio_LIBS"])
fi

# check for libcfitsio headers
if test "${cfitsio}" = "true"; then
  lal_pre_cfitsio_CPPFLAGS="$CPPFLAGS"
  LALSUITE_ADD_FLAGS([C],[${CFITSIO_CFLAGS}],[])
  AC_CHECK_HEADERS([fitsio.h],,
                   [cfitsio="false"]
                   [CPPFLAGS="$lal_pre_cfitsio_CPPFLAGS"])
fi
LALSUITE_ENABLE_MODULE([CFITSIO])

# check for pss
if test "${pss}" = "true"; then
  PKG_CHECK_MODULES([PSS],[pss],[true],[false])
  lal_pre_pss_LIBS="$LIBS"
  LALSUITE_ADD_FLAGS([C],[],[${PSS_LIBS}])
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
  LALSUITE_ADD_FLAGS([C],[${PSS_CFLAGS}],[])
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
LALSUITE_ENABLE_MODULE([PSS])

# check for libbambi
# Note:  if done in AS_IF() then AC_REQUIRE'd checks are done when not wanted, causing problems
BAMBI_ENABLE_VAL="DISABLED"
hbf=false
if test "$bambi" = "true"; then
  LALAPPS_CHECK_BAMBI
fi
AM_CONDITIONAL([HAVE_BAMBI],[${hbf}])
AM_CONDITIONAL([HAVE_BAMBI_MPI],[${bambimpi}])

# use c++ mode for boinc build
if test "${boinc}" = "true" ; then
   AC_LANG_PUSH([C++])
fi

# check for lal libararies and headers
LALSUITE_CHECK_LIB([LAL],[6.12.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.12.0],[XLALCacheImport],[lal/LALCache.h])
LALSUITE_CHECK_OPT_LIB([LALFrame],[1.1.0],[LALFrOpen],[lal/LALFrameL.h])
LALSUITE_CHECK_OPT_LIB([LALMetaIO],[1.2.0],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_OPT_LIB([LALXML],[1.1.3],[XLALXMLFilePrintElements],[lal/LALXML.h])
LALSUITE_CHECK_OPT_LIB([LALSimulation],[1.1.0],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_OPT_LIB([LALBurst],[1.2.0],[XLALEPSearch],[lal/EPSearch.h])
LALSUITE_CHECK_OPT_LIB([LALInspiral],[1.5.0],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_OPT_LIB([LALPulsar],[1.9.0],[XLALComputeFstat],[lal/ComputeFstat.h])
LALSUITE_CHECK_OPT_LIB([LALInference],[1.3.0],[LALInferenceAddVariable],[lal/LALInference.h])
LALSUITE_CHECK_OPT_LIB([LALStochastic],[1.1.4],[LALStochasticOptimalFilter],[lal/StochasticCrossCorrelation.h])

# check that lal was built with BOINC support, if required
LALSUITE_CHECK_BOINC

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
  LALSUITE_ADD_FLAGS([],[],[-all-static])
  AC_DISABLE_SHARED
  AC_DEFINE([LALAPPS_CONDOR],[1],[LALApps is condor compiled])
fi

# add -all-static to AM_LDFLAGS for --enable-static-binaries
if test "${static_binaries}" = "true"; then
  LALSUITE_ADD_FLAGS([],[],[-all-static])
  AC_DISABLE_SHARED
fi

# force all python modules into the exec dir to prevent mixed C and
# pure-python library code from being split across lib and lib64 on
# RedHat-derived systems.
AC_SUBST([pkgpythondir], ["\${pkgpyexecdir}"])

# output
AC_OUTPUT
AC_MSG_NOTICE([
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
* BAMBI library support is $BAMBI_ENABLE_VAL
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
])
