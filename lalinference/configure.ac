AC_PREREQ([2.63])
AC_INIT([lalinference],[1.1.0.1],[lal-discuss@gravity.phys.uwm.edu])
AC_CONFIG_HEADERS([include/config.h src/LALInferenceConfig.h])
AC_CONFIG_SRCDIR([src/LALInference.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalinference.spec \
  lalinference.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  python/Makefile \
  python/__init__.py \
  swig/Makefile \
  swig/octave/Makefile \
  swig/python/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([1.11 foreign color-tests parallel-tests])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR
LALSUITE_DISTCHECK_CONFIGURE_FLAGS

# set library api info
# 1. start with version information of "0:0:0" for each libtool library
# 2. update the version information only immediately before a public
#    release of your software. More frequent updates are unnecessary,
#    and only guarantee that the current interface number gets larger
#    faster.
# 3. if the library source code has changed at all since the last
#    update, then increment revision ("c:r:a" becomes "c:r+1:a").
# 4. if any interfaces have been added, removed, or changed since the
#    last update, increment current, and set revision to 0.
# 5. if any interfaces have been added since the last public release,
#    then increment age.
# 6. if any interfaces have been removed since the last public release,
#    then set age to 0.
AC_SUBST([LIBCURRENT],[5])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# check for git repo
LALSUITE_CHECK_GIT_REPO

# nightly build
LALSUITE_ENABLE_NIGHTLY

# debug
LALSUITE_ENABLE_DEBUG

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# optional lal libraries
LALSUITE_ENABLE_LALFRAME
LALSUITE_ENABLE_LALMETAIO
LALSUITE_ENABLE_LALINSPIRAL
LALSUITE_ENABLE_LALPULSAR
LALSUITE_ENABLE_LALXML

# enable SWIG wrapping modules
LALSUITE_ENABLE_SWIG

# put version/configure info in config header
lalinference_version_major=`echo "$VERSION" | cut -d. -f1`
lalinference_version_minor=`echo "$VERSION" | cut -d. -f2`
lalinference_version_micro=`echo "$VERSION" | cut -d. -f3`
lalinference_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalinference_version_micro" && lalinference_version_micro=0
test -z "$lalinference_version_devel" && lalinference_version_devel=0
lalinspiral_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALINFERENCE_VERSION],["$VERSION"],[LALInference Version])
AC_DEFINE_UNQUOTED([LALINFERENCE_VERSION_MAJOR],[$lalinference_version_major],
                   [LALInference Version Major Number])
AC_DEFINE_UNQUOTED([LALINFERENCE_VERSION_MINOR],[$lalinference_version_minor],
                   [LALInference Version Minor Number])
AC_DEFINE_UNQUOTED([LALINFERENCE_VERSION_MICRO],[$lalinference_version_micro],
                   [LALInference Version Micro Number])
AC_DEFINE_UNQUOTED([LALINFERENCE_VERSION_DEVEL],[$lalinference_version_devel],
                   [LALInference Version Devel Number])
AC_DEFINE_UNQUOTED([LALINFERENCE_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALInference Configure Arguments])
AC_DEFINE_UNQUOTED([LALINFERENCE_CONFIGURE_DATE],["$lalinspiral_configure_date"],
                   [LALInference Configure Date])

# workaround broken libtool multilib support on sl6
LALSUITE_MULTILIB_LIBTOOL_HACK

# use silent build rules
AM_SILENT_RULES([yes])

# check for C/C++ compilers
LALSUITE_PROG_CC_CXX

# checks for programs
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED
AM_PATH_PYTHON([2.6],,[AC_MSG_ERROR([Python-2.6, or higher, is required])])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for gsl
PKG_CHECK_MODULES([GSL],[gsl],[true],[false])
CPPFLAGS="$CPPFLAGS $GSL_CFLAGS"
LIBS="$LIBS $GSL_LIBS"
AC_CHECK_LIB([gslcblas],[main],,[AC_MSG_ERROR([could not find the gsl library])])
AC_CHECK_LIB([gsl],[gsl_strerror],,[AC_MSG_ERROR([could not find the gsl library])])
LALSUITE_ENABLE_FAST_GSL

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
LIBS="$LIBS $METAIO_LIBS"
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])

# check for libmetaio headers
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

if test "$lalxml" = "true"; then
  # check for libxml2
  PKG_CHECK_MODULES([LIBXML2],[libxml-2.0 >= 2.6],[true],[false])
  CPPFLAGS="$CPPFLAGS $LIBXML2_CFLAGS"
  LIBS="$LIBS $LIBXML2_LIBS"
  AC_CHECK_LIB([xml2],[xmlInitParser],,[AC_MSG_WARN([cannot find the xml2 library])])
fi

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.9.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.9.0],[LALOpenDataFile],[lal/FileIO.h])
LALSUITE_CHECK_LIB([LALMetaIO],[1.0.5],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_LIB([LALFrame],[1.0.7],[LALFrOpen],[lal/LALFrameL.h])
LALSUITE_CHECK_LIB([LALSimulation],[0.6.0],[XLALSimDetectorStrainREAL8TimeSeries],[lal/LALSimulation.h])
LALSUITE_CHECK_LIB([LALInspiral],[1.4.1],[LALInspiralParameterCalc],[lal/LALInspiral.h])
LALSUITE_CHECK_LIB([LALPulsar],[1.6.0],[XLALComputeFaFb],[lal/ComputeFstat.h])

# check for lalxml
LALSUITE_CHECK_OPT_LIB([LALXML],[1.1.5],[XLALXMLFilePrintElements],[lal/LALXML.h])

# configure SWIG wrapping modules
LALSUITE_USE_SWIG([XLALInference LALInference XLAL LAL])
LALSUITE_SWIG_DEPENDS([lal],[true])

# add gcc specific flags
LALSUITE_ADD_GCC_FLAGS

# output
LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE([LAL_HEADER_LINK])
AC_PROG_MAKE_SET
AC_OUTPUT
echo "
==================================================
LALInference has now been successfully configured:

* LALXML library support is $LALXML_ENABLE_VAL
* SWIG wrapping module for Octave is $SWIG_OCTAVE_ENABLE_VAL
* SWIG wrapping module for Python is $SWIG_PYTHON_ENABLE_VAL

and will be installed under the directory:

${prefix}

Now run 'make' to build LALInference,
and run 'make install' to install LALInference.
==================================================
"
