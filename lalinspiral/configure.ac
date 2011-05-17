AC_INIT([lalinspiral],[1.1.0.2],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h src/LALInspiralConfig.h])
AC_CONFIG_SRCDIR([src/LALInspiral.h])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalinspiral.spec \
  lalinspiral.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  test/Makefile
])
AM_INIT_AUTOMAKE([foreign])

AC_CANONICAL_HOST

LALSUITE_ARG_VAR

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
AC_SUBST([LIBCURRENT],[1])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# default install path
AC_PREFIX_DEFAULT([/opt/lscsoft/lalinspiral])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# debug
LALSUITE_ENABLE_DEBUG

# put version/configure info in config header
lalinspiral_version_major=`echo "$VERSION" | cut -d. -f1`
lalinspiral_version_minor=`echo "$VERSION" | cut -d. -f2`
lalinspiral_version_micro=`echo "$VERSION" | cut -d. -f3`
lalinspiral_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalinspiral_version_micro" && lalinspiral_version_micro=0
test -z "$lalinspiral_version_devel" && lalinspiral_version_devel=0
lalinspiral_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALINSPIRAL_VERSION],["$VERSION"],[LALInspiral Version])
AC_DEFINE_UNQUOTED([LALINSPIRAL_VERSION_MAJOR],[$lalinspiral_version_major],
                   [LALInspiral Version Major Number])
AC_DEFINE_UNQUOTED([LALINSPIRAL_VERSION_MINOR],[$lalinspiral_version_minor],
                   [LALInspiral Version Minor Number])
AC_DEFINE_UNQUOTED([LALINSPIRAL_VERSION_MICRO],[$lalinspiral_version_micro],
                   [LALInspiral Version Micro Number])
AC_DEFINE_UNQUOTED([LALINSPIRAL_VERSION_DEVEL],[$lalinspiral_version_devel],
                   [LALInspiral Version Devel Number])
AC_DEFINE_UNQUOTED([LALINSPIRAL_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALInspiral Configure Arguments])
AC_DEFINE_UNQUOTED([LALINSPIRAL_CONFIGURE_DATE],["$lalinspiral_configure_date"],
                   [LALInspiral Configure Date])

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

# optional lal libraries
LALSUITE_ENABLE_LALFRAME

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99],[AC_PROG_CC_C99],[LALSUITE_AC_PROG_CC_C99])

# use silent build rules if available
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

# checks for programs
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
LT_INIT
AM_PATH_PYTHON([2.4],,[AC_MSG_ERROR([Python-2.4, or higher, is required])])

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

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# metaio
PKG_CHECK_MODULES([METAIO],[libmetaio],[true],[false])
CPPFLAGS="$CPPFLAGS $METAIO_CFLAGS"
LIBS="$LIBS $METAIO_LIBS"
AC_CHECK_LIB([metaio],[MetaioOpen],,[AC_MSG_ERROR([could not find the metaio library])])

# check for libmetaio headers
AC_CHECK_HEADERS([metaio.h],,[AC_MSG_ERROR([could not find the metaio.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.4.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.4.0],[LALOpenDataFile],[lal/FileIO.h])
LALSUITE_CHECK_LIB([LALMetaIO],[1.0.1],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])
LALSUITE_CHECK_OPT_LIB([LALFrame],[1.0.1],[LALFrOpen],[lal/LALFrameL.h])

# add gcc specific flags
if test "$GCC" = yes; then
  CFLAGS="$CFLAGS $lal_gcc_flags"

  # add mac os x specific flags
  if test "x$MACOSX_VERSION" != "x"; then
    CFLAGS="$CFLAGS -mmacosx-version-min=10.4"
  fi
fi

# output
LAL_HEADER_LINK=${ac_aux_dir}/header-link.mk
AC_SUBST_FILE([LAL_HEADER_LINK])
AC_PROG_MAKE_SET
AC_OUTPUT
echo "
================================================================

        LALInspiral has now been successfully configured
        and will be installed under the directory:

          ${prefix}

        Now run 'make' to build LALInspiral
        and run 'make install' to install LALInspiral

================================================================="
