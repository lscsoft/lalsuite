AC_INIT([lalburst],[1.0],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h src/LALBurstConfig.h])
AC_CONFIG_SRCDIR([include/config.h.in])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalburst.spec \
  lalburst.pc \
  Makefile \
  gnuscripts/Makefile \
  debian/Makefile \
  include/Makefile \
  include/lal/Makefile \
  src/Makefile \
  test/Makefile
])
AC_CONFIG_FILES([src/git_version],[chmod +x src/git_version])
AM_INIT_AUTOMAKE([foreign])

AC_CANONICAL_HOST

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
AC_SUBST([LIBCURRENT],[0])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# default install path
AC_PREFIX_DEFAULT([/opt/lscsoft/lalburst])

# put version/configure info in config header
lalburst_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALBURST_VERSION],["$VERSION"],[LALBurst Version])
AC_DEFINE_UNQUOTED([LALBURST_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALBurst Configure Arguments])
AC_DEFINE_UNQUOTED([LALBURST_CONFIGURE_DATE],["$lalburst_configure_date"],
                   [LALBurst Configure Date])

# enable strict gcc flags
LAL_ENABLE_GCC_FLAGS

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

# check for lal
LALSUITE_CHECK_LIB([LAL],[lal],[6.3.2.2],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSUPPORT],[lalsupport],[6.3.2.2],[LALOpenDataFile],[lal/FileIO.h])
LALSUITE_CHECK_LIB([LALMETAIO],[lalmetaio],[1.0],[XLALCreateMetaTableDir],[lal/LIGOLwXMLRead.h])

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
AC_SUBST_FILE(LAL_HEADER_LINK)
AC_PROG_MAKE_SET
AC_OUTPUT
echo "
================================================================

        LAL Burst has now been successfully configured
        and will be installed under the directory:

          ${prefix}

        Now run 'make' to build LAL Burst
        and run 'make install' to install LAL Burst

================================================================="
