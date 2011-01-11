AC_INIT([lalframe],[1.0.1.1],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h src/LALFrameConfig.h])
AC_CONFIG_SRCDIR([src/LALFrameIO.c])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalframe.spec \
  lalframe.pc \
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
AC_SUBST([LIBAGE],[1])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# default install path
AC_PREFIX_DEFAULT([/opt/lscsoft/lalframe])

# nightly build
LALSUITE_ENABLE_NIGHTLY

# put version/configure info in config header
lalframe_version_major=`echo "$VERSION" | cut -d. -f1`
lalframe_version_minor=`echo "$VERSION" | cut -d. -f2`
lalframe_version_micro=`echo "$VERSION" | cut -d. -f3`
lalframe_version_devel=`echo "$VERSION" | cut -d. -f4-`
test -z "$lalframe_version_micro" && lalframe_version_micro=0
test -z "$lalframe_version_devel" && lalframe_version_devel=0
lalframe_configure_date=`date +"%Y-%m-%dT%H:%M:%S%z"`
AC_DEFINE_UNQUOTED([LALFRAME_VERSION],["$VERSION"],[LALFrame Version])
AC_DEFINE_UNQUOTED([LALFRAME_VERSION_MAJOR],[$lalframe_version_major],
                   [LALFrame Version Major Number])
AC_DEFINE_UNQUOTED([LALFRAME_VERSION_MINOR],[$lalframe_version_minor],
                   [LALFrame Version Minor Number])
AC_DEFINE_UNQUOTED([LALFRAME_VERSION_MICRO],[$lalframe_version_micro],
                   [LALFrame Version Micro Number])
AC_DEFINE_UNQUOTED([LALFRAME_VERSION_DEVEL],[$lalframe_version_devel],
                   [LALFrame Version Devel Number])
AC_DEFINE_UNQUOTED([LALFRAME_CONFIGURE_ARGS],["$ac_configure_args"],
                   [LALFrame Configure Arguments])
AC_DEFINE_UNQUOTED([LALFRAME_CONFIGURE_DATE],["$lalframe_configure_date"],
                   [LALFrame Configure Date])

# enable strict gcc flags
LALSUITE_ENABLE_GCC_FLAGS

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
AC_CHECK_PROGS([M4],[gm4 m4],[m4])

# check for supported mac os x version
LALSUITE_OSX_VERSION_CHECK

# link tests using libtool
LALSUITE_USE_LIBTOOL

# check for header files
AC_HEADER_STDC

# check for gethostname in unistd.h
AC_MSG_CHECKING([for gethostname prototype in unistd.h])
AC_EGREP_HEADER([gethostname],[unistd.h],[AC_MSG_RESULT([yes])]
AC_DEFINE([HAVE_GETHOSTNAME_PROTOTYPE],[1],
          [Define if gethostname prototype is in unistd.h]),
          [AC_MSG_RESULT]([no]))

# check for required libraries
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([could not find the math library])])

# check for frame library
PKG_CHECK_MODULES([LIBFRAME],[libframe],[true],[false])
CPPFLAGS="$CPPFLAGS $LIBFRAME_CFLAGS"
LIBS="$LIBS $LIBFRAME_LIBS"
AC_CHECK_LIB([Frame],[FrLibIni],,[AC_MSG_ERROR([could not find the frame library])])

# check for frame library headers
AC_CHECK_HEADERS([FrameL.h],,[AC_MSG_ERROR([could not find the FrameL.h header])])

# check for lal libraries and headers
LALSUITE_CHECK_LIB([LAL],[6.5.0],[LALVersion],[lal/LALStdio.h])
LALSUITE_CHECK_LIB([LALSupport],[6.5.0],[LALOpenDataFile],[lal/FileIO.h])

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

        LALFrame has now been successfully configured
        and will be installed under the directory:

          ${prefix}

        Now run 'make' to build LALFrame
        and run 'make install' to install LALFrame

================================================================="
