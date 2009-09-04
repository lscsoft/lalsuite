AC_INIT([lalstochastic],[1.0],[lal-discuss@gravity.phys.uwm.edu])
AM_CONFIG_HEADER([include/config.h])
AC_CONFIG_SRCDIR([include/config.h.in])
AC_CONFIG_AUX_DIR([gnuscripts])
AC_CONFIG_MACRO_DIR([gnuscripts])
AC_CONFIG_FILES([ \
  lalstochastic.spec \
  lalstochastic.pc \
  Makefile \
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
AC_SUBST([LIBCURRENT],[0])
AC_SUBST([LIBREVISION],[0])
AC_SUBST([LIBAGE],[0])
AC_SUBST([LIBVERSION],[${LIBCURRENT}:${LIBREVISION}:${LIBAGE}])

# default install path
AC_PREFIX_DEFAULT([/opt/lscsoft/lalstochastic])

# enable strict gcc flags
LAL_ENABLE_GCC_FLAGS

# check for c compiler
m4_pattern_allow([AC_PROG_CC_C99])
m4_ifdef([AC_PROG_CC_C99],[AC_PROG_CC_C99],[LAL_AC_PROG_CC_C99])

# use silent build rules if available
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])],)

# checks for programs
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL

# check for supported mac os x version
if test "x$build_vendor" = "xapple"; then
  AC_CHECK_PROGS(SW_VERS, sw_vers)
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
AC_CHECK_LIB([m],[main],,[AC_MSG_ERROR([cannot find the math library])])

# check for lal
PKG_CHECK_MODULES([LAL],[lal >= 6.3],[true],[false])
CPPFLAGS="$CPPFLAGS $LAL_CFLAGS"
LIBS="$LIBS $LAL_LIBS"

# check for lalsupport, needed for test suite
PKG_CHECK_MODULES([LALSUPPORT],[lalsupport],[true],[false])

# check libraries if not a top-level build
CPPFLAGS="$CPPFLAGS $LAL_CFLAGS $LALSUPPORT_CFLAGS"
LIBS="$LALSUPPORT_LIBS $LAL_LIBS $LIBS"
if test "$LALSUITE_BUILD" = "true" ; then
	test -n "$LALSUPPORT_LIBS" && lalsupport="true"
else
	AC_CHECK_LIB([lal],[main],,[AC_MSG_ERROR(cannot find lal library)])
	AC_CHECK_LIB([lalsupport],[LALOpenDataFile],[lalsupport=true],[lalsupport=false])
fi
AM_CONDITIONAL([LALSUPPORT], [test x$lalsupport = xtrue])

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

        LAL Stochastic has now been successfully configured
        and will be installed under the directory:

          ${prefix}

        Now run 'make' to build LAL support
        and run 'make install' to install LAL Stochastic

================================================================="
