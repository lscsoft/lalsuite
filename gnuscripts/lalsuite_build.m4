AC_DEFUN([LALSUITE_ENABLE_MODULE],[
AM_CONDITIONAL([$1],[test x$$2 = xtrue])
eval $1_ENABLE_VAL="`eval test "$$2" = "true" && echo "ENABLED" || echo "DISABLED"`"
])

AC_DEFUN([LALSUITE_CHECK_LIB],[
PKG_CHECK_MODULES([$1],[$2 >= $3],[eval $2="true"],[eval $2="false"])
if test "$$2" = "true"; then
  CPPFLAGS="$CPPFLAGS $$1_CFLAGS"
  LIBS="$LIBS $$1_LIBS"
  if test "$LALSUITE_BUILD" = "true"; then
    AC_DEFINE([HAVE_LIB$1],[1],[Define to 1 if you have the $1 library])
    eval $2="true"
  else
    AC_CHECK_LIB([$2],[$4],[eval $2="true"],[AC_MSG_ERROR([could not find the $1 library])])
    AC_DEFINE([HAVE_LIB$1],[1],[Define to 1 if you have the $1 library])
  fi
else
  AC_MSG_ERROR([could not find the $1 library])
fi
LALSUITE_ENABLE_MODULE([$1],[$2])
])

AC_DEFUN([LALSUITE_CHECK_OPT_LIB],[
PKG_CHECK_MODULES([$1],[$2 >= $3],[eval $2="true"],[eval $2="false"])
if test "$$2" = "true"; then
  if test "$LALSUITE_BUILD" = "true"; then
    AC_DEFINE([HAVE_LIB$1],[1],[Define to 1 if you have the $1 library])
    eval $2="true"
    CPPFLAGS="$CPPFLAGS $$1_CFLAGS"
    LIBS="$LIBS $$1_LIBS"
  else
    CPPFLAGS="$CPPFLAGS $$1_CFLAGS"
    LIBS="$LIBS $$1_LIBS"
    AC_CHECK_LIB([$2],[$4],[eval $2="true"],[eval $2=false
      AC_MSG_WARN([could not find the $1 library])])
    if test "$$2" = true; then
      AC_DEFINE([HAVE_LIB$1],[1],[Define to 1 if you have the $1 library])
    fi
  fi
fi
LALSUITE_ENABLE_MODULE([$1],[$2])
])
