# lalsuite_gccflags.m4 - macros to set strict gcc flags
#
# serial 9

AC_DEFUN([LALSUITE_ENABLE_GCC_FLAGS],
[AC_ARG_ENABLE([gcc_flags],
  AC_HELP_STRING([--enable-gcc-flags],[turn on strict gcc warning flags (default=yes)]),
  [case "${enableval}" in
     yes) DO_ENABLE_LALSUITE_GCC_FLAGS;;
     no) ;;
     *) DO_ENABLE_LALSUITE_GCC_FLAGS;;
   esac ],
   [ DO_ENABLE_LALSUITE_GCC_FLAGS ] )
])

AC_DEFUN([DO_ENABLE_LALSUITE_GCC_FLAGS],
[
  lal_gcc_flags="-g3 -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -Wnested-externs -Wno-format-zero-length -fno-strict-aliasing"

  # check if compiler support -Wno-unused-result
  my_save_cflags="$CFLAGS"
  CFLAGS=-Wno-unused-result
  AC_MSG_CHECKING([whether CC supports -Wno-unused-result])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])]
    [lal_gcc_flags="${lal_gcc_flags} -Wno-unused-result"],
    [AC_MSG_RESULT([no])]
  )
  CFLAGS="$my_save_cflags"

  # don't use -Werror in LALApps
  case ${PACKAGE} in
    lalapps) ;;
    *) lal_gcc_flags="${lal_gcc_flags} -Werror" ;;
  esac

# comment out usage of -pedantic flag due to gcc bug 7263
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=7263
#
#  case $host_cpu-$host_os in
#    *i386-darwin*) lal_gcc_flags="${lal_gcc_flags} -pedantic" ;;
#    *x86_64-darwin*) lal_gcc_flags="${lal_gcc_flags} -pedantic" ;;
#    *) lal_gcc_flags="${lal_gcc_flags} -pedantic-errors" ;;
#  esac
])

AC_DEFUN([LALSUITE_ADD_GCC_FLAGS],
[
  if test "x${GCC}" = "xyes"; then
    # don't use gcc flags when cuda is enabled
    if test "x${cuda}" = "xtrue"; then
      AC_MSG_NOTICE([CUDA support is enabled, disabling GCC flags])
    else
      CFLAGS="${CFLAGS} ${lal_gcc_flags}"
    fi

    # add mac os x specific flags
    if test "x${MACOSX_VERSION}" != "x"; then
      CFLAGS="${CFLAGS} -mmacosx-version-min=10.4"
    fi

    # ignore unused flags with clang/clang++
    AS_IF([test -n "${CLANG_CC}"], [CFLAGS="${CFLAGS} -Xcompiler -Qunused-arguments"])
    AS_IF([test -n "${CLANG_CXX}"], [CXXFLAGS="${CXXFLAGS} -Xcompiler -Qunused-arguments"])
  fi
])
