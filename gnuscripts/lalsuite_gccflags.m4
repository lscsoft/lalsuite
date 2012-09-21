# lalsuite_gccflags.m4 - macros to set strict gcc flags
#
# serial 7

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
  lal_gcc_flags="-g3 -Wall -W -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -Wnested-externs -Wno-format-zero-length -fno-strict-aliasing -Wno-unused-result"

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
