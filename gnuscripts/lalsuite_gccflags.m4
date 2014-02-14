# -*- mode: autoconf; -*-
# lalsuite_gccflags.m4 - macros to set strict gcc flags
#
# serial 19

AC_DEFUN([LALSUITE_ENABLE_GCC_FLAGS],[
  # $0: enable GCC warning flags
  AC_ARG_ENABLE([gcc_flags],
    AC_HELP_STRING([--enable-gcc-flags],[turn on strict GCC warning flags [default=yes]]),
    AS_CASE(["${enableval}"],
      [yes|no],[:],
      AC_MSG_ERROR([bad value '${enableval}' for --enable-gcc-flags])
    ),[
      enable_gcc_flags=yes
    ]
  )
  # end $0
])

AC_DEFUN([LALSUITE_ADD_GCC_FLAGS],[
  # $0: add GCC warning flags
  AS_IF([test "x${GCC}" = xyes && test "x${enable_gcc_flags}" = xyes],[

    gcc_flags="-g3 -Wall -W -Werror -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -fno-common -Wnested-externs -Wno-format-zero-length -fno-strict-aliasing"

    # check if compiler supports -Wno-unused-result
    LALSUITE_PUSH_UVARS
    LALSUITE_CLEAR_UVARS
    CFLAGS="-Wno-unused-result"
    AC_MSG_CHECKING([whether compiler supports -Wno-unused-result])
    AC_COMPILE_IFELSE([
      AC_LANG_PROGRAM([])
    ],[
      AC_MSG_RESULT([yes])
      gcc_flags="${gcc_flags} -Wno-unused-result"
    ],[
      AC_MSG_RESULT([no])
    ])
    LALSUITE_POP_UVARS

# comment out usage of -pedantic flag due to gcc bug 7263
# http://gcc.gnu.org/bugzilla/show_bug.cgi?id=7263
#
#  AS_CASE(["${host_cpu}-${host_os}"],
#      [*i386-darwin*], [gcc_flags="${gcc_flags} -pedantic"],
#      [*x86_64-darwin*], [gcc_flags="${gcc_flags} -pedantic"],
#      [gcc_flags="${gcc_flags} -pedantic-errors"])

    # add mac os x specific flags
    AS_IF([test "x${MACOSX_VERSION}" != "x"],[
      gcc_flags="${gcc_flags} -mmacosx-version-min=10.4"
      # check if compiler supports -Wl,-no_compact_unwind
      LALSUITE_PUSH_UVARS
      LALSUITE_CLEAR_UVARS
      LDFLAGS="-Wl,-no_compact_unwind"
      AC_MSG_CHECKING([whether linker supports -Wl,-no_compact_unwind])
      AC_LINK_IFELSE([
        AC_LANG_PROGRAM([])
      ],[
        AC_MSG_RESULT([yes])
        gcc_ldflags="-Wl,-no_compact_unwind"
      ],[
        AC_MSG_RESULT([no])
      ])
      LALSUITE_POP_UVARS
    ])

    gcc_cflags="${gcc_flags}"
    gcc_cxxflags="${gcc_cxxflags}"

    # don't warn about unknown warning options for clang/clang++
    AS_IF([test -n "${CLANG_CC}"], [gcc_cflags="-Wno-unknown-warning-option ${gcc_cflags}"])
    AS_IF([test -n "${CLANG_CXX}"], [gcc_cxxflags="-Wno-unknown-warning-option ${gcc_cxxflags}"])

    # ignore unused flags with clang/clang++
    AS_IF([test -n "${CLANG_CC}"], [gcc_cflags="${gcc_cflags} -Xcompiler -Qunused-arguments"])
    AS_IF([test -n "${CLANG_CXX}"], [gcc_cxxflags="${gcc_cxxflags} -Xcompiler -Qunused-arguments"])

    # don't use gcc flags when cuda is enabled
    AS_IF([test "x${cuda}" = "xtrue"],[
      AC_MSG_NOTICE([CUDA support is enabled, disabling GCC flags])
    ],[
      LALSUITE_ADD_FLAGS([C],[${gcc_cflags}],[])
      LALSUITE_ADD_FLAGS([CXX],[${gcc_cxxflags}],[])
      LALSUITE_ADD_FLAGS([LD],[${gcc_ldflags}],[])
    ])

  ])
  # end $0
])
