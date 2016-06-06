# -*- mode: autoconf; -*-
# lalsuite_gccflags.m4 - macros to set strict gcc flags
#
# serial 26

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
    AC_LANG_PUSH([C])
    gcc_cflags=
    gcc_ldflags=

    # check for GCC warning flags
    # FIXME: Remove -Wno-misleading-indentation when indentation issues are fixed
    LALSUITE_CHECK_COMPILE_FLAGS([
      -g3
      -Wall
      -W
      -Werror
      -Wmissing-prototypes
      -Wstrict-prototypes
      -Wshadow
      -Wpointer-arith
      -Wcast-qual
      -Wcast-align
      -Wwrite-strings
      -fno-common
      -Wnested-externs
      -Wno-format-zero-length
      -fno-strict-aliasing
      -Wno-unused-result
      -Wno-unknown-pragmas
      -Wno-misleading-indentation
      -Qunused-arguments
      ],[gcc_cflags="${gcc_cflags} ${flag}"]
    )

    # check for Mac OS X specific flags
    AS_IF([test "x${MACOSX_VERSION}" != "x"],[
      LALSUITE_CHECK_LINK_FLAGS([
        -Wl[,]-no_compact_unwind
        ],[gcc_ldflags="${gcc_ldflags} ${flag}"]
      )
    ])

    # add flags
    LALSUITE_ADD_FLAGS([C],[${gcc_cflags}],[${gcc_ldflags}])

    AC_LANG_POP([C])
  ])
  # end $0
])
