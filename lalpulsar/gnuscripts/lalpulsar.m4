# -*- mode: autoconf; -*-
# lalpulsar.m4 - lalpulsar specific macros
#
# serial 1

AC_DEFUN([LALPULSAR_WITH_SIMD],[
  # $0: check for SIMD extensions
  AC_ARG_WITH(
    [simd],
    AC_HELP_STRING([--with-simd],[use SIMD extensions @<:@default: yes@:>@]),
    [],[with_simd=yes]
  )
  LALSUITE_PUSH_UVARS
  LALSUITE_CLEAR_UVARS
  SIMD_FLAGS=
  AS_CASE([${with_simd}],
    [yes],[
      AX_EXT
      AX_CHECK_COMPILE_FLAG([-mfpmath=sse],[SIMD_FLAGS="${SIMD_FLAGS} -mfpmath=sse"],[:],[-Werror])
      AX_GCC_ARCHFLAG([yes],[SIMD_FLAGS="${SIMD_FLAGS} ${ax_cv_gcc_archflag}"])
    ],
    [no],[:],
    [AC_MSG_ERROR([bad value '${with_simd}' for --with-simd])]
  )
  AS_IF([test "x${SIMD_FLAGS}" != x],[
    CFLAGS="${SIMD_FLAGS} -O0"
    AC_MSG_CHECKING([whether C compiler assembles basic math with SIMD extensions])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[[
double volatile a = 1.2;
double volatile b = 3.4;
double volatile c = a * b;
]])],[
      AC_MSG_RESULT([yes])
    ],[
      AC_MSG_RESULT([no; disabling SIMD extensions])
      SIMD_FLAGS=
    ])
  ])
  LALSUITE_POP_UVARS
  AS_IF([test "x${SIMD_FLAGS}" != x],[
    LALSUITE_ADD_FLAGS([C],[${SIMD_FLAGS}],[])
    SIMD_ENABLE_VAL=ENABLED
  ],[
    SIMD_ENABLE_VAL=DISABLED
  ])
  # end $0
])
