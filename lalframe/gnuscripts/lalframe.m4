# lalframe.m4 - lalframe specific macros
#
# serial 7

AC_DEFUN([LALFRAME_ENABLE_FRAMEC],
[AC_ARG_ENABLE(
    [framec],
    AC_HELP_STRING([--enable-framec],[enable support for FrameC library]),[
        AS_CASE([${enableval}],
            [yes],[framec="true"],
            [no],[framec="false"],
            [AC_MSG_ERROR([bad value ${enableval} for --enable-framec])])
    ],[framec="false"])
])

AC_DEFUN([LALFRAME_ENABLE_FRAMEL],
[AC_ARG_ENABLE(
    [framel],
    AC_HELP_STRING([--enable-framel],[enable support for FrameL library]),[
        AS_CASE([${enableval}],
            [yes],[framel="true"],
            [no],[framel="false"],
            [AC_MSG_ERROR([bad value ${enableval} for --enable-framel])])
    ],[framel="true"])
])
