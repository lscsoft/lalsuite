# lalframe.m4 - lalframe specific macros
#
# serial 4

AC_DEFUN([LALFRAME_ENABLE_FRAMEC],[
AC_ARG_ENABLE(
    [framec],
    AC_HELP_STRING([--enable-framec],[use the FrameC frame library]),[
        AS_CASE([${enableval}],
            [yes],[framec=true],
            [no],[framec=false],
            [AC_MSG_ERROR([bad value ${enableval} for --enable-framec])])
    ],[framec=false])
])

AC_DEFUN([LALFRAME_ENABLE_INOTIFY],
[AC_ARG_ENABLE(
  [inotify],
  AC_HELP_STRING([--enable-inotify],[use Linux filesystem change notification for low latency frame input @<:@default=check@:>@]),[
    AS_CASE([${enableval}],
      [yes],[inotify="true"],
      [no],[inotify="false"],
      [AC_MSG_ERROR(bad value ${enableval} for --enable-inotify)])
  ],[inotify="true"])
])
