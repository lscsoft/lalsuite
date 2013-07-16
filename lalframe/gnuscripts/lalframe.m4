# lalframe.m4 - lalframe specific macros
#
# serial 3

AC_DEFUN([LALFRAME_WITH_FRAME_LIBRARY],[
AC_ARG_WITH(
    [frame_library],
    AC_HELP_STRING([--with-frame-library],[specify frame library to use [framel]]),[
        AS_CASE([${with_frame_library}],
            [framec],[framec="true"],
            [framel],[framel="true"],
            [AC_MSG_ERROR([Unknown Frame Library ${with_frame_library}])])
    ],[framec="false"; framel="true"])
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
