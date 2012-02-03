# lalframe.m4 - lalframe specific macros
#
# serial 1

AC_DEFUN([LALFRAME_ENABLE_INOTIFY],
[AC_ARG_ENABLE(
  [inotify],
  AC_HELP_STRING([--enable-inotify],[use Linux filesystem change notification for low latency frame input @<:@default=check@:>@]),
  [ case "${enableval}" in
      yes) inotify=true;;
      no) inotify=false;;
      *) AC_MSG_ERROR(bad value ${enableval} for --enable-inotify);;
    esac
  ], [ inotify=true ] )
])
