#!/bin/sh

# make sure configure has been run
if ! test -f configure ; then
  echo "!!! no file configure: run 00boot and configure first"
  exit 1
fi

# make sure that configure has been run
if ! ( test -f config.status && test -f lal-config && test -f Makefile ) ; then
  echo "!!! configure needs to be re-run"
  exit 1
fi

LALVERSION="`./lal-config --version`"
RCSTAG=`echo '$Name$' | sed 's/\\$//g' | sed 's/Name: //'`
LOGFILE=${LOGFILE:-testscript.log}
test -f $LOGFILE && mv -f $LOGFILE $LOGFILE.old

cleanup () {
  test -f lal-$LALVERSION.tar.gz && rm -f lal-$LALVERSION.tar.gz
  test -d lal-$LALVERSION && chmod -R +w lal-$LALVERSION && rm -rf lal-$LALVERSION
}

fail () {
  message="!!! `date`: Failed.  Consult file $LOGFILE."
  echo $message >&3
  echo $message
  # cleanup
  exit 1
}

trap 'echo "!!! `date`: Interrupt" >&4 ; cleanup ; exit 1' 1 2 15

# Reassign stdout and stderr to logfile
# Use descriptors 3 and 4 as original stdout and stderr
exec 3>&1
exec 4>&2
exec 1>$LOGFILE
exec 2>&1
set -x

message=">>> `date`: Testing LAL $LALVERSION $RCSTAG"
echo $message >&3
echo $message

make distcheck DISTCHECK_CONFIGURE_FLAGS='--with-cc="$(CC)" --with-extra-cflags="$(CFLAGS)" --with-extra-cppflags="$(CPPFLAGS)" --with-extra-ldflags="$(LDFLAGS)" --with-extra-libs="$(LIBS)"' || fail

message=">>> `date`: Finished"
echo $message >&3
echo $message

echo "******************************************************" >&3
echo "*                                                    *" >&3
echo "* Congratulations: LAL has been successfully tested! *" >&3
echo "*                                                    *" >&3
echo "******************************************************" >&3

cleanup
exit
