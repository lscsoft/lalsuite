#!/bin/sh

# make sure that 00boot and configure have run
if ! test -f configure.ac ; then
  if test -f configure.in ; then
    echo "!!! autoconf-2.13 has been used: this script may not work"
  else
    echo "!!! no file configure.ac: run 00boot and configure first"
    exit 1
  fi
fi
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




















# How to fail
fail () {
  echo "!!! `date`: Failed.  Consult file `basename $LOG`."
  exit 1;
}

# Help message
helpmsg="Usage $0 [options]
Options: [defaults in brackets after description]"
helpmsg="$helpmsg
  --help                print this message"
helpmsg="$helpmsg
  --fresh               fresh start and stop [no]"
helpmsg="$helpmsg
  --include-deps        pass option --include-deps to automake on 00boot [no]"
helpmsg="$helpmsg
  --use-tar=TAR         use the program TAR [tar]"
helpmsg="$helpmsg
  --use-make=MAKE       use the program MAKE [make]"
helpmsg="$helpmsg
  --prefix=PREFIX       set prefix to PREFIX [/tmp/`logname`]"
helpmsg="$helpmsg
  --extra-config-args=ARGS configure with extra arguments ARGS [none]"
helpmsg="$helpmsg
  --without-gcc-flags   do not add --with-gcc-flags to configure arguments"
helpmsg="$helpmsg
  --disable-frame       do not add --enable-frame to configure arguments"
helpmsg="$helpmsg
  --disable-mpi         do not add --enable-mpi to configure arguments"
helpmsg="$helpmsg
  --disable-static      do not build static library"
helpmsg="$helpmsg
  --disable-shared      do not build shared library"
helpmsg="$helpmsg
  --production          test build in LDAS production mode"


# Setable variables

FRESH="no"
TAR="tar"
MAKE="make"
PREFIX="/tmp/`logname`"
CONFIGARGS=""
INCLUDEDEPS=""
GCCFLAGS="yes"
FRAME="yes"
MPI="yes"
STATIC="yes"
SHARED="yes"
PRODUCTION="no"

# Process args
while test $# -gt 0 ; do
  option=$1
  case "$option" in
    -*=*) optarg=`echo "$option" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
  esac
  case $option in
    -h | -help | --help ) echo "$helpmsg"; exit 0;;
    -f | -fresh | --fresh ) FRESH="yes";;
    -include-deps | --include-deps ) INCLUDEDEPS=" --include-deps";;
    -use-tar=* | --use-tar=* ) TAR="$optarg";;
    -use-make=* | --use-make=* ) MAKE="$optarg";;
    -prefix=* | --prefix=* ) PREFIX="$optarg";;
    -extra-config-args=* | --extra-config-args=* ) CONFIGARGS="$CONFIGARGS $optarg";;
    -without-gcc-flags | --without-gcc-flags ) GCCFLAGS="no";;
    -disable-frame | --disable-frame ) FRAME="no";;
    -disable-mpi | --disable-mpi ) MPI="no";;
    -disable-static | --disable-static ) STATIC="no";;
    -disable-shared | --disable-shared ) SHARED="no";;
    -production | --production ) PRODUCTION="yes";;
    *) echo "unrecognized option $option"; exit 1;;
  esac
  shift
done

if test "x$GCCFLAGS" = "xyes" ; then
  CONFIGARGS="$CONFIGARGS --with-gcc-flags"
fi
if test "x$FRAME" = "xyes" ; then
  CONFIGARGS="$CONFIGARGS --enable-frame"
fi
if test "x$MPI" = "xyes" ; then
  CONFIGARGS="$CONFIGARGS --enable-mpi"
fi
if test "x$STATIC" = "xno" ; then
  CONFIGARGS="$CONFIGARGS --disable-static"
fi
if test "x$SHARED" = "xno" ; then
  CONFIGARGS="$CONFIGARGS --disable-shared"
fi
if test "x$PRODUCTION" = "xyes" ; then
  CONFIGARGS="--enable-mpi --disable-static --disable-debug"
  STATIC="no"
fi

echo ">>> Using configure arguments: $CONFIGARGS"

TOPDIR="`pwd`"
LOG="$TOPDIR/testscript.log"

# Start fresh

lal=`sed -n 's/package=//p' 00boot`
ver=`sed -n 's/version=//p' 00boot`
lalver="$lal-$ver"

echo ">>> Testing $lalver [sending output to file `basename $LOG`]"

test -f $LOG && rm -f $LOG
test -d $lalver               && rm -rf $lalver 
rm -f $PREFIX/lib/lib$lal.*
test -d $PREFIX/include/$lal  && rm -rf $PREFIX/include/$lal
test -d $PREFIX/doc/$lalver   && rm -rf $PREFIX/doc/$lalver

if test -f Makefile -a "x$FRESH" = "xyes" ; then
  MESSAGE=">>> `date`: Start fresh"
  echo "$MESSAGE" > $LOG
  echo "$MESSAGE"
  $MAKE cvs-clean >> $LOG 2>&1 || fail
fi

MESSAGE=">>> `date`: Make distribution"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

./00boot $INCLUDEDEPS >> $LOG 2>&1 || fail
./configure >> $LOG 2>&1 || fail
$MAKE dist >> $LOG 2>&1 || fail
tarfile="$lalver.tar.gz"
test -f $tarfile || fail
$TAR zxvf $tarfile >> $LOG 2>&1 || fail
rm -f $tarfile
test -d $lalver || fail
cd $lalver || fail
mkdir tmp || fail
cd tmp || fail

MESSAGE=">>> `date`: Test build"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

../configure --prefix=$PREFIX $CONFIGARGS >> $LOG 2>&1 || fail
$MAKE >> $LOG 2>&1 || fail
$MAKE check >> $LOG 2>&1 || fail
$MAKE dvi >> $LOG 2>&1 || fail
$MAKE install >> $LOG 2>&1 || fail
test -f $PREFIX/include/$lal/LALStdlib.h || fail
if test "x$STATIC" = "xyes" ; then
  test -f $PREFIX/lib/lib$lal.a || fail
fi
if test "x$SHARED" = "xyes" ; then
  test -f $PREFIX/lib/lib$lal.so || fail
fi
$MAKE uninstall >> $LOG 2>&1 || fail
test ! -f $PREFIX/include/$lal/LALStdlib.h || fail
test ! -f $PREFIX/lib/lib$lal.a || fail
test ! -f $PREFIX/lib/lib$lal.so || fail
$MAKE distclean >> $LOG 2>&1 || fail

MESSAGE=">>> `date`: Cleanup"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

cd .. || fail
rm -rf tmp || fail
cd .. || fail
rm -rf $lalver || fail

if test "x$FRESH" = "xyes" ; then
  $MAKE cvs-clean >> $LOG 2>&1 || fail
fi

MESSAGE=">>> `date`: Finished"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

echo "******************************************************"
echo "*                                                    *"
echo "* Congratulations: LAL has been successfully tested! *"
echo "*                                                    *"
echo "******************************************************"
exit 0
