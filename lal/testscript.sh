#!/bin/sh +x

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
  --config-args=ARGS    configure with arguments ARGS
                        [\"--enable-frame --enable-mpi --with-gcc-flags\"]"
helpmsg="$helpmsg
  --extra-config-args=ARGS configure with extra arguments ARGS [none]"

# Setable variables

FRESH="no"
TAR="tar"
MAKE="make"
PREFIX="/tmp/`logname`"
CONFIGARGS="--with-gcc-flags --enable-frame --enable-mpi"
INCLUDEDEPS=""

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
    -config-args=* | --config-args=* ) CONFIGARGS="$optarg";;
    -extra-config-args=* | --extra-config-args=* ) CONFIGARGS="$CONFIGARGS $optarg";;
    *) echo "unrecognized option $option"; exit 1;;
  esac
  shift
done

TOPDIR="`pwd`"
LOG="$TOPDIR/testscript.log"

# Start fresh

lal=`sed -n 's/package=//p' 00boot`
ver=`sed -n 's/version=//p' 00boot`
lalver="$lal-$ver"

echo ">>> Testing $lalver [sending output to file `basename $LOG`]"

test -f $LOG && rm -f $LOG
test -d $lalver               && rm -rf $lalver 
test -f $PREFIX/lib/lib$lal.* && rm -f  $PREFIX/lib/lib$lal.*
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
test -f $PREFIX/lib/lib$lal.a || fail
test -f $PREFIX/lib/lib$lal.so || fail
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


###
### OLD VERSION
###

## Be sure to set all important paths in the environment variables:
##
##   C_INCLUDE_PATH
##   LIBRARY_PATH
##   LD_LIBRARY_PATH
##   PATH

TOPDIR="`pwd`"
LOG="$TOPDIR/testscript.log"

# Configure arguments for various tests:

TEST0CONFARGS="--prefix=/tmp"
TEST1CONFARGS="$TEST0CONFARGS --enable-frame"
TEST2CONFARGS="$TEST1CONFARGS --enable-mpi"
TEST3CONFARGS="$TEST2CONFARGS --with-gcc-flags"
TEST4CONFARGS="$TEST3CONFARGS --disable-debug"

# The function that performs the tests:

test0() {
  $BUILDDIR/configure $TEST0CONFARGS &&         \
  $MAKE &&                                      \
  $MAKE check &&                                \
  $MAKE dvi &&                                  \
  $MAKE install &&                              \
  test -f /tmp/include/lal/LALStdlib.h &&       \
  test -f /tmp/lib/liblal.a &&                  \
  test -f /tmp/lib/liblal.so &&                 \
  $MAKE uninstall &&                            \
  test ! -f /tmp/include/lal/LALStdlib.h &&     \
  test ! -f /tmp/lib/liblal.a &&                \
  test ! -f /tmp/lib/liblal.so &&               \
  $MAKE distclean;
}

test1() {
  $BUILDDIR/configure $TEST1CONFARGS &&         \
  $MAKE &&                                      \
  $MAKE check &&                                \
  $MAKE install &&                              \
  test -f /tmp/include/lal/LALStdlib.h &&       \
  test -f /tmp/lib/liblal.a &&                  \
  test -f /tmp/lib/liblal.so &&                 \
  $MAKE uninstall &&                            \
  test ! -f /tmp/include/lal/LALStdlib.h &&     \
  test ! -f /tmp/lib/liblal.a &&                \
  test ! -f /tmp/lib/liblal.so &&               \
  $MAKE distclean;
}

test2() {
  $BUILDDIR/configure $TEST2CONFARGS &&         \
  $MAKE &&                                      \
  $MAKE check &&                                \
  $MAKE install &&                              \
  test -f /tmp/include/lal/LALStdlib.h &&       \
  test -f /tmp/lib/liblal.a &&                  \
  test -f /tmp/lib/liblal.so &&                 \
  $MAKE uninstall &&                            \
  test ! -f /tmp/include/lal/LALStdlib.h &&     \
  test ! -f /tmp/lib/liblal.a &&                \
  test ! -f /tmp/lib/liblal.so &&               \
  $MAKE distclean;
}

test3() {
  $BUILDDIR/configure $TEST3CONFARGS &&         \
  $MAKE &&                                      \
  $MAKE check &&                                \
  $MAKE install &&                              \
  test -f /tmp/include/lal/LALStdlib.h &&       \
  test -f /tmp/lib/liblal.a &&                  \
  test -f /tmp/lib/liblal.so &&                 \
  $MAKE uninstall &&                            \
  test ! -f /tmp/include/lal/LALStdlib.h &&     \
  test ! -f /tmp/lib/liblal.a &&                \
  test ! -f /tmp/lib/liblal.so &&               \
  $MAKE distclean;
}

test4() {
  $BUILDDIR/configure $TEST4CONFARGS &&         \
  $MAKE &&                                      \
  $MAKE install &&                              \
  test -f /tmp/include/lal/LALStdlib.h &&       \
  test -f /tmp/lib/liblal.a &&                  \
  test -f /tmp/lib/liblal.so &&                 \
  $MAKE uninstall &&                            \
  test ! -f /tmp/include/lal/LALStdlib.h &&     \
  test ! -f /tmp/lib/liblal.a &&                \
  test ! -f /tmp/lib/liblal.so &&               \
  $MAKE distclean;
}

# Start fresh

rm -f $LOG

MESSAGE="*** `date`: Start fresh"
echo "$MESSAGE" > $LOG
echo "$MESSAGE"

rm -rf tmp
rm -f  /tmp/lib/liblal.*
rm -rf /tmp/include/lal
rm -rf /tmp/doc/lal-*.*

if test -f Makefile
then
  make cvs-clean >> $LOG 2>&1 || exit 1
fi
./00boot >> $LOG 2>&1

# Perform tests of CVS archive

MESSAGE="*** `date`: Perform tests of CVS archive"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

BUILDDIR='.'
echo "*** `date`:   test0" ; test0 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test1" ; test1 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test2" ; test2 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test3" ; test3 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test4" ; test4 >> $LOG 2>&1 || exit 1

# Make a distribution extract it

MESSAGE="*** `date`: Create and extract distribution"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

./configure >> $LOG 2>&1 || exit 1
make dist >> $LOG 2>&1 || exit 1
$TAR zxvf lal-*.*.tar.gz >> $LOG 2>&1 || exit 1
rm -f lal-*.*.tar.gz >> $LOG 2>&1 || exit 1
mv lal-*.* tmp >> $LOG 2>&1 || exit 1
cd tmp >> $LOG 2>&1 || exit 1

# Perform tests of distribution

MESSAGE="*** `date`: Perform tests of distribution"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

BUILDDIR='.'
echo "*** `date`:   test0" ; test0 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test1" ; test1 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test2" ; test2 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test3" ; test3 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test4" ; test4 >> $LOG 2>&1 || exit 1

# Perform tests in sub-directory

MESSAGE="*** `date`: Perform tests of distribution in subdirectory"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

mkdir tmp && cd tmp || exit 1

BUILDDIR='..'
echo "*** `date`:   test0" ; test0 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test1" ; test1 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test2" ; test2 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test3" ; test3 >> $LOG 2>&1 || exit 1
echo "*** `date`:   test4" ; test4 >> $LOG 2>&1 || exit 1

# Clean up

MESSAGE="*** `date`: Clean-up"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

cd .. || exit 1
rm -rf tmp || exit 1
cd .. || exit 1
rm -rf tmp || exit 1

rm -f  /tmp/lib/liblal.*
rm -rf /tmp/include/lal
rm -rf /tmp/doc/lal-*.*

# Successful test

MESSAGE="*** `date`: Finished"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

echo "******************************************************"
echo "*                                                    *"
echo "* Congratulations: LAL has been successfully tested! *"
echo "*                                                    *"
echo "******************************************************"

exit 0
