#!/bin/sh

# Help message
helpmsg="Usage $0 [options]
Options: [defaults in brackets after description]"
helpmsg="$helpmsg
  --help                print this message"
helpmsg="$helpmsg
  --use-tar=TAR         use the program TAR [tar]"
helpmsg="$helpmsg
  --use-make=MAKE       use the program MAKE [make]"
helpmsg="$helpmsg
  --prefix=PREFIX       set prefix to PREFIX [/tmp]"
helpmsg="$helpmsg
  --config-args=ARGS    configure with arguments ARGS
                        [\"--enable-frame --enable-mpi --with-gcc-flags\"]"

# Setable variables
TAR="tar"
MAKE="make"
PREFIX="/tmp"
CONFIGARGS="--with-gcc-flags --enable-frame --enable-mpi"

# Process args
while test $# -gt 0 ; do
  option=$1
  case "$option" in
    -*=*) optarg=`echo "$option" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
    *) optarg= ;;
  esac
  case $option in
    -h | -help | --help ) echo "$helpmsg"; exit 0;;
    -use-tar=* | --use-tar=* ) TAR="$optarg";;
    -use-make=* | --use-make=* ) MAKE="$optarg";;
    -prefix=* | --prefix=* ) PREFIX="$optarg";;
    -config-args=* | --config-args=* ) CONFIGARGS="$optarg";;
    *) echo "unrecognized option $option"; exit 1;;
  esac
  shift
done

TOPDIR="`pwd`"
LOG="$TOPDIR/testscript.log"

# Start fresh

lal=`grep AM_INIT_AUTOMAKE configure.in | sed 's/[^a-z]//g'`
ver=`grep AM_INIT_AUTOMAKE configure.in | sed 's/[^\.0-9]//g'`
lalver="$lal-$ver"

echo ">>> Testing $lalver [sending output to file `basename $LOG`]"

test -f $LOG && rm -f $LOG

MESSAGE=">>> `date`: Start fresh"
echo "$MESSAGE" > $LOG
echo "$MESSAGE"

test -d $lalver               && rm -rf $lalver 
test -f $PREFIX/lib/lib$lal.* && rm -f  $PREFIX/lib/lib$lal.*
test -d $PREFIX/include/$lal  && rm -rf $PREFIX/include/$lal
test -d $PREFIX/doc/$lalver   && rm -rf $PREFIX/doc/$lalver

if test -f Makefile ; then
  $MAKE cvs-clean >> $LOG 2>&1 || exit 1
fi

MESSAGE=">>> `date`: Make distribution"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

./00boot >> $LOG 2>&1 || exit 1
./configure >> $LOG 2>&1 || exit 1
$MAKE dist >> $LOG 2>&1 || exit 1
tarfile="$lalver.tar.gz"
test -f $tarfile || exit 1
$TAR zxvf $tarfile >> $LOG 2>&1 || exit 1
rm -f $tarfile
test -d $lalver || exit 1
cd $lalver || exit 1
mkdir tmp || exit 1
cd tmp || exit 1

MESSAGE=">>> `date`: Test build"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

../configure --prefix=$PREFIX $CONFIGARGS >> $LOG 2>&1 || exit 1
$MAKE >> $LOG 2>&1 || exit 1
$MAKE check >> $LOG 2>&1 || exit 1
$MAKE dvi >> $LOG 2>&1 || exit 1
$MAKE install >> $LOG 2>&1 || exit 1
test -f $PREFIX/include/$lal/LALStdlib.h || exit 1
test -f $PREFIX/lib/lib$lal.a || exit 1
test -f $PREFIX/lib/lib$lal.so || exit 1
$MAKE uninstall >> $LOG 2>&1 || exit 1
test ! -f $PREFIX/include/$lal/LALStdlib.h || exit 1
test ! -f $PREFIX/lib/lib$lal.a || exit 1
test ! -f $PREFIX/lib/lib$lal.so || exit 1
$MAKE distclean >> $LOG 2>&1 || exit 1

MESSAGE=">>> `date`: Cleanup"
echo "$MESSAGE" >> $LOG
echo "$MESSAGE"

cd .. || exit 1
rm -rf tmp || exit 1
cd .. || exit 1
rm -rf $lalver || exit 1
$MAKE cvs-clean >> $LOG 2>&1 || exit 1

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
