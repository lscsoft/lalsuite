#!/bin/sh

## Be sure to set all important paths in the environment variables:
##
##   C_INCLUDE_PATH
##   LIBRARY_PATH
##   LD_LIBRARY_PATH
##   PATH

TAR="gtar"
MAKE="gmake"
TOPDIR=`pwd`
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

echo "*** `date`: Start fresh" | tee $LOG

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

echo "*** `date`: Perform tests of CVS archive"

BUILDDIR='.'
test0 >> $LOG 2>&1 || exit 1
test1 >> $LOG 2>&1 || exit 1
test2 >> $LOG 2>&1 || exit 1
test3 >> $LOG 2>&1 || exit 1
test4 >> $LOG 2>&1 || exit 1

# Make a distribution extract it

echo "*** `date`: Create and extract distribution"

./configure >> $LOG 2>&1 || exit 1
make dist >> $LOG 2>&1 || exit 1
$TAR zxvf lal-*.*.tar.gz >> $LOG 2>&1 || exit 1
rm -f lal-*.*.tar.gz >> $LOG 2>&1 || exit 1
mv lal-*.* tmp >> $LOG 2>&1 || exit 1
cd tmp >> $LOG 2>&1 || exit 1

# Perform tests of distribution

echo "*** `date`: Perform tests of distribution"

BUILDDIR='.'
test0 >> $LOG 2>&1 || exit 1
test1 >> $LOG 2>&1 || exit 1
test2 >> $LOG 2>&1 || exit 1
test3 >> $LOG 2>&1 || exit 1
test4 >> $LOG 2>&1 || exit 1

# Perform tests in sub-directory

echo "*** `date`: Perform tests of distribution in subdirectory"

mkdir tmp && cd tmp || exit 1

BUILDDIR='..'
test0 >> $LOG 2>&1 || exit 1
test1 >> $LOG 2>&1 || exit 1
test2 >> $LOG 2>&1 || exit 1
test3 >> $LOG 2>&1 || exit 1
test4 >> $LOG 2>&1 || exit 1

# Clean up

echo "*** `date`: Clean-up"

cd .. || exit 1
rm -rf tmp || exit 1
cd .. || exit 1
rm -rf tmp || exit 1

rm -f  /tmp/lib/liblal.*
rm -rf /tmp/include/lal
rm -rf /tmp/doc/lal-*.*

# Successful test

echo "******************************************************"
echo "*                                                    *"
echo "* Congratulations: LAL has been successfully tested! *"
echo "*                                                    *"
echo "******************************************************"
echo `date`

exit 0

