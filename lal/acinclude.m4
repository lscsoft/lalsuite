dnl acinclude.m4

AC_DEFUN(AC_LIB_FFTW,
[AC_SEARCH_LIBS(fftw_one, sfftw fftw, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW library libfftw.a         *"
 echo "* You must install FFTW (v >= 2.0) on your system.           *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If libfftw.a is in, e.g., /home/bob/lib, you must include  *"
 echo "* this in the LIBRARY_PATH environment variable.             *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
AC_MSG_ERROR([FFTW must be properly installed.])
], -lm)
])

AC_DEFUN(AC_LIB_RFFTW,
[AC_SEARCH_LIBS(rfftw_one, srfftw rfftw, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW library librfftw.a        *"
 echo "* You must install FFTW (v >= 2.0) on your system.           *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If librfftw.a is in, e.g., /home/bob/lib, you must include *"
 echo "* this in the LIBRARY_PATH environment variable.             *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
AC_MSG_ERROR(FFTW must be properly installed.)
], -lm)
])

AC_DEFUN(AC_HEADER_FFTW,
[AC_CHECK_HEADER(fftw.h, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW header fftw.h             *"
 echo "* You must install FFTW on your system.                      *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "**************************************************************"
])])

AC_DEFUN(AC_HEADER_RFFTW,
[AC_CHECK_HEADER(rfftw.h, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW header rfftw.h            *"
 echo "* You must install FFTW on your system.                      *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "**************************************************************"
])])


AC_DEFUN(AC_FFTW_WORKS,
[AC_MSG_CHECKING(whether FFTW works)
AC_TRY_RUN([
#include <stdio.h>
#include <sfftw.h>
int main() { return sizeof(fftw_real)-4; }
],
AC_MSG_RESULT(yes)
AC_DEFINE(FFTW_PREFIX),
AC_TRY_RUN([
#include <stdio.h>
#include <fftw.h>
int main() { return sizeof(fftw_real)-4; }
],
AC_MSG_RESULT(yes),
AC_MSG_RESULT(no)
echo "**************************************************************"
echo "* WARNING: FFTW does not seem to be working.                 *"
echo "* Possible problems:                                         *"
echo "*   - FFTW version < 2.0                                     *"
echo "*   - Compiler could not find header sfftw.h or fftw.h       *"
echo "*   - FFTW was configured with the --enable-float option     *"
echo "*                                                            *"
echo "* Remove the file config.cache before re-running configure.  *"
echo "*                                                            *"
echo "* Please see the instructions in the file INSTALL.           *"
echo "**************************************************************"
AC_MSG_ERROR([FFTW must be properly installed.])
))])

