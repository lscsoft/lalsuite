dnl acinclude.m4

AC_DEFUN(AC_LIB_FFTW,
[AC_SEARCH_LIBS(fftw_one, sfftw fftw, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW library libsfftw.a        *"
 echo "* or libfftw.a                                               *"
 echo "*                                                            *"
 echo "* You must install FFTW (v >= 2.0) on your system.           *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If libfftw.a or libsfftw.a is in, e.g., /home/bob/lib, you *"
 echo "* must set the LIBS environment variable:                    *"
 echo "*                                                            *"
 echo "*   LIBS=-L/home/bob/lib                                     *"
 echo "*                                                            *"
 echo "* or:                                                        *"
 echo "*                                                            *"
 echo "*   setenv LIBS -L/home/bob/lib                              *"
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
 echo "* WARNING: Could not find the FFTW library libsrfftw.a       *"
 echo "* or librfftw.a                                              *"
 echo "*                                                            *"
 echo "* You must install FFTW (v >= 2.0) on your system.           *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If librfftw.a or libsrfftw.a is in, e.g., /home/bob/lib,   *"
 echo "* you must set the LIBS environment variable:                *"
 echo "*                                                            *"
 echo "*   LIBS=-L/home/bob/lib                                     *"
 echo "*                                                            *"
 echo "* or:                                                        *"
 echo "*                                                            *"
 echo "*   setenv LIBS -L/home/bob/lib                              *"
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
 echo "* WARNING: Could not find the FFTW header fftw.h or sfftw.h  *"
 echo "*                                                            *"
 echo "* You must install FFTW on your system.                      *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If fftw.h or sfftw.h is in, e.g., /home/bob/include, you   *"
 echo "* must set the CPPFLAGS environment variable:                *"
 echo "*                                                            *"
 echo "*   CPPFLAGS=-I/home/bob/include                             *"
 echo "*                                                            *"
 echo "* or:                                                        *"
 echo "*                                                            *"
 echo "*   setenv CPPFLAGS -I/home/bob/include                      *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
 echo "**************************************************************"
])])

AC_DEFUN(AC_HEADER_RFFTW,
[AC_CHECK_HEADER(rfftw.h, ,
[echo "**************************************************************"
 echo "* WARNING: Could not find the FFTW header srfftw.h or        *"
 echo "* rfftw.h                                                    *"
 echo "*                                                            *"
 echo "* You must install FFTW on your system.                      *"
 echo "* FFTW is avaliable from http://www.fftw.org                 *"
 echo "* FFTW must be configured with the --enable-float argument.  *"
 echo "* Install FFTW on your system using the commands:            *"
 echo "*                                                            *"
 echo "*   ./configure --enable-float                               *"
 echo "*   make                                                     *"
 echo "*   make install                                             *"
 echo "*                                                            *"
 echo "* If rfftw.h or srfftw.h is in, e.g., /home/bob/include, you *"
 echo "* must set the CPPFLAGS environment variable:                *"
 echo "*                                                            *"
 echo "*   CPPFLAGS=-I/home/bob/include                             *"
 echo "*                                                            *"
 echo "* or:                                                        *"
 echo "*                                                            *"
 echo "*   setenv CPPFLAGS -I/home/bob/include                      *"
 echo "*                                                            *"
 echo "* Remove the file config.cache before re-running configure.  *"
 echo "*                                                            *"
 echo "* Please see the instructions in the file INSTALL.           *"
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
echo "*   - FFTW was not configured with the --enable-float option *"
echo "*                                                            *"
echo "* Remove the file config.cache before re-running configure.  *"
echo "*                                                            *"
echo "* Please see the instructions in the file INSTALL.           *"
echo "**************************************************************"
AC_MSG_ERROR([FFTW must be properly installed.])
))])

