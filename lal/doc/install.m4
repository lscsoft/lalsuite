dnl $Id$
dnl process this file with m4 to get installation instructions
INSTALLATION INSTRUCTIONS
changequote(@@,@@)
This file describes how to perform a basic install of LAL.  This is not
quite a minimal install -- support for the frame library is included -- but it
is probably the best basic installation.  If you are following the instructions
from this document, please note that the instructions are written for use with
a Bourne-like shell (e.g., bash) rather than a C-shell (e.g., tcsh).  Really
the only difference is that you need to change statements like

        MYENV=value
        export MYENV

(Bourne shell syntax) to

        setenv MYENV value

The instructions in this file are extracted from the shell scripts:

        doc/lal-preinstall-x.sh
        doc/lal-install-x.sh

Also extracted are nicely-formatted versions of these shell scripts:

        doc/lal-preinstall.sh
        doc/lal-install.sh

If this file (README.install) or these scripts are not present in the version
of LAL that you obtained from CVS, it is because they have not yet been
auto-generated.  Run ./00boot to create them.

You can edit these scripts (as appropriate) and run them to install LAL.
(But note that these scripts run in a Bourne shell, so you don't need to
change how the environment variables are set in these scripts.)
Or just follow the instructions below.


1. PRE-INSTALLATION

Building LAL requires some software to be pre-installed.  This section
describes how to download and install this software if it is missing from
your system.  If you wish to install pre-compiled RPMs you can skip to that
section [Section 1.2].  The next section describes how to build the required
software from source.

1.1 BUILDING REQUIRED SOFTWARE FROM SOURCE

include(@@lal-preinstall.txt@@)

1.2 INSTALLING REQUIRED SOFTWARE FROM BINARY RPMS

WARNING: These instructions are intended to be helpful for users of RPM who
wish to install the software required by LAL system-wide.  However, these
instructions and the provided RPMs are not officially supported, and may not
work on all systems.

The commands listed below are appropriate for a Bourne-shell (e.g., bash);
they will need to be modified appropriately for C-shells (e.g., tcsh).
This is where to get sources:
        
        LALRPMURL=http://www.lsc-group.phys.uwm.edu/lal/sources/rpms

get required autoconf, automake, fftw3, gsl, libframe, and libmetaio
you can use "lynx -dump" or "wget -O-" instead of "curl"

        curl $LALRPMURL/fftw3-3.0.1-2.i386.rpm > fftw-3.0.1-2.i386.rpm
        curl $LALRPMURL/fftw3-devel-3.0.1-2.i386.rpm > fftw-devel-3.0.1-2.i386.rpm
        curl $LALRPMURL/gsl-1.5-1.i386.rpm > gsl-1.5-1.i386.rpm
        curl $LALRPMURL/gsl-devel-1.5-1.i386.rpm > gsl-devel-1.5-1.i386.rpm
        curl $LALRPMURL/libframe-6.14-1.i386.rpm > libframe-6.14-1.i386.rpm
        curl $LALRPMURL/libframe-devel-6.14-1.i386.rpm > libframe-devel-6.14-1.i386.rpm
        curl $LALRPMURL/libframe-utils-6.14-1.i386.rpm > libframe-utils-6.14-1.i386.rpm
        curl $LALRPMURL/libmetaio-5.4-3.i386.rpm > libmetaio-5.4-3.i386.rpm
        curl $LALRPMURL/libmetaio-devel-5.4-3.i386.rpm > libmetaio-devel-5.4-3.i386.rpm

now login as root and install the RPMs

        rpm -Uvh fftw-3.0.1-2.i386.rpm fftw-devel-3.0.1-2.i386.rpm
        rpm -Uvh gsl-1.5-1.i386.rpm gsl-devel-1.5-1.i386.rpm
        rpm -Uvh libframe-6.14-1.i386.rpm libframe-devel-6.14-1.i386.rpm libframe-utils-6.14-1.i386.rpm
        rpm -Uvh libmetaio-5.4-3.i386.rpm libmetaio-devel-5.4-3.i386.rpm


1.3 BUILDING REQUIRED SOFTWARE FROM SOURCE RPMS

If you cannot use these binary RPMs and you wish to build from source RPMs,
follow these instructions.

first download the source RPMs

        curl $LALRPMURL/fftw3-3.0.1-2.src.rpm > fftw3-3.0.1-2.src.rpm
        curl $LALRPMURL/gsl-1.5-1.src.rpm > gsl-1.5-1.src.rpm
        curl $LALRPMURL/libframe-6.14-1.src.rpm > libframe-6.14-1.src.rpm
        curl $LALRPMURL/libmetaio-5.4-3.src.rpm > libmetaio-5.4-3.src.rpm

now login as root and build binary RPMs from the source RPMs (this will
take a very long time)

        rpmbuild --rebuild fftw3-3.0.1-2.src.rpm
        rpmbuild --rebuild gsl-1.5-1.src.rpm
        rpmbuild --rebuild libframe-6.14-1.src.rpm
        rpmbuild --rebuild libmetaio-5.4-3.src.rpm

finally install the new binary RPMs

        cd /usr/src/redhat/RPMS
        rpm -Uvh i386/fftw3-3.0.1-2.i386.rpm
        rpm -Uvh i386/gsl-1.5-1.i386.rpm
        rpm -Uvh i386/gsl-devel-1.5-1.i386.rpm
        rpm -Uvh i386/libframe-6.14-1.i386.rpm
        rpm -Uvh i386/libframe-devel-6.14-1.i386.rpm
        rpm -Uvh i386/libmetaio-5.4-3.i386.rpm
        rpm -Uvh i386/libmetaio-devel-5.4-3.i386.rpm


2. INSTALLING LAL

2.1 BASIC BUILD OF LAL

include(@@lal-install.txt@@)

2.2 MORE DETAILS

Other useful make targets are:

        make dvi                # make documentation
        make check              # run the basic test suite
        make uninstall          # uninstall the library and header files
        make clean              # clean up compiled code (as before "make")
        make distclean          # clean up distribution (as before "configure")
        make cvs-clean          # clean up to cvs files (as before "00boot")

see the file INSTALL for additional details on configuring LAL.


3. SYSTEM-SPECIFIC INSTALLATION INSTRUCTIONS

SGI running IRIX 6.5 with gcc:

  * Configure with the option --with-cc="gcc -n32".

  * If you put shared objects (e.g., of the frame library) in non-standard
    places and are hoping to use LD_LIBRARY_PATH to locate them, you may need
    to set the environment variable LD_LIBRARYN32_PATH too.

  * If you have command-lines that are too long, you'll need to change
    the length of the lines allowed.  To do this use systune -i
    (or perhaps systune -r):

        $ systune -r
        systune-> ncargs 204800
        systune-> quit
    
    This increases the command line length maximum until reboot.
    Change it permanently with systune -b.

Alpha running Linux with ccc --- making shared libraries:

  * Problem: libtool doesn't know that ccc makes position-independent code
    by default (I hope ccc makes PIC by default...).

  * Solution: trick the configure script into thinking that you are using
    OSF/1 by using the --host option:

      ./configure --host=alpha-dec-osf3 --enable-shared

  * Note: use the same trick to make shared libraries for fftw!


HP-UX 10.20 may need

   --with-extra-cppflags="-D_HPUX_SOURCE"


Mac OS X (10.2.x, possibly 10.1.x, but NOT 10.3.x) with bundled cc/gcc:

  * Configure with:  --with-extra-cflags="-D_ANSI_SOURCE -no-cpp-precomp"

  * Note: I (Jolien) don't need these with 10.2 ... perhaps it depends on the
    version of the developer tools.  Also, do NOT use these flags with 10.3.


TROUBLESHOOTING

* If you need to re-run configure after it has failed while checking for a
  working FFTW, FrameL, or MPI, make sure to remove the file config.cache.

* The configure script assumes that ranlib is necessary unless it cannot find
  the program in your path.  If ranlib is on your path and you don't need
  ranlib, set the environment RANLIB to echo.

* "make dvi" must be run after "make" since make dvi requires the program
  laldoc must be compiled for "make dvi" to work.

* If you want to use a different latex program than the one chosen by the
  configure script, set it in the environment variable LATEX.  Also, if you
  want to have different tex flags (e.g., you want to disable the batchmode
  that the configure script uses) set the TEXFLAGS environment variable
  to the flags you want (or a space if you don't want any flags used).

* If you want to make a shared library version (default) of LAL with frame
  library and/or MPI interface, you need to make a shared library version of
  fftw, the frame library, and mpi too.  To make a static LAL library only,
  use the --disable-shared option when configuring LAL.

