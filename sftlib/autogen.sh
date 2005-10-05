#!/bin/sh
# $Id$
#
# Generate ./configure from config.in and Makefile.in from Makefile.am.
# This also adds files like missing,depcomp,install-sh to the source
# direcory. To update these files at a later date use:
#	autoreconf -f -i -v


# Cygwin?
test -x /usr/bin/uname && /usr/bin/uname | grep -i CYGWIN >/dev/null &&
{
    # Enable strict case checking
    # (to avoid e.g "DIST_COMMON = ... ChangeLog ..." in Makefile.in)
    export CYGWIN="${CYGWIN}${CYGWIN:+ }check_case:strict"

    # Check for Unix text file type
    echo > dostest.tmp
    test "`wc -c < dostest.tmp`" -eq 1 ||
        echo "Warning: DOS text file type set, 'make dist' and related targets will not work."
    rm -f dostest.tmp
}

echo "Bootstrapping configure script and makefiles:"

ACLOCAL=aclocal
AUTOHEADER=autoheader
AUTOMAKE=automake
AUTOCONF=autoconf

if $ACLOCAL && $AUTOHEADER && $AUTOMAKE --foreign && $AUTOCONF; then

    echo "Done, now you can run ./configure"
    exit 0
else
    echo "Something failed .... please check error-message and re-run when fixed."
    echo "exiting..."
    exit 1
fi
