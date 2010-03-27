#!/bin/sh

# Set to location where NDS2 was installed
INSTALLDIR=

# Change as needed
PYTHONINCDIR=/usr/include/python2.4

# Add -DDEBUG if desired
gcc -fPIC -c _nds2.c -I${PYTHONINCDIR} -I${INSTALLDIR}/include

ld -shared -lsasl2 -o _nds2.so _nds2.o ${INSTALLDIR}/lib/libndsclient.a

# For OS X, use this:
# gcc -dynamiclib -lsasl2 -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,1.0,-current_version,1.0 -o _nds2.so _nds2.o ${INSTALLDIR}/lib/libndsclient.a
