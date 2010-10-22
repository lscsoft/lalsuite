#!/bin/bash

cwd=`pwd`;
topdir="${cwd}/../../..";
doxydir="${topdir}/lal/doc/doxygen"
cmdline="find $topdir -path $doxydir -prune -o \( -regex '.*[.]eps$' -print -o -regex '.*[.]png$' -print -o -regex '.*[.]pdf$' -print \)"
allfigs=`eval $cmdline`;

figdir="./figures";
latexdir="./latex";

if [ ! -d "$latexdir" ]; then
    mkdir -p $latexdir;
fi

if [ -d "$figdir" ]; then
    rm -rf $figdir;
fi
mkdir -p $figdir;

for i in $allfigs; do
    ln -s $i $figdir &> /dev/null;
done

