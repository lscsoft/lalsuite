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

if [ ! -d "$figdir" ]; then
    mkdir -p $figdir;
fi

for i in $allfigs; do
    ln -f -s $i $figdir &> /dev/null;
    ln -f -s $i $latexdir &> /dev/null;
done

