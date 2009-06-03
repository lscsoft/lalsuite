#!/bin/sh

cwd=`pwd`;
pkgdir="$cwd/../../packages";

epsfigures=`find $pkgdir -name "*.eps" -print | grep -e"[^.].*[.]eps$"`;
pngfigures=`find $pkgdir -name "*.png" -print`;
figdir="./figures";
latexdir="./latex";

if [ ! -d "$latexdir" ]; then
    mkdir -p $latexdir;
fi

if [ ! -d "$figdir" ]; then
    rm -f $figdir;
    mkdir -p $figdir;
fi

for i in $epsfigures $cwd/*.eps; do
    ln -sf $i $figdir &> /dev/null;
    # hack to get around doxygen missing some eps-figures for latex
    ln -sf $i $latexdir &> /dev/null;
done

for i in $pngfigures $cwd/*.png; do
    ln -sf $i $figdir &> /dev/null;
done
echo
