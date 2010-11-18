#!/bin/bash

printerror() {
    echo ERROR: "$@" >& 2
    exit 1
}
export -f printerror

[ -z "$abs_builddir"   ] && printerror $0 needs '$abs_builddir'
[ -z "$abs_srcdir"     ] && printerror $0 needs '$abs_srcdir'
[ -z "$abs_top_srcdir" ] && printerror $0 needs '$abs_top_srcdir'
[ -z "$MKDIR_P"        ] && printerror $0 needs '$MKDIR_P'
[ -z "$LN_S"           ] && printerror $0 needs '$LN_S'

figdir="$abs_builddir/figures"
latexdir="$abs_builddir/latex"

cmd=""
for basedir in "$figdir" "$latexdir"; do
    if [ ! -d "$basedir" ] && ! ( $MKDIR_P "$basedir" ); then
        printerror $0 could not make directory "'$basedir'"
    fi
    cmd="$cmd ( cd '$basedir' && rm -f '%f' && $LN_S '%p' '%f' ) ||"
    cmd="$cmd printerror $0 failed to link '%f' in '$basedir';\n"
done

find $abs_top_srcdir -path $abs_srcdir -prune -o \
     -regex '.*\.\(eps\|png\|pdf\)$' -printf "$cmd" | /bin/bash
