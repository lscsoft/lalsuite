#!/bin/sh

## ----- allow local testing outside of 'make'
if test -n "$1"; then
    srcdir="$1";
else
    srcdir=".";
fi
if test -n "$2"; then
    prefix="$2";
else
    prefix="";
fi

pwd0=`pwd`;
cd ${srcdir};

outfile="misc/${prefix}GitID.h"
tmpfile="misc/__${prefix}GitID.h"

## ---------- read out git-log of last commit --------------------
fmt="format:${prefix}CommitID: %H %n${prefix}CommitDate: %aD %n${prefix}CommitAuthor: %ae %n${prefix}CommitTitle: %s";
logcmd="git log -1 --pretty='$fmt'";

## first check if the git log works at all: could be a) no git installed or b) no git-repository
if eval "$logcmd > /dev/null 2>&1"; then
    git_log_ok="true";
else
    git_log_ok="false";
fi

if test "$git_log_ok" = "true"; then
    git_log=`eval $logcmd`;
else
    git_log="${prefix}CommitID: unknown.";
fi

## ---------- check for modified/added git-content [ignores untracked files!] ----------
statuscmd="git status";
## first check if git status is actually available:
if test "$git_log_ok" = "true"; then
    if eval `$statuscmd -a 2>&1 1> /dev/null` ; then
	git_status_ok="true";
    fi
else
    git_status_ok="false";
fi

if test "$git_status_ok" = "true"; then
    if eval "$statuscmd -a > /dev/null 2>&1"; then
	git_status="UNCLEAN: some modifications were not commited!";
    else
	git_status="CLEAN. All modifications commited.";
    fi
else	## no git or git-repository
    git_status="unknown.";
fi

git_log="${git_log}
${prefix}GitStatus: ${git_status}";

## ---------- parse output and generate ${prefix}GitID.h header file --------------------

## make sure the 'git_log' string doesn't contain any double-quotes or $-signs, from commit-messages,
## to avoid messing up the C-string and ident-keywords
git_log_safe=`echo "$git_log" | sed -e"s/\"/''/g" | sed -e"s/[$]/_/g" `;

## put proper $ quotations for each line to form proper 'ident' keyword strings:
git_log_ident=`echo "$git_log_safe" | sed -e"s/\(.*\)/\"$\1 $ \\\\\\n\"/g"`;

cat > $tmpfile <<EOF
#ifndef ${prefix}GITID_H
#define ${prefix}GITID_H
#include <lal/LALRCSID.h>
NRCSID (${prefix}GitID, $git_log_ident);
#endif
EOF

if diff $tmpfile $outfile  >/dev/null 2>&1; then
    ## files are identical, no update required
    exit 0;
else
    ## files diff: update
    cp $tmpfile $outfile;
fi

cd ${pwd0}
