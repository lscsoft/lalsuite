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
fmt_id="format:%H";
logcmd_id="git log -1 --pretty='$fmt_id'";
fmt_date="format:%ai";
logcmd_date="git log -1 --pretty='$fmt_date'";
fmt_author="format:%ae";
logcmd_author="git log -1 --pretty='$fmt_author'";
fmt_title="format:%s";
logcmd_title="git log -1 --pretty='$fmt_title'";

## first check if the git log works at all: could be a) no git installed or b) no git-repository
if eval "$logcmd_id > /dev/null 2>&1"; then
    git_log_ok="true";
else
    git_log_ok="false";
fi

if test "$git_log_ok" = "true"; then
    git_id=`eval $logcmd_id`;
    git_date=`eval $logcmd_date`;
    git_author=`eval $logcmd_author`;
    git_title=`eval $logcmd_title`;
    git_date_utc=`date -ud "$git_date" +"%Y-%m-%d %H:%M:%S %z"`
else
    git_id="unknown.";
    git_date="1980-01-06 00:00:00 +0000";
    git_author="unknown.";
    git_title="unknown.";
    git_date_utc=$git_date;
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

## make sure the 'git_title' string doesn't contain any double-quotes or $-signs, from commit-messages,
## to avoid messing up the C-string and ident-keywords
git_title_safe=`echo "$git_title" | sed -e"s/\"/''/g" | sed -e"s/[$]/_/g" `;

## put proper $ quotations for each line to form proper 'ident' keyword strings:
cat > $tmpfile <<EOF
#ifndef ${prefix}GITID_H
#define ${prefix}GITID_H
#include <lal/LALRCSID.h>
NRCSID (${prefix}GitID, "\$${prefix}CommitID: $git_id  \$ \n"
"\$${prefix}CommitDate: $git_date_utc  \$ \n"
"\$${prefix}CommitAuthor: $git_author  \$ \n"
"\$${prefix}CommitTitle: $git_title_safe \$ \n"
"\$${prefix}GitStatus: $git_status \$ ");
NRCSID (${prefix}GitCommitID, "\$${prefix}CommitID: $git_id \$ ");
NRCSID (${prefix}GitCommitDate, "\$${prefix}CommitDate: $git_date_utc \$ ");
NRCSID (${prefix}GitCommitAuthor, "\$${prefix}CommitAuthor: $git_author \$ ");
NRCSID (${prefix}GitCommitTitle, "\$${prefix}CommitTitle: $git_title_safe \$ ");
NRCSID (${prefix}GitGitStatus, "\$${prefix}GitStatus: $git_status \$ ");
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
