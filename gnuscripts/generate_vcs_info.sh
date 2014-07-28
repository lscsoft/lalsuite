#!/bin/sh

default_git_path=git
default_am_v_gen=
git_path=$default_git_path
am_v_gen=$default_am_v_gen

# print usage message and exit
usage() {
    cat <<EOF
usage: $0 [options] package_name output_file source_file
options [defaults in brackets]:
  -h, --help              display this help and exit
  --git-path=GIT          path to git program [$default_git_path]
  --am-v-gen=\$(AM_V_GEN)  generation message [${default_am_v_gen:-"none"}]
EOF
    exit $1
}

# parse options
prev=
for option; do
    if test -n "$prev" ; then
        eval $prev=\$option
        prev=
        shift
        continue
    fi

    case $option in
    *=?*) optarg=`expr "X$option" : '[^=]*=\(.*\)'` ;;
    *=)   optarg= ;;
    *)    optarg=yes ;;
    esac
    
    case $option in
    -help | --help | --hel | --he | -h) usage 0;;
    --git-path)   prev=git_path;;
    --git-path=*) git_path=$optarg;;
    --am-v-gen)   prev=am_v_gen;;
    --am-v-gen=*) am_v_gen=$optarg;;
    -*) echo "unrecognized option \`$option'" 1>&2; usage 1;;
    *) break;; # a non-option: done with option parsing
    esac

    shift
done

if test -n "$prev"; then
    option=--`echo $prev | sed 's/_/-/g'`
    echo "missing argument to $option" 1>&2
    exit 1
fi

# make sure there the right number of remaining arguments
if test $# != 3; then
    echo "wrong number of arguments" 1>&2
    usage 1
fi

# set argument variables
am_v_gen="`expr "X$am_v_gen" : 'X\([ ]*GEN[ ]*\)'`"
package="${1?"package_name not set"}"
output="${2?"output_file not set"}"
source="${3?"input_file not set"}"
srcdir="`dirname "$source"`"

# use git log to get the important fields
git_id="`cd "$srcdir" && $git_path log -1 --pretty=format:%H`"
git_date="`cd "$srcdir" && env LC_TIME=C TZ=GMT0 $git_path log -1 --date=local --pretty=format:%cd`"
git_author="`cd "$srcdir" && $git_path log -1 --pretty=format:"%an <%ae>"`"
git_committer="`cd "$srcdir" && $git_path log -1 --pretty=format:"%cn <%ce>"`"

# extract relevant fields of %c date format
git_date_mon=`echo "$git_date" | cut -d " " -f2`
git_date_day=`echo "$git_date" | cut -d " " -f3`
git_date_hms=`echo "$git_date" | cut -d " " -f4`
git_date_year=`echo "$git_date" | cut -d " " -f5`

# convert month name to number
case "$git_date_mon" in
    Jan) git_date_mon=01;; Feb) git_date_mon=02;; Mar) git_date_mon=03;;
    Apr) git_date_mon=04;; May) git_date_mon=05;; Jun) git_date_mon=06;;
    Jul) git_date_mon=07;; Aug) git_date_mon=08;; Sep) git_date_mon=09;;
    Oct) git_date_mon=10;; Nov) git_date_mon=11;; Dec) git_date_mon=12;;
esac

# ISO-8601 representation of date in UTC
git_date="$git_date_year-$git_date_mon-$git_date_day $git_date_hms +0000"

# determine git tag, if present
git_tag="`cd "$srcdir" && $git_path describe --exact-match --tags $git_id 2>/dev/null || echo "None"`"

# determine git branch
git_branch="`cd "$srcdir" && $git_path rev-parse --symbolic-full-name HEAD 2>/dev/null`"
case "$git_branch" in
    HEAD) git_branch="None";;
    *) git_branch="`basename "$git_branch"`";;
esac

# update index
(cd "$srcdir" && $git_path update-index -q --refresh)

# determine state of git repository
if ! eval "(cd \"$srcdir\" && $git_path diff-files --quiet)"; then
    git_status="UNCLEAN: Modified working tree"
elif ! eval "(cd \"$srcdir\" && $git_path diff-index --cached --quiet HEAD)"; then
    git_status="UNCLEAN: Modified index"
else
    git_status="CLEAN: All modifications committed"
fi

# determine builder
git_builder_name="`cd "$srcdir" && $git_path config user.name 2>/dev/null || echo "Unknown User"`"
git_builder_email="`cd "$srcdir" && $git_path config user.email 2>/dev/null`"
git_builder="$git_builder_name <$git_builder_email>"

# transliterate package name
LALPackage=$package
LALPACKAGE=`echo $package | sed 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
lalpackage=`echo $package | sed 'y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/'`
package=`echo $package | sed 's/^LAL//'`

# determine appropriate header files
if test "$lalpackage" = "lalapps"; then
    pkg_vcs_info_header="${LALPackage}VCSInfo.h"
    pkg_config_header="config.h"
else
    pkg_vcs_info_header="lal/${LALPackage}VCSInfo.h"
    pkg_config_header="lal/${LALPackage}Config.h"
fi

# sed command to apply to source file
sedcmd=":t
s%@PACKAGE_NAME@%$LALPackage%;t t
s%@PACKAGE_NAME_UCASE@%$LALPACKAGE%;t t
s%@PACKAGE_NAME_LCASE@%$lalpackage%;t t
s%@PACKAGE_NAME_NOLAL@%$package%;t t
s%@PACKAGE_VCS_INFO_HEADER@%$pkg_vcs_info_header%;t t
s%@PACKAGE_CONFIG_HEADER@%$pkg_config_header%;t t
s%@ID@%$git_id%;t t
s%@DATE@%$git_date%;t t
s%@BRANCH@%$git_branch%;t t
s%@TAG@%$git_tag%;t t
s%@AUTHOR@%$git_author%;t t
s%@COMMITTER@%$git_committer%;t t
s%@STATUS@%$git_status%;t t
s%@BUILDER@%$git_builder%;t t
"

# use sed to perform substitutions; store in a temporary file
test -f "$output.tmp" && rm -f "$output.tmp"
sed "$sedcmd" "$source" > "$output.tmp"

# determine if output file needs to be modified
if test -r "$output" && cmp "$output" "$output.tmp" >/dev/null; then
    rm -f "$output.tmp"
else
    mv -f "$output.tmp" "$output"
    test -n "$am_v_gen" && echo "$am_v_gen $output"
fi
