#!/bin/sh

fail() {
  echo "FAIL: $0"
  exit 1
}

srcdir=${2:-.}
dirname=`cd "${srcdir}" && pwd`
test -d "${dirname:-/dev/null/nowhere}" || fail
destdir=`test -d "${3:-"${srcdir}"}" && cd "${3:-"${srcdir}"}" && pwd || echo "${3}"`
outfile="${1:-"-"}"
if test "${outfile}" = "-"; then
  out=''
else
  test -f "${outfile}" && rm -f "${outfile}"
  out=' >> "${outfile}"'
fi
frfiles=`ls "${dirname}"/*.gwf`
test -n "${frfiles}" || fail

IFS='
'
for file in ${frfiles}; do
  basename=`basename "${file}" .gwf`
  IFS_save="$IFS"
  IFS=-
  command='printf "%s %s %s %s "'
  eval $command $basename $out || fail
  IFS="$IFS_save"
  command='echo "file://localhost${destdir}/${basename}.gwf"'
  eval $command $out || fail
done
