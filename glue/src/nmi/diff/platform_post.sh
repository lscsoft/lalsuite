#!/bin/sh

# exit immediately if any command exits with a non-zero status.
set -e
# treat unset variables as an error when performing parameter expansion.
set -u
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

# unpack XML diff results for review
tar zxvf results.tar.gz \*.xml.diff
