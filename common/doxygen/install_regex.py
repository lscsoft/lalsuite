# install_regex.py - build regex for installing Doxygen documentation

__author__ = "Karl Wette <karl.wette@ligo.org>"
__copyright__ = "Copyright (C) 2015 Karl Wette"

import os
import sys


# print error message and exit
def fail(msg):
    sys.stderr.write("%s: %s\n" % (sys.argv[0], msg))
    sys.exit(1)


# get input arguments
install_dir, install_dirmap = sys.argv[1:]

# install regex built from install directory map
# - make 'to_dir' relative to install_dir to ensure
#   Doxygen documentation does not contain absolute paths
#   and hence is relocatable
for elem in install_dirmap.split():
    (from_dir, to_dir) = elem.split(":")
    if len(from_dir) == 0:
        fail("from-directory in install directory map is empty")
    if len(to_dir) == 0:
        fail("to-directory in install directory map is empty")
    print("s|%s|%s|g" % (from_dir, os.path.relpath(to_dir, install_dir)))

# install regex for top-level navigation tabs
print(r"s|\.\./\.\./\.\./\(lal[a-z]*\)/doxygen/out/index\.html|../\1/index.html|g")
print(r"s|\.\./\.\./\.\./doxygen/out/index\.html|../lalsuite/index.html|g")
print(r"s|\.\./\.\./\(lal[a-z]*\)/doxygen/out/index\.html|../\1/index.html|g")

# output
print("p")
