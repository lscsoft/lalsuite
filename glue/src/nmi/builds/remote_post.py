#!/usr/bin/env python
#
# This script is responsible for packing up a completed build's
# installation dir or, if the build failed, any build artifacts
# needed for debugging.
#
# "Packing up" means simply creating a results.tar.gz file, which
# Metronome will subsequently look for and transfer back to the submit
# host automatically.

import os
import tarfile
from optparse import OptionParser

def OptionParserInit():
    parser = OptionParser(version="%prog $Id$")
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="print verbose status messages to stdout")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose",
                      help="don't print status messages to stdout")
    return parser

# set up basic -v -q options
parser = OptionParserInit()
(options, args) = parser.parse_args()

# if any part of the build failed, we want to pack up src/ for debugging
debugging_files = []
if os.getenv("_NMI_STEP_FAILED") is not None:
    debugging_files = [ "src" ]

tar = tarfile.open("results.tar.gz", "w:gz")
for name in [ "opt" ] + debugging_files:
    tar.add(name)

if options.verbose:
    tar.list(verbose=False)

tar.close()
