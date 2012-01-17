# check SWIG Python module wrapping the LALFrame library
# Author: Karl Wette, 2011

import os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# check module load
import swiglalframe
from swiglalframe import *
msg("passed module load")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
