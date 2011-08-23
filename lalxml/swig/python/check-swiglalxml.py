# check SWIG Python module wrapping the LALXML library
# Author: Karl Wette, 2011

import os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# check module load
import swiglalxml
from swiglalxml import *
msg("passed module load")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
