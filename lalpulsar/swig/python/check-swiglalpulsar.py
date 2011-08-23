# check SWIG Python module wrapping the LALPulsar library
# Author: Karl Wette, 2011

import os

def msg(str):
    print(os.path.basename(__file__) + ": " + str)

# check module load
import swiglalpulsar
from swiglalpulsar import *
msg("passed module load")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
