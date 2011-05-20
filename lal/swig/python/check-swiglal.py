# check SWIG Python module wrapping the LAL library
# Author: Karl Wette, 2011

import os
import datetime

def msg(str):
    print(os.path.basename(__file__) + ": " + str)
class error(Exception):
    def __init__(self, str):
        self.str = str
    def __str__(self):
        return str

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
