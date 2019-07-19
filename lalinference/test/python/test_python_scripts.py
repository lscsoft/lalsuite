#!/usr/bin/env python
#
# Test --help for all python scripts

import os
import sys
from pathlib import Path
from subprocess import check_call

import pytest

# get paths to all python module files
HERE = Path(__file__).parent.absolute()
PYTHON_DIR = (HERE.parent.parent / "python").absolute()
SCRIPTS = set()
for pyf in PYTHON_DIR.glob("*.py"):
    # build system creates a shell wrapper around each .py script
    # so we want to actually execute that, this also allows us to
    # only pick up scripts that are to be installed
    shf = Path(str(pyf)[:-3])
    if shf.is_file():
        SCRIPTS.add(str(shf.name))

# list of excluded scripts:
EXCLUDE = {
    "lalinference_burst_pp_pipe",  # requires lalapps
    "lalinference_cpnest",  # requires cpnest
    "lalinference_pipe",  # requires lalapps
    "lalinference_pp_pipe",  # requires lalapps
    "lalinference_tiger_pipe",  # requires lalapps
    "rapidpe_calculate_overlap",  # requires sklearn
    "rapidpe_compute_intrinsic_grid",  # requires sklearn
}


@pytest.mark.parametrize('script', sorted(SCRIPTS-EXCLUDE))
def test_help(script):
    os.chdir(str(PYTHON_DIR))
    check_call("./{} --help".format(script), shell=True)


if __name__ == "__main__":
    if "-v" not in " ".join(sys.argv[1:]):  # default to verbose
        sys.argv.append("-v")
    sys.exit(pytest.main(args=[__file__] + sys.argv[1:]))
