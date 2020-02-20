#!/usr/bin/env python
#
# Test --help for all python scripts

import os
import sys
import warnings
try:
    from pathlib import Path
except ImportError as exc:
    warnings.warn(str(exc))
    sys.exit(77)
from subprocess import check_call
try:
    from unittest import mock
except ImportError:  # python < 3
    import mock
    FileNotFoundError = IOError

import pytest

HERE = Path(__file__).parent.absolute()
TOPDIR = HERE.parent.parent
PYTHONDIR = TOPDIR / "python"


def _ignore_symlinks(self):
    """Patch function for py._path.local.LocalPath.realpath

    LALSuite symlinks _this_ script into multiple subpackage directories,
    this function stops pytest from automatically resolving the symlink
    back to the absolute path, since we need this script to sit beside
    the relevant exclude files.
    """
    return self


def read_excludes(source):
    """Read all excluded file paths from the given source file
    """
    excludes = set()
    with open(str(source), "r") as fobj:
        for line in fobj:
            if isinstance(line, bytes):
                line = line.decode("utf-8")
            content = line.strip().split("#", 1)[0].strip()
            if content:
                excludes.add(content)
    return excludes


def find_scripts(path):
    scripts = set()
    for pyf in path.glob("*.py"):
        # build system creates a shell wrapper around each .py script
        # so we want to actually execute that, this also allows us to
        # only pick up scripts that are to be installed
        shf = Path(str(pyf)[:-3])
        if shf.is_file():
            scripts.add(str(shf.name))
    return scripts


SCRIPTS = find_scripts(PYTHONDIR)
try:
    EXCLUDE = read_excludes(HERE / "exclude-scripts.txt")
except FileNotFoundError:  # no exclusion file
    EXCLUDE = set()


@pytest.mark.parametrize('script', sorted(SCRIPTS))
def test_help(script):
    if script in EXCLUDE:
        pytest.skip("excluded {}".format(str(script)))
    os.chdir(str(PYTHONDIR))
    check_call("./{} --help".format(script), shell=True)


# run from command line
if __name__ == "__main__":
    args = sys.argv[1:] or [
        "-v",
        "-rs",
        "--junit-xml=junit-scripts.xml",
    ]
    sys.exit(pytest.main(args=[__file__] + args))
    # run pytest with patch to not resolve symlinks
    with mock.patch("py._path.local.LocalPath.realpath", _ignore_symlinks):
        sys.exit(pytest.main(args=[__file__] + args))
