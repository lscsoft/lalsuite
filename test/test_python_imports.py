#!/usr/bin/env python
#
# Test modules import

import sys
import warnings
from importlib import import_module
try:
    from pathlib import Path
except ImportError as exc:
    warnings.warn(str(exc))
    sys.exit(77)
try:
    from unittest import mock
except ImportError:  # python < 3
    import mock
    FileNotFoundError = IOError

import pytest

HERE = Path(__file__).parent.absolute()
TOPDIR = HERE.parent.parent
SUBPACKAGE = TOPDIR.name
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


def find_modules(path):
    """Returns the paths to all python module files
    """
    return {str(x.relative_to(path.parent)) for x in path.rglob("*.py")}


def path_to_name(filepath):
    name = filepath.replace(".py", "").replace(r"/", ".")
    if name.endswith(".__init__"):
        name = name[:-9]
    return name


FILENAMES = find_modules(PYTHONDIR / SUBPACKAGE)
try:
    EXCLUDE = read_excludes(HERE / "exclude-modules.txt")
except FileNotFoundError:  # no exclusion file
    EXCLUDE = set()
MODULES = [path_to_name(x) for x in sorted(FILENAMES)]
EXCLUDE = set(map(path_to_name, EXCLUDE))


@pytest.mark.parametrize("module", MODULES)
def test_import(module):
    if module in EXCLUDE:
        pytest.skip("excluded {}".format(str(module)))
    import_module(module)


# run from command-line
if __name__ == "__main__":
    if "-v" not in " ".join(sys.argv[1:]):  # default to verbose
        sys.argv.append("-v")
    sys.argv.append("-rs")
    # run pytest with patch to not resolve symlinks
    with mock.patch("py._path.local.LocalPath.realpath", _ignore_symlinks):
        sys.exit(pytest.main(args=[__file__] + sys.argv[1:]))
