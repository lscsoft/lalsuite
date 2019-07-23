# Test utilities for LALSuite

This directory contains generic testing utilities for LALSuite.

## `test_python_imports.py`

This script discovers all python modules in the `/<subpackage>/python/<subpackage>/`
tree and sets up a unit test to run `import <module>`.

This script must be symlinked to `/<subpackage>/test/python/` and
referenced from the `Makefile.am` in that directory to run properly.

It can be accompanied by `exclude-modules.txt` (in the same directory) that lists out
the _file paths_ of modules to exclude from the test relative to the `/<subpackage>/python/`
directory, e.g. to exclude `lal.git_version` from the tests you would put

```
lal/git_version.py
```

into `exclude-modules.txt`.

## `test_python_scripts.py`

This script discovers all python scripts in the `/<subpackage>/python` directory
and sets up a unit test to run `<script> --help` in a shell.

This script must be symlinked to `/<subpackage>/test/python/` and
referenced from the `Makefile.am` in that directory to run properly.

It can be accompanied by `exclude-scripts.txt` (in the same directory) that lists out
the name of scripts (no directory, no extension) to exclude from the test,
e.g. to exclude `lal_my_awesome_tool` (built from `/lal/python/lal_my_awesome_tool.py`)
from the tests you would put

```
lal_my_awesome_tool
```

into `exclude-scripts.txt`.
