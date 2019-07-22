#!/usr/bin/env python
#
# Test modules import

import sys
from importlib import import_module
from pathlib import Path

import pytest

# get paths to all python module files
HERE = Path(__file__).parent.absolute()
PYTHON_DIR = HERE.parent.parent / "python"
MODULES = {str(x.relative_to(PYTHON_DIR)) for x in (PYTHON_DIR / "lalinference").rglob("*.py")}

# manually add exclusions here
EXCLUDE = set([
   "lalinference/plot/cylon.py",  # not distributed
   "lalinference/popprior/__init__.py",  # not distributed
   "lalinference/popprior/calculations_rho_voronoi_cheb_edit.py",  # not distributed
   "lalinference/popprior/lagrangian_fit.py",  # not distributed
   "lalinference/popprior/uberbank_database.py",  # not distributed
   "lalinference/tiger/combine_posteriors.py",  # not distributed
])

# build list of module names
INCLUDE = []
for _module in sorted(MODULES - EXCLUDE):
    _name = _module.replace(".py", "").replace(r"/", ".")
    if _name.endswith(".__init__"):
        _name = _name[:-9]
    INCLUDE.append(_name)


@pytest.mark.parametrize("module", INCLUDE)
def test_import(module):
    import_module(module)


if __name__ == "__main__":
    if "-v" not in " ".join(sys.argv[1:]):  # default to verbose
        sys.argv.append("-v")
    sys.exit(pytest.main(args=[__file__] + sys.argv[1:]))
