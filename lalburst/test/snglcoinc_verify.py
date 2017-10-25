#!/usr/bin/env python

import doctest, sys
from lalburst import snglcoinc

if __name__ == '__main__':
	sys.exit(bool(doctest.testmod(snglcoinc).failed))
