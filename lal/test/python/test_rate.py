import doctest

try:
	from igwn_ligolw import ligolw
except ImportError as exc:
	import warnings
	import sys
	warnings.warn(str(exc))
	sys.exit(77)

from lal import rate

if __name__ == '__main__':
	doctest.testmod(rate)
