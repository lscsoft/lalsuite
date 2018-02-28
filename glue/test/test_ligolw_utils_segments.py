#!/usr/bin/env python

import doctest
from glue.ligolw.utils import segments as ligolw_segments

if __name__ == '__main__':
	doctest.testmod(ligolw_segments)
