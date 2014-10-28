#!/usr/bin/env python
#
# Copyright (C) 2014  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import sys
try:
	import lal
except ImportError:
	print "WARNING:  lal module not found.  test disabled"
	sys.exit(0)
try:
	lal.MarcumQmodified
except AttributeError:
	print "WARNING:  lal module does not correspond to this source tree.  test disabled"
	sys.exit(0)


#
# =============================================================================
#
#                                  Test Data
#
# =============================================================================
#


#
# test data copied from testmarq.f90 in reference implementation by Gil et
# al. published as part of
#
# A. Gil, J. Segura, and N. M. Temme.  Algorithm 939:  Computation of the
# Marcum Q-Function.  ACM Transactions on Mathematical Software (TOMS),
# Volume 40 Issue 3, April 2014, Article No. 20. arXiv:1311.0681
#

mu = [1.0, 3.0, 4.0, 6.0, 8.0, 10.0, 20.0, 22.0, 25.0, 27.0, 30.0, 32.0, 40.0, 50.0, 200.0, 350.0, 570.0, 1000.0]

x = [0.3, 2.0, 8.0, 25.0, 13.0, 45.0, 47.0, 100.0, 85.0, 120.0, 130.0, 140.0, 30.0, 40.0, 0.01, 100.0, 1.0, 0.08]

y = [0.01, 0.1, 50.0, 10.0, 15.0, 25.0, 30.0, 150.0, 60.0, 205.0, 90.0, 100.0, 120.0, 150.0, 190.0, 320.0, 480.0, 799.0]

p = [.7382308441994e-02, .2199222796783e-04, .9999999768807, .1746869995977e-03, .1483130042637, .1748328323235e-03, .1340769184710e-04, .9646591215441, .1783991673043e-04, .9994542406431, .1220231636641e-05, .1757487653307e-05, .9999894753719, .9999968347378, .2431297758708, .3851423018735e-09, .2984493152360e-04, .4191472999694e-11]

q = [.9926176915580, .9999780077720, .2311934913546e-07, .9998253130004, .8516869957363, .9998251671677, .9999865923082, .3534087845586e-01, .9999821600833, .5457593568564e-03, .9999987797684, .9999982425123, .1052462813144e-04, .3165262228904e-05, .7568702241292, .9999999996149, .9999701550685, .9999999999958]

for mu, x, y, p, q in zip(mu, x, y, p, q):
	lal_q = lal.MarcumQmodified(mu, x, y)
	rel_err_q = abs(q - lal_q) / q
	print "Q(%g,%g,%g) = %.16g\t%.16g\t%.16g" % (mu, x, y, q, lal_q, rel_err_q)
	assert rel_err_q < 1e-12


#
# =============================================================================
#
#                                 Extra Tests
#
# =============================================================================
#


#
# exercise specific code paths
#


# series expansion in section 3, x < 30
print
print lal.MarcumQmodified(1., 18, 8.)

# asymptotic expansion in section 4.1, xsi > 30 && M*M < 2 * xsi
# for x < y
print
print lal.MarcumQmodified(1., 32., 40.5)
# for x == y
print lal.MarcumQmodified(1., 32., 32.)
# for y > x
print lal.MarcumQmodified(1., 32., 8.)

# recurrence relation in (14), f1 < y < f2 && M < 135.
print
print lal.MarcumQmodified(50., 32., 80.)

# asymptotic expansion in section 4.2, f1 < y < f2 && M >= 135.
#print
#print lal.MarcumQmodified(150., 32., 180.)

# quadrature integration in section 5
print
print lal.MarcumQmodified(50., 32., 98.)

# more tests
print
print lal.MarcumQmodified(1., 31.99999, 32.00001), lal.MarcumQmodified(1., 32.00001, 31.99999)
print lal.MarcumQmodified(1., 31.9999999999, 32.0000000001), lal.MarcumQmodified(1., 32.0000000001, 31.9999999999)
print lal.MarcumQmodified(1., 32., 32.), lal.MarcumQmodified(1., 32., 32.)
