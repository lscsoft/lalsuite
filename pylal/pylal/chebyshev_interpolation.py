# Copyright (C) 2011  Drew Keppel
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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

import scipy
import numpy

def U_cheby_nodes_2D(N1, N2):
	"""
	A function to compute the positions (x,y) of the nodes of the
	2-dimensional Chebyshev polynomials of degree (N1,N2).
	"""

	x_cheby = numpy.array([numpy.cos(numpy.pi*(2*idx+1)/N1/2) for idx in range(N1)])
	y_cheby = numpy.array([numpy.cos(numpy.pi*(2*idx+1)/N2/2) for idx in range(N2)])
	return scipy.meshgrid(x_cheby, x_cheby)


def factorial(x):
	"""
	A wrapper for scipy's factorial.
	"""

	return scipy.factorial(x,exact=1)

def T_cheby_2D(x, y, n1, N1, n2, N2):
	"""
	A function that returns the (n,m)th 2-dimensional Chebyshev
	polynomials of degree (N,M) at the locations given in x and y.
	"""

	ux = numpy.zeros(scipy.shape(x))
	uy = numpy.zeros(scipy.shape(y))
	for thisn1 in range(n1/2+1):
		ux += (x**2.-1.)**(thisn1) * x**(n1-2*thisn1) * factorial(n1)/factorial(2*thisn1)/factorial(n1-2*thisn1)
	for thisn2 in range(n2/2+1):
		uy += (y**2.-1.)**(thisn2) * y**(n2-2*thisn2) * factorial(n2)/factorial(2*thisn2)/factorial(n2-2*thisn2)

	w = 1.
	if n1:
		w *= N1/2.
	else:
		w *= N1
	if n2:
		w *= N2/2.
	else:
		w *= N2
	return ux*uy/w**.5

def interpolation_kernel_values(x, y, N1, N2):
	"""
	A function that creates all of the 2-dimensional Chebyshev polynomials
	of degree (N1,N2) evalulated at positions given by (x,y). 
	"""

	Us = []
	for idx in range(N1):
		Us.append([])
		for jdx in range(N2):
			Us[idx].append(T_cheby_2D(x, y, idx, N1, jdx, N2))
	Us = numpy.array(Us)
	return Us

def interpolation_coefficients(func_measurements, Us, N1, N2):
	"""
	A function that takes function measurements at the interpolation
	locations and computes the interpolation coefficients associated with
	the 2-dimensional Chebyshev polynomials of degree (N1,N2) (Us).
	"""

	coeffs = []
	for idx in range(N1):
		coeffs.append([])
		for jdx in range(N2):
			coeffs[idx].append(sum((func_measurements*Us[idx][jdx]).flatten()))
	return coeffs

def reconstruct_interpolation(coeffs, Vs, N1, N2):
	"""
	A function that takes interpolation coefficients associated with
	the 2-dimensional Chebyshev polynomials of degree (N1,N2) and
	recontructs the interpolated function at locations where the
	2-dimensional Chebyshev polynomials of degree (N1,N2) (Vs) have been sampled.
	"""

	func_interpolated = numpy.zeros(numpy.shape(Vs[0,0]))
	for idx in range(N1):
		for jdx in range(N2):
			func_interpolated += coeffs[idx][jdx]*Vs[idx,jdx]
	return func_interpolated
