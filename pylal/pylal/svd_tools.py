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

def project_waveforms_onto_bases(T, U, S):
	"""
	A function to project waveforms given as rows of T onto the basis
	vectors given as rows of U and normalized by the singular values S.
	This returns the reconstruction coefficients associated with these
	waveforms, basis vectors, and singular values.
	"""
	V = []
	for idx,t in enumerate(T):
		print >> sys.stderr, "\r\ttemplate %i/%i"%(idx+1,len(T)),
		template_copy = numpy.array([t for x in range(len(S))])
		V.append((U*template_copy).sum(axis=1)/S)
	print >> sys.stderr, ''
	return numpy.array(V)

def reconstruct_waveforms_from_SVD(V, S, U, n, normalize=False):
	"""
	A function to reconstruct waveforms from the reconstruction
	coefficients V, singular values S, and basis vectors U. The top n
	basis vectors are used for the reconstruction. If the normalize flag
	is set to True, the waveforms are normalized when before being
	returned.
	"""
	U = numpy.dot(numpy.diag(S), U)[:n,:]
	T = numpy.tensordot(V[:,:n], U, axes=1)
	if normalize:
		for idx,t in enumerate(T):
			T[idx] /= sum(t*t)**.5
	return T
