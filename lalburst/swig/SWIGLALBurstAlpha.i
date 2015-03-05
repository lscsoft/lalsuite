//
// Copyright (C) 2011--2014 Karl Wette
// Copyright (C) 2015 Kipp Cannon
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

///
/// \defgroup SWIGLALBurstAlpha_i Interface SWIGLALBurstAlpha.i
/// \ingroup lalburst_swig
/// \brief SWIG code which must appear \e before the LALBurst headers.
/// \author Karl Wette
///

/// for XLALEPGetTimingParameters()
%apply int *INOUT { int *psd_length };

// Local Variables:
// mode: c
// End:
