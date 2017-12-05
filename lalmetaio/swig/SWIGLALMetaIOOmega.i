//
// Copyright (C) 2011--2014 Karl Wette
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
/// \defgroup SWIGLALMetaIOOmega_i Interface SWIGLALMetaIOOmega.i
/// \ingroup lalmetaio_swig
/// \brief SWIG code which must appear \e after the LALMetaIO headers.
/// \author Karl Wette
///

///
/// # Specialised wrapping of ::SnglInspiralTable
///

%{
int tagSnglInspiralTable_end_time_get(SnglInspiralTable *self) {
  return self->end.gpsSeconds;
}
void tagSnglInspiralTable_end_time_set(SnglInspiralTable *self, int val) {
  self->end.gpsSeconds = val;
}
int tagSnglInspiralTable_end_time_ns_get(SnglInspiralTable *self) {
  return self->end.gpsNanoSeconds;
}
void tagSnglInspiralTable_end_time_ns_set(SnglInspiralTable *self, int val) {
  self->end.gpsNanoSeconds = val;
}
%}

///
/// Extend the ::SnglInspiralTable class.
%extend tagSnglInspiralTable {
  /// <ul><li>

  /// Export .end integer and nanosecond parts as .end_time and
  /// .end_time_ns for compatibility with glue
  int end_time;
  int end_time_ns;

  /// </li></ul>
}
///

// Local Variables:
// mode: c
// End:
