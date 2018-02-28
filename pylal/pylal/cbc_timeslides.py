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

from glue import offsetvector
from pylal import git_version


__version__ = "git id %s" % git_version.id
__date__ = git_version.date


def parse_lalapps_thinca_slidespec(slidespec):
  """
  Accepts a string in the format
  count:instrument=offset[,instrument=offset...] and returns the
  tuple (count, {instrument: offset, ...})
  Example:

  >>> parse_inspiral_num_slides_slidespec("3:H1=0,H2=5,L1=10")
  (3, {'H2': 5.0, 'H1': 0.0, 'L1': 10.0})
  """
  count, offsets = slidespec.strip().split(":")
  tokens = offsets.strip().split(",")
  offsetvect = offsetvector.offsetvector( (instrument.strip(), float(offset)) \
      for instrument, offset in (token.strip().split("=") for token in tokens) )
  return int(count), offsetvect


def Inspiral_Num_Slides_Iter(count, offsets):
  '''
  This generator yields a sequence of time slide dictionaries in the
  style of lalapps_thinca's time slides.  Each resulting dictionary
  maps instrument to offset.  The input is a count of time slides (an
  integer), and a dictionary mapping instrument to offset.  The
  output dictionaries describe time slides that are integer multiples
  of the input time shifts.
  Example (formatted for clarity):

  >>> list(Inspiral_Num_Slides_Iter(3, {"H1": 0.0, "H2": 5.0, "L1": 10.0}))
  [{'H2': -15.0, 'H1': -0.0, 'L1': -30.0},
   {'H2': -10.0, 'H1': -0.0, 'L1': -20.0},
   {'H2': -5.0, 'H1': -0.0, 'L1': -10.0},
   {'H2': 0.0, 'H1': 0.0, 'L1': 0.0},
   {'H2': 5.0, 'H1': 0.0, 'L1': 10.0},
   {'H2': 10.0, 'H1': 0.0, 'L1': 20.0},
   {'H2': 15.0, 'H1': 0.0, 'L1': 30.0}]

  Output time slides are integer multiples of the input time shift 
  vector in the range [-count, +count], including zero, and are 
  returned in increasing order of multiplier.
  '''
  for n in range(-count, +count + 1):
    yield offsetvector.offsetvector( (instrument, offset * n) \
        for instrument, offset in offsets.items() )

