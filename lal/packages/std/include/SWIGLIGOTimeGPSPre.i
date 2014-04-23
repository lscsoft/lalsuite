//
//  Copyright (C) 2011--2014 Karl Wette
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with with program; see the file COPYING. If not, write to the
//  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//  MA  02111-1307  USA
//

// SWIG typemaps, methods and operators for LIGOTimeGPS.
// Author: Karl Wette

// This file contains SWIG code which must appear *before*
// the definition of the LIGOTimeGPS struct.

// Only in SWIG interface.
#if defined(SWIG) && !defined(SWIGXML)

// Specialised input typemaps for LIGOTimeGPS structs.
// Accepts a SWIG-wrapped LIGOTimeGPS or a double as input.
%typemap(in, noblock=1, fragment=SWIG_AsVal_frag(double))
  struct tagLIGOTimeGPS (void *argp = 0, int res = 0),
  const struct tagLIGOTimeGPS (void *argp = 0, int res = 0)
{
  res = SWIG_ConvertPtr($input, &argp, $&descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    double val = 0;
    res = SWIG_AsVal(double)($input, &val);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    } else {
      XLALGPSSetREAL8(&$1, val);
    }
  } else {
    if (!argp) {
      %argument_nullref("$type", $symname, $argnum);
    } else {
      $&ltype temp = %reinterpret_cast(argp, $&ltype);
      $1 = *temp;
      if (SWIG_IsNewObj(res)) {
        %delete(temp);
      }
    }
  }
}
%typemap(freearg) struct tagLIGOTimeGPS, const struct tagLIGOTimeGPS "";

// Specialised input typemaps for pointers to LIGOTimeGPS.
// Accepts a SWIG-wrapped LIGOTimeGPS or a double as input.
%typemap(in, noblock=1, fragment=SWIG_AsVal_frag(double))
  struct tagLIGOTimeGPS* (LIGOTimeGPS tmp, void *argp = 0, int res = 0),
  const struct tagLIGOTimeGPS* (LIGOTimeGPS tmp, void *argp = 0, int res = 0)
{
  res = SWIG_ConvertPtr($input, &argp, $descriptor, $disown | %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    double val = 0;
    res = SWIG_AsVal(double)($input, &val);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    } else {
      $1 = %reinterpret_cast(XLALGPSSetREAL8(&tmp, val), $ltype);
    }
  } else {
    $1 = %reinterpret_cast(argp, $ltype);
  }
}
%typemap(freearg) struct tagLIGOTimeGPS*, const struct tagLIGOTimeGPS* "";

#endif // SWIG && !SWIGXML
