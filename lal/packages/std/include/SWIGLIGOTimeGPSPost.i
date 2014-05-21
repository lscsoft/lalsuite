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

// This file contains SWIG code which must appear *after*
// the definition of the LIGOTimeGPS struct.

// Only in SWIG interface.
#if defined(SWIG) && !defined(SWIGXML)

// Allocate a new LIGOTimeGPS.
#define %swiglal_new_LIGOTimeGPS() (LIGOTimeGPS*)(XLALCalloc(1, sizeof(LIGOTimeGPS)))

// Extend the LIGOTimeGPS class.
%extend tagLIGOTimeGPS {

  // Construct a new LIGOTimeGPS from a real number.
  tagLIGOTimeGPS(REAL8 t) {
    return XLALGPSSetREAL8(%swiglal_new_LIGOTimeGPS(), t);
  }

  // Construct a new LIGOTimeGPS from integer seconds and nanoseconds.
  tagLIGOTimeGPS(INT4 gpssec) {
    return XLALGPSSet(%swiglal_new_LIGOTimeGPS(), gpssec, 0);
  }
  tagLIGOTimeGPS(INT4 gpssec, INT8 gpsnan) {
    return XLALGPSSet(%swiglal_new_LIGOTimeGPS(), gpssec, gpsnan);
  }

  // Construct a new LIGOTimeGPS from a string
  tagLIGOTimeGPS(const char *str) {
    LIGOTimeGPS *gps = %swiglal_new_LIGOTimeGPS();
    char *end = NULL;
    if (XLALStrToGPS(gps, str, &end) < 0 || *end != '\0') {
      XLALFree(gps);
      xlalErrno = XLAL_EFUNC;   // Silently signal an error to constructor
      return NULL;
    }
    return gps;
  }

  // Operators are implemented by defining Python-style __operator__
  // methods (since LAL is C99, we don't have C++ operators available).
  // Many SWIG language modules will automatically map these functions
  // to scripting-language operations in their runtime code. In some
  // cases (ironically, Python itself when using -builtin!) additional
  // directives are needed in the scripting-language-specific interface.

  // Return new LIGOTimeGPSs which are the positive and negative values of $self.
  LIGOTimeGPS* __pos__() {
    return XLALINT8NSToGPS(%swiglal_new_LIGOTimeGPS(), +XLALGPSToINT8NS($self));
  }
  LIGOTimeGPS* __neg__() {
    return XLALINT8NSToGPS(%swiglal_new_LIGOTimeGPS(), -XLALGPSToINT8NS($self));
  }

  // Return a new LIGOTimeGPS which is the absolute value of $self.
  LIGOTimeGPS* __abs__() {
    return XLALINT8NSToGPS(%swiglal_new_LIGOTimeGPS(), llabs(XLALGPSToINT8NS($self)));
  }

  // Return whether a LIGOTimeGPS is non-zero.
  bool __nonzero__() {
    return $self->gpsSeconds || $self->gpsNanoSeconds;
  }

  // Return integer representations of the LIGOTimeGPS seconds.
  int __int__() {
    return $self->gpsSeconds;
  }
  long __long__() {
    return $self->gpsSeconds;
  }

  // Return a floating-point representation of a LIGOTimeGPS.
  double __float__() {
    return XLALGPSGetREAL8($self);
  }

  // Return a string representation of a LIGOTimeGPS. Because
  // XLALGPSToStr() allocates a new string using LAL memory,
  // %newobject is used to make SWIG use a 'newfree' typemap,
  // where the string is freed; SWIG will have already copied it
  // to a native scripting-language string to return as output.
  %newobject __str__;
  %typemap(newfree) char* __str__ "XLALFree($1);";
  char* __str__() {
    return XLALGPSToStr(NULL, $self);
  }
  %newobject __repr__;
  %typemap(newfree) char* __repr__ "XLALFree($1);";
  char* __repr__() {
    return XLALGPSToStr(NULL, $self);
  }

  // Return the hash value of a LIGOTimeGPS.
  long __hash__() {
    long hash = (long)$self->gpsSeconds ^ (long)$self->gpsNanoSeconds;
    return hash == -1 ? -2 : hash;
  }

  // Binary operators need only be implemented for two LIGOTimeGPS
  // arguments; the specialised input typemaps defined above will
  // then allow other suitable types to be substituted as arguments.

  // Return the addition of two LIGOTimeGPSs
  LIGOTimeGPS* __add__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSAddGPS(retn, gps);
  }
  LIGOTimeGPS* __radd__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSAddGPS(retn, $self);
  }

  // Return the subtraction of two LIGOTimeGPSs
  LIGOTimeGPS* __sub__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSSetREAL8(retn, XLALGPSDiff(retn, gps));
  }
  LIGOTimeGPS* __rsub__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSSetREAL8(retn, XLALGPSDiff(retn, $self));
  }

  // Return the multiplication of two LIGOTimeGPSs
  LIGOTimeGPS* __mul__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSMultiply(retn, XLALGPSGetREAL8(gps));
  }
  LIGOTimeGPS* __rmul__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSMultiply(retn, XLALGPSGetREAL8($self));
  }

  // Return the floating-point division of two LIGOTimeGPSs
  LIGOTimeGPS* __div__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSDivide(retn, XLALGPSGetREAL8(gps));
  }
  LIGOTimeGPS* __rdiv__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSDivide(retn, XLALGPSGetREAL8($self));
  }

  // Return the integer division of two LIGOTimeGPSs
  LIGOTimeGPS* __floordiv__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSSetREAL8(retn, floor(XLALGPSGetREAL8(XLALGPSDivide(retn, XLALGPSGetREAL8(gps)))));
  }
  LIGOTimeGPS* __rfloordiv__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSSetREAL8(retn, floor(XLALGPSGetREAL8(XLALGPSDivide(retn, XLALGPSGetREAL8($self)))));
  }

  // Return the modulus of two LIGOTimeGPSs
  LIGOTimeGPS* __mod__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *$self;
    return XLALGPSSetREAL8(retn, fmod(XLALGPSGetREAL8(retn), XLALGPSGetREAL8(gps)));
  }
  LIGOTimeGPS* __rmod__(LIGOTimeGPS* gps) {
    LIGOTimeGPS* retn = %swiglal_new_LIGOTimeGPS();
    *retn = *gps;
    return XLALGPSSetREAL8(retn, fmod(XLALGPSGetREAL8(retn), XLALGPSGetREAL8($self)));
  }

  // Comparison operators between two LIGOTimeGPSs are generated by
  // the following SWIG macro. NAME is the name of the Python operator
  // and OP is the C operator. The correct comparison NAME is obtained
  // by comparing the result of XLALGPSCmp() against zero using OP.
  %define %swiglal_LIGOTimeGPS_comparison_operator(NAME, OP)
    bool __##NAME##__(LIGOTimeGPS *gps) {
      return XLALGPSCmp($self, gps) OP 0;
    }
  %enddef
  %swiglal_LIGOTimeGPS_comparison_operator(lt, < );
  %swiglal_LIGOTimeGPS_comparison_operator(le, <=);
  %swiglal_LIGOTimeGPS_comparison_operator(eq, ==);
  %swiglal_LIGOTimeGPS_comparison_operator(ne, !=);
  %swiglal_LIGOTimeGPS_comparison_operator(gt, > );
  %swiglal_LIGOTimeGPS_comparison_operator(ge, >=);

  // Return the number of nanoseconds in a LIGOTimeGPS
  INT8 ns() {
    return XLALGPSToINT8NS($self);
  }

} // %extend tagLIGOTimeGPS

#endif // SWIG && !SWIGXML
