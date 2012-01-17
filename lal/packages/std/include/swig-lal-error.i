//
//  Copyright (C) 2011 Karl Wette
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

// XLAL/LAL error handlers which generate native scripting-language exceptions
// Author: Karl Wette, 2011

// only if included as a SWIG interface file
#ifdef SWIG

// XLAL error handler which generates native scripting-language exceptions
%header %{

  void swiglal_XLALErrorHandler(const char *func, const char *file, int line, int errnum) {
    // print an error message describing the XLAL error
    XLALPerror(func, file, line, errnum);
    // raise an exception through SWIG, with
    // a message string describing the XLAL error
    (void) SWIG_Error(SWIG_RuntimeError, XLALErrorString(errnum));
  }

%}

// replace default XLAL error handler when loading SWIG module
%init %{
  XLALErrorHandler = swiglal_XLALErrorHandler;
%}

// LAL raise/abort hooks which generate native scripting-language exceptions
%header %{

  #include <signal.h>

  int swiglal_LALRaise(int sig, const char *fmt, ...) {
    // print the error message
    va_list ap;
    va_start(ap, fmt);
    (void) vfprintf(stderr, fmt, ap);
    va_end(ap);
    // raise an exception through SWIG, with
    // a message string describing the signal
    (void) SWIG_Error(SWIG_RuntimeError, strsignal(sig));
    return 0;
  }

  void swiglal_LALAbort(const char *fmt, ...) {
    // print the error message
    va_list ap;
    va_start(ap, fmt);
    (void) vfprintf(stderr, fmt, ap);
    va_end(ap);
    // raise an exception through SWIG
    (void) SWIG_Error(SWIG_RuntimeError, "LAL aborted");
  }

%}

// replace default LAL raise/abort hooks when loading SWIG module
%init %{
  lalRaiseHook = swiglal_LALRaise;
  lalAbortHook = swiglal_LALAbort;
%}

#endif // SWIG
