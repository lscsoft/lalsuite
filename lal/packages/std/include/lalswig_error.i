//
//  Copyright (C) 2011, 2012 Karl Wette
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

// SWIG interface LAL error handling
// Author: Karl Wette

// Only in SWIG interface.
#ifdef SWIG

// Replace default LAL raise/abort hooks with custom handlers which
// translate them into an XLAL error, so that they will be caught
// by the SWIG %exception handler.
%header %{

  #include <signal.h>

  // Print the supplied error message, then raise an XLAL error.
  static int lalswig_LALRaise(int sig, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    (void) vfprintf(stderr, fmt, ap);
    va_end(ap);
    (void) fprintf(stderr, "LALRaise: %s\n", strsignal(sig));
    XLALSetErrno(XLAL_EFAILED);
    return 0;
  }

  // Print the supplied error message, then raise an XLAL error.
  static void lalswig_LALAbort(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    (void) vfprintf(stderr, fmt, ap);
    va_end(ap);
    XLALSetErrno(XLAL_EFAILED);
  }

%}
%init %{
  lalRaiseHook = lalswig_LALRaise;
  lalAbortHook = lalswig_LALAbort;
%}

#endif // SWIG
