/*
*  Copyright (C) 2004-2008 Bernd Machenschalk
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* MS-Windows specific functions */

static volatile const char *rcsid_win_lib_cpp = "$Id$";

/* This needs to be compiled in C++ mode ! */

#include "win_lib.h"
#ifdef _MSC_VER
#include <atlbase.h>
#else
#include <math.h>
#include <stdlib.h>
#endif
#include <float.h>
#include <limits>
#include <string.h>

using namespace std;

void sleep(unsigned int s) {
#ifdef _MSC_VER
  Sleep(s*1000L);
#else
  _sleep(s);
#endif
}

int finite(double x) {
  return(_finite(x));
}

float get_float_snan(void) {
  char*idummy;
  return (numeric_limits<float>::signaling_NaN());
  idummy = (char*)rcsid_win_lib_cpp;
  idummy = (char*)rcsid_win_lib_h;
}

/*
  replacement for the asctime_r and gmtime_r functions
  M$ made the asctime and gmtime functions "thread-safe"
  (whatever this means in Windows) instead of adding
  asctime_r and gmtime_r . So for our purpose we implement
  them the trivial way:
*/
struct tm *gmtime_r(const time_t *t, struct tm *s) {
  memcpy(s,gmtime(t),sizeof(struct tm));
  return(s);
}
char *asctime_r(const struct tm *t, char *s) {
  strcpy(s,asctime(t));
  return(s);
}
