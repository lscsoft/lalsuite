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
#include <stdio.h>
#include <windows.h> // don't move this earlier

using namespace std;

/* from BOINC */
#define FILE_RETRY_INTERVAL 5
/* links to BOINC's util.o */
extern double dtime();
extern void boinc_sleep(double);
/* from BOINC's util.h */
static inline double drand() {
    return (double)rand()/(double)RAND_MAX;
}

/* weird re-implementation of sleep based on Win32 APIs */
void sleep(unsigned int s) {
#ifdef _MSC_VER
  Sleep(s*1000L);
#else
  _sleep(s);
#endif
}

/*
  provided finite() as a function so that linking LAL
  with win_lib.o works
*/ 
int finite(double x) {
  return(_finite(x));
}

/* for testing FPE */
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

/*
  more atomic replacement functions for the non-atomic boinc_rename()
  eah_rename() is just a copy of boinc_rename() in boinc/lib/filesys.cpp
  with the only difference that it calls eah_rename_aux() instead of
  boinc_rename_aux(). eah_rename_aux() uses MoveFileEx() where
  boinc_rename_aux() uses the boinc_delete(newf); MoveFile(old,newf)
  sequence (not even checking the return code of boinc_delete()).
*/

static int eah_rename_aux(const char* old, const char* newf) {
  int err = 0;
  static OSVERSIONINFO osv = {0};

  /* don't know how expensive GetVersionEx() is, better call it only once */
  if (osv.dwOSVersionInfoSize == 0) {
    osv.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    if(GetVersionEx(&osv) == 0) {
      /* this should never happen. If it does, at least reset the tested MajorVersion */
      fprintf(stderr,"Warning: GetVersionEx() failed (%d)\n",GetLastError());
      osv.dwMajorVersion = 0;
    }
  }

  if(osv.dwMajorVersion >= 5) {

    /* Windows >= Win2k supports MoveFileEx(), which is really atomic */
    if (MoveFileEx(old, newf, MOVEFILE_REPLACE_EXISTING))
      return(0);
    err = GetLastError();

  } else {

    /* copy the new file and then delete the old one should be better than what
       boinc_rename() does, however only for small files like our checkpoint file */
    CopyFile(old,newf,false);
    err = GetLastError();
    if(!err) {
      DeleteFile(old);
      err = GetLastError();
    }
  }

  return(err);
}

/* this is just a copy of boinc_rename() which calls
   eah_rename_aux() instead of boinc_rename_aux() */
int eah_rename(const char* oldf, const char* newf) {
  int retval=0;
  retval = eah_rename_aux(oldf, newf);
  if (retval) {
    double start = dtime();
    do {
      boinc_sleep(drand()*2); // avoid lockstep
      retval = eah_rename_aux(oldf, newf);
      if (!retval) break;
    } while (dtime() < start + FILE_RETRY_INTERVAL);
  }
  return retval;
}
