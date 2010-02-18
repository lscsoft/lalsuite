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
#include <windows.h> // don't move this further up

using namespace std;

/*
  from BOINC
*/
/* from BOINC's filesys.h */
#define FILE_RETRY_INTERVAL 5
extern "C" int boinc_delete_file(const char*);
/* links to BOINC's util.o */
extern double dtime();
extern void boinc_sleep(double);
/* from BOINC's util.h */
static inline double drand() {
    return (double)rand()/(double)RAND_MAX;
}



/*
  sleep() based on Win32 APIs
*/
void sleep(unsigned int s) {
  Sleep(s*1000L);
}


/*
  provides index() as a function for linking LAL with win_lib.o
*/ 
char *index(const char *s, int c) {
  return strchr((char*)s,c);
}


/*
  provides finite() as a function for linking LAL with win_lib.o
*/ 
int finite(double x) {
  return(_finite(x));
}


/*
  for testing FPE
*/
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
  atomic replacement functions for the (currently) non-atomic boinc_rename()
*/

static int eah_rename_aux(const char* old, const char* newf, const bool w2k) {
#ifdef _WIN32
  if(w2k) {
    MoveFileEx(old, newf, MOVEFILE_REPLACE_EXISTING);
  } else {
    CopyFile(old, newf, false);
  }
  return(GetLastError());
#else
  return(rename(old, newf));
#endif
}


int eah_rename(const char* old, const char* newf) {
  int retval = 0;
  static bool w2k = false;

#ifdef _WIN32
  static OSVERSIONINFO osv = {0};

  /* don't know how expensive GetVersionEx() is, better call it only once */
  if (osv.dwOSVersionInfoSize == 0) {
    osv.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    if(GetVersionEx(&osv) == 0) {
      /* this should never happen */
      fprintf(stderr,"WARNING: GetVersionEx() failed (%d)\n", GetLastError());
    } else {
      fprintf(stderr,"INFO: Major Windows version: %d\n", osv.dwMajorVersion);
      w2k = (osv.dwMajorVersion >= 5);
    }
  }
#endif

  /* copied from the original boinc_rename() */
  retval = eah_rename_aux(old,newf,w2k);
  if (retval) {
    double start = dtime();
    do {
      boinc_sleep(drand()*2); // avoid lockstep
      retval = eah_rename_aux(old,newf,w2k);
      if (!retval) break;
    } while (dtime() < start + FILE_RETRY_INTERVAL);
  }

  /* return if there was an error, keeping the old file */
  if(retval)
    return(retval);

#ifdef _WIN32
  /* if we used CopyFile() we still have to delete the old file */
  if(!w2k) {
    retval = boinc_delete_file(old);
    if(retval)
      fprintf(stderr,"WARNING: boinc_delete(%s) failed (%d)\n", old, retval);
  }
#endif

  /* an error while deleting is not fatal */
  return(0);
}
