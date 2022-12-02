/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix, Holger Pletsch
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include "config.h"

#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include <lal/StringInput.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>
#include <lal/HeapToplist.h>
#include <lal/LALPulsarVCSInfo.h>

#include "GCTtoplist.h"

/* Windows specifics */
#ifdef _WIN32

/* errno */
#include <errno.h>
#ifndef _doserrno
extern int _doserrno;
#endif

/* snprintf */
#define snprintf _snprintf

/* fsync */
#define fsync _commit
#define fileno _fileno

/* finite */
#include <float.h>
#define finite _finite

#else /* WIN32 */

/* errno */
#include <errno.h>

/* this is defined in C99 and *should* be in math.h. Long term
   protect this with a HAVE_FINITE */
int finite(double);

#endif /* WIN32 */

#ifdef DEBUG_SORTING
static FILE*debugfp = NULL;
#endif

/* define min macro if not already defined */
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* maximum number of succesive failures before switching off syncing */
#define SYNC_FAIL_LIMIT 5

/* output file column headings string, globally defined from HSGCT program */
extern char *global_column_headings_stringp;

/* local prototypes */
static void reduce_gctFstat_toplist_precision(toplist_t *l);
static int _atomic_write_gctFstat_toplist_to_file(toplist_t *l, const char *filename, UINT4*checksum, int write_done);
static int print_gctFstatline_to_str(GCTtopOutputEntry fline, char* buf, int buflen);
static int print_single_detector_quantities_to_str(char *outstr, size_t outstrlen, const REAL4 *quantities, const UINT4 numDetectors);
static int print_single_detector_intval_to_str( char *outstr, size_t outstrlen, const INT4 *quantities, const UINT4 numDetectors );
static int write_gctFstat_toplist_item_to_fp(GCTtopOutputEntry fline, FILE*fp, UINT4*checksum);

/* ordering function for sorting the list */
static int gctFstat_result_order(const void *a, const void *b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20lf  %20lf\n%20lf  %20lf\n%20lf  %20lf\n\n",
	    ((const GCTtopOutputEntry*)a)->Freq,  ((const GCTtopOutputEntry*)b)->Freq,
	    ((const GCTtopOutputEntry*)a)->Alpha, ((const GCTtopOutputEntry*)b)->Alpha,
	    ((const GCTtopOutputEntry*)a)->Delta, ((const GCTtopOutputEntry*)b)->Delta,
	    ((const GCTtopOutputEntry*)a)->F1dot, ((const GCTtopOutputEntry*)b)->F1dot);
#endif
  if      (((const GCTtopOutputEntry*)a)->Freq  < ((const GCTtopOutputEntry*)b)->Freq)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->Freq  > ((const GCTtopOutputEntry*)b)->Freq)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->Alpha < ((const GCTtopOutputEntry*)b)->Alpha)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->Alpha > ((const GCTtopOutputEntry*)b)->Alpha)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->Delta < ((const GCTtopOutputEntry*)b)->Delta)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->Delta > ((const GCTtopOutputEntry*)b)->Delta)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->F1dot < ((const GCTtopOutputEntry*)b)->F1dot)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->F1dot > ((const GCTtopOutputEntry*)b)->F1dot)
    return 1;
  else
    return 0;
}

/* ordering function defining the toplist: SORT BY avTwoF */
static int gctFstat_smaller(const void*a, const void*b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20u  %20u\n\n",
	    ((const GCTtopOutputEntry*)a)->avTwoF,  ((const GCTtopOutputEntry*)b)->avTwoF,
	    ((const GCTtopOutputEntry*)a)->nc, ((const GCTtopOutputEntry*)b)->nc);
#endif
  if      (((const GCTtopOutputEntry*)a)->avTwoF < ((const GCTtopOutputEntry*)b)->avTwoF)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF > ((const GCTtopOutputEntry*)b)->avTwoF)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->nc < ((const GCTtopOutputEntry*)b)->nc)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->nc > ((const GCTtopOutputEntry*)b)->nc)
    return -1;
  else
    return(gctFstat_result_order(a,b));
}


/* ordering function defining the toplist: SORT BY Numbercount */
static int gctNC_smaller(const void*a, const void*b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20u  %20u\n\n",
	    ((const GCTtopOutputEntry*)a)->avTwoF,  ((const GCTtopOutputEntry*)b)->avTwoF,
	    ((const GCTtopOutputEntry*)a)->nc, ((const GCTtopOutputEntry*)b)->nc);
#endif
  if      (((const GCTtopOutputEntry*)a)->nc < ((const GCTtopOutputEntry*)b)->nc)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->nc > ((const GCTtopOutputEntry*)b)->nc)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF < ((const GCTtopOutputEntry*)b)->avTwoF)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF > ((const GCTtopOutputEntry*)b)->avTwoF)
    return -1;
  else
    return(gctFstat_result_order(a,b));
}


/* ordering function defining the toplist: SORT BY BSGL */
static int gctBSGL_smaller(const void*a, const void*b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20lf  %20lf\n\n",
	    ((const GCTtopOutputEntry*)a)->avTwoF,  ((const GCTtopOutputEntry*)b)->avTwoF,
	    ((const GCTtopOutputEntry*)a)->log10BSGL, ((const GCTtopOutputEntry*)b)->log10BSGL);
#endif
  if      (((const GCTtopOutputEntry*)a)->log10BSGL < ((const GCTtopOutputEntry*)b)->log10BSGL)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->log10BSGL > ((const GCTtopOutputEntry*)b)->log10BSGL)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF < ((const GCTtopOutputEntry*)b)->avTwoF)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF > ((const GCTtopOutputEntry*)b)->avTwoF)
    return -1;
  else
    return(gctFstat_result_order(a,b));
}
/* ordering function defining the toplist: SORT BY BSGLtL */
static int gctBSGLtL_smaller(const void*a, const void*b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20lf  %20lf\n\n",
	    ((const GCTtopOutputEntry*)a)->avTwoF,  ((const GCTtopOutputEntry*)b)->avTwoF,
	    ((const GCTtopOutputEntry*)a)->log10BSGLtL, ((const GCTtopOutputEntry*)b)->log10BSGLtL);
#endif
  if      (((const GCTtopOutputEntry*)a)->log10BSGLtL < ((const GCTtopOutputEntry*)b)->log10BSGLtL)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->log10BSGLtL > ((const GCTtopOutputEntry*)b)->log10BSGLtL)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF < ((const GCTtopOutputEntry*)b)->avTwoF)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF > ((const GCTtopOutputEntry*)b)->avTwoF)
    return -1;
  else
    return(gctFstat_result_order(a,b));
}
/* ordering function defining the toplist: SORT BY BtSGLtL */
static int gctBtSGLtL_smaller(const void*a, const void*b) {
#ifdef DEBUG_SORTING
  if(debugfp)
    fprintf(debugfp,"%20lf  %20lf\n%20lf  %20lf\n\n",
	    ((const GCTtopOutputEntry*)a)->avTwoF,  ((const GCTtopOutputEntry*)b)->avTwoF,
	    ((const GCTtopOutputEntry*)a)->log10BtSGLtL, ((const GCTtopOutputEntry*)b)->log10BtSGLtL);
#endif
  if      (((const GCTtopOutputEntry*)a)->log10BtSGLtL < ((const GCTtopOutputEntry*)b)->log10BtSGLtL)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->log10BtSGLtL > ((const GCTtopOutputEntry*)b)->log10BtSGLtL)
    return -1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF < ((const GCTtopOutputEntry*)b)->avTwoF)
    return 1;
  else if (((const GCTtopOutputEntry*)a)->avTwoF > ((const GCTtopOutputEntry*)b)->avTwoF)
    return -1;
  else
    return(gctFstat_result_order(a,b));
}


/* functions for qsort based on the above ordering functions */
static int gctFstat_final_qsort(const void*a, const void*b) {
  void const* const* pa = (void const* const*)a;
  void const* const* pb = (void const* const*)b;
  return(gctFstat_result_order(*pa,*pb));
}


/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
int create_gctFstat_toplist(toplist_t**tl, UINT8 length, UINT4 whatToSortBy) {
#ifdef DEBUG_SORTING
  if(!debugfp)
    debugfp=fopen("debug_sort","w");
#endif

  if (whatToSortBy==SORTBY_NC) {
    return( create_toplist(tl, length, sizeof(GCTtopOutputEntry), gctNC_smaller) );
  }
  else if (whatToSortBy==SORTBY_BSGL) {
    return( create_toplist(tl, length, sizeof(GCTtopOutputEntry), gctBSGL_smaller) );
  }
  else if (whatToSortBy==SORTBY_BSGLtL) {
    return( create_toplist(tl, length, sizeof(GCTtopOutputEntry), gctBSGLtL_smaller) );
  }
  else if (whatToSortBy==SORTBY_BtSGLtL) {
    return( create_toplist(tl, length, sizeof(GCTtopOutputEntry), gctBtSGLtL_smaller) );
  }
  else {
    return( create_toplist(tl, length, sizeof(GCTtopOutputEntry), gctFstat_smaller) );
  }

}

/* frees the space occupied by the toplist
   NOTE: toplist must not contain any allocated structs */
void free_gctFstat_toplist(toplist_t**l) {
  free_toplist(l);
} /* free_gctFstat_toplist() */


/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist and
   look for the now smallest one.
   Returns 1 if the element was actually inserted, 0 if not. */
int insert_into_gctFstat_toplist(toplist_t*tl, GCTtopOutputEntry * elem) {
  if ( !tl )
    return 0;
  else
    return(insert_into_toplist(tl, (void*) elem));
}

/* (q)sort the toplist according to the sorting function. */
void sort_gctFstat_toplist(toplist_t*l) {
  qsort(l->heap,l->elems,sizeof(char*),gctFstat_final_qsort);
}

/* Prints a Toplist line to a string buffer.
   Separate function to assure consistency of output and reduced precision for sorting */
static int print_gctFstatline_to_str(GCTtopOutputEntry fline, char* buf, int buflen) {

  /* add extra output field for line-robust statistic BSGL */
  char BSGLstr[256] = "";	/* defaults to empty */
  if ( fline.log10BSGL > -LAL_REAL4_MAX*0.2 ) /* if --computeBSGL=FALSE, the log10BSGL field was initialised to -LAL_REAL4_MAX; if --computeBSGL=TRUE, it is at least -LAL_REAL4_MAX*0.1 */
    {
      snprintf ( BSGLstr, sizeof(BSGLstr), " %.6f", fline.log10BSGL );
      print_single_detector_quantities_to_str ( BSGLstr, sizeof(BSGLstr), fline.avTwoFX, fline.numDetectors );
    } /* if fline.log10BSGL */

  /* add extra output fields for line-robust statistic BSGLtL and BtSGLtL */
  char BSGLtLstr[256] = "";	/* defaults to empty */
  if ( fline.log10BSGLtL > -LAL_REAL4_MAX*0.2 ) {
    snprintf ( BSGLtLstr, sizeof(BSGLtLstr), " %.6f", fline.log10BSGLtL );
  }
  char BtSGLtLstr[256] = "";	/* defaults to empty */
  if ( fline.log10BtSGLtL > -LAL_REAL4_MAX*0.2 ) {
    snprintf ( BtSGLtLstr, sizeof(BtSGLtLstr), " %.6f", fline.log10BtSGLtL );
  }

  /* add extra output fields for max2F over segments */
  char maxTwoFstr[256] = "";	/* defaults to empty */
  if ( fline.maxTwoFl >= 0.0 ) /* this was initialised to -1.0 and is only >= 0.0 if actually computed in inner loop */
    {
      snprintf ( maxTwoFstr, sizeof(maxTwoFstr), " %.6f %d", fline.maxTwoFl, fline.maxTwoFlSeg );
      print_single_detector_quantities_to_str ( maxTwoFstr, sizeof(maxTwoFstr), fline.maxTwoFXl, fline.numDetectors );
      print_single_detector_intval_to_str ( maxTwoFstr, sizeof(maxTwoFstr), fline.maxTwoFXlSeg, fline.numDetectors );
    }

  /* add extra output fields for recalculated statistics */
  char recalcStr[256] = "";	/* defaults to empty */
  if ( fline.avTwoFrecalc >= 0.0 ) /* this was initialised to -1.0 and is only >= 0.0 if actually recomputed in recalcToplistStats step */
    {
      char buf0[256];
      snprintf ( recalcStr, sizeof(recalcStr), " %.6f", fline.avTwoFrecalc );
      if ( fline.log10BSGLrecalc > -LAL_REAL4_MAX*0.2 ) /* if --computeBSGL=FALSE, the log10BSGLrecalc field was initialised to -LAL_REAL4_MAX; if --computeBSGL=TRUE, it is at least -LAL_REAL4_MAX*0.1 */
        {
          snprintf ( buf0, sizeof(buf0), " %.6f", fline.log10BSGLrecalc );
          strcat ( recalcStr, buf0 );
        } /* if ( fline.log10BSGL > -LAL_REAL4_MAX*0.2 ) */
      print_single_detector_quantities_to_str ( recalcStr, sizeof(recalcStr), fline.avTwoFXrecalc, fline.numDetectors );

      if ( fline.twoFloudestSeg >= 0.0 ) /* this was initialised to -1.0 and is only >= 0.0 if actually recomputed in recalcToplistStats step */
      {
        snprintf ( buf0, sizeof(buf0), " %d %.6f", fline.loudestSeg, fline.twoFloudestSeg );
        strcat ( recalcStr, buf0 );
        print_single_detector_quantities_to_str ( recalcStr, sizeof(recalcStr), fline.twoFXloudestSeg, fline.numDetectors );
      } /* if ( fline.twoFloudestSeg >= 0.0 ) */

      if ( fline.log10BSGLtLrecalc > -LAL_REAL4_MAX*0.2 ) /* if --computeBSGL=FALSE, the log10BSGLtLrecalc field was initialised to -LAL_REAL4_MAX; if --computeBSGL=TRUE, it is at least -LAL_REAL4_MAX*0.1 */
        {
          snprintf ( buf0, sizeof(buf0), " %.6f", fline.log10BSGLtLrecalc );
          strcat ( recalcStr, buf0 );
        } /* if ( fline.log10BSGLtL > -LAL_REAL4_MAX*0.2 ) */

    } /* if avTwoFrecalc */

  int len;
  if (fline.have_f3dot){
  len = snprintf(buf, buflen,
                     "%.16g %.13g %.13g %.13g %.13g %.13g %d %.6f%s%s%s%s%s\n",
                     fline.Freq,
                     fline.Alpha,
                     fline.Delta,
                     fline.F1dot,
                     fline.F2dot,
                     fline.F3dot,
                     fline.nc,
                     fline.avTwoF,
                     BSGLstr,
                     BSGLtLstr,
                     BtSGLtLstr,
                     maxTwoFstr,
                     recalcStr
                 );
  }
  else {
  len = snprintf(buf, buflen,
                     "%.16g %.13g %.13g %.13g %.13g %d %.6f%s%s%s%s%s\n",
                     fline.Freq,
                     fline.Alpha,
                     fline.Delta,
                     fline.F1dot,
                     fline.F2dot,
                     fline.nc,
                     fline.avTwoF,
                     BSGLstr,
                     BSGLtLstr,
                     BtSGLtLstr,
                     maxTwoFstr,
                     recalcStr
                 );
}
  return len;

} /* print_gctFstatline_to_str() */

static int print_single_detector_quantities_to_str ( char *outstr, size_t outstrlen, const REAL4 *quantities, const UINT4 numDetectors ) {
  const char *fn = __func__;

  UINT4 len = strlen ( outstr );

  char buf0[outstrlen];
  for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
    snprintf ( buf0, sizeof(buf0), " %.6f", quantities[X] );
    len = strlen ( outstr ) + strlen ( buf0 ) + 1;
    if ( len > outstrlen ) {
      XLALPrintError ("%s: assembled output string too long! (%d > %zu)\n", fn, len, outstrlen);
      break;	/* we can't really terminate with error in this function, but at least we avoid crashing */
    }
    strcat ( outstr, buf0 );
  } /* for X < numDet */

  return len;

} /* print_single_detector_quantities_to_str() */


static int print_single_detector_intval_to_str ( char *outstr, size_t outstrlen, const INT4 *quantities, const UINT4 numDetectors ) {
  const char *fn = __func__;

  UINT4 len = strlen ( outstr );

  char buf0[outstrlen];
  for ( UINT4 X = 0; X < numDetectors ; X ++ ) {
    snprintf ( buf0, sizeof(buf0), " %d", quantities[X] );
    len = strlen ( outstr ) + strlen ( buf0 ) + 1;
    if ( len > outstrlen ) {
      XLALPrintError ("%s: assembled output string too long! (%d > %zu)\n", fn, len, outstrlen);
      break;    /* we can't really terminate with error in this function, but at least we avoid crashing */
    }
    strcat ( outstr, buf0 );
  } /* for X < numDet */

  return len;

} /* print_single_detector_quantities_to_str() */





/* writes an GCTtopOutputEntry line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given */
static int write_gctFstat_toplist_item_to_fp(GCTtopOutputEntry fline, FILE*fp, UINT4*checksum) {
  char linebuf[512];
  UINT4 i;

  UINT4 length = print_gctFstatline_to_str(fline, linebuf, sizeof(linebuf)-1);

  if(length>sizeof(linebuf)-1) {
    return -1;
  }

  if (checksum)
    for(i=0;i<length;i++)
      *checksum += linebuf[i];

  XLAL_LAST_ELEM(linebuf) = '\0';

  return(fprintf(fp,"%s",linebuf));
}



/* Reduces the precision of all elements in the toplist which influence the sorting order.
   To be called before sorting and finally writing out the list */
static void reduce_gctFstatline_precision(void*line) {
  char linebuf[256];
  print_gctFstatline_to_str((*(GCTtopOutputEntry*)line), linebuf, sizeof(linebuf));
  sscanf(linebuf,
         "%" LAL_REAL8_FORMAT
         " %" LAL_REAL8_FORMAT
         " %" LAL_REAL8_FORMAT
         " %" LAL_REAL8_FORMAT
         "%*s\n",
         &((*(GCTtopOutputEntry*)line).Freq),
         &((*(GCTtopOutputEntry*)line).Alpha),
         &((*(GCTtopOutputEntry*)line).Delta),
         &((*(GCTtopOutputEntry*)line).F1dot));
}

static void reduce_gctFstat_toplist_precision(toplist_t *l) {
  go_through_toplist(l,reduce_gctFstatline_precision);
}


/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   Returns something <0 on error */
int write_gctFstat_toplist_to_fp(toplist_t*tl, FILE*fp, UINT4*checksum) {
  UINT8 c=0,i;
  INT8 r;
  if(checksum)
    *checksum = 0;
  for(i=0;i<tl->elems;i++)
    if ((r = write_gctFstat_toplist_item_to_fp(*((GCTtopOutputEntry*)(void*)(tl->heap[i])), fp, checksum)) < 0) {
      LogPrintf (LOG_CRITICAL, "Failed to write toplistitem to output fp: %d: %s %" LAL_UINT8_FORMAT "\n",
		 errno,strerror(errno),i);
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      return(r);
    } else
      c += r;
  return(c);
}


/* function that does the actual work of atomic_write_gctFstat_toplist_to_file(),
   appending a %DONE marker if specified (not when called from atomic_write_gctFstat_toplist_to_file().
   NOTE that the checksum will be a little wrong when %DOME is appended, as this line is not counted */
static int _atomic_write_gctFstat_toplist_to_file(toplist_t *l, const char *filename, UINT4*checksum, int write_done) {
  char* tempname;
  INT4 length=0;
  FILE * fpnew;
  UINT4 s;
  int ret;

#define TEMP_EXT ".tmp"
  s = strlen(filename)+strlen(TEMP_EXT)+1;
  tempname = (char*)malloc(s);
  if(!tempname) {
    LogPrintf (LOG_CRITICAL, "Could not allocate new filename\n");
    return(-1);
  }
  strcpy(tempname,filename);
  strcat(tempname,TEMP_EXT);

  fpnew=fopen(tempname, "wb");
  if(!fpnew) {
    LogPrintf (LOG_CRITICAL, "Failed to open temp gctFstat file \"%s\" for writing: %d: %s\n",
	       tempname,errno,strerror(errno));
#ifdef _MSC_VER
    LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
    free(tempname);
    return -1;
  }

  /* when done, write code version and command line as comment in the result file */
  if (write_done) {
    int a;
    CHAR *VCSInfoString;

    /* write the version string */
    if ( (VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ")) == NULL ) {
      LogPrintf (LOG_CRITICAL, "XLALVCSInfoString failed.\n");
      length = -1;
    } else {
      ret = fprintf(fpnew,"%s", VCSInfoString);
      XLALFree(VCSInfoString);
      if (ret < 0)
	length = ret;
      else
	length += ret;
    }

    /* write the command-line */
    if (length >= 0) {
      for(a=0;a<global_argc;a++) {
	ret = fprintf(fpnew,"%%%% argv[%d]: '%s'\n", a, global_argv[a]);
	if (ret < 0) {
	  length = ret;
	  break;
	} else
	  length += ret;
      }
    }

    /* write internal toplist sorting as header line
       NOTE this does not necessarily correspond to the actual output-file row sorting, which by default is done by frequency */
    if (length >= 0) {
      const CHAR *sortstat = NULL;
      if ( l->smaller == gctFstat_smaller )
        sortstat = "<2F>";
      else if ( l->smaller == gctNC_smaller )
        sortstat = "nc";
      else if ( l->smaller == gctBSGL_smaller )
        sortstat = "BSGL";
      else if ( l->smaller == gctBSGLtL_smaller )
        sortstat = "BSGLtL";
      else if ( l->smaller == gctBtSGLtL_smaller )
        sortstat = "BtSGLtL";
      else {
        LogPrintf (LOG_CRITICAL, "Failed to write toplist sorting line, toplist is sorted by unknown statistic.\n");
        length = -1;
      }
      if (length >= 0) {
        ret = fprintf(fpnew,"%%%% candidates selected by %s as toplist statistic\n", sortstat);
        if (ret < 0)
          length = ret;
        else
          length += ret;
      }
    }

    /* write column headings line */
    if (length >= 0) {
      ret = fprintf(fpnew,"%%%% columns:\n%%%% %s\n", global_column_headings_stringp);
      if (ret < 0)
        length = ret;
      else
        length += ret;
    }

  }

  /* write the actual toplist */
  if (length >= 0) {
    ret = write_gctFstat_toplist_to_fp(l,fpnew,checksum);
      if (ret < 0)
	length = ret;
      else
	length += ret;
  }

  /* write the done marker if told to */
  if ((write_done) && (length >= 0)) {
    ret = fprintf(fpnew,"%%DONE\n");
    if (ret < 0)
      length = ret;
    else
      length += ret;
  }

  fclose(fpnew);

  if (length < 0) {
    LogPrintf (LOG_CRITICAL, "Failed to write temp gctFstat file \"%s\": %d: %s\n",
	       tempname,errno,strerror(errno));
#ifdef _MSC_VER
    LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
    free(tempname);
    return(length);
  }

  if((ret = rename(tempname, filename))) {
    LogPrintf (LOG_CRITICAL, "Failed to rename gctFstat file to \"%s\": %d: %s(%d)\n",
	       filename,ret,strerror(errno),errno);
#ifdef _MSC_VER
    LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
    free(tempname);
    return -1;
  }

  free(tempname);
  return length;
}


/* New easier checkpointing - simply dump the whole toplist (plus a counter and
   a checksum) into a binary file.
   The heap structure array is hard to dump because it's based on pointers, so it
   is restored after reading the data back in by sorting the list once.
*/

/** log an I/O error, i.e. source code line no., ferror, errno and strerror, and doserrno on Windows, too */
#ifndef __func__
#ifdef __FUNCTION__
#define __func__  __FUNCTION__
#else
#define __func__  ""
#endif
#endif

#ifdef _WIN32
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: %s (%s:%d): doserr:%d, errno:%d: %s\n",\
	      mess,filename,__func__,__FILE__,__LINE__,_doserrno,errno,strerror(errno))
#else
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: %s (%s:%d): errno:%d: %s\n",\
	      mess,filename,__func__,__FILE__,__LINE__,errno,strerror(errno))
#endif

/* dumps toplist to a temporary file, then renames the file to filename */
int write_gct_checkpoint(const char*filename, toplist_t*tl, toplist_t*t2, toplist_t*t3, UINT4 counter, BOOLEAN do_sync) {
#define TMP_EXT ".tmp"
  char*tmpfilename;
  FILE*fp;
  UINT4 len;
  UINT4 checksum = 0;
  static UINT4 sync_fail_counter = 0;

  /* do nothing with an empty filename */
  if (!filename) {
    LogPrintf(LOG_DETAIL, "%s(): checkpoint filename NULL\n", __func__);
    return(0);
  }

  /* construct temporary filename */
  len = strlen(filename)+strlen(TMP_EXT)+1;
  tmpfilename=LALCalloc(len,sizeof(char));
  if(!tmpfilename){
    LogPrintf(LOG_CRITICAL,"Couldn't allocate tmpfilename\n");
    return(-2);
  }
  strcpy(tmpfilename,filename);
  strcat(tmpfilename,TMP_EXT);

  /* open tempfile */
  fp=fopen(tmpfilename,"wb");
  if(!fp) {
    LOGIOERROR("Couldn't open",tmpfilename);
    LALFree(tmpfilename);
    return(-1);
  }

  /* write number of elements */
  len = fwrite(&(tl->elems), sizeof(tl->elems), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write elems to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
    LALFree(tmpfilename);
    return(-1);
  }

  /* write data */
  len = fwrite(tl->data, tl->size, tl->elems, fp);
  if(len != tl->elems) {
    LOGIOERROR("Couldn't write data to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %zu\n", len, tl->elems);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
    LALFree(tmpfilename);
    return(-1);
  }

  /* dump heap order */
  for(UINT4 i = 0; i < tl->elems; i++) {
    UINT4 idx = (tl->heap[i] - tl->data) / tl->size;
    len = fwrite(&idx, sizeof(idx), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't write idx to", tmpfilename);
      LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
      return(-1);
    }
    for(len = 0; len < sizeof(idx); len++)
      checksum += *(((char*)&idx) + len);
  }

  if (t2) {
    /* write number of elements */
    len = fwrite(&(t2->elems), sizeof(t2->elems), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't write elems to", tmpfilename);
      LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
      return(-1);
    }

    /* write data */
    len = fwrite(t2->data, t2->size, t2->elems, fp);
    if(len != t2->elems) {
      LOGIOERROR("Couldn't write data to", tmpfilename);
      LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %zu\n", len, t2->elems);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
      return(-1);
    }

    /* dump heap order */
    for(UINT4 i = 0; i < t2->elems; i++) {
      UINT4 idx = (t2->heap[i] - t2->data) / t2->size;
      len = fwrite(&idx, sizeof(idx), 1, fp);
      if(len != 1) {
	LOGIOERROR("Couldn't write idx to", tmpfilename);
	LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
	if(fclose(fp))
	  LOGIOERROR("In addition: couldn't close", tmpfilename);
	LALFree(tmpfilename);
	return(-1);
      }
      for(len = 0; len < sizeof(idx); len++)
	checksum += *(((char*)&idx) + len);
    }
  } /* if t2 */

  /* FIXME: this is getting silly, put this into a function and call it for the 
     different toplists */

  if (t3) {
    /* write number of elements */
    len = fwrite(&(t3->elems), sizeof(t3->elems), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't write elems to", tmpfilename);
      LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
      if(fclose(fp))
        LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
      return(-1);
    }

    /* write data */
    len = fwrite(t3->data, t3->size, t3->elems, fp);
    if(len != t3->elems) {
      LOGIOERROR("Couldn't write data to", tmpfilename);
      LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %zu\n", len, t3->elems);
      if(fclose(fp))
        LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
      return(-1);
    }

    /* dump heap order */
    for(UINT4 i = 0; i < t3->elems; i++) {
      UINT4 idx = (t3->heap[i] - t3->data) / t3->size;
      len = fwrite(&idx, sizeof(idx), 1, fp);
      if(len != 1) {
        LOGIOERROR("Couldn't write idx to", tmpfilename);
        LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
        if(fclose(fp))
          LOGIOERROR("In addition: couldn't close", tmpfilename);
        LALFree(tmpfilename);
        return(-1);
      }
      for(len = 0; len < sizeof(idx); len++)
        checksum += *(((char*)&idx) + len);
    }
  } /* if t3 */


  /* write counter */
  len = fwrite(&counter, sizeof(counter), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write counter to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
    return(-1);
  }

  /* calculate checksum */
  for(len = 0; len < sizeof(tl->elems); len++)
    checksum += *(((char*)&(tl->elems)) + len);
  for(len = 0; len < (tl->elems * tl->size); len++)
    checksum += *(((char*)tl->data) + len);
  if (t2) {
    for(len = 0; len < sizeof(t2->elems); len++)
      checksum += *(((char*)&(t2->elems)) + len);
    for(len = 0; len < (t2->elems * t2->size); len++)
      checksum += *(((char*)t2->data) + len);
  }

  if (t3) {
    for(len = 0; len < sizeof(t3->elems); len++)
      checksum += *(((char*)&(t3->elems)) + len);
    for(len = 0; len < (t3->elems * t3->size); len++)
      checksum += *(((char*)t3->data) + len);
  }

  for(len = 0; len < sizeof(counter); len++)
    checksum += *(((char*)&counter) + len);

  /* write checksum */
  len = fwrite(&checksum, sizeof(checksum), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write checksum to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
      LALFree(tmpfilename);
    return(-1);
  }

  if ( do_sync && (sync_fail_counter < SYNC_FAIL_LIMIT) ) {
    /* make sure the data ends up on disk */
    if(fsync(fileno(fp))) {
      LOGIOERROR("Couldn't sync", tmpfilename);
      sync_fail_counter++;
      if (sync_fail_counter >= SYNC_FAIL_LIMIT)
        LogPrintf(LOG_NORMAL,"WARNING: syncing disabled\n");
      } else {
        sync_fail_counter = 0;
    }
  }

  /* close tempfile */
  if(fclose(fp)) {
    LOGIOERROR("Couldn't close", tmpfilename);
    LALFree(tmpfilename);
    return(-1);
  }

  /* rename to filename */
  if(rename(tmpfilename,filename)) {
    LOGIOERROR("Couldn't rename\n", tmpfilename);
    LALFree(tmpfilename);
    return(-1);
  }

  /* all went well */
  LALFree(tmpfilename);
  return(0);
} /* write_gct_checkpoint() */


int read_gct_checkpoint(const char*filename, toplist_t*tl, toplist_t*t2, toplist_t*t3, UINT4*counter) {
  FILE*fp;
  UINT4 len;
  UINT4 checksum, indexsum = 0;

  /* counter should be 0 if we couldn't read a checkpoint */
  *counter = 0;

  /* do nothing with an empty filename */
  if (!filename) {
    LogPrintf(LOG_DETAIL, "%s(): checkpoint filename NULL\n", __func__);
    return(0);
  }

  /* try to open file */
  fp = fopen(filename, "rb");
  if(!fp) {
    if(errno == ENOENT) {
      LogPrintf(LOG_NORMAL,"INFO: No checkpoint %s found - starting from scratch\n", filename);
      clear_toplist(tl);
      if (t2) clear_toplist(t2);
      return(1);
    } else {
      LOGIOERROR("Checkpoint found but couldn't open", filename);
      clear_toplist(tl);
      if (t2) clear_toplist(t2);
      return(-1);
    }
  }

  /* read number of elements */
  len = fread(&(tl->elems), sizeof(tl->elems), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read elems from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    return(-1);
  }
  /* sanity check */
  if (tl->elems > tl->length) {
    LogPrintf(LOG_CRITICAL,
	      "Number of elements read larger than length of toplist: %zu > %zu\n",
	      tl->elems, tl->length);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    return(-2);
  }

  /* read data */
  len = fread(tl->data, tl->size, tl->elems, fp);
  if(len != tl->elems) {
    LOGIOERROR("Couldn't read data from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %zu\n", len, tl->elems);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    clear_toplist(tl);
    if (t2) clear_toplist(t2);
    return(-1);
  }

  /* read heap order */
  for(UINT4 i = 0; i < tl->elems; i++) {
    UINT4 idx;
    len = fread(&idx, sizeof(idx), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't read idx from", filename);
      LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n",len,1);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", filename);
      return(-1);
    }
    tl->heap[i] = (char*)(tl->data + idx * tl->size);
    for(len = 0; len < sizeof(idx); len++)
      indexsum += *(((char*)&idx) + len);
  }

  if (t2) {
    /* read number of elements */
    len = fread(&(t2->elems), sizeof(t2->elems), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't read elems from", filename);
      LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", filename);
      return(-1);
    }
    /* sanity check */
    if (t2->elems > t2->length) {
      LogPrintf(LOG_CRITICAL,
		"Number of elements read larger than length of toplist: %zu > %zu\n",
		t2->elems, t2->length);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", filename);
      return(-2);
    }

    /* read data */
    len = fread(t2->data, t2->size, t2->elems, fp);
    if(len != t2->elems) {
      LOGIOERROR("Couldn't read data from", filename);
      LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %zu\n", len, t2->elems);
      if(fclose(fp))
	LOGIOERROR("In addition: couldn't close", filename);
      clear_toplist(tl);
      clear_toplist(t2);
      return(-1);
    }

    /* read heap order */
    for(UINT4 i = 0; i < t2->elems; i++) {
      UINT4 idx;
      len = fread(&idx, sizeof(idx), 1, fp);
      if(len != 1) {
	LOGIOERROR("Couldn't read idx from", filename);
	LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n",len,1);
	if(fclose(fp))
	  LOGIOERROR("In addition: couldn't close", filename);
	return(-1);
      }
      t2->heap[i] = (char*)(t2->data + idx * t2->size);
      for(len = 0; len < sizeof(idx); len++)
	indexsum += *(((char*)&idx) + len);
    }
  } /* if (t2) */

  if (t3) {
    /* read number of elements */
    len = fread(&(t3->elems), sizeof(t3->elems), 1, fp);
    if(len != 1) {
      LOGIOERROR("Couldn't read elems from", filename);
      LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
      if(fclose(fp))
        LOGIOERROR("In addition: couldn't close", filename);
      return(-1);
    }
    /* sanity check */
    if (t3->elems > t3->length) {
      LogPrintf(LOG_CRITICAL,
                "Number of elements read larger than length of toplist: %zu > %zu\n",
                t3->elems, t3->length);
      if(fclose(fp))
        LOGIOERROR("In addition: couldn't close", filename);
      return(-2);
    }

    /* read data */
    len = fread(t3->data, t3->size, t3->elems, fp);
    if(len != t3->elems) {
      LOGIOERROR("Couldn't read data from", filename);
      LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %zu\n", len, t3->elems);
      if(fclose(fp))
        LOGIOERROR("In addition: couldn't close", filename);
      clear_toplist(tl);
      clear_toplist(t2);
      clear_toplist(t3);
      return(-1);
    }

    /* read heap order */
    for(UINT4 i = 0; i < t3->elems; i++) {
      UINT4 idx;
      len = fread(&idx, sizeof(idx), 1, fp);
      if(len != 1) {
        LOGIOERROR("Couldn't read idx from", filename);
        LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n",len,1);
        if(fclose(fp))
          LOGIOERROR("In addition: couldn't close", filename);
        return(-1);
      }
      t3->heap[i] = (char*)(t3->data + idx * t3->size);
      for(len = 0; len < sizeof(idx); len++)
        indexsum += *(((char*)&idx) + len);
    }
  } /* if (t3) */


  /* read counter */
  len = fread(counter, sizeof(*counter), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read counter from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    clear_toplist(tl);
    if (t2) clear_toplist(t2);
    if (t3) clear_toplist(t3);
    return(-1);
  }

  /* read checksum */
  len = fread(&checksum, sizeof(checksum), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read checksum from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    if (t2) clear_toplist(t2);
    if (t3) clear_toplist(t3);

    clear_toplist(tl);
    return(-1);
  }

  /* close file */
  if(fclose(fp)) {
    LOGIOERROR("Couldn't close", filename);
    clear_toplist(tl);
    if (t2) clear_toplist(t2);
    if (t3) clear_toplist(t3);

    return(-1);
  }

  /* verify checksum */
  checksum -= indexsum;
  for(len = 0; len < sizeof(tl->elems); len++)
    checksum -= *(((char*)&(tl->elems)) + len);
  for(len = 0; len < (tl->elems * tl->size); len++)
    checksum -= *(((char*)tl->data) + len);
  if (t2) {
    for(len = 0; len < sizeof(t2->elems); len++)
      checksum -= *(((char*)&(t2->elems)) + len);
    for(len = 0; len < (t2->elems * t2->size); len++)
      checksum -= *(((char*)t2->data) + len);
  }
  if (t3) {
    for(len = 0; len < sizeof(t3->elems); len++)
      checksum -= *(((char*)&(t3->elems)) + len);
    for(len = 0; len < (t3->elems * t3->size); len++)
      checksum -= *(((char*)t3->data) + len);
  }
  for(len = 0; len < sizeof(*counter); len++)
    checksum -= *(((char*)counter) + len);
  if(checksum) {
    LogPrintf(LOG_CRITICAL,"Checksum error: %d\n", checksum);
    clear_toplist(tl);
    if (t2) clear_toplist(t2);
    if (t3) clear_toplist(t3);

    return(-2);
  }

  /* all went well */
  LogPrintf(LOG_DEBUG,"Successfully read checkpoint:%d\n", *counter);

  return(0);
} /* read_gct_checkpoint() */


/**
 * removes a checkpoint
 * returns 0 on success, errno on failure
 */
int clear_gct_checkpoint(const char*filename) {
  /* do nothing with an empty filename */
  if (!filename) {
    LogPrintf(LOG_DETAIL, "%s(): checkpoint filename NULL\n", __func__);
    return(0);
  }

  if (unlink(filename)) {
    LOGIOERROR("Couldn't delete checkpoint",filename);
    return (errno);
  }
  return(0);
} /* clear_gct_checkpoint() */


#ifdef DEBUG_SORTING
static void dump_heap_order(const toplist_t*tl, const char*name) {
  unsigned int i;
  FILE*fp;
  if((fp=fopen(name,"w"))) {
    for(i = 0; i < tl->elems; i++) {
      fprintf(fp,"%u\n",(unsigned int)((tl->heap[i] - tl->data) / sizeof(GCTtopOutputEntry)));
    }
    fclose(fp);
  }
}

static void sort_gctFstat_toplist_debug(toplist_t*l) {
  if(!debugfp)
    debugfp=fopen("debug_sort","w");
  sort_gctFstat_toplist(l);
  if(debugfp) {
    fclose(debugfp);
    debugfp=NULL;
  }
}
#endif

int write_hfs_oputput(const char*filename, toplist_t*tl) {
  /* reduce the precision of the calculated values before doing the sort to
     the precision we will write the result with. This should ensure a sorting
     order that looks right to the validator, too */
  reduce_gctFstat_toplist_precision(tl);
#ifdef DEBUG_SORTING
  dump_heap_order(tl,"heap_before.dump");
  sort_gctFstat_toplist_debug(tl);
  dump_heap_order(tl,"heap_after.dump");
#else
  sort_gctFstat_toplist(tl);
#endif
  return(_atomic_write_gctFstat_toplist_to_file(tl, filename, NULL, 1));
}
