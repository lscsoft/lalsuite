#include "FstatToplist.h"
#include "HeapToplist.h"
#include <lal/StringInput.h> /* for LAL_REAL8_FORMAT etc. */

#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>

#ifdef USE_BOINC
#include "filesys.h"
#define fopen boinc_fopen
#define rename boinc_rename
#endif

#include <lal/LogPrintf.h>

RCSID("$Id$");



/* MSC specifics */
#ifdef _MSC_VER

/* errno */
extern int errno;
extern int _doserrno;

/* snprintf */
#define LALSnprintf _snprintf

/* finite */
#include <float.h>
#define finite _finite

#else /* MSC */

/* errno */
#include <errno.h>

/* this is defined in C99 and *should* be in math.h. Long term
   protect this with a HAVE_FINITE */
int finite(double);

#endif  /* MSC */



/* define min macro if not already defined */
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif



/* local prototypes */
static void reduce_fstat_toplist_precision(toplist_t *l);
static int print_fstatline_to_str(FstatOutputEntry fline, char* buf, int buflen);


/* ordering function for sorting the list */
static int fstat_toplist_qsort_function(const void *a, const void *b) {
    if      (((FstatOutputEntry*)a)->Freq  < ((FstatOutputEntry*)b)->Freq)
	return -1;
    else if (((FstatOutputEntry*)a)->Freq  > ((FstatOutputEntry*)b)->Freq)
	return 1;
    else if (((FstatOutputEntry*)a)->Alpha < ((FstatOutputEntry*)b)->Alpha)
	return -1;
    else if (((FstatOutputEntry*)a)->Alpha > ((FstatOutputEntry*)b)->Alpha)
	return 1;
    else if (((FstatOutputEntry*)a)->Delta < ((FstatOutputEntry*)b)->Delta)
	return -1;
    else if (((FstatOutputEntry*)a)->Delta > ((FstatOutputEntry*)b)->Delta)
	return 1;
    else if (((FstatOutputEntry*)a)->f1dot < ((FstatOutputEntry*)b)->f1dot)
	return -1;
    else if (((FstatOutputEntry*)a)->f1dot > ((FstatOutputEntry*)b)->f1dot)
	return 1;
    else
	return 0;
}

/* ordering function defining the toplist */
static int fstat_smaller(const void*a, const void*b) {
  if     (((FstatOutputEntry*)a)->Fstat < ((FstatOutputEntry*)b)->Fstat)
    return(1);
  else if(((FstatOutputEntry*)a)->Fstat > ((FstatOutputEntry*)b)->Fstat)
    return(-1);
  else
    return(fstat_toplist_qsort_function(a,b));
}


/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
int create_fstat_toplist(toplist_t**tl, UINT8 length) {
  return(create_toplist(tl, length, sizeof(FstatOutputEntry), fstat_smaller));
}


/* frees the space occupied by the toplist */
void free_fstat_toplist(toplist_t**l) {
  free_toplist(l);
}


/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist and
   look for the now smallest one.
   Returns 1 if the element was actually inserted, 0 if not. */
int insert_into_fstat_toplist(toplist_t*tl, FstatOutputEntry elem) {
  if ( !tl )
    return 0;
  else
    return(insert_into_toplist(tl, (void*)&elem));
}


/* (q)sort the toplist according to the sorting function. */
void sort_fstat_toplist(toplist_t*l) {
  qsort_toplist(l,fstat_toplist_qsort_function);
}


/* reads a (created!) toplist from an open filepointer
   returns the number of bytes read,
   -1 if the file contained a syntax error,
   -2 if given an improper toplist */
int read_fstat_toplist_from_fp(toplist_t*l, FILE*fp, UINT4*checksum, UINT4 maxbytes) {
    CHAR line[256];       /* buffer for reading a line */
    UINT4 items, lines;   /* number of items read from a line, linecounter */
    UINT4 len, chars = 0; /* length of a line, total characters read from the file */
    UINT4 i;              /* loop counter */
    CHAR lastchar;        /* last character of a line read, should be newline */
    FstatOutputEntry FstatLine;
    REAL8 epsilon=1e-5;

    /* basic check that the list argument is valid */
    if(!l)
	return -2;

    /* make sure the line buffer is terminated correctly */
    line[sizeof(line)-1]='\0';

    /* init the checksum if given */
    if(checksum)
	*checksum = 0;

    /* set maxbytes to maximum if zero */
    if (maxbytes == 0)
	maxbytes--;

    lines=1;
    while(fgets(line,sizeof(line)-1, fp)) {

	len = strlen(line);
	chars += len;

	if (len==0) {
	  LogPrintf (LOG_CRITICAL, "Line %d is empty.\n", lines);
	    return -1;
	}
	else if (line[len-1] != '\n') {
	  LogPrintf (LOG_CRITICAL, 
		     "Line %d is too long or has no NEWLINE. First %d chars are:\n'%s'\n",
		     lines,sizeof(line)-1, line);
	  return -1;
	}
      
	items = sscanf (line,
			 "%" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT "%c",
			&FstatLine.Freq,
			&FstatLine.Alpha,
			&FstatLine.Delta,
			&FstatLine.f1dot,
			&FstatLine.Fstat,
			&lastchar);

	/* check the values scanned */
	if (
	    items != 6 ||

	    !finite(FstatLine.Freq)	||
	    !finite(FstatLine.f1dot)	||
	    !finite(FstatLine.Alpha)	||
	    !finite(FstatLine.Delta)	||
	    !finite(FstatLine.Fstat)	||

	    FstatLine.Freq  < 0.0                    ||
	    FstatLine.Alpha <         0.0 - epsilon  ||
	    FstatLine.Alpha >   LAL_TWOPI + epsilon  ||
	    FstatLine.Delta < -0.5*LAL_PI - epsilon  ||
	    FstatLine.Delta >  0.5*LAL_PI + epsilon  ||
            FstatLine.Fstat < 0.0                    ||

	    lastchar != '\n'
	    ) {
	    LogPrintf (LOG_CRITICAL, 
		       "Line %d has invalid values.\n"
		       "First %d chars are:\n"
		       "%s\n"
		       "All fields should be finite\n"
		       "1st and 2nd field should be positive.\n" 
		       "3rd field should lie between 0 and %1.15f.\n" 
		       "4th field should lie between %1.15f and %1.15f.\n",
		       lines, sizeof(line)-1, line,
		       (double)LAL_TWOPI, (double)-LAL_PI/2.0, (double)LAL_PI/2.0);
	    return -1;
        }

	if (checksum)
	    for(i=0;i<len;i++)
		*checksum += line[i];
	
	insert_into_toplist(l, &FstatLine);
	lines++;

	/* NOTE: it *CAN* happen (and on Linux it DOES) that the fully buffered Fstat stream
	 * gets written to the File at program termination.
	 * This does not seem to happen on Mac though, most likely due to different
	 * exit()-calls used (_exit() vs exit() etc.....)
	 *
	 * The bottom-line is: the File-contents CAN legally extend beyond maxbytes,
	 * which is why we'll ensure here that we don't actually read more than 
	 * maxbytes.
	 */
	if ( chars == maxbytes )
	  {
	    LogPrintf (LOG_DEBUG, "Read exactly %d == maxbytes from Fstat-file, that's enough.\n", 
		       chars);
	    break;
	  }
	/* however, if we've read more than maxbytes, something is gone wrong */
	if ( chars > maxbytes )
	  {
	    LogPrintf (LOG_CRITICAL, "Read %d bytes > maxbytes %d from Fstat-file ... corrupted.\n",
		       chars, maxbytes );
	    return -1;
	  }

    } /* while (fgets() ) */

    return 0;

} /* read_fstat_toplist_from_fp() */


/* Prints a Tooplist line to a string buffer.
   Separate function to force consistency of output and reduced precision for sorting */
static int print_fstatline_to_str(FstatOutputEntry fline, char* buf, int buflen) {
      return(LALSnprintf(buf, buflen,
		  /* output precision: choose by following (generous!) significant-digit constraints:
		   * Freq:1e-13 
		   * Alpha,Delta:1e-7 
		   * f1dot:1e-5
		   * F:1e-6 
		   */
		  "%.13g %.7g %.7g %.5g %.6g\n",
		  fline.Freq,
		  fline.Alpha,
		  fline.Delta,
		  fline.f1dot,
		  fline.Fstat));
}


/* writes an FstatOutputEntry line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given */
int write_fstat_toplist_item_to_fp(FstatOutputEntry fline, FILE*fp, UINT4*checksum) {
    char linebuf[256];
    UINT4 i;

    UINT4 length = print_fstatline_to_str(fline, linebuf, sizeof(linebuf));
    
    if(length>sizeof(linebuf)) {
       return -1;
    }

    if (checksum)
	for(i=0;i<length;i++)
	    *checksum += linebuf[i];

    return(fprintf(fp,"%s",linebuf));
}


/* Reduces the precision of all elements in the toplist which influence the sorting order.
   To be called before sorting and finally writing out the list */
static void reduce_fstatline_precision(const void*line) {
  char linebuf[256];
  print_fstatline_to_str((*(FstatOutputEntry*)line), linebuf, sizeof(linebuf));
  sscanf(linebuf,
	 "%" LAL_REAL8_FORMAT
	 " %" LAL_REAL8_FORMAT
	 " %" LAL_REAL8_FORMAT
	 " %" LAL_REAL8_FORMAT
	 "%*s\n",
	 &((*(FstatOutputEntry*)line).Freq),
	 &((*(FstatOutputEntry*)line).Alpha),
	 &((*(FstatOutputEntry*)line).Delta),
	 &((*(FstatOutputEntry*)line).f1dot));
}

static void reduce_fstat_toplist_precision(toplist_t *l) {
  go_through_toplist(l,reduce_fstatline_precision);
}


/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   Returns something <0 on error */
int write_fstat_toplist_to_fp(toplist_t*tl, FILE*fp, UINT4*checksum) {
   UINT8 c=0,i;
   INT8 r;
   if(checksum)
       *checksum = 0;
   for(i=0;i<tl->elems;i++)
     if ((r = write_fstat_toplist_item_to_fp(*((FstatOutputEntry*)(tl->heap[i])), fp, checksum)) < 0) {
       LogPrintf (LOG_CRITICAL, "Failed to write toplistitem to output fp: %d: %s\n",
		  errno,strerror(errno));
#ifdef _MSC_VER
       LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      return(r);
     } else
       c += r;
   return(c);
}


/* writes the given toplitst to a temporary file, then renames the
   temporary file to filename. The name of the temporary file is
   derived from the filename by appending ".tmp". Returns the number
   of chars written or -1 if the temp file could not be opened. */
int atomic_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum) {
    char tempname[MAXFILENAMELENGTH];
    INT4 length;
    FILE * fpnew;

    strncpy(tempname,filename,sizeof(tempname)-4);
    strcat(tempname,".tmp");
    fpnew=fopen(tempname, "wb");
    if(!fpnew) {
      LogPrintf (LOG_CRITICAL, "Failed to open temp Fstat file \"%s\" for writing: %d: %s\n",
		 tempname,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      return -1;
    }
    length = write_fstat_toplist_to_fp(l,fpnew,checksum);
    fclose(fpnew);
    if (length < 0) {
      LogPrintf (LOG_CRITICAL, "Failed to write temp Fstat file \"%s\": %d: %s\n",
		 tempname,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      return(length);
    }
    if(rename(tempname, filename)) {
      LogPrintf (LOG_CRITICAL, "Failed to rename Fstat file to \"%s\": %d: %s\n",
		 filename,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      return -1;
    } else
      return length;
}


/* meant for the final writing of the toplist
   - reduces toplist precision
   - sorts the toplist
   - the calls atomic_write_fstat_toplist_to_file() */
int final_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum) {
  reduce_fstat_toplist_precision(l);
  sort_fstat_toplist(l);
  return(atomic_write_fstat_toplist_to_file(l,filename,checksum));
}
