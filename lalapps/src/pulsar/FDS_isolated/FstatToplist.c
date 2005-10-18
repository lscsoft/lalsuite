#include "FstatToplist.h"
#include <lal/StringInput.h> /* for LAL_REAL8_FORMAT etc. */

#include <lal/LALStdio.h>

#ifdef USE_BOINC
#include "filesys.h"
#define fopen boinc_fopen
#define rename boinc_rename
#endif

#include "LogPrintf.h"

RCSID("$Id$");

/* MSC specifics */
#ifdef _MSC_VER

/* snprintf */
#define LALSnprintf _snprintf

/* finite */
#include <float.h>
#define finite _finite

#else /* MSC */

/* this is defined in C99 and *should* be in math.h. Long term
   protect this with a HAVE_FINITE */
int finite(double);

#endif  /* MSC */

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
int create_toplist(toplist_t**tl, UINT8 length) {
    UINT8 i;

    toplist_t *thetoplist = malloc(sizeof(toplist_t));
    if(!thetoplist)
	return(-1);

    thetoplist->length = length;
    thetoplist->elems = 0;
    thetoplist->smallest = 0;
    
    thetoplist->data = malloc(length * sizeof(TOPLISTLINE));
    if(!thetoplist->data) {
	free(thetoplist);
	return(-1);
    }

    thetoplist->sorted = malloc(length * sizeof(TOPLISTLINE*));
    if(!thetoplist->sorted) {
	free(thetoplist->data);
	free(thetoplist);
	return(-1);
    }

    for(i=0; i<length; i++)
	thetoplist->sorted[i] = &(thetoplist->data[i]);

    *tl=thetoplist;
    return(0);
}


/* frees the space occupied by the toplist */
void free_toplist(toplist_t**l) {
    if(l)
	if(*l){
	    if((**l).data)
		free((**l).data);
	    free(*l);
	    *l=NULL;
	}
}


/* Inserts an element in to the toplist either if there is space left
   or the element is larger than the smallest element in the toplist.
   In the latter case, remove the smallest element from the toplist and
   look for the now smallest one.
   Returns 1 if the element was actually inserted, 0 if not. */
int insert_into_toplist(toplist_t*tl, TOPLISTLINE elem) {
 
  if ( !tl )
    return 0;

    /* check if the toplist is full, if not, just add the new element */
    if (tl->elems < tl->length) {
	/* just add the new element if the toplist is not full yet */
	tl->data[tl->elems] = elem;
	/* update the smallest if necessary */
	if ((tl->elems==0) || (elem.Fstat < tl->data[tl->smallest].Fstat))
	    tl->smallest = tl->elems;
	tl->elems++;
	return(1); /* element was inserted */

    /* if the toplist is full, we don't need to do anything if the new
       element is even smaller than the smallest one in the toplist */
    } else if (elem.Fstat > tl->data[tl->smallest].Fstat) {
	UINT8 i; /* loop counter */
	/* replace the smaleest element so far with the new one */
	tl->data[tl->smallest] = elem;
	/* find smallest element again */
	tl->smallest = 0;
	for(i=1; i<tl->elems; i++)
	   if(tl->data[i].Fstat < tl->data[tl->smallest].Fstat)
	      tl->smallest = i;
	return(1); /* element was inserted */
    }

    /* if we end up here, the element was not inserted */
    return(0);
}


/* Writes the toplist to an (already open) filepointer
   Returns the number of written charactes
   Returns something <0 on error */
int write_toplist_to_fp(toplist_t*tl, FILE*fp, UINT4*checksum) {
   UINT8 i,c=0;
   if(checksum)
       *checksum = 0;
   for(i=0;i<tl->elems;i++)
       c += write_toplist_item_to_fp(*(tl->sorted[i]), fp, checksum);
   return(c);
}


/* ordering function for sorting the list */
static int _toplist_qsort_function(const void *ppa, const void *ppb) {
    const TOPLISTLINE**pa = ppa;
    const TOPLISTLINE**pb = ppb;
    const TOPLISTLINE*a = *pa;
    const TOPLISTLINE*b = *pb;

    if (a->Freq < b->Freq)
	return -1;
    else if (a->Freq > b->Freq)
	return 1;
    else if (a->Alpha < b->Alpha)
	return -1;
    else if (a->Alpha > b->Alpha)
	return 1;
    else if (a->Delta < b->Delta)
	return -1;
    else if (a->Delta > b->Delta)
	return 1;
    else
	return 0;
}

/* (q)sort the toplist according to the sorting function.
   This actually only updates the "sorted" pointers */
void sort_toplist(toplist_t*l) {
    qsort(l->sorted,
	  l->elems,
	  sizeof(l->sorted),
	  _toplist_qsort_function);
}

/* reads a (created!) toplist from an open filepointer
   returns the number of bytes read,
   -1 if the file contained a syntax error, -2 if given an improper toplist
*/
int read_toplist_from_fp(toplist_t*l, FILE*fp, UINT4*checksum, UINT4 maxbytes) {
    CHAR line[256];       /* buffer for reading a line */
    UINT4 items, lines;   /* number of items read from a line, linecounter */
    UINT4 len, chars = 0; /* length of a line, total characters read from the file */
    UINT4 i;              /* loop counter */
    CHAR lastchar;        /* last character of a line read, should be newline */
    TOPLISTLINE FstatLine;
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
	
	insert_into_toplist(l, FstatLine);
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

} /* read_toplist_from_fp() */

/* writes an TOPLISTLINE line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given */
int write_toplist_item_to_fp(TOPLISTLINE fline, FILE*fp, UINT4*checksum) {
    char linebuf[256];
    UINT4 i;

    UINT4 length =
      LALSnprintf(linebuf, sizeof(linebuf),
		  /* output precision: choose by following (generous!) relative-error constraints:
		   * Freq:1e-9 
		   * Alpha,Delta:1e-9 
		   * f1dot:1e-9 
		   * F:1e-6 
		   */
		  "%10f %10f %10f %10g %7f\n", 
		  fline.Freq,
		  fline.Alpha,
		  fline.Delta,
		  fline.f1dot,
		  fline.Fstat);
    
    if(length>sizeof(linebuf))
	return -1;

    if (checksum)
	for(i=0;i<length;i++)
	    *checksum += linebuf[i];

    return(fprintf(fp,"%s",linebuf));
}

/* writes the given toplitst to a temporary file, then renames the
   temporary file to filename. The name of the temporary file is
   derived from the filename by appending ".tmp". Returns the number
   of chars written or -1 if the temp file could not be opened. */
int atomic_write_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum) {
    char tempname[MAXFILENAMELENGTH];
    UINT4 length;
    FILE * fpnew;

    strncpy(tempname,filename,sizeof(tempname)-4);
    strcat(tempname,".tmp");
    fpnew=fopen(tempname, "wb");
    if(!fpnew)
	return -1;
    length = write_toplist_to_fp(l,fpnew,checksum);
    fclose(fpnew);
    if(rename(tempname, filename))
      return -1;
    else
      return length;
}
