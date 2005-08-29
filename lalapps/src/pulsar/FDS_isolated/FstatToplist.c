#include "FstatToplist.h"
#include <lal/StringInput.h> /* for LAL_REAL8_FORMAT etc. */
#include <stdio.h> /* for rename */

RCSID("$Id$");

/* this is defined in C99 and *should* be in math.h. Long term
   protect this with a HAVE_FINITE */
#ifdef _MSC_VER
#include <float.h>
#define finite _finite
#else
int finite(double);
#endif

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
    
    thetoplist->data = malloc(length * sizeof());
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
 
    /* check if the toplist is full, if not, just add the new element */
    if (tl->elems < tl->length) {
	/* just add the new element if the toplist is not full yet */
	tl->data[tl->elems] = elem;
	/* update the smallest if necessary */
	if ((tl->elems==0) || (elem.max < tl->data[tl->smallest].max))
	    tl->smallest = tl->elems;
	tl->elems++;
	return(1); /* element was inserted */

    /* if the toplist is full, we don't need to do anything if the new
       element is even smaller than the smallest one in the toplist */
    } else if (elem.max > tl->data[tl->smallest].max) {
	UINT8 i; /* loop counter */
	/* replace the smaleest element so far with the new one */
	tl->data[tl->smallest] = elem;
	/* find smallest element again */
	tl->smallest = 0;
	for(i=1; i<tl->elems; i++)
	   if(tl->data[i].max < tl->data[tl->smallest].max)
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
    char inline[256];     /* buffer for reading a line */
    UINT4 items, lines;   /* number of items read from a line, linecounter */
    UINT4 len, chars = 0; /* length of a line, total characters read from the file */
    UINT4 i;              /* loop counter */
    char lastchar;        /* last character of a line read, should be newline */
    TOPLISTLINE FstatLine;
    REAL8 epsilon=1e-5;

    /* basic check that the list argument is valid */
    if(!l)
	return -2;

    /* make sure the inline buffer is terminated correctly */
    inline[sizeof(inline)-1]='\0';

    /* init the checksum if given */
    if(checksum)
	*checksum = 0;

    /* set maxbytes to maximum if zero */
    if (maxbytes == 0)
	maxbytes--;

    lines=1;
    while(fgets(inline,min(sizeof(inline)-1,maxbytes-chars),fp)) {

	len = strlen(inline);
	chars += len;

	if (len==0 || inline[len-1] != '\n') {
	    LALPrintError(
                "Line too long or has no NEWLINE. First %d chars are:\n%s\n",
                sizeof(inline)-1, inline);
	    return -1;
	}
      
	items = sscanf (inline,
			 "%" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_INT4_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			" %" LAL_REAL8_FORMAT
			"%c",
			&FstatLine.Freq,
			&FstatLine.f1dot,
			&FstatLine.Alpha,
			&FstatLine.Delta,
			&FstatLine.Nbins,
			&FstatLine.mean,
			&FstatLine.std,
			&FstatLine.max,
			&lastchar);

	/* check the values scanned */
	if (
	    items != 9 ||

	    !finite(FstatLine.Freq)	||
	    !finite(FstatLine.f1dot)	||
	    !finite(FstatLine.Alpha)	||
	    !finite(FstatLine.Delta)	||
	    !finite(FstatLine.Nbins)	||
	    !finite(FstatLine.mean)	||
	    !finite(FstatLine.std)	||
	    !finite(FstatLine.max)	||

	    FstatLine.Freq  < 0.0                    ||
	    FstatLine.f1dot < 0.0                    ||
	    FstatLine.Alpha <         0.0 - epsilon  ||
	    FstatLine.Alpha >   LAL_TWOPI + epsilon  ||
	    FstatLine.Delta < -0.5*LAL_PI - epsilon  ||
	    FstatLine.Delta >  0.5*LAL_PI + epsilon  ||

	    lastchar != '\n'
	    ) {
	    LALPrintError(
		"Line %d has invalid values.\n"
		"First %d chars are:\n"
		"%s\n"
		"All fields should be finite\n",
		"1st and 2nd field should be positive.\n" 
		"3rd field should lie between 0 and %1.15f.\n" 
		"4th field should lie between %1.15f and %1.15f.\n",
		lines, sizeof(inline)-1, inline,
		(double)LAL_TWOPI, (double)-LAL_PI/2.0, (double)LAL_PI/2.0);
	    return -1;
        }

	if (checksum)
	    for(i=0;i<len;i++)
		*checksum += inline[i];
	
	insert_into_toplist(l, FstatLine);
	lines++;
    }
    return 0;
}

/* writes an TOPLISTLINE line to an open filepointer.
   Returns the number of chars written, -1 if in error
   Updates checksum if given */
int write_toplist_item_to_fp(TOPLISTLINE fline, FILE*fp, UINT4*checksum) {
    char linebuf[256];
    UINT4 i;

    UINT4 length =
	snprintf(linebuf, sizeof(linebuf),
		 "%1.13e %e %e %e %d %e %e %1.15e\n",
		 fline.Freq,
		 fline.f1dot,
		 fline.Alpha,
		 fline.Delta,
		 fline.Nbins,
		 fline.mean,
		 fline.std,
		 fline.max);

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
    char tempname[256];
    UINT4 length;

    strncpy(tempname,filename,sizeof(tempname)-4);
    strcat(tempname,".tmp");
    FILE *fpnew=fopen(tempname, "w");
    if(!fpnew)
	return -1;
    length = write_toplist_to_fp(l,fpnew,checksum);
    fclose(fpnew);
    rename(tempname, filename);
    return length;
}
