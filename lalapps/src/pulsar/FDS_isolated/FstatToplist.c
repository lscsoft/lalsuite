#include "ComputeFStatistic.h"

/* this is defined in C99 and *should* be in math.h. Long term
   protect this with a HAVE_FINITE */
#ifdef _MSC_VER
#include <float.h>
#define finite _finite
#else
int finite(double);
#endif


/* toplist structure based on FstatsClusterOutput */
typedef struct {
    UINT8 length;   /* the length (maximal number of entries) of the toplist */
    UINT8 elems;    /* number of elements currently in the toplist */
    UINT8 smallest; /* index of the smallest element in the toplist */
    FstatsClusterOutput *data; /* points to the actual data */
} toplist;


/* creates a toplist with length elements,
   returns -1 on error (usually out of memory), else 0 */
int create_toplist(toplist**tl, UINT8 length) {

    toplist *thetoplist = malloc(sizeof(toplist));
    if(!thetoplist)
	return(-1);

    thetoplist->length = length;
    thetoplist->elems = 0;
    thetoplist->smallest = 0;
    
    thetoplist->data = malloc(length * sizeof(FstatsClusterOutput));
    if(!thetoplist->data) {
	free(thetoplist);
	return(-1);
    }

    *tl=thetoplist;
    return(0);
}


/* frees the space occupied by the toplist */
void free_toplist(toplist**l) {
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
int insert_into_toplist(toplist*tl, FstatsClusterOutput elem) {
 
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
int write_toplist_to_fp(toplist*tl, FILE*fp) {
   UINT8 i,c=0;
   for(i=0;i<tl->elems;i++)
       c += fprintf(fp,"%e %e %e %e %d %e %e %e\n",
		    tl->data[i].Freq,
		    tl->data[i].f1dot,
		    tl->data[i].Alpha,
		    tl->data[i].Delta,
		    tl->data[i].Nbins,
		    tl->data[i].mean,
		    tl->data[i].std,
		    tl->data[i].max);
   return(c);
}

/* reads a (created!) toplist from an open filepointer
   returns -1 if the file contained a syntax error, -2 if given an improper toplist
*/
int read_toplist_from_fp(toplist*l, FILE*fp) {
    char inline[255]; /* buffer for reading a line */
    UINT4 items, lines; /* number of items read from a line, linecounter */
    char lastchar;    /* last character of a line read, should be newline */
    FstatsClusterOutput FstatLine;
    REAL8 epsilon=1e-5;

    /* basic check that the list argument is valid */
    if(!l)
	return -2;

    lines=1;
    while(fgets(inline,sizeof(inline),fp)) {

	if (strlen(inline)==0 || inline[strlen(inline)-1] != '\n') {
	    LALPrintError(
                "Line too long or has no NEWLINE.  First %d chars are:\n%s\n",
                sizeof(inline), inline);
	    return -1;
	}
      

	items = sscanf (inline, 
			"%" LAL_REAL8_FORMAT
			"%" LAL_REAL8_FORMAT
			"%" LAL_REAL8_FORMAT
			"%" LAL_REAL8_FORMAT
			"%" LAL_REAL8_FORMAT
			"%" LAL_INT4_FORMAT
			"%" LAL_REAL8_FORMAT
			"%" LAL_REAL8_FORMAT
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

	    lastchar != '\0'
	    ) {
	    LALPrintError(
		"Line %d has invalid values.\n"
		"First 255 chars are:\n"
		"%s\n"
		"All fields should be finite\n",
		"1st and 2nd field should be positive.\n" 
		"3rd field should lie between 0 and %1.15f.\n" 
		"4th field should lie between %1.15f and %1.15f.\n",
		lines, inline, (double)LAL_TWOPI, (double)-LAL_PI/2.0, (double)LAL_PI/2.0);
	    return -1;
        }
	
	insert_into_toplist(l, FstatLine);
	lines++;
    }
    return 0;
}
