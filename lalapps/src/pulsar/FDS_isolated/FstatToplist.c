/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix
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

#include <unistd.h>
#include "FstatToplist.h"
#include "HeapToplist.h"
#include <lal/StringInput.h> /* for LAL_REAL8_FORMAT etc. */

#include <lal/LALStdio.h>
#include <lal/LogPrintf.h>

#if defined(USE_BOINC) || defined(EAH_BOINC)
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
#define snprintf _snprintf

/* fsync */
#define fsync _commit
#define fileno _fileno

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

/* maximum number of succesive failures before switching off syncing */
#define SYNC_FAIL_LIMIT 5

/* local prototypes */
static void reduce_fstat_toplist_precision(toplist_t *l);
static int _atomic_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum, int write_done);
static int print_fstatline_to_str(FstatOutputEntry fline, char* buf, int buflen);


/* ordering function for sorting the list */
static int fstat_toplist_qsort_function(const void *a, const void *b) {
    if      (((const FstatOutputEntry*)a)->Freq  < ((const FstatOutputEntry*)b)->Freq)
	return -1;
    else if (((const FstatOutputEntry*)a)->Freq  > ((const FstatOutputEntry*)b)->Freq)
	return 1;
    else if (((const FstatOutputEntry*)a)->Alpha < ((const FstatOutputEntry*)b)->Alpha)
	return -1;
    else if (((const FstatOutputEntry*)a)->Alpha > ((const FstatOutputEntry*)b)->Alpha)
	return 1;
    else if (((const FstatOutputEntry*)a)->Delta < ((const FstatOutputEntry*)b)->Delta)
	return -1;
    else if (((const FstatOutputEntry*)a)->Delta > ((const FstatOutputEntry*)b)->Delta)
	return 1;
    else if (((const FstatOutputEntry*)a)->f1dot < ((const FstatOutputEntry*)b)->f1dot)
	return -1;
    else if (((const FstatOutputEntry*)a)->f1dot > ((const FstatOutputEntry*)b)->f1dot)
	return 1;
    else
	return 0;
}

/* ordering function defining the toplist */
static int fstat_smaller(const void*a, const void*b) {
  if     (((const FstatOutputEntry*)a)->Fstat < ((const FstatOutputEntry*)b)->Fstat)
    return(1);
  else if(((const FstatOutputEntry*)a)->Fstat > ((const FstatOutputEntry*)b)->Fstat)
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
    0 if we found a %DONE marker at the end,
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

        if (!strncmp(line,"%DONE\n",strlen("%DONE\n"))) {
	  LogPrintf(LOG_NORMAL,"WARNING: found end marker - the task was already finished\n");
	  return(0);
        }

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
		       "1st field should be positive.\n" 
		       "2nd field should lie between 0 and %1.15f.\n" 
		       "3rd field should lie between %1.15f and %1.15f.\n",
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

    return chars;

} /* read_fstat_toplist_from_fp() */


/* Prints a Tooplist line to a string buffer.
   Separate function to assure consistency of output and reduced precision for sorting */
static int print_fstatline_to_str(FstatOutputEntry fline, char* buf, int buflen) {
  return(snprintf(buf, buflen,
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

    UINT4 length = print_fstatline_to_str(fline, linebuf, sizeof(linebuf)-1);
    
    if(length>sizeof(linebuf)-1) {
       return -1;
    }

    if (checksum)
	for(i=0;i<length;i++)
	    *checksum += linebuf[i];

    linebuf[sizeof(linebuf)-1] = '\0';

    return(fprintf(fp,"%s",linebuf));
}


/* Reduces the precision of all elements in the toplist which influence the sorting order.
   To be called before sorting and finally writing out the list */
static void reduce_fstatline_precision(void*line) {
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


/* writes the given toplitst to a temporary file, then renames the temporary file to filename.
   The name of the temporary file is derived from the filename by appending ".tmp". Returns the
   number of chars written or -1 if the temp file could not be opened.
   This just calls _atomic_write_fstat_toplist_to_file() telling it not to write a %DONE marker*/
int atomic_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum) {
  return(_atomic_write_fstat_toplist_to_file(l, filename, checksum, 0));
}

/* function that does the actual work of atomic_write_fstat_toplist_to_file(),
   appending a %DONE marker if specified (not when called from atomic_write_fstat_toplist_to_file().
   NOTE that the checksum will be a little wrong when %DOME is appended, as this line is not counted */
static int _atomic_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum, int write_done) {
    char* tempname;
    INT4 length;
    FILE * fpnew;
    UINT4 s;

#define TEMP_EXT ".tmp"
    s = strlen(filename)+strlen(TEMP_EXT)+1;
    tempname = (char*)malloc(s);
    if(!tempname) {
      LogPrintf (LOG_CRITICAL, "Could not allocate new filename\n");
      return(-1);
    }
    strncpy(tempname,filename,s);
    strncat(tempname,TEMP_EXT,s);

    fpnew=fopen(tempname, "wb");
    if(!fpnew) {
      LogPrintf (LOG_CRITICAL, "Failed to open temp Fstat file \"%s\" for writing: %d: %s\n",
		 tempname,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      free(tempname);
      return -1;
    }
    length = write_fstat_toplist_to_fp(l,fpnew,checksum);

    if ((write_done) && (length >= 0)) {
      int ret;
      ret = fprintf(fpnew,"%%DONE\n");
      if (ret < 0)
	length = ret;
      else
	length += ret;
    }

    fclose(fpnew);
    
    if (length < 0) {
      LogPrintf (LOG_CRITICAL, "Failed to write temp Fstat file \"%s\": %d: %s\n",
		 tempname,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif
      free(tempname);
      return(length);
    }

    if(rename(tempname, filename)) {
      LogPrintf (LOG_CRITICAL, "Failed to rename Fstat file to \"%s\": %d: %s\n",
		 filename,errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif 
      free(tempname);
      return -1;
    }

    free(tempname);
    return length;
}


/* meant for the final writing of the toplist
   - reduces toplist precision
   - sorts the toplist
   - then calls atomic_write_fstat_toplist_to_file() */
int final_write_fstat_toplist_to_file(toplist_t *l, char *filename, UINT4*checksum) {
  reduce_fstat_toplist_precision(l);
  sort_fstat_toplist(l);
  return(atomic_write_fstat_toplist_to_file(l,filename,checksum));
}




/* sets up a FStatCheckpointFile from parameters */
int fstat_cpt_file_create (FStatCheckpointFile **cptf,
			   CHAR  *filename,
			   UINT4 bufsize,
			   UINT4 maxsize,
			   toplist_t*list) {

  /* input sanity checks */
  if ( (cptf == NULL) ||
       (*cptf != NULL) ||
       (list == NULL) ||
       (filename == NULL) ||
       (strlen (filename) == 0) ) {
    LogPrintf (LOG_CRITICAL, "ERROR: error in input parameters (fstat_cpt_file_create)\n");
    return(-1);
  }

  /* allocation */
  *cptf = LALMalloc(sizeof(FStatCheckpointFile));
  if (!(*cptf)) {
    LogPrintf (LOG_CRITICAL, "ERROR: out of memeory (fstat_cpt_file_create)\n");
    return(-1);
  }

  (*cptf)->filename = LALMalloc(strlen(filename)+1);
  if (!((*cptf)->filename)) {
    LogPrintf (LOG_CRITICAL, "ERROR: out of memeory (fstat_cpt_file_create)\n");
    LALFree(*cptf);
    *cptf = NULL;
    return(-1);
  }

  if (bufsize > 0) {
    (*cptf)->buffer = LALMalloc(bufsize);
    if (!((*cptf)->buffer)) {
      LogPrintf (LOG_CRITICAL, "ERROR: out of memeory (fstat_cpt_file_create)\n");
      LALFree(*cptf);
      LALFree((*cptf)->filename);
      *cptf = NULL;
      return(-1);
    }
  }

  /* initialization */
  strncpy((*cptf)->filename,filename,strlen(filename)+1);

  (*cptf)->bytes = 0;
  (*cptf)->bufsize = bufsize;
  (*cptf)->maxsize = maxsize;
  (*cptf)->checksum = 0;
  (*cptf)->fp = NULL;
  (*cptf)->list = list;
  return(0);
}



/* destroys a FStatCheckpointFile structure */
int fstat_cpt_file_destroy (FStatCheckpointFile **cptf) {
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }
  if((*cptf)->filename)
    LALFree((*cptf)->filename);
  if((*cptf)->buffer)
    LALFree((*cptf)->buffer);
  LALFree(*cptf);
  *cptf = NULL;
  return(0);
}



/* opens the file named in the structure (for writing) and attaches an output buffer if specified */
int fstat_cpt_file_open (FStatCheckpointFile *cptf) {
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }
  /* try open for appending if the file already exists */
  cptf->fp = fopen(cptf->filename, "rb+");
  if (!(cptf->fp)) {
    LogPrintf (LOG_NORMAL, "WARNING: Couldn't open existing checkpointing toplist file %s\n",cptf->filename);
    /* that didn't work,'s try opening a new file for writing */
    cptf->fp = fopen(cptf->filename, "wb+");
  }
  /* leave here if we couldn't open a checkpointed file at all */
  if (!(cptf->fp)) {
    LogPrintf (LOG_CRITICAL, "ERROR: Couldn't open new checkpointing toplist file %s\n",cptf->filename);
    return(-1);
  }
  /* attach a buffer large enough that no output is written to disk
   * unless we fflush(). Policy is fully buffered. */
  if (cptf->bufsize > 0)
    setvbuf(cptf->fp, cptf->buffer, _IOFBF, cptf->bufsize);
  return(0);
}



/* flushes the checkpoint file (only useful if buffered) */
int fstat_cpt_file_flush (FStatCheckpointFile *cptf) {
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }
  if (!(cptf->fp)) {
    LogPrintf (LOG_CRITICAL, "ERROR: invalid checkpointing toplist file pointer\n");
    return(-1);
  }
  return(fflush(cptf->fp));
}



/* returns information for checkpointing */
extern int fstat_cpt_file_info (FStatCheckpointFile *cptf,
				CHAR**filename, UINT4*bytes, UINT4*checksum) {
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }
  if (filename)
    *filename = cptf->filename;
  if (bytes)
    *bytes = cptf->bytes;
  if (checksum)
    *checksum = cptf->checksum;
  return(0);
}


/* closes and compacts the file */
int fstat_cpt_file_close(FStatCheckpointFile*cptf) {
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }
  fclose(cptf->fp);

  /* reduce the precision of the calculated values to the ones we may have
     read in from a previous output (checkpoint) to achieve a proper sorting order */
  reduce_fstat_toplist_precision(cptf->list);

  /* completely sort the list before writing it (a heap is only partially sorted) */
  sort_fstat_toplist(cptf->list);

  /* write out the list a final time, this time with a %DONE marker at the end */
  return(_atomic_write_fstat_toplist_to_file
	 (cptf->list, cptf->filename, &(cptf->checksum), 1));
}


/* adds an item to the toplist and keeps the file consistent, i.e.
   adds the entry to the file if it was really inserted
   returns 1 if the item was actually inserted, 0 if not,
   -1 in case of an error
 */
int fstat_cpt_file_add (FStatCheckpointFile*cptf, FstatOutputEntry line) {
  int ret, bytes;

  ret = insert_into_toplist(cptf->list, &line);
  if (ret) {
    bytes = write_fstat_toplist_item_to_fp(line, cptf->fp, &(cptf->checksum));
    if (bytes < 0) {
      LogPrintf(LOG_CRITICAL, "Failed to write toplist item to file: %d: %s\n",
		errno,strerror(errno));
#ifdef _MSC_VER
      LogPrintf (LOG_CRITICAL, "Windows system call returned: %d\n", _doserrno);
#endif 
      return(-1);
    }
    cptf->bytes += bytes;

    if ((int)cptf->bytes != ftell(cptf->fp)) 
      LogPrintf(LOG_DEBUG,"ERROR: bytecount mismatch: returned: %d, bytes: %u, file: %ld\n",
		bytes, cptf->bytes, ftell(cptf->fp));
  }

  return(ret);
}


/* reads a written checkpointed toplist back into a toplist
   returns 0 if successfully read a list,
    1 if we found a %DONE marker at the end,
   -1 if the file contained a syntax error,
   -2 if given an improper toplist */
int fstat_cpt_file_read (FStatCheckpointFile*cptf, UINT4 checksum_should, UINT4 maxbytes) {
  INT4  bytes;
  UINT4 checksum_read;
  if (!cptf) {
    LogPrintf (LOG_CRITICAL, "ERROR: FStatCheckpointFile is NULL\n");
    return(-1);
  }

  /* the only actual call here, everything else is error handling */
  bytes = read_fstat_toplist_from_fp(cptf->list, cptf->fp, &checksum_read, maxbytes);

  LogPrintf (LOG_DEBUG, "DEBUG: read_fstat_toplist_from_fp() returned %d\n", bytes);

  cptf->bytes = 0;
  cptf->checksum = 0;

  if (bytes == 0) {
    LogPrintf (LOG_DEBUG, "DEBUG: found %DONE end marker\n");
    return(1);
  } else if (bytes == -2) {
    LogPrintf (LOG_CRITICAL, "ERROR: invalid toplist\n");
    return(bytes);
  } else if (bytes == -1) {
    LogPrintf (LOG_CRITICAL, "ERROR: format error in toplist\n");
    rewind(cptf->fp);
    clear_toplist(cptf->list);
    return(bytes);
  } else if (checksum_read != checksum_should) {
    LogPrintf (LOG_CRITICAL, "ERROR: checksum error in toplist %d / %d\n",
	       checksum_should, checksum_read);
    rewind(cptf->fp);
    clear_toplist(cptf->list);
    return(-1);
  }

  if (ftell(cptf->fp) > bytes) {
    LogPrintf(LOG_DEBUG,"ERROR: File length mismatch: bytes: %u, file: %d, fixing\n",
	      bytes, ftell(cptf->fp));
    fseek(cptf->fp, bytes, SEEK_SET);  
  }
  cptf->bytes = bytes;
  cptf->checksum = checksum_read;

  if ((int)(cptf->bytes) != ftell(cptf->fp)) 
    LogPrintf(LOG_DEBUG,"ERROR: read: %u, file: %ld\n", cptf->bytes, ftell(cptf->fp));

  return(0);
}


/* compacts a checkpointed toplist, i.e. re-writes it from scratch from memory */
int fstat_cpt_file_compact(FStatCheckpointFile*cptf) {
  INT4  bytes;
  UINT4 checksum;

  LogPrintf(LOG_NORMAL, "Compacting toplist file\n");
  if (cptf->fp)
    fclose(cptf->fp);
  else
    LogPrintf(LOG_CRITICAL, "ERROR: Toplist filepointer is NULL\n");

  bytes = atomic_write_fstat_toplist_to_file(cptf->list, cptf->filename, &checksum);
  if (bytes < 0) {
    return(bytes);
  }

  fstat_cpt_file_open(cptf);
  fseek(cptf->fp, 0L, SEEK_END);
  cptf->bytes = bytes;
  cptf->checksum = checksum;

  return(bytes);
}



/* New easier checkpointing - simply dump the whole toplist (plus a counter and
   a checksum) into a binary file.
   The heap structure array is hard to dump because it's based on pointers, so it
   is restored after reding the data back in by sorting the list once.
*/

/** log an I/O error, i.e. source code line no., ferror, errno and strerror, and doserrno on Windows, too */
#ifdef _MSC_VER
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, doserr:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,_doserrno,((fp)?(ferror(fp)):0),errno,strerror(errno))
#else
#define LOGIOERROR(mess,filename) \
    LogPrintf(LOG_CRITICAL, "ERROR: %s %s: line:%d, ferr:%d, errno:%d: %s\n",\
	      mess,filename,__LINE__,((fp)?(ferror(fp)):0),errno,strerror(errno))
#endif

/* dumps toplist to a temporary file, then renames the file to filename */
int write_hs_checkpoint(const char*filename, toplist_t*tl, UINT4 counter, BOOLEAN do_sync) {
#define TMP_EXT ".tmp"
  char*tmpfilename;
  FILE*fp;
  UINT4 len;
  UINT4 checksum;
  static UINT4 sync_fail_counter = 0;

  /* construct temporary filename */
  len = strlen(filename)+strlen(TMP_EXT)+1;
  tmpfilename=LALCalloc(len,sizeof(char));
  if(!tmpfilename){
    LogPrintf(LOG_CRITICAL,"Couldn't allocate tmpfilename\n");
    return(-2);
  }
  strncpy(tmpfilename,filename,len);
  strncat(tmpfilename,TMP_EXT,len);

  /* calculate checksum */
  checksum = 0;
  for(len = 0; len < sizeof(tl->elems); len++)
    checksum += *(((char*)&(tl->elems)) + len);
  for(len = 0; len < (tl->elems * tl->size); len++)
    checksum += *(((char*)tl->data) + len);
  for(len = 0; len < sizeof(counter); len++)
    checksum += *(((char*)&counter) + len);

  /* open tempfile */
  fp=fopen(tmpfilename,"wb");
  if(!fp) {
    LOGIOERROR("Couldn't open",tmpfilename);
    return(-1);
  }

  /* write number of elements */
  len = fwrite(&(tl->elems), sizeof(tl->elems), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write elems to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
    return(-1);
  }

  /* write data */
  len = fwrite(tl->data, tl->size, tl->elems, fp);
  if(len != tl->elems) {
    LOGIOERROR("Couldn't write data to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n", len, tl->elems);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
    return(-1);
  }

  /* write counter */
  len = fwrite(&counter, sizeof(counter), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write counter to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
    return(-1);
  }

  /* write checksum */
  len = fwrite(&checksum, sizeof(checksum), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't write checksum to", tmpfilename);
    LogPrintf(LOG_CRITICAL,"fwrite() returned %d, length was %d\n",len,1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", tmpfilename);
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
    return(-1);
  }

  /* rename to filename */
  if(rename(tmpfilename,filename)) {
    LOGIOERROR("Couldn't rename\n", tmpfilename);
    return(-1);
  }

  /* all went well */
  return(0);
}


int read_hs_checkpoint(const char*filename, toplist_t*tl, UINT4*counter) {
  FILE*fp;
  UINT4 len;
  UINT4 checksum;

  /* counter should be 0 if we couldn't read a checkpoint */
  *counter = 0;

  /* try to open file */
  fp = fopen(filename, "rb");
  if(!fp) {
    LogPrintf(LOG_NORMAL,"INFO: Couldn't open checkpoint %s\n", filename);
    clear_toplist(tl);
    return(1);
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
	      "Number of elements read larger than length of toplist: %d, > %d\n",
	      tl->elems, tl->length);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    return(-2);
  }

  /* read data */
  len = fread(tl->data, tl->size, tl->elems, fp);
  if(len != tl->elems) {
    LOGIOERROR("Couldn't read data from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, tl->elems);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    clear_toplist(tl);
    return(-1);
  }

  /* read counter */
  len = fread(counter, sizeof(*counter), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read counter from", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    clear_toplist(tl);
    return(-1);
  }

  /* read checksum */
  len = fread(&checksum, sizeof(checksum), 1, fp);
  if(len != 1) {
    LOGIOERROR("Couldn't read checksum to", filename);
    LogPrintf(LOG_CRITICAL,"fread() returned %d, length was %d\n", len, 1);
    if(fclose(fp))
      LOGIOERROR("In addition: couldn't close", filename);
    clear_toplist(tl);
    return(-1);
  }

  /* close file */
  if(fclose(fp)) {
    LOGIOERROR("Couldn't close", filename);
    clear_toplist(tl);
    return(-1);
  }

  /* verify checksum */
  for(len = 0; len < sizeof(tl->elems); len++)
    checksum -= *(((char*)&(tl->elems)) + len);
  for(len = 0; len < (tl->elems * tl->size); len++)
    checksum -= *(((char*)tl->data) + len);
  for(len = 0; len < sizeof(*counter); len++)
    checksum -= *(((char*)counter) + len);
  if(checksum) {
    LogPrintf(LOG_CRITICAL,"Checksum error: %d\n", checksum);
    clear_toplist(tl);
    return(-2);
  }

  /* restore Heap structure by sorting */
  for(len = 0; len < tl->elems; len++)
    tl->heap[len] = tl->data + len * tl->size;
  qsort_toplist_r(tl,fstat_smaller);

  /* all went well */
  LogPrintf(LOG_DEBUG,"Successfully read checkpoint\n");

  return(0);
}


int write_hs_oputput(const char*filename, toplist_t*tl) {
  /* reduce the precision of the calculated values before doing the sort to
     the precision we will write the result with. This should ensure a sorting
     order that looks right to the validator, too */
  reduce_fstat_toplist_precision(tl);
  sort_fstat_toplist(tl);
  return(_atomic_write_fstat_toplist_to_file(tl, filename, NULL, 1));
}
