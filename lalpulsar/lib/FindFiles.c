/*
 * Copyright (C) 2010, 2012, 2014, 2016, 2021, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2009, 2011 Adam Mercer
 * Copyright (C) 2004--2006, 2008, 2013 Reinhard Prix
 * Copyright (C) 2004--2008, 2010 Bernd Machenschalk
 * Copyright (C) 2004, 2005 Alicia Sintes
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/*---------- includes ----------*/

#include <string.h>
#include <strings.h>
#include <ctype.h>

#ifndef _MSC_VER
#include <dirent.h>
#else
#include <io.h>
#endif

#include <lal/ConfigFile.h>

#include "SFTinternal.h"

/*---------- internal prototypes ----------*/

static BOOLEAN is_pattern(const char*c); /* filename string is a glob-style pattern */
static int amatch(char *str, char *p);	/* glob pattern-matcher (public domain)*/

/*========== function definitions ==========*/

/**
 * Returns a list of filenames matching the input argument, which may be one of
 * the following:
 *   - <tt>\<file\>[;\<file\>;...]</tt>: a list of filenames.
 *   - <tt>\<glob\>[;\<glob\>;...]</tt>: a list of glob-like pattern(s) such
 *     as <tt>*.sft</tt>, <tt>./dir/\htmlonly\endhtmlonly*.sft</tt>, etc.
 *   - <tt>list:\<filelist\></tt>: a file containing a list of filenames.
 *     Prefixes of the form <tt>file:\htmlonly\endhtmlonly//localhost/</tt>
 *     or <tt>file:\htmlonly\endhtmlonly///</tt> are removed.
 *
 * Note: the list of filenames is returned sorted alphabetically.
 */
LALStringVector *
XLALFindFiles (const CHAR *globstring)
{
#ifndef _MSC_VER
  DIR *dir;
  struct dirent *entry;
#else
  intptr_t dir;
  struct _finddata_t entry;
  CHAR* ptr3;
#endif
  CHAR *dname;
  const CHAR *ptr1, *ptr2;
  CHAR *fpattern;
  size_t dirlen;
  CHAR **filelist = NULL;
  UINT4 numFiles = 0, newNumFiles = 0;
  LALStringVector *ret = NULL;
  UINT4 j;
  UINT4 namelen;
  CHAR *thisFname = NULL;

  XLAL_CHECK_NULL ( globstring != NULL, XLAL_EINVAL );

#define FILE_SEPARATOR ';'
  if ( (ptr2 = strchr (globstring, FILE_SEPARATOR)) )
    { /* globstring is multi-pattern ("pattern1;pattern2;pattern3") */
      /* call XLALFindFiles() with every pattern found in globstring */

      ptr1 = (const CHAR*)globstring;
      while ( (ptr2 = strchr (ptr1, FILE_SEPARATOR)) )
	{
	  /* ptr1 points to the beginning of a pattern, ptr2 to the end */

	  /* copy the current name to thisFname */
	  namelen = ptr2 - ptr1;
	  if ((thisFname = LALRealloc(thisFname, (namelen+1)*sizeof(CHAR))) == NULL) {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    if(filelist)
	      LALFree (filelist);
	    XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  }
	  strncpy(thisFname,ptr1,namelen);
	  thisFname[namelen] = '\0';

	  /* call XLALFindFiles(thisFname) */
	  ret = XLALFindFiles(thisFname);

	  /* append the output (if any) to the existing filelist */
	  if (ret) {
	    newNumFiles = numFiles + ret->length;

	    if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	      XLALDestroyStringVector(ret);
	      LALFree(thisFname);
	      XLAL_ERROR_NULL ( XLAL_ENOMEM );
	    }

	    for(j=0; j < ret->length; j++)
	      filelist[numFiles+j] = ret->data[j];
	    LALFree(ret->data);
	    LALFree(ret);
	    numFiles = newNumFiles;
	  } else {
	    for (j=0; j < numFiles; j++)
	      LALFree (filelist[j]);
	    if(filelist)
	      LALFree (filelist);
	    LALFree(thisFname);
	    XLAL_ERROR_NULL ( XLAL_EFUNC);
	  }

	  /* skip the separator */
	  ptr1 = ptr2 + 1;
	} /* while */

      LALFree(thisFname);

      ret = XLALFindFiles(ptr1);
      if (ret) {
	newNumFiles = numFiles + ret->length;

	if ((filelist = LALRealloc (filelist, (newNumFiles) * sizeof(CHAR*))) == NULL) {
	  XLALDestroyStringVector(ret);
	  XLAL_ERROR_NULL ( XLAL_ENOMEM );
	}

	for(j=0; j < ret->length; j++)
	  filelist[numFiles+j] = ret->data[j];
	LALFree(ret->data);
	LALFree(ret);
	numFiles = newNumFiles;
      }

    } /* if multi-pattern */

  /* read list of file names from a "list file" */
#define LIST_PREFIX "list:"
  else if (strncmp(globstring, LIST_PREFIX, strlen(LIST_PREFIX)) == 0) {
    LALParsedDataFile *list = NULL;
    CHAR* listfname = NULL;

    /* extract list file name */
    if ((listfname = XLALStringDuplicate(globstring + strlen(LIST_PREFIX))) == NULL) {
      XLAL_ERROR_NULL ( XLAL_ENOMEM ) ;
    }
#undef LIST_PREFIX

    /* read list of file names from file */
    if (XLALParseDataFile(&list, listfname) != XLAL_SUCCESS) {
      XLAL_ERROR_NULL ( XLAL_EFUNC, "Could not parse list file '%s'\n",listfname );
    }

    /* allocate "filelist" */
    numFiles = list->lines->nTokens;
    if (numFiles == 0) {
      XLALPrintWarning("\n%s: List file '%s' contains no file names\n", __func__, listfname);
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      XLAL_ERROR_NULL ( XLAL_EINVAL );
    }
    if ((filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
      LALFree(listfname);
      XLALDestroyParsedDataFile(list);
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
    }

    /* copy file names from "list" to "filelist" */
    for (j = 0; j < numFiles; ++j) {
      ptr1 = list->lines->tokens[j];

      /* these prefixes are added to file names by e.g. ligo_data_find */
#define FILE_PREFIX "file://localhost/"
      if (strncmp(ptr1, FILE_PREFIX, strlen(FILE_PREFIX)) == 0) {
	ptr1 += strlen(FILE_PREFIX) - 1;
      }
#undef FILE_PREFIX
      else
#define FILE_PREFIX "file:///"
      if (strncmp(ptr1, FILE_PREFIX, strlen(FILE_PREFIX)) == 0) {
	ptr1 += strlen(FILE_PREFIX) - 1;
      }
#undef FILE_PREFIX

      /* allocate "filelist", and cleanup if it fails  */
      if ((filelist[j] = LALCalloc(1, strlen(ptr1) + 1)) == NULL) {
	while (j-- > 0)
	  LALFree(filelist[j]);
	LALFree(filelist);
	LALFree(listfname);
	XLALDestroyParsedDataFile(list);
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      /* copy string */
      strcpy(filelist[j], ptr1);

    }

    /* cleanup */
    LALFree(listfname);
    XLALDestroyParsedDataFile(list);

  } /* if list file */

  else if (is_pattern(globstring))

    { /* globstring is a single glob-style pattern */

      /* First we separate the globstring into directory-path and file-pattern */

#ifndef _WIN32
#define DIR_SEPARATOR '/'
#else
#define DIR_SEPARATOR '\\'
#endif

      /* any path specified or not ? */
      ptr1 = strrchr (globstring, DIR_SEPARATOR);
      if (ptr1)
	{ /* yes, copy directory-path */
	  dirlen = (size_t)(ptr1 - globstring) + 1;
	  if ( (dname = LALCalloc (1, dirlen)) == NULL)
	    XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  strncpy (dname, globstring, dirlen);
	  dname[dirlen-1] = '\0';

	  ptr1 ++; /* skip dir-separator */
	  /* copy the rest as a glob-pattern for matching */
	  if ( (fpattern = LALCalloc (1, strlen(ptr1) + 1)) == NULL )
	    {
	      LALFree (dname);
	      XLAL_ERROR_NULL ( XLAL_ENOMEM );
	    }
	  strcpy (fpattern, ptr1);

	} /* if ptr1 */
      else /* no pathname given, assume "." */
	{
	  if ( (dname = LALCalloc(1, 2)) == NULL)
            XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  strcpy (dname, ".");

	  if ( (fpattern = LALCalloc(1, strlen(globstring)+1)) == NULL)
	    {
	      LALFree (dname);
              XLAL_ERROR_NULL ( XLAL_ENOMEM );
	    }
	  strcpy (fpattern, globstring);	/* just file-pattern given */
	} /* if !ptr */


#ifndef _MSC_VER
      /* now go through the file-list in this directory */
      if ( (dir = opendir(dname)) == NULL) {
	XLALPrintError ("Can't open data-directory `%s`\n", dname);
	LALFree (dname);
        XLAL_ERROR_NULL ( XLAL_EIO );
      }
#else
      if ((ptr3 = (CHAR*)LALMalloc(strlen(dname)+3)) == NULL)
	return(NULL);
      sprintf(ptr3,"%s\\*",dname);
      dir = _findfirst(ptr3,&entry);
      LALFree(ptr3);
      if (dir == -1) {
	XLALPrintError ("Can't find file for pattern `%s`\n", ptr3);
	LALFree (dname);
        XLAL_ERROR_NULL ( XLAL_EIO );
      }
#endif

#ifndef _MSC_VER
      while ( (entry = readdir (dir)) != NULL )
#else
      do
#endif
	{
#ifndef _MSC_VER
	  thisFname = entry->d_name;
#else
	  thisFname = entry.name;
#endif

	  /* now check if glob-pattern fpattern matches the current filename */
	  if ( amatch(thisFname, fpattern)
	       /* and check if we didnt' match some obvious garbage like "." or ".." : */
	       && strcmp( thisFname, ".") && strcmp( thisFname, "..") )
	    {

	      numFiles ++;
	      if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
		LALFree (dname);
		LALFree (fpattern);
                XLAL_ERROR_NULL ( XLAL_ENOMEM );
	      }

	      namelen = strlen(thisFname) + strlen(dname) + 2 ;

	      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
		for (j=0; j < numFiles; j++)
		  LALFree (filelist[j]);
		LALFree (filelist);
		LALFree (dname);
		LALFree (fpattern);
                XLAL_ERROR_NULL ( XLAL_ENOMEM );
	      }

	      sprintf(filelist[numFiles-1], "%s%c%s", dname, DIR_SEPARATOR, thisFname);

	    } /* if filename matched pattern */

	} /* while more directory entries */
#ifdef _MSC_VER
      while ( _findnext (dir,&entry) == 0 );
#endif

#ifndef _MSC_VER
      closedir (dir);
#else
      _findclose(dir);
#endif

      LALFree (dname);
      LALFree (fpattern);

    } /* if is_pattern */

  else

    { /* globstring is a single simple filename */
      /* add it to the list of filenames as it is */

      numFiles++;
      if ( (filelist = LALRealloc (filelist, numFiles * sizeof(CHAR*))) == NULL) {
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      namelen = strlen(globstring) + 1;
      if ( (filelist[ numFiles - 1 ] = LALCalloc (1, namelen)) == NULL) {
	LALFree (filelist);
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      strcpy(filelist[numFiles-1], globstring );
    }

  /* ok, did we find anything? */
  if (numFiles == 0)
    XLAL_ERROR_NULL ( XLAL_EINVAL );


  /* make a LALStringVector from the list of filenames */
  if ( (ret = LALCalloc (1, sizeof (LALStringVector) )) == NULL)
    {
      for (j=0; j<numFiles; j++)
	LALFree (filelist[j]);
      LALFree (filelist);
      XLAL_ERROR_NULL ( XLAL_ENOMEM );
    }
  ret->length = numFiles;
  ret->data = filelist;

  /* sort this alphabetically (in-place) */
  if(numFiles>1)
    XLALSortStringVector (ret);

  return (ret);

} /* XLALFindFiles() */


/* filename string is a glob-style pattern, i.e. it contains '*' or '?' or '[' */
static BOOLEAN is_pattern(const char*c) {
  while((*c != '\0') && (*c != '*') && (*c != '?') && (*c != '['))
    c++;
  return(*c != '\0');
}


/*======================================================================*/
/*
 * robust glob pattern matcher
 * ozan s. yigit/dec 1994
 * public domain
 *
 * glob patterns:
 *	*	matches zero or more characters
 *	?	matches any single character
 *	[set]	matches any character in the set
 *	[^set]	matches any character NOT in the set
 *		where a set is a group of characters or ranges. a range
 *		is written as two characters seperated with a hyphen: a-z denotes
 *		all characters between a to z inclusive.
 *	[-set]	set matches a literal hypen and any character in the set
 *	[]set]	matches a literal close bracket and any character in the set
 *
 *	char	matches itself except where char is '*' or '?' or '['
 *	\char	matches char, including any pattern character
 *
 * examples:
 *	a*c		ac abc abbc ...
 *	a?c		acc abc aXc ...
 *	a[a-z]c		aac abc acc ...
 *	a[-a-z]c	a-c aac abc ...
 *
 */

#ifndef NEGATE
#define NEGATE	'^'			/* std cset negation char */
#endif

static int
amatch(char *str, char *p)
{
	int negate;
	int match;
	int c;

	while (*p) {
		if (!*str && *p != '*')
			return FALSE;

		switch (c = *p++) {

		case '*':
			while (*p == '*')
				p++;

			if (!*p)
				return TRUE;

			if (*p != '?' && *p != '[' && *p != '\\')
				while (*str && *p != *str)
					str++;

			while (*str) {
				if (amatch(str, p))
					return TRUE;
				str++;
			}
			return FALSE;

		case '?':
			if (*str)
				break;
			return FALSE;
/*
 * set specification is inclusive, that is [a-z] is a, z and
 * everything in between. this means [z-a] may be interpreted
 * as a set that contains z, a and nothing in between.
 */
		case '[':
			if (*p != NEGATE)
				negate = FALSE;
			else {
				negate = TRUE;
				p++;
			}

			match = FALSE;

			while (!match && (c = *p++)) {
				if (!*p)
					return FALSE;
				if (*p == '-') {	/* c-c */
					if (!*++p)
						return FALSE;
					if (*p != ']') {
						if (*str == c || *str == *p ||
						    (*str > c && *str < *p))
							match = TRUE;
					}
					else {		/* c-] */
						if (*str >= c)
							match = TRUE;
						break;
					}
				}
				else {			/* cc or c] */
					if (c == *str)
						match = TRUE;
					if (*p != ']') {
						if (*p == *str)
							match = TRUE;
					}
					else
						break;
				}
			}

			if (negate == match)
				return FALSE;
/*
 * if there is a match, skip past the cset and continue on
 */
			while (*p && *p != ']')
				p++;
			if (!*p++)	/* oops! */
				return FALSE;
			break;

		case '\\':
			if (*p)
				c = *p++;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
			__attribute__ ((fallthrough));
#endif
		default:
			if (c != *str)
				return FALSE;
			break;

		}
		str++;
	}

	return !*str;
}
