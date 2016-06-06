/*
 * Copyright (C) 2015 Reinhard Prix
 * Copyright (C) 2010 Larne Pekowsky
 * Copyright (C) 2004, 2005 Reinhard Prix
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/FileIO.h>
#include <lal/StreamInput.h>
#include <lal/AVFactories.h>
#include <lal/StringVector.h>

#include <lal/UserInputParse.h>

#include <lal/ConfigFile.h>

// ---------- local defines ----------
#define WHITESPACE " \t"

#define TRUE   (1==1)
#define FALSE  (1==0)

// ---------- local prototypes ----------
static void cleanConfig ( char *text );

// ==================== function definitions ==========
/**
 * Parse an ASCII data-file into a pre-cleaned array of lines.
 * The cleaning gets rid of comments ('\#', '\%'), empty lines,
 * and performs line-continuation if '\\' is found at EOL
 *
 * NOTE: This function can transparently detect and read gzip-compressed
 * data-files, independently of filename-extension
 *
 * NOTE2: allows passing of *file-contents* directly instead of filename to read
 * by passing a string like "{ file-contents }" instead of a file-path,
 * ie if first character == '{' and last character == '}'
 * This is useful to allow ascii-file user inputs to be transparently
 * passed as filenames or direct contents
 */
int
XLALParseDataFile (LALParsedDataFile **cfgdata, /**< [out] pre-parsed data-file lines */
                   const CHAR *path             /**< [in] file-path of config-file to be read */
                   )
{
  XLAL_CHECK ( (cfgdata != NULL) && (*cfgdata == NULL), XLAL_EINVAL );
  XLAL_CHECK ( path != NULL, XLAL_EINVAL );

  char *dataBuffer;
  if ( (path[0] == '{') && (path[strlen(path)-1] == '}') )
    {
      XLAL_CHECK ( (dataBuffer = XLALStringDuplicate ( path + 1 )) != NULL, XLAL_EFUNC );
      dataBuffer[strlen(dataBuffer)-1] = 0;
    }
  else
    {
      XLAL_CHECK ( (dataBuffer = XLALFileLoad ( path )) != NULL, XLAL_EFUNC );
    }

  if ( XLALParseDataFileContent ( cfgdata, dataBuffer ) != XLAL_SUCCESS ) {
    XLALFree ( dataBuffer );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  XLALFree ( dataBuffer );

  return XLAL_SUCCESS;

} /* XLALParseDataFile() */

int
XLALParseDataFileContent (LALParsedDataFile **cfgdata,  /**< [out] pre-parsed data-file lines */
                          const CHAR *string            /**< [in] string-contents of config-file: can get modified! */
                          )
{
  XLAL_CHECK ( (cfgdata != NULL) && (*cfgdata == NULL), XLAL_EINVAL );
  XLAL_CHECK ( string != NULL, XLAL_EINVAL );

  char *rawdata;
  XLAL_CHECK ( (rawdata = XLALMalloc ( strlen(string) + 1 )) != NULL, XLAL_ENOMEM );
  strcpy ( rawdata, string );   // keep local copy for modifying

  /* get rid of comments and do line-continuation */
  cleanConfig ( rawdata );

  LALParsedDataFile *cfg;
  XLAL_CHECK ( (cfg = XLALCalloc (1, sizeof(*cfg))) != NULL, XLAL_ENOMEM );

  /* parse this into individual lines */
  int err = XLALCreateTokenList ( &(cfg->lines), rawdata, "\n");
  if (err) {
    XLALFree (cfg);
    XLALFree (rawdata);
    XLAL_ERROR ( XLAL_EFUNC, "XLALCreateTokenList() failed.\n" );
  }
  XLALFree (rawdata);

  /* initialize the 'wasRead' flags for the lines */
  if ( cfg->lines->nTokens )
    {
      int len = cfg->lines->nTokens * sizeof(cfg->wasRead[0]);
      if ( (cfg->wasRead = XLALCalloc(1, len )) == NULL )
        {
          XLALFree (cfg->lines);
          XLALFree (cfg);
          XLAL_ERROR ( XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );
        }
    }
  else {
    cfg->wasRead = NULL;
  }

  (*cfgdata) = cfg;

  return XLAL_SUCCESS;

} // XLALParseDataFileContent()


/**
 * Free memory associated with a LALParsedDataFile structure.
 */
void
XLALDestroyParsedDataFile (LALParsedDataFile *cfgdata)  /**< [in] config-file data */
{
  if ( cfgdata == NULL ) {
    return;
  }

  XLALDestroyTokenList ( cfgdata->lines );
  XLALFree ( cfgdata->wasRead );
  XLALFree ( cfgdata );

  return;
} /* XLALDestroyParsedDataFile() */


/**
 * Function to determine whether a given section secName exists in the parsed
 * config-file contents cfgdata.
 *
 * \note: this function tolerates NULL input as secName, cfgdata, or cfgdata->lines,
 * in which case the answer is simply 'FALSE'.
 */
int
XLALConfigSectionExists ( const LALParsedDataFile *cfgdata,     /**< [in] pre-parsed config-data */
                          const CHAR *secName)                  /**< [in] section-name to read */
{
  UINT4 i;
  size_t sec_searchlen = 0;

  /* If there's no config file, or no section, then   */
  /* the section isn;t in the config file, return 0   */
  if ( secName == NULL || cfgdata == NULL || cfgdata->lines == NULL )
    {
      return FALSE;
    }

  sec_searchlen = strlen(secName);

  for (i = 0; i < cfgdata->lines->nTokens; i++)
    {
      /* Is this the start of a new section? */
      if (cfgdata->lines->tokens[i][0] == '[')
        {
          /* If we're looking for a particular section, is this it? */
          if (strncmp(cfgdata->lines->tokens[i] + 1, secName, sec_searchlen) == 0)
            {
              return TRUE;
            }
         }
    }

  return FALSE;

} /* XLALConfigSectionExists() */

/**
 * Function to find all sections in given config-file contents cfgdata.
 *
 * A section start is defined by a string "[ section-name ]" found at the beginning of a line
 * The first non-section part of a config-file is referred to as the "default" section,
 * which is included in the returned list of section-names provided it is not empty.
 *
 */
LALStringVector *
XLALListConfigFileSections ( const LALParsedDataFile *cfgdata )    /**< [in] pre-parsed config-data */
{
  XLAL_CHECK_NULL ( (cfgdata != NULL) && ( cfgdata->lines != NULL), XLAL_EINVAL );

  const TokenList *lines = cfgdata->lines;

  LALStringVector *sections;
  XLAL_CHECK_NULL ( (sections = XLALCalloc ( 1, sizeof(*sections) ) ) != NULL, XLAL_ENOMEM );   // empty string vector

  if ( lines->tokens[0][0] != '[' ) // there is a non-empty 'default' section
    {
      XLAL_CHECK_NULL ( (sections = XLALAppendString2Vector ( sections, "default" )) != NULL, XLAL_EFUNC );
    } // if non-empty default section

  for ( UINT4 i = 0; i < lines->nTokens; i++ )
    {
      const CHAR *thisLine = lines->tokens[i];
      XLAL_CHECK_NULL ( thisLine != NULL, XLAL_EINVAL );
      /* Is this the start of a new section? */
      if ( thisLine[0] == '[' )
        {
          UINT4 len = strlen ( thisLine );
          XLAL_CHECK_NULL ( thisLine[len-1] == ']', XLAL_EINVAL, "Invalid section start '%s'\n", thisLine );

          const char *secName0 = thisLine + 1;  // skip '['
          char *secName;
          XLAL_CHECK_NULL ( (secName = XLALDeblankString ( secName0, len - 2 )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (sections = XLALAppendString2Vector ( sections, secName )) != NULL, XLAL_EFUNC );
          XLALFree ( secName );
        } // if section found

    } // for i < numLines

  return sections;

} // XLALListConfigFileSections()



/**
 * String parser for config-file: can read config-variables of the form VARIABLE [=:] VALUE.
 * Input is a TokenList containing the 'logical' lines of the cleaned config-file
 *
 * \note Opening and closing quotes (\' or \") are removed from the returned string.
 */
int
XLALReadConfigSTRINGVariable ( CHAR **varp,                     //!< [out] return string value if found (NULL otherwise)
                               LALParsedDataFile *cfgdata,      //!< [in,out] pre-parsed config-data
                               const CHAR * secName,            //!< [in] section name in which to find variable 'varName', NULL='default' section
                               const CHAR * varName,            //!< [in] variable name to be read
                               BOOLEAN *wasRead                 //!< [out] did we succeed in reading? */
                               )
{
  // check input consistency
  XLAL_CHECK ( (varp != NULL) && (*varp == NULL), XLAL_EINVAL );
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata->lines != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata->lines->nTokens > 0, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata->wasRead != NULL, XLAL_EINVAL );

  BOOLEAN inRightSection = FALSE;

  (*wasRead) = FALSE;

  // If we haven't been asked for a section then we want the
  // "default" section, which starts at the top of the file without any section heading
  if ( secName == NULL )
    {
      inRightSection = TRUE;
    }

  /* find the variable-name in the token-list (and in the right section, if given) */
  size_t searchlen = strlen ( varName );

  for ( UINT4 i = 0; i < cfgdata->lines->nTokens; i++ )
    {
      if (cfgdata->lines->tokens[i][0] == '[')       /* Is this the start of a new section? */
        {
          if ( inRightSection ) {
            // if we previously were in the right section, it means we've now left it
            // and therefore didn't find the variable we were looking for, so we return,
            // but this is not an error!
            return XLAL_SUCCESS;
          }

          if ( (secName != NULL) && ( strncmp ( cfgdata->lines->tokens[i] + 1, secName, strlen(secName)) == 0) )
            {
              inRightSection = TRUE;
            }
        } // end: if start of new section found
      else
        {
          if ( !inRightSection ) {
            continue;
          }

          UINT4 varlen = strcspn ( cfgdata->lines->tokens[i], WHITESPACE "=:" );        /* get length of variable-name */
          XLAL_CHECK ( varlen > 0, XLAL_EDOM, "Parsing error: nonexistent variable name in '%s'\n", cfgdata->lines->tokens[i] );

          // pre-select based on length of variable-name
          if ( varlen != searchlen ) {
            continue;
          }

          // same len, but are they identical ?
          if ( strncmp ( varName, cfgdata->lines->tokens[i], varlen ) == 0 )
            {
              char *strVal = cfgdata->lines->tokens[i] + varlen;
              strVal += strspn ( strVal, WHITESPACE "=:" );     // skip all whitespace and define-chars
              XLAL_CHECK ( XLALParseStringValueAsSTRING ( varp, strVal ) == XLAL_SUCCESS, XLAL_EFUNC ); // copy and remove quotes (if any)
              cfgdata->wasRead[i] = 1;
              (*wasRead) = TRUE;
              break; // exit loop, we've found it
            } // end: if found variable

        } // end: if in right section

    } // end: for i < num_lines

  return XLAL_SUCCESS;

} // XLALReadConfigSTRINGVariable()

// ------------------------------------------------------------
// define type-specific wrappers to the generic XLALReadConfigSTRINGVariable() function,
// using a template macro:
#define DEFINE_XLALREADCONFIGVARIABLE(TYPE,CTYPE)                       \
DECLARE_XLALREADCONFIGVARIABLE(TYPE,CTYPE)                              \
{                                                                       \
 /* first read the value as a string */                                 \
 CHAR *valString = NULL;                                                \
 XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead ) == XLAL_SUCCESS, XLAL_EFUNC ); \
 if ( ! (*wasRead ) ) {                                                 \
   return XLAL_SUCCESS;                                                 \
 }                                                                      \
 XLAL_CHECK ( XLALParseStringValueAs ##TYPE ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC ); \
 XLALFree (valString);                                                  \
 return XLAL_SUCCESS;                                                   \
}

DEFINE_XLALREADCONFIGVARIABLE(BOOLEAN,BOOLEAN);
DEFINE_XLALREADCONFIGVARIABLE(INT4,INT4);
DEFINE_XLALREADCONFIGVARIABLE(INT8,INT8);
DEFINE_XLALREADCONFIGVARIABLE(REAL8,REAL8);
DEFINE_XLALREADCONFIGVARIABLE(STRINGVector,LALStringVector*);
DEFINE_XLALREADCONFIGVARIABLE(EPOCH,LIGOTimeGPS);
DEFINE_XLALREADCONFIGVARIABLE(RAJ,REAL8);
DEFINE_XLALREADCONFIGVARIABLE(DECJ,REAL8);
// ------------------------------------------------------------


/**
 * Return a list of unread config-file entries, NULL if none found (without error).
 */
UINT4Vector *
XLALConfigFileGetUnreadEntries ( const LALParsedDataFile *cfgdata	///< [in] config-file data
                                 )
{
  XLAL_CHECK_NULL ( cfgdata != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( cfgdata->lines != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( cfgdata->lines->nTokens > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL ( cfgdata->wasRead != NULL, XLAL_EINVAL );

  UINT4Vector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCalloc ( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );

  for (UINT4 i=0; i < cfgdata->lines->nTokens; i++)
    {
      // We don't require section headers to be marked as read
      // ie we consider the config-file to be fully read/parsed if
      // if every value in every section has been parsed.
      if (cfgdata->lines->tokens[i][0] == '[') {
        continue;
      }

      if ( ! cfgdata->wasRead[i] ) {
        ret->length ++;
        XLAL_CHECK_NULL ( (ret->data = XLALRealloc ( ret->data, ret->length * sizeof(ret->data[0]) )) != NULL, XLAL_ENOMEM );
        ret->data[ret->length-1] = i;
      }

    } // for i < numLines

  if ( ret->length == 0 )
    {
      XLALFree ( ret );
      ret = NULL;
    }

  return ret;

} // XLALConfigFileGetUnreadEntries()

/*----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------
 * cleanConfig(): do some preprocessing on the config-file, namely 'erase'
 * all comments by '\n', and glue '\'-continued lines
 *----------------------------------------------------------------------*/
void
cleanConfig ( char *text )
{
  if ( text == NULL ) {
    return;
  }

  size_t len;
  CHAR *ptr, *ptr2, *eol;
  BOOLEAN inQuotes = 0;
  INT4 inBracesCount = 0;
  /*----------------------------------------------------------------------
   * RUN 1: clean out comments, by replacing them by '\n'
   */
  ptr = text;
  while ( *ptr )
    {
      if ( (*ptr) == '\"' ) {
        inQuotes = !inQuotes;   /* flip state */
      }
      if ( (*ptr) == '{' ) {
        inBracesCount ++;
      }
      if ( (*ptr) == '}' ) {
        inBracesCount --;
      }

      if ( ((*ptr) == '#') || ( (*ptr) == '%') ) {
        if ( !inQuotes )        /* only consider as comments if not quoted */
          {
            len = strcspn (ptr, "\n");
            memset ( (void*)ptr, '\n', len);
          }
      }

      // replace un-quoted ';' {iff outside of any braces} by '\n' to allow semi-colons to separate assignments
      if ( (!inQuotes) && (inBracesCount == 0) && ((*ptr) == ';') ) {
        (*ptr) = '\n';
      }
      // replace DOS-style '\r' EOL characters by '\n'
      if ( (*ptr) == '\r' ) {
        (*ptr) = '\n';
      }

      ptr ++;

    } /* while *ptr */

  /*----------------------------------------------------------------------
   * RUN 2: do line-gluing when '\' is found at end-of-line
   */
  ptr = text;
  while ( (ptr = strchr(ptr, '\\')) != NULL )
    {
      if ( ptr[1] == '\n' )
        {
          /* ok, now it gets a bit tricky: to avoid getting spurious spaces from
           * the line-continuation, we shift the rest of the file forward by 2 positions
           * to nicely fit to the previous line...
           */
          len = strlen (ptr+2);
          memmove(ptr, ptr+2, len+1);   /* move the whole rest (add +1 for '\0') */
        }
      else
        {
          ptr ++;
        }
    } /* while '\' found in text */

  /*----------------------------------------------------------------------
   * RUN 3: turn all tabs into single spaces..
   */
  ptr = text;
  while ( (ptr = strchr(ptr, '\t')) != NULL ) {
    *ptr = ' ';
  }

  /*----------------------------------------------------------------------
   * RUN 4: get rid of initial and trailing whitespace (replace it by '\n')
   */
  ptr = text;
  char *endptr = text + strlen(text);   // points to closing '\0' character in input-string
  while (ptr < endptr )
    {
      eol = strchr (ptr, '\n'); /* point to end-of-line */

      len = strspn (ptr, WHITESPACE);
      if (len) { memset ( (void*)ptr, '\n', len); }

      if (eol != NULL) {
        ptr = eol;
      }
      else {
        ptr = strchr (ptr, '\0'); /* or end of file */
      }

      /* clean away all trailing whitespace of last line*/
      ptr2 = ptr - 1;
      while ( ptr2 >= text && ( strspn ( ptr2, WHITESPACE ) != 0 ) ) {
        *ptr2-- = '\n';
      }

      /* step to next line */
      ptr += 1;
    } // while ptr < end

  return;

} /* cleanConfig() */
