/*
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
#include <lal/ParseStringValue.h>

#include <lal/ConfigFile.h>

// ---------- local defines ----------
#define FMT_STRING "string"    /* reading in quoted strings needs some special treatment */
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
 * Parser for config-file: can read config-variables of the form
 * VARIABLE [=:] VALUE.
 * Input is a TokenList containing the 'logical' lines of the cleaned config-file
 *
 * - <tt>param->varName</tt> is the name of the config-variable to read
 * - <tt>param->fmt</tt>     is the format string to use for reading
 * - <tt>param->secName</tt> is the section name within which to find varName. NULL means 'default' section.
 * - <tt>param->strictness</tt>   what to do if variable not found: ignore, warn, error
 *
 * \note a special format-string is FMT_STRING, which means read the whole remaining line (without initial whitespace!)
 * which is different from \"\%s\"! (which would read only one word)
 * In this case, this also does the memory-allocation!
 *
 */
int
XLALReadConfigVariable ( void *varp,                      /**< [out] result gets written here! */
                         const LALParsedDataFile *cfgdata, /**< [in] pre-parsed config-data */
                         const LALConfigVar *param,       /**< [in]  var-name, section-name, fmt-string, strictness */
                         BOOLEAN *wasRead)                /**< [out] did we succeed in reading? */
{
  // check input consistency
  XLAL_CHECK ( varp != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata->lines != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (cfgdata->lines->nTokens == 0) || (cfgdata->wasRead != NULL), XLAL_EINVAL );
  XLAL_CHECK ( param->varName != NULL, XLAL_EINVAL );
  XLAL_CHECK ( param->fmt != NULL, XLAL_EINVAL );

  CHAR *found    = NULL;
  INT2 ret = 0;

  UINT4 i;
  INT4 linefound     = -1;
  INT4 section_found = -1;

  size_t len;
  size_t searchlen = strlen (param->varName);
  size_t sec_searchlen = 0;

  if (param->secName == NULL)
    {
      /* If we haven't been asked for a section then we want the
         "default" section, which starts at the top of the file */
      sec_searchlen = 0;
      section_found = 1;
    } else {
      sec_searchlen = strlen(param->secName);
    }

  *wasRead = FALSE;

  /* let's look for the variable-name in the token-list (has to at beginning of line!) */
  for (i = 0; i < cfgdata->lines->nTokens; i++)
    {
      /* Is this the start of a new section? */
      if (cfgdata->lines->tokens[i][0] == '[')
        {
          /* If we're looking for a particular section, is this it? */
          if (sec_searchlen > 0 && strncmp(cfgdata->lines->tokens[i] + 1, param->secName, sec_searchlen) == 0)
            {
              section_found = i;
            }
          else
            section_found = -1;  /* We might have moved out of the right section */
         }
      else if (section_found > -1) /* Not section start, are we in the one we want? */
        {
          len = strcspn (cfgdata->lines->tokens[i], WHITESPACE "=:");	/* get length of variable-name */
          if (len == 0)
                {			/* malformed token-list */
                  XLAL_ERROR ( XLAL_EDOM, "Parsing error: nonexistent variable name in '%s'\n", cfgdata->lines->tokens[i] );
                }

          /* pre-select based on length of variable-name */
          if (len != searchlen)
            continue;

          /* same len, but are they identical ? */
          if (strncmp (param->varName, cfgdata->lines->tokens[i], len) == 0)
                {
                  found = cfgdata->lines->tokens[i] + len;
                  found += strspn (found, WHITESPACE "=:");	/* skip all whitespace and define-chars */
                  linefound = i;
                  break;		/* ok, we've found it */
                }
        }           /* if section found */
    }				/* for lines */

  if (!found)
    {
      switch (param->strictness)
        {
        case CONFIGFILE_IGNORE:
          return 0;
          break;
        case CONFIGFILE_WARN:
          if (lalDebugLevel & LALWARNING)
            {
              if (sec_searchlen > 0)
                XLALPrintError ("%s: Warning: Config-file variable '%s' in section '%s' was not found!\n",
                                __func__, param->varName, param->secName);
              else
                XLALPrintError ("%s: Warning: Config-file variable '%s' was not found!\n", __func__, param->varName );
            }
          return 0;
          break;
        case CONFIGFILE_ERROR:
        default:
          if (sec_searchlen > 0) {
            XLALPrintError ( "%s: Error: Config-file variable '%s' in section '%s' was not found!\n", __func__, param->varName, param->secName);
          } else {
            XLALPrintError ( "%s: Error: Config-file variable '%s' was not found!\n", __func__, param->varName);
          }
          XLAL_ERROR ( XLAL_EDOM );
          break;
        } /* switch (strictness) */

    } /* if not found */

  /* now read the value into the variable */

  /* reading a quoted string needs some special treatment: */
  if (!strcmp (param->fmt, FMT_STRING))
    {
      /* NOTE: varp here is supposed to be a pointer to CHAR* !! */
      CHAR **cstr = (CHAR **) varp;

      XLAL_CHECK ( (*cstr) == NULL, XLAL_EINVAL );
      (*cstr) = (CHAR *) XLALMalloc (strlen (found) + 1);
      strcpy ((*cstr), found);
      ret = 1;
    }
  else				/* but the default case is just sscanf... */
    ret = sscanf (found, param->fmt, varp);

  if ((ret == 0) || (ret == EOF))
    {
      XLALPrintError ("%s: ERROR: Config-file variable %s was not readable using the format %s\n\n", __func__, param->varName, param->fmt);
      XLAL_ERROR ( XLAL_EINVAL );
    }

  /* ok, we have successfully read in the config-variable: let's make a note of it */
  cfgdata->wasRead[linefound] = 1;

  *wasRead = TRUE;

  return XLAL_SUCCESS;

} /* XLALReadConfigVariable() */


/**
 * Type-specialization of generic reading-function XLALReadConfigVariable() to BOOLEAN variables.
 */
int
XLALReadConfigBOOLVariable (BOOLEAN *varp,              /**< [out] variable to store result */
                            const LALParsedDataFile *cfgdata,   /**< [in] pre-parsed config-data */
                            const CHAR *secName,        /**< [in] section-name to read */
                            const CHAR *varName,        /**< [in] variable-name to read */
                            BOOLEAN *wasRead            /**< [out] did we succeed in reading? */
                            )
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  /* first read the value as a string */
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsBOOLEAN ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree (valString);

  return XLAL_SUCCESS;

} /* XLALReadConfigBOOLVariable() */


/**
 * Read a signed integer value.
 *
 * Note: Rather than using the sscanf() for the string-to-integer conversion,
 * we do our own parsing, in order to check for spurious extra characters at the end.
 * This allows us to catch user-input mistakes like specifying "1e3" for an int-variable, which
 * would silently be converted as '1' by sscanf("%ld").
 */
int
XLALReadConfigINT4Variable (INT4 *varp,
                            const LALParsedDataFile *cfgdata,
                            const CHAR *secName,
                            const CHAR *varName,
                            BOOLEAN *wasRead)
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsINT4 ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree ( valString );

  return XLAL_SUCCESS;

} /* XLALReadConfigINT4Variable() */


/**
 * Type-specialization of generic reading-function XLALReadConfigVariable() to REAL8 variables.
 */
int
XLALReadConfigREAL8Variable (REAL8 *varp,
                             const LALParsedDataFile *cfgdata,
                             const CHAR *secName,
                             const CHAR *varName,
                             BOOLEAN *wasRead)
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsREAL8 ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree ( valString );

  return XLAL_SUCCESS;

} /* XLALReadConfigREAL8Variable() */


/**
 * Type-specialization of generic reading-function XLALReadConfigVariable() to
 * STRING variables
 * \note this means the rest of the line (skipping initial whitespace and excluding trailing comments), NOT "%s"!
 * \par Note2: if string is quoted by ", everything within quotes is read, and the quotes are removed here
 *
 */
int
XLALReadConfigSTRINGVariable (CHAR ** varp,		/**< [out] string, allocated here! */
                              const LALParsedDataFile * cfgdata, /**< [in] pre-parsed config-data */
                              const CHAR * secName,	/**< [in] section-name to be read */
                              const CHAR * varName,	/**< [in] variable-name to be read */
                              BOOLEAN * wasRead)	/**< [out] did we succeed in reading? */
{
  XLAL_CHECK ( (*varp) == NULL, XLAL_EINVAL );

  LALConfigVar param = { 0, 0, 0, CONFIGFILE_IGNORE };
  CHAR *str = NULL;
  CHAR *ret = NULL;

  param.secName = secName;
  param.varName = varName;
  param.fmt = FMT_STRING;
  param.strictness = CONFIGFILE_IGNORE;

  if ( XLALReadConfigVariable ((void *) &str, cfgdata, &param, wasRead) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALReadConfigVariable() failed\n", __func__ );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  if (*wasRead && (str != NULL))
    {
      INT2 numQuotes = 0;
      CHAR *ptr = str;
      /* count number of quotation marks */
      while ((ptr = strchr (ptr, '"')))
        {
          numQuotes++;
          ptr++;
        }			/* while quotes found */

      /* check balanced quotes (don't allow escaping for now) */
      if ((numQuotes != 0) && (numQuotes != 2))
        {
          XLAL_ERROR ( XLAL_EDOM, "Parsing error: unmatched quotes\n");
        }
      if (numQuotes == 2)
        {
          /* allowed only at end and beginning */
          if ((str[0] != '"') || (str[strlen (str) - 1] != '"'))
            {
              XLAL_ERROR ( XLAL_EDOM, "Parsing error: quotes only allowed beginning and end of string '%s'\n", str );
            }
          /* quotes ok, now remove them */
          XLAL_CHECK ( (ret = XLALMalloc (strlen (str) - 2 + 1)) != NULL, XLAL_ENOMEM );
          str[strlen (str) - 1] = 0;
          strcpy (ret, str + 1);
          XLALFree (str);
        }			/* if 2 quotation marks */
      else
        ret = str;		/* no quotes, just return string */

      *varp = ret;

    }				/* if wasRead */
  else
    *varp = NULL;

  return XLAL_SUCCESS;

} /* XLALReadConfigSTRINGVariable() */

/**
 * Read a GPS 'epoch' value, specified either as GPS or MJD(TT) string.
 *
 */
int
XLALReadConfigEPOCHVariable ( LIGOTimeGPS *varp,
                              const LALParsedDataFile *cfgdata,
                              const CHAR *secName,
                              const CHAR *varName,
                              BOOLEAN *wasRead)
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( varp, valString ) != NULL, XLAL_EFUNC );

  XLALFree ( valString );

  return XLAL_SUCCESS;

} /* XLALReadConfigEPOCHVariable() */


/**
 * Type-specialization of generic reading-function XLALReadConfigVariable() to RAJ variables
 * (allowing for either radians or "hours:minutes:seconds" input format).
 */
int
XLALReadConfigRAJVariable ( REAL8 *varp,
                            const LALParsedDataFile *cfgdata,
                            const CHAR *secName,
                            const CHAR *varName,
                            BOOLEAN *wasRead)
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsRAJ ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree ( valString );

  return XLAL_SUCCESS;

} /* XLALReadConfigRAJVariable() */

/**
 * Type-specialization of generic reading-function XLALReadConfigVariable() to DECJ variables
 * (allowing for either radians or "degrees:minutes:seconds" input format).
 */
int
XLALReadConfigDECJVariable ( REAL8 *varp,
                             const LALParsedDataFile *cfgdata,
                             const CHAR *secName,
                             const CHAR *varName,
                             BOOLEAN *wasRead)
{
  (*wasRead) = FALSE;
  /* first read the value as a string */
  CHAR *valString = NULL;
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfgdata, secName, varName, wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  if ( ! (*wasRead ) ) { // if nothing was read and XLALReadConfigSTRINGVariable() didn't throw an error, then we're ok and return
    return XLAL_SUCCESS;
  }
  XLAL_CHECK ( valString != NULL, XLAL_EFAILED, "Got NULL string after reading config-variable '%s' in section '%s'\n", varName, secName ? secName: "default" );

  XLAL_CHECK ( XLALParseStringValueAsDECJ ( varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALFree ( valString );

  return XLAL_SUCCESS;

} /* XLALReadConfigDECJVariable() */


/**
 * Check if all lines of config-file have been successfully read in
 * and issue a warning or error (depending on strictness) if not.
 */
int
XLALCheckConfigReadComplete (const LALParsedDataFile *cfgdata,  /**< [in] config-file data */
                             ConfigStrictness strict)           /**< [in] what to do if unparsed lines */
{
  XLAL_CHECK ( cfgdata != NULL, XLAL_EINVAL );
  XLAL_CHECK ( cfgdata->lines != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (cfgdata->lines->nTokens == 0) || (cfgdata->wasRead != NULL), XLAL_EINVAL );

  for (UINT4 i=0; i < cfgdata->lines->nTokens; i++)
    {
      /* Don't require section headers to be marked read */
      /* This has the effect of considering a file to have */
      /* been read completely if every value in every section */
      /* has been read. */
      if (cfgdata->lines->tokens[i][0] == '[')
        continue;

      if (cfgdata->wasRead[i] == 0)
        {
          switch (strict)
            {
            case CONFIGFILE_IGNORE:
              continue;
            case CONFIGFILE_WARN:
              XLALPrintError ( "%s: Warning: Ignoring unknown config-file entry '%s'.\n", __func__, cfgdata->lines->tokens[i] );
              continue;
            case CONFIGFILE_ERROR:
            default:
              XLALPrintError ( "%s: ERROR: config-file entry #%d has not been read!\n", __func__, i);
              XLALPrintError ( "Line was: '%s'\n", cfgdata->lines->tokens[i]);
              XLAL_ERROR ( XLAL_EDOM );
            } /* switch strict */
        } /* if some line not read */

    } /* for i < lines */

  return XLAL_SUCCESS;

} /* XLALCheckConfigReadComplete() */

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
