/*
 * Copyright (C) 2004 Reinhard Prix
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

#ifndef _CONFIGFILE_H  /* Double-include protection. */
#define _CONFIGFILE_H

#include <lal/LALDatatypes.h>
#include <lal/StringInput.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \addtogroup ConfigFile_h
 * \author Reinhard Prix
 * \brief Module for general parsing of simple ASCII-based config-files.
 *
 * ### Description ###
 *
 * This module provides routines for reading formatted
 * config-files containing definitions of the form <tt>variable = value</tt>.
 * The general syntax is somewhat similar to the one provided by the
 * perl-module <tt>ConfigParser</tt> (cf.
 * http://www.python.org/doc/current/lib/module-ConfigParser.html)
 *
 * Comments are allowed using either '<tt>\#</tt>' or <tt>\%</tt>.
 * You can also use line-continuation  using a '<tt>\\</tt>' at the end of the line.
 * Also note that comment-signs '<tt>\#\%</tt>' within double-quotes &quot;...&quot;
 * are <em>not</em> treated as comment-characters.
 * Semi-colons <tt>;</tt> are ignored, but can be used to separate several assignments on the same line.
 * The general syntax is best illustrated
 * using a simple example:
 * \code
 * # comment line
 * var1 = 1.0; var2 = 3.1;   ## several assignments on a line, separated by ';'
 * somevar = some text.\
 * You can also use\
 * line-continuation
 * var3 = 4      # whatever that means
 * note = "this is also possible, and # here does nothing"
 * a_switch = true  #possible values: 0,1,true,false,yes,no, case insensitive
 * ...
 * \endcode
 *
 * Note that TABS generally get replaced by a single space, which can be
 * useful in the case of line-continuation (see example). All leading and
 * trailing spaces in are ignore (except within double-quotes).
 *
 * The general approach of reading from such a config-file, is to first
 * call XLALParseDataFile() which loads and pre-parses the contents of the
 * config-file into the structure LALParsedDataFile. Then one can read in
 * config-variables either using one of the type-strict custom-wrappers
 * <tt>XLALReadConfig<TYPE>Variable()</tt> or the general-purpose reading function
 * XLALReadConfigVariable().
 *
 * A boolean variable read by XLALReadConfigBOOLVariable() can have any of the values
 * <tt>{1, 0, yes, no, true, false}</tt>, where the comparison is done
 * <em>case-insensitively</em>, i.e. you can also use 'True' or 'FALSE'....
 *
 * If one wishes a tight sytnax for the config-file, one can check
 * that there are no illegal entries in the config-file. This is done
 * by checking at the end that all config-file entries have been
 * successfully parsed, using:
 * XLALCheckConfigReadComplete(), where \a strictness is either
 * CONFIGFILE_WARN or CONFIGFILE_ERROR.
 * In the first case only a warning is issued, while in the second it is
 * treated as a LAL-error if some config-file entries have not been
 * read-in. (The use of this function is optional).
 *
 * The configfile-data should be freed at the end using
 * XLALDestroyParsedDataFile().
 *
 * \par Notes
 *
 * XLALReadConfigSTRINGVariable() and XLALReadConfigSTRINGVariable() are not
 * the same as using <tt>%quot;\%s&quot;</tt> as a format string, as they read the
 * <em>rest</em> of the logical line (excluding comments) as a string.
 *
 * In the case of XLALReadConfigSTRINGVariable(), the required
 * memory is allocated and has to be freed by the caller, while for
 * XLALReadConfigSTRINGVariable() the caller has to provide a
 * CHARVector of length N, which defines the maximum length of
 * string to be read.
 *
 * \note instead of using these functions directly, it might be
 * more convenient to use the \ref UserInput_h.
 *
 */
/*@{*/

/** Levels of strictness for config-file parsing. */
typedef enum {
  CONFIGFILE_IGNORE = 0,	/**< ignore missing config-variable or unparsed config-entries */
  CONFIGFILE_WARN,		/**< issue a warning but don't report an error. */
  CONFIGFILE_ERROR,		/**< issue an error-message and report a LAL-error */
  CONFIGFILE_LAST
} ConfigStrictness;


/**
 * This structure defines a config-variable to be read in using the
 * general-purpose reading function LALReadConfigVariable().
 */
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagLALConfigVar, secName, varName, fmt));
#endif /* SWIG */
typedef struct tagLALConfigVar {
  const CHAR *secName;          /**< Section name within which to find varName.  May be NULL */
  const CHAR *varName;		/**< Variable-name to be read in the config-file */
  const CHAR *fmt;		/**< Format string for reading (<tt>sscanf()</tt>-style) */
  ConfigStrictness strictness;	/**< what to do if variable not found: ignore, warn, error */
} LALConfigVar;


/**
 * This structure is returned by LALParseDataFile() and holds the contents of an
 * ASCII data-file in a pre-parsed form, namely stripped from all comments ('\#', '\%'),
 * spurious whitespaces, and separated into lines (taking into account line-continuation
 * by '\\' at the end of lines).
 * This is used as the input structure in the config-variable reading routines.
 */
typedef struct tagLALParsedDataFile {
  TokenList *lines;	/**< list of pre-parsed data-file lines */
  BOOLEAN *wasRead;	/**< keep track of successfully read lines for strictness-checking */
} LALParsedDataFile;


/* Function prototypes */
int XLALParseDataFile (LALParsedDataFile **cfgdata, const CHAR *fname);
int XLALParseDataFileContent (LALParsedDataFile **cfgdata, const CHAR *string );

void XLALDestroyParsedDataFile (LALParsedDataFile *cfgdata);

int XLALConfigSectionExists(const LALParsedDataFile *, const CHAR *);
LALStringVector *XLALListConfigFileSections ( const LALParsedDataFile *cfgdata );

int
XLALReadConfigBOOLVariable (BOOLEAN *varp,
                            const LALParsedDataFile *cfgdata,
                            const CHAR *secName,
                            const CHAR *varName,
                            BOOLEAN *wasRead);
int
XLALReadConfigINT4Variable (INT4 *varp,
                           const LALParsedDataFile *cfgdata,
                           const CHAR *secName,
                           const CHAR *varName,
                           BOOLEAN *wasRead);

int
XLALReadConfigREAL8Variable (REAL8 *varp,
                            const LALParsedDataFile *cfgdata,
                            const CHAR *secName,
                            const CHAR *varName,
                            BOOLEAN *wasRead);

int
XLALReadConfigSTRINGVariable (CHAR **varp,
                             const LALParsedDataFile *cfgdata,
                             const CHAR *secName,
                             const CHAR *varName,
                             BOOLEAN *wasRead);

int
XLALReadConfigSTRINGNVariable (CHARVector *varp,
                              const LALParsedDataFile *cfgdata,
                              const CHAR *secName,
                              const CHAR *varName,
                              BOOLEAN *wasRead);

int
XLALReadConfigVariable (void *varp,
                       const LALParsedDataFile *cfgdata,
                       const LALConfigVar *param,
                       BOOLEAN *wasRead);

int XLALCheckConfigReadComplete (const LALParsedDataFile *cfgdata, ConfigStrictness strict);

/* ========== DEPRECATED LAL INTERFACE FUNCTIONS, which have been replaced by XLAL functions,
 * These functions are just wrappers around the XLAL functions
 */


/** \name Error codes */
/*@{*/
#define CONFIGFILEH_ENULL               1       /**< Arguments contained an unexpected null pointer. */
#define CONFIGFILEH_EFILE               2       /**< File error. */
#define CONFIGFILEH_EVAR                3       /**< Config variable not found. */
#define CONFIGFILEH_EFMT                4       /**< Config variable not readable using given format-string. */
#define CONFIGFILEH_ETOKENS             5       /**< The input ConfigData seems corrupted. */
#define CONFIGFILEH_ENONULL             6       /**< Output pointer is not NULL */
#define CONFIGFILEH_EUNKNOWN            8       /**< Unknown config-file entry found */
#define CONFIGFILEH_EMEM                9       /**< Out of memory */
#define CONFIGFILEH_EBOOL               10      /**< Illegal BOOLEAN entry */
#define CONFIGFILEH_ESTRING             11      /**< Malformed quoted string */
#define CONFIGFILEH_EXLAL               12      /**< Failure in XLAL function */
/*@}*/

/** \cond DONT_DOXYGEN */
#define CONFIGFILEH_MSGENULL            "Arguments contained an unexpected null pointer."
#define CONFIGFILEH_MSGEFILE		"File error."
#define CONFIGFILEH_MSGEVAR		"Config variable not found."
#define CONFIGFILEH_MSGEFMT		"Config variable not readable using given format-string."
#define CONFIGFILEH_MSGETOKENS		"The input ConfigData seems corrupted."
#define CONFIGFILEH_MSGENONULL		"Output pointer is not NULL"
#define CONFIGFILEH_MSGEUNKNOWN		"Unknown config-file entry found"
#define CONFIGFILEH_MSGEMEM		"Out of memory"
#define CONFIGFILEH_MSGEBOOL		"Illegal BOOLEAN entry"
#define CONFIGFILEH_MSGESTRING		"Malformed quoted string"
#define CONFIGFILEH_MSGEXLAL		"Failure in XLAL function"
/** \endcond */

/**
 * \name Deprecated LAL-interface
 * These functions are deprecated, and you should user their XLAL-equivalents instead.
 */
/*@{*/
void LALParseDataFile (LALStatus *, LALParsedDataFile **cfgdata, const CHAR *fname);
void LALDestroyParsedDataFile (LALStatus *, LALParsedDataFile **cfgdata);

void LALReadConfigBOOLVariable (LALStatus *, BOOLEAN *varp, const LALParsedDataFile *cfgdata, const CHAR *varName, BOOLEAN *wasRead);
void LALReadConfigINT4Variable (LALStatus *, INT4 *varp, const LALParsedDataFile *cfgdata, const CHAR *varName, BOOLEAN *wasRead);
void LALReadConfigREAL8Variable (LALStatus *, REAL8 *varp, const LALParsedDataFile *cfgdata, const CHAR *varName, BOOLEAN *wasRead);
void LALReadConfigSTRINGVariable (LALStatus *, CHAR **varp, const LALParsedDataFile *cfgdata, const CHAR *varName, BOOLEAN *wasRead);
void LALReadConfigSTRINGNVariable (LALStatus *, CHARVector *varp, const LALParsedDataFile *cfgdata, const CHAR *varName, BOOLEAN *wasRead);
void LALReadConfigVariable (LALStatus *, void *varp, const LALParsedDataFile *cfgdata, const LALConfigVar *param, BOOLEAN *wasRead);
void LALCheckConfigReadComplete (LALStatus *, const LALParsedDataFile *cfgdata, ConfigStrictness strict);
/*@}*/

/*@}*/


/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
