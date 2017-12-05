/*
 * Copyright (C) 2015 Reinhard Prix
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
 * \defgroup ConfigFile_h Header ConfigFile.h
 * \ingroup lal_support
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
 * A boolean variable read by XLALReadConfigBOOLEANVariable() can have any of the values
 * <tt>{1, 0, yes, no, true, false}</tt>, where the comparison is done
 * <em>case-insensitively</em>, i.e. you can also use 'True' or 'FALSE'....
 *
 * If one wishes a tight sytnax for the config-file, one can check
 * that there are no illegal entries in the config-file. This is done
 * by checking at the end that all config-file entries have been
 * successfully parsed, using: XLALConfigFileGetUnreadEntries()
 *
 * The configfile-data should be freed at the end using
 * XLALDestroyParsedDataFile().
 *
 * \par Notes
 *
 * XLALReadConfigSTRINGVariable() read the <em>rest</em> of the logical line (excluding comments) as a string,
 * and removes any surrounding quotes \' or \".
 *
 *
 * \note instead of using these functions directly, it might be
 * more convenient to use the \ref UserInput_h.
 *
 */
/*@{*/

/**
 * This structure is returned by XLALParseDataFile() and holds the contents of an
 * ASCII data-file in a pre-parsed form, namely stripped from all comments ('\#', '\%'),
 * spurious whitespaces, and separated into lines (taking into account line-continuation
 * by '\\' at the end of lines).
 * This is used as the input structure in the config-variable reading routines.
 */
typedef struct tagLALParsedDataFile {
  TokenList *lines;     /**< list of pre-parsed data-file lines */
  BOOLEAN *wasRead;     /**< keep track of successfully read lines */
} LALParsedDataFile;


/* Function prototypes */
int XLALParseDataFile (LALParsedDataFile **cfgdata, const CHAR *fname);
int XLALParseDataFileContent (LALParsedDataFile **cfgdata, const CHAR *string );

void XLALDestroyParsedDataFile (LALParsedDataFile *cfgdata);

int XLALConfigSectionExists(const LALParsedDataFile *, const CHAR *);
LALStringVector *XLALListConfigFileSections ( const LALParsedDataFile *cfgdata );
UINT4Vector *XLALConfigFileGetUnreadEntries ( const LALParsedDataFile *cfgdata );

// ---------- type-specific parser prototypes generated via template
#define DECLARE_XLALREADCONFIGVARIABLE(TYPE,CTYPE)                      \
  int XLALReadConfig ##TYPE## Variable ( CTYPE *varp, LALParsedDataFile *cfgdata, const CHAR *secName, const CHAR *varName, BOOLEAN *wasRead )

DECLARE_XLALREADCONFIGVARIABLE(STRING,CHAR*);
DECLARE_XLALREADCONFIGVARIABLE(BOOLEAN,BOOLEAN);
DECLARE_XLALREADCONFIGVARIABLE(INT4,INT4);
DECLARE_XLALREADCONFIGVARIABLE(INT8,INT8);
DECLARE_XLALREADCONFIGVARIABLE(REAL8,REAL8);
DECLARE_XLALREADCONFIGVARIABLE(STRINGVector,LALStringVector*);
DECLARE_XLALREADCONFIGVARIABLE(EPOCH,LIGOTimeGPS);
DECLARE_XLALREADCONFIGVARIABLE(RAJ,REAL8);
DECLARE_XLALREADCONFIGVARIABLE(DECJ,REAL8);
// ------------------------------------------------------------



/*@}*/


/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
