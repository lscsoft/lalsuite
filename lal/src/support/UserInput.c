/*
 * Copyright (C) 2015 Reinhard Prix
 * Copyright (C) 2010 Reinhard Prix (xlalified)
 * Copyright (C) 2004, 2005, 2015 Reinhard Prix
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

// ---------- local includes ----------
#include <lal/LALStdio.h>
#include <lal/LALgetopt.h>
#include <lal/LogPrintf.h>
#include <lal/LALString.h>
#include <lal/Date.h>
#include <lal/StringVector.h>
#include <lal/AVFactories.h>

#include <lal/UserInputParse.h>
#include <lal/UserInputPrint.h>

#include <lal/UserInput.h>
// ---------- local defines ----------
#define TRUE  (1==1)
#define FALSE (1==0)


// ---------- local Macro definitions ----------

// ----- macro template for defining registration functions for UserInput variables
#define DEFN_REGISTER_UVAR(UTYPE,CTYPE)                     \
DECL_REGISTER_UVAR(UTYPE,CTYPE)                             \
{                                                                       \
  return XLALRegisterUserVar (name, UVAR_TYPE_ ## UTYPE, optchar, category, helpstr, cvar); \
}

// ---------- local type definitions ----------

// Define the type of a "user-variable": bool, int, real or string, ...
typedef enum {
  UVAR_TYPE_START=0, 		// internal start marker for range checking

  UVAR_TYPE_BOOLEAN, 		// boolean
  UVAR_TYPE_INT4,    		// integer
  UVAR_TYPE_REAL8,   		// float
  UVAR_TYPE_EPOCH,   		// time 'epoch', specified in either GPS or MJD(TT) format, translated into GPS
  UVAR_TYPE_RAJ,     		// sky equatorial longitude (aka right-ascencion or RA), in either radians or hours:minutes:seconds format, translated into radians
  UVAR_TYPE_DECJ,    		// sky equatorial latitude (aka declination or DEC), in either radians or degrees:minutes:seconds format, translated into radians

  UVAR_TYPE_STRING, 		// normal string
  UVAR_TYPE_STRINGVector,	// list of comma-separated strings
  UVAR_TYPE_REAL8Vector,	// list of comma-separated REAL8's
  UVAR_TYPE_INT4Vector,		// list of comma-separated INT4's

  UVAR_TYPE_END      	// internal end marker for range checking
} UserVarType;

//
// Linked list to hold the complete information about the user-variables.
//
typedef struct tagLALUserVariable {
  const CHAR *name;			// full name
  UserVarType type;			// variable type: BOOLEAN, INT4, REAL8, ...
  CHAR optchar;				// cmd-line character
  const CHAR *help;			// help-string
  void *varp;				// pointer to the actual C-variable
  UserVarCategory category;		// category (optional, required, developer, ... )
  BOOLEAN was_set;			// was this set by the user in any way? (ie vie cmdline or cfg-file)
  struct tagLALUserVariable *next; 	// linked list
} LALUserVariable;

// ---------- local prototypes ----------
int XLALRegisterUserVar (const CHAR *name, UserVarType type, CHAR optchar, UserVarCategory category, const CHAR *helpstr, void *cvar);
void check_and_mark_as_set ( LALUserVariable *varp );

// ----- define templated registration functions for all supported UVAR_TYPE_ 'UTYPES'
DEFN_REGISTER_UVAR(BOOLEAN,BOOLEAN);
DEFN_REGISTER_UVAR(INT4,INT4);
DEFN_REGISTER_UVAR(REAL8,REAL8);
DEFN_REGISTER_UVAR(RAJ,REAL8);
DEFN_REGISTER_UVAR(DECJ,REAL8);
DEFN_REGISTER_UVAR(EPOCH,LIGOTimeGPS);
DEFN_REGISTER_UVAR(STRING,CHAR*);
DEFN_REGISTER_UVAR(STRINGVector,LALStringVector*);
DEFN_REGISTER_UVAR(REAL8Vector,REAL8Vector*);
DEFN_REGISTER_UVAR(INT4Vector,INT4Vector*);

// ----- define helper types for casting
typedef int (*parserT)(void*, const char*);
typedef void (*destructorT)(void*);
typedef char *(*printerT)(const void*);

// ----- handy macro to simplify adding 'regular' entries for new UTYPES into UserVarTypeMap
#define REGULAR_MAP_ENTRY(UTYPE,DESTRUCTOR)                             \
  [UVAR_TYPE_##UTYPE] = { #UTYPE, (parserT)XLALParseStringValueAs##UTYPE,	(printerT)XLALPrintStringValueOf##UTYPE, (destructorT)DESTRUCTOR }

// ---------- HOWTO add new UserInput variable types ----------
// In order to add a new type \<UTYPE\> to be handled by the UserInput module, you just need to
// 1) add an entry 'UVAR_TYPE_\<UTYPE\>' in the UserVarType enum
// 2) provide
//   a)  a parser function XLALParseStringValueAs\<UTYPE\>() (recommended to be added in \ref UserInputParse_h)
//   b)  a printer function XLALPrintStringValueOf\<UTYPE\>() (recommended to be added in \ref UserInputPrint_h)
//   c)  a unit test for the new parser+printer, ideally checking identity of print(parse(x)) or parse(print(x))
// 3) generate a corresponding registration function declaration + definition using the macro-templates
//    DECL_REGISTER_UVAR_AS<VALUE|POINTER>() and DEFN_REGISTER_UVAR_AS<VALUE|POINTER>(),
// 4) add an entry in the master map 'UserInputTypeMap', specifying the parser, printer and (if required) a destructor
//    If these follow the standard naming and API, the template macro REGULAR_MAP_ENTRY() can be used for that.
//
// ---------- Master 'map' defining all UserInput types and how to handle them ----------
// in particular, lists their name, and how to parse and print them, and (if required) how to destroy them
static const struct
{
  const char *const name;			///< type name
  int (*parser)(void*, const char*);		///< parser function to parse string as this type
  char *(*printer)(const void *);		///< 'printer' function returning string value for given type
  void (*destructor)(void*);			///< destructor for this variable type, NULL if none required
} UserVarTypeMap[UVAR_TYPE_END]
= {
  // either use 'manual' entries of the form
  // [UVAR_TYPE_\<UTYPE\>] = { "\<UTYPE\>",	(parserT)XLALParseStringValueAs\<UTYPE\>, (printerT)XLALPrintStringValueOf\<UTYPE\>, (destructorT)XLALDestroy\<UTYPE\> },
  // or the convenience macro for cases using 'standard' function names and API
  // REGULAR_MAP_ENTRY ( \<UTYPE\>, XLALDestroy\<UTYPE\> ),
  REGULAR_MAP_ENTRY ( BOOLEAN, NULL ),
  REGULAR_MAP_ENTRY ( INT4, NULL ),
  REGULAR_MAP_ENTRY ( REAL8, NULL ),
  REGULAR_MAP_ENTRY ( STRING, XLALFree ),
  REGULAR_MAP_ENTRY ( STRINGVector, XLALDestroyStringVector ),
  REGULAR_MAP_ENTRY ( EPOCH, NULL ),
  REGULAR_MAP_ENTRY ( RAJ, NULL ),
  REGULAR_MAP_ENTRY ( DECJ, NULL ),
  REGULAR_MAP_ENTRY ( REAL8Vector, XLALDestroyREAL8Vector ),
  REGULAR_MAP_ENTRY ( INT4Vector, XLALDestroyINT4Vector )
};


// ---------- The module-local linked list to hold the user-variables
static LALUserVariable UVAR_vars;	// empty head
static const CHAR *program_name;	// keep a pointer to the program name


// ==================== Function definitions ====================

/**
 * \ingroup UserInput_h
 * Internal function: Register a user-variable with the module.
 * Effectively put an appropriate entry into UVAR_vars
 *
 * Checks that long- and short-options are unique, an error is returned
 * if a previous option name collides.
 *
 * \note don't use this function directly, as it is not type-safe!!
 * ==> use the type-safe macro XLALRegisterUvarMember(name,type,option,category,help) instead!
 */
int
XLALRegisterUserVar ( const CHAR *name,		/**< name of user-variable to register */
                      UserVarType type,		/**< variable type (int,bool,string,real) */
                      CHAR optchar,		/**< optional short-option character */
                      UserVarCategory category,		/**< sets category to this */
                      const CHAR *helpstr,	/**< help-string explaining this input-variable */
                      void *cvar		/**< pointer to the actual C-variabe to link to this user-variable */
                      )
{
  XLAL_CHECK ( name != NULL, XLAL_EINVAL );
  XLAL_CHECK ( (category > UVAR_CATEGORY_START) && (category < UVAR_CATEGORY_END), XLAL_EINVAL );
  XLAL_CHECK ( cvar != NULL, XLAL_EINVAL );

  // find end of uvar-list && check that neither short- nor long-option are taken already
  LALUserVariable *ptr = &UVAR_vars;
  while ( ptr->next != NULL )
    {
      ptr = ptr->next;

      // long-option name taken already?
      XLAL_CHECK ( (name == NULL) || (ptr->name == NULL) || (strcmp ( name, ptr->name ) != 0), XLAL_EINVAL, "Long-option name '--%s' already taken!\n", name );
      // short-option character taken already?
      XLAL_CHECK ( (optchar == 0) || (ptr->optchar == 0) || (optchar != ptr->optchar), XLAL_EINVAL, "Short-option '-%c' already taken (by '--%s')!\n", ptr->optchar, ptr->name );

    } // while ptr->next

  // append new entry at the end
  XLAL_CHECK ( (ptr->next = XLALCalloc (1, sizeof(LALUserVariable))) != NULL, XLAL_ENOMEM );

  // set pointer to newly created entry and fill in values
  ptr = ptr->next;

  ptr->name 	= name;
  ptr->type 	= type;
  ptr->optchar 	= optchar;
  ptr->help 	= helpstr;
  ptr->varp 	= cvar;
  ptr->category = category;

  return XLAL_SUCCESS;

} // XLALRegisterUserVar()

/**
 * Free all memory associated with user-variable linked list
 */
void
XLALDestroyUserVars ( void )
{
  LALUserVariable *ptr = &(UVAR_vars);
  LALUserVariable *lastptr = NULL;

  // step through user-variables: free list-entries and all allocated strings
  while ( (ptr=ptr->next) != NULL )
    {
      XLAL_CHECK_VOID ( (ptr->type > UVAR_TYPE_START) && (ptr->type < UVAR_TYPE_END), XLAL_EFAILED, "Invalid UVAR_TYPE '%d' outside of [%d,%d]\n", ptr->type, UVAR_TYPE_START+1, UVAR_TYPE_END-1 );

      // is there a destructor function registered for this type?
      if ( UserVarTypeMap [ ptr->type ].destructor != NULL )
        {
          UserVarTypeMap [ ptr->type ].destructor ( *(CHAR**)ptr->varp );
          *(CHAR**)ptr->varp = NULL;
        }

      /* free list-entry behind us (except for the head) */
      if ( lastptr != NULL ) {
        XLALFree ( lastptr );
      }

      lastptr = ptr;

    } // while ptr->next

  if ( lastptr != NULL ) {
    XLALFree ( lastptr );
  }

  // clean head
  memset (&UVAR_vars, 0, sizeof(UVAR_vars));

  return;

} /* XLALDestroyUserVars() */


/**
 * Parse command-line into UserVariable array
 */
int
XLALUserVarReadCmdline ( int argc, char *argv[] )
{
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Input error, NULL argv[] pointer passed.\n" );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "Internal error, no UVAR memory allocated. Did you register any user-variables?" );

  LALUserVariable *ptr;
  UINT4 pos;

  // ---------- build optstring of short-options
  UINT4 numvars = 0;
  char optstring[512] = "\0";	// string of short-options
  ptr = &UVAR_vars;	// set to empty head
  pos = 0;
  while ( (ptr = ptr->next) != NULL )
    {
      numvars ++;			/* counter number of user-variables */
      if (ptr->optchar == 0) {		/* if no short-option given, ignore */
	continue;
      }
      optstring[pos++] = ptr->optchar;
      optstring[pos++] = ':';		/* everything but bool takes an argument */
      if (ptr->type == UVAR_TYPE_BOOLEAN) {	/* but for BOOL its optional */
	optstring[pos++] = ':';
      }
    } // while ptr->next
  optstring[pos] = '\0';

  // ---------- fill option-struct for long-options
  struct LALoption *long_options = LALCalloc (1, (numvars+1) * sizeof(struct LALoption));
  ptr = &UVAR_vars;	// start again from beginning: empty head
  pos = 0;
  while ( (ptr= ptr->next) != NULL)
    {
      if (ptr->name == NULL) {		/* if no long-option given, ignore */
	continue;
      }
      long_options[pos].name 	= ptr->name;
      long_options[pos].has_arg = (ptr->type == UVAR_TYPE_BOOLEAN) ? optional_argument : required_argument;
      long_options[pos].flag = NULL;	// get val returned from LALgetopt_long()
      long_options[pos].val 	= 0;	// we use longindex to find long-options
      pos ++;
    } // while ptr->next

  // null-terminate array
  long_options[pos].name = 0;
  long_options[pos].has_arg = 0;
  long_options[pos].flag = 0;
  long_options[pos].val = 0;

  /* NOTE: in case we get called several times, we have to make sure here that getopt() gets
   * properly reset/initialized. We do this using the (undocumented) feature of GNU getopt
   * of setting optind to 0. As we're linking our private version of GNU getopt, this should be
   * guaranteed to work.
   *
   * Bruce's notes: read LALgetopt_long() source code, and in particular
   * _getopt_internal() to see what is initialized.
   */
  LALoptind = 0; 	// reset our local LALgetopt(), LALgetopt_long()

  // ---------- parse the command-line
  int longindex = -1;
  INT4 c;
  while ( (c = LALgetopt_long(argc, argv, optstring, long_options, &longindex)) != -1 )
    {
      XLAL_CHECK (c != '?', XLAL_EINVAL, "Unkown command-line option encountered, see '%s --help' for usage-help\n\n", argv[0] );
      if (c != 0) 	// find short-option character
	{
	  ptr = &UVAR_vars;
	  do {
	    if (c == ptr->optchar) {
	      break;
            }
	  } while ( (ptr=ptr->next) != NULL);
	} // end: if short-option given
      else	// find long-option: returned in longindex
	{
	  ptr = &UVAR_vars;
	  while ( (ptr=ptr->next) != NULL) {
	    if ( !strcmp (long_options[longindex].name, ptr->name) ) {
	      break;
            }
          }
	} // end: if long-option

      XLAL_CHECK ( ptr != NULL, XLAL_EFAILED, "ERROR: failed to find matching option ... this points to a coding-error!\n" );
      XLAL_CHECK ( (ptr->type > UVAR_TYPE_START) && (ptr->type < UVAR_TYPE_END), XLAL_EFAILED, "Invalid UVAR_TYPE '%d' outside of [%d,%d]\n", ptr->type, UVAR_TYPE_START+1, UVAR_TYPE_END-1 );

      switch (ptr->type)
	{
	case UVAR_TYPE_BOOLEAN:
	  // subtlety with optional arguments: it's not necessarily found in the *same* argv-entry
          // eg, if no '=' was used, so we have to check for that case by hand:
	  // if the next entry is not an option, take it as an argument
	  if ( (LALoptarg == NULL) && (LALoptind < argc) && (argv[LALoptind][0] != '-') && (argv[LALoptind][0] != '@') )
            {
              LALoptarg = argv[LALoptind];
              LALoptind ++;
            }

	  if ( LALoptarg == NULL ) { // if no argument given, defaults to TRUE
            *(BOOLEAN*)(ptr->varp) = TRUE;
          } else {
            XLAL_CHECK ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
          }
	  break;

	default:
          // all other UVAR_TYPE_ types can be handled canonically: first destroy previous value, the parse new one
          if ( UserVarTypeMap [ ptr->type ].destructor != NULL )
            {
              UserVarTypeMap [ ptr->type ].destructor( *(char**)ptr->varp );
              *(char**)ptr->varp = NULL;
            } // if a destructor was registered
          XLAL_CHECK ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
	  break;

	} // switch ptr->type

      check_and_mark_as_set ( ptr );

    } // while LALgetopt_long()

  // ---------- check if there's any non-option strings left (except for a config-file specification '@file')
  if ( (LALoptind == argc - 1) && (argv[LALoptind][0] == '@' ) ) {
    LALoptind ++;	// advance counter in case of one config-file specification (only one allowed)
  }
  if ( LALoptind < argc ) // still stuff left? ==> error
    {
      XLALPrintError ( "\nGot non-option ARGV-elements: [ ");
      while (LALoptind < argc) {
        if ( argv[LALoptind][0] == '@' ) { LALoptind ++; continue; }	// don't list config-file entries here
        XLALPrintError ("%s ", argv[LALoptind++]);
      }
      XLALPrintError(" ]\n");
      XLAL_ERROR ( XLAL_EINVAL );
    } // trailing non-option arguments found

  XLALFree (long_options);
  long_options=NULL;

  return XLAL_SUCCESS;

} // XLALUserVarReadCmdline()


/**
 * Read config-variables from cfgfile and parse into input-structure.
 * An error is reported if the config-file reading fails, but the
 * individual variable-reads are treated as optional
 */
int
XLALUserVarReadCfgfile ( const CHAR *cfgfile ) 	   /**< [in] name of config-file */
{
  XLAL_CHECK ( cfgfile != NULL, XLAL_EINVAL );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No memory allocated in UVAR_vars.next, did you register any user-variables?\n" );

  LALParsedDataFile *cfg = NULL;
  XLAL_CHECK ( XLALParseDataFile ( &cfg, cfgfile ) == XLAL_SUCCESS, XLAL_EFUNC );

  // step through all user-variable: read those with names from config-file
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->name == NULL) {	// ignore name-less user-variable
	continue;
      }

      XLAL_CHECK ( (ptr->type > UVAR_TYPE_START) && (ptr->type < UVAR_TYPE_END), XLAL_EFAILED, "Invalid UVAR_TYPE '%d' outside of [%d,%d]\n", ptr->type, UVAR_TYPE_START+1, UVAR_TYPE_END-1 );

      BOOLEAN wasRead;
      CHAR *valString = NULL;       // first read the value as a string
      XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( wasRead ) {	// if successful, parse this as the desired type
        XLAL_CHECK ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLALFree (valString);
        check_and_mark_as_set ( ptr );
      }

    } // while ptr->next

  // ok, that should be it: check if there were more definitions we did not read
  UINT4Vector *unread = XLALConfigFileGetUnreadEntries ( cfg );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALConfigFileGetUnreadEntries() failed\n");
  if ( unread != NULL )
    {
      XLALPrintWarning ("The following entries in config-file '%s' have not been parsed:\n", cfgfile );
      for ( UINT4 i = 0; i < unread->length; i ++ ) {
        XLALPrintWarning ("%s\n", cfg->lines->tokens[ unread->data[i] ] );
      }
      XLALDestroyUINT4Vector ( unread );
    }

  XLALDestroyParsedDataFile ( cfg );

  return XLAL_SUCCESS;

} // XLALUserVarReadCfgfile()

/**
 * Assemble all help-info from uvars into a help-string.
 */
CHAR *
XLALUserVarHelpString ( const CHAR *progname )
{
  XLAL_CHECK_NULL ( progname != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?\n" );

  BOOLEAN showDeveloperOptions  = (lalDebugLevel >= LALWARNING);	// only output for lalDebugLevel >= warning
  BOOLEAN showDeprecatedOptions  = (lalDebugLevel >= LALINFO);		// only output for lalDebugLevel >= info

  // ---------- ZEROTH PASS: find longest long-option and type names, for proper output formatting
  LALUserVariable *ptr = &UVAR_vars;
  UINT4 nameFieldLen = 0;
  UINT4 typeFieldLen = 0;
  BOOLEAN haveDeveloperOptions = 0;
  BOOLEAN haveDeprecatedOptions = 0;

  while ( (ptr=ptr->next) != NULL )
    {
      if (ptr->category == UVAR_CATEGORY_DEVELOPER) {
        haveDeveloperOptions = 1;
      }
      if ( ptr->category == UVAR_CATEGORY_DEPRECATED ) {
        haveDeprecatedOptions = 1;
      }
      if ( (ptr->category == UVAR_CATEGORY_DEVELOPER) && !showDeveloperOptions ) {
	continue;	// skip developer options if not requested
      }
      if ( (ptr->category == UVAR_CATEGORY_DEPRECATED) && !showDeprecatedOptions ) {
	continue;	// skip deprecated options if not requested
      }
      if ( ptr->category == UVAR_CATEGORY_OBSOLETE ) {
	continue;	// always skip obsolete options to hide them completely from help string
      }

      UINT4 len;
      // max variable name length
      if ( ptr->name != NULL )
        {
          len = strlen ( ptr->name );
          nameFieldLen = (len > nameFieldLen) ? len : nameFieldLen;
        }

      // max type name length
      len = strlen ( UserVarTypeMap[ptr->type].name );
      typeFieldLen = (len > typeFieldLen) ? len : typeFieldLen;

    } // while ptr=ptr->next

  CHAR fmtStr[256];		// for building a dynamic format-string
  snprintf ( fmtStr, sizeof(fmtStr), "  %%s --%%-%ds   %%-%ds  %%s [%%s]\n", nameFieldLen, typeFieldLen );
  XLAL_LAST_ELEM(fmtStr)=0;

  CHAR defaultstr[256]; 	// for display of default-value
  CHAR strbuf[512];

  CHAR *helpstr_regular    = NULL;
  CHAR *helpstr_developer  = NULL;
  CHAR *helpstr_deprecated = NULL;
  // ---------- provide header line: info about config-file reading

  snprintf (strbuf, sizeof(strbuf), "Usage: %s [@ConfigFile] [options], where options are:\n\n", progname);
  XLAL_LAST_ELEM(strbuf) = 0;
  XLAL_CHECK_NULL ( (helpstr_regular = XLALStringDuplicate ( strbuf )) != NULL, XLAL_EFUNC );

  // ---------- MAIN LOOP: step through all user variables and add entry to appropriate help string
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL )	// header always empty
    {
      if ( ptr->category == UVAR_CATEGORY_REQUIRED ) {
	strcpy (defaultstr, "REQUIRED");
      }
      else if ( ptr->category == UVAR_CATEGORY_HELP ) {
        strcpy ( defaultstr, "");
      }
      else // write the current default-value into a string
	{
	  CHAR *valstr;
	  XLAL_CHECK_NULL ( (valstr = UserVarTypeMap [ ptr->type ].printer( ptr->varp )) != NULL, XLAL_EFUNC );
	  strncpy ( defaultstr, valstr, sizeof(defaultstr) );	// cut short for default-entry
	  XLAL_LAST_ELEM(defaultstr) = 0;
	  XLALFree (valstr);
	}

      CHAR optstr[10];
      if (ptr->optchar != 0) {
	sprintf (optstr, "-%c,", ptr->optchar);
      } else {
	strcpy (optstr, "   ");
      }

      snprintf ( strbuf, sizeof(strbuf),  fmtStr,
                 optstr,
                 ptr->name ? ptr->name : "-NONE-",
                 UserVarTypeMap[ptr->type].name,
                 ptr->help ? ptr->help : "-NONE-",
                 defaultstr
                 );
      XLAL_LAST_ELEM(strbuf) = 0;

      // now append new line to the appropriate helpstring
      switch ( ptr->category )
        {
        case UVAR_CATEGORY_DEVELOPER:
          if ( showDeveloperOptions ) {
            helpstr_developer = XLALStringAppend ( helpstr_developer, strbuf );
          }
          break;

        case UVAR_CATEGORY_DEPRECATED:
          if ( showDeprecatedOptions ) {
            helpstr_deprecated = XLALStringAppend ( helpstr_deprecated, strbuf );
          }
          break;

        default:
          helpstr_regular = XLALStringAppend ( helpstr_regular, strbuf );
          break;
        } // switch category

    } // while ptr=ptr->next

  CHAR *helpstr = NULL;
  XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, helpstr_regular )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, "\n" )) != NULL, XLAL_EFUNC );

  // handle output of developer options, if requested
  if ( haveDeveloperOptions )
    {
      if ( !showDeveloperOptions )
        {
          const char *str = " ---------- Use help with lalDebugLevel >= warning to also see all 'developer' options ----------\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
        }
      else
        {
          const char *str = " ---------- The following are 'developer'-options not useful for most users:----------\n\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, helpstr_developer )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, "\n" )) != NULL, XLAL_EFUNC );
        }
    } // if haveDeveloperOptions

  // handle output of deprecated options, if requested
  if ( haveDeprecatedOptions )
    {
      if ( !showDeprecatedOptions )
        {
          const char *str = " ---------- Use help with lalDebugLevel >= info to also see deprecated options ----------\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
        }
      else
        {
          const char *str = " ---------- The following are *DEPRECATED* options that shouldn't be used any more:----------\n\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, helpstr_deprecated )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, "\n" )) != NULL, XLAL_EFUNC );
        }
    } // if haveDeprecatedOptions

  XLALFree ( helpstr_regular );
  XLALFree ( helpstr_developer );
  XLALFree ( helpstr_deprecated );

  return helpstr;

} // XLALUserVarHelpString()


/**
 * Put all the pieces together, and basically does everything:
 * get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 */
int
XLALUserVarReadAllInput ( int argc, char *argv[] )
{
  XLAL_CHECK ( argc > 0, XLAL_EINVAL );
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );
  XLAL_CHECK ( argv[0] != NULL, XLAL_EINVAL );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  program_name = argv[0];	// keep a modul-local pointer to the executable name

  // ---------- pre-process command-line: have we got a config-file ?
  CHAR* cfgfile_name = NULL;
  for ( INT4 i = 1; i < argc; i++ )
    {
      char *argi = argv[i];
      XLAL_CHECK ( argi != NULL, XLAL_EINVAL, "argc = %d, but argv[%d] == NULL!\n", argc, i );

      if ( argi[0] == '@' )
	{
	  XLAL_CHECK ( cfgfile_name == NULL, XLAL_EINVAL, "Can only handle *one* config-file passed on commandline!\n" );
	  argi ++;
          XLAL_CHECK ( (cfgfile_name = XLALStringDuplicate ( argi )) != NULL, XLAL_EFUNC );
	} // if argument starts with '@' -> config-file

    } // for i < argc

  // ---------- if config-file specified, read from that first
  if ( cfgfile_name != NULL )
    {
      XLAL_CHECK ( XLALUserVarReadCfgfile ( cfgfile_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALFree (cfgfile_name);
    }

  // ---------- now parse cmdline: overloads previous config-file settings
  XLAL_CHECK ( XLALUserVarReadCmdline ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- handle special options that need some action ----------
  BOOLEAN skipCheckRequired = FALSE;
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL )
    {
      if ( (ptr->category == UVAR_CATEGORY_HELP) && ( *((BOOLEAN*)ptr->varp) ) )
	{
	  CHAR *helpstring;
	  XLAL_CHECK ( ( helpstring = XLALUserVarHelpString(argv[0])) != NULL, XLAL_EFUNC );
          printf ("\n%s\n", helpstring);
	  XLALFree (helpstring);
	  return XLAL_SUCCESS;
	} // if help requested

      // check 'special' category, which suppresses the CheckRequired test
      if ( (ptr->category == UVAR_CATEGORY_SPECIAL) && ptr->was_set ) {
	skipCheckRequired = TRUE;
      }

      // handle DEPRECATED options by outputting a warning:
      if ( ptr->category == UVAR_CATEGORY_DEPRECATED && ptr->was_set ) {
        XLALPrintError ("Option '%s' is DEPRECATED: %s\n", ptr->name, ptr->help );	// we output warning on error-level to make this very noticeable!
      }

      // handle DEPREC_ERROR options by throwing an error:
      if ( ptr->category == UVAR_CATEGORY_OBSOLETE && ptr->was_set ) {
        XLAL_ERROR ( XLAL_EINVAL, "Option '%s' is OBSOLETE: %s\n", ptr->name, ptr->help );
      }

    } // while ptr = ptr->next

  // check that all required input-variables have been specified
  if ( !skipCheckRequired ) {
    XLAL_CHECK ( XLALUserVarCheckRequired() == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALUserVarReadAllInput()


/**
 * Has this user-variable been set by the user?
 * returns 1 (=TRUE) or 0 (=FALSE) on success, error-code otherwise
 */
int
XLALUserVarWasSet ( const void *cvar )
{
  XLAL_CHECK ( cvar != NULL, XLAL_EINVAL );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  // find this variable name in the list of registered user-variables
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL )
    {
      if ( ptr->varp == cvar) {
        break;
      }
    } // while ptr = ptr->next

  XLAL_CHECK ( ptr != NULL, XLAL_EINVAL, "Variable pointer passed UVARwasSet is not a registered User-variable\n" );

  // we found it: has it been set by user?
  if ( ptr->was_set ) {
    return 1;
  } else {
    return 0;
  }

} // XLALUserVarWasSet()


/**
 * Check that all required user-variables have been set successfully.
 * Print error if not
 */
int
XLALUserVarCheckRequired (void)
{
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  // go through list of uvars
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL )
    {
      XLAL_CHECK ( ! ( ptr->category == UVAR_CATEGORY_REQUIRED && !ptr->was_set), XLAL_EFAILED, "Required user-variable '%s' has not been specified!\n\n", ptr->name );
    }

  return XLAL_SUCCESS;

} // XLALUserVarCheckRequired()


/**
 * Return a log-string representing the <em>complete</em> user-input.
 * <em>NOTE:</em> we only record user-variables that have been set
 * by the user.
 */
CHAR *
XLALUserVarGetLog ( UserVarLogFormat format 	/**< output format: return as config-file or command-line */
                    )
{
  XLAL_CHECK_NULL ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );
  XLAL_CHECK_NULL ( format < UVAR_LOGFMT_LAST, XLAL_EINVAL );

  CHAR *record = NULL;

  if ( format == UVAR_LOGFMT_CMDLINE ) {
    XLAL_CHECK_NULL ( (record = XLALStringAppend ( record, program_name)) != NULL, XLAL_EFUNC );
  }

  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr = ptr->next) )
    {
      if ( ! ptr->was_set ) { // skip unset variables
	continue;
      }

      CHAR *valstr;
      XLAL_CHECK_NULL ( (valstr = UserVarTypeMap [ ptr->type ].printer( ptr->varp )) != NULL, XLAL_EFUNC );

      char append[256];
      switch (format)
	{
	case UVAR_LOGFMT_CFGFILE:
	  snprintf (append, sizeof(append), "%s = %s;\n", ptr->name, valstr);
	  break;

	case UVAR_LOGFMT_CMDLINE:
	  snprintf (append, sizeof(append), " --%s=%s", ptr->name, valstr);
	  break;

	case UVAR_LOGFMT_PROCPARAMS:
	  snprintf (append, sizeof(append), "--%s = %s :%s;", ptr->name, valstr, UserVarTypeMap[ptr->type].name );
	  break;

	default:
          XLAL_ERROR_NULL ( XLAL_EINVAL, "Unknown format for recording user-input: '%i'\n", format );
	  break;
	} // switch (format)
      XLAL_LAST_ELEM(append) = 0;

      XLAL_CHECK_NULL ( (record = XLALStringAppend (record, append)) != NULL, XLAL_EFUNC );

      XLALFree (valstr);
    } // while ptr=ptr->next

  return record;

} // XLALUserVarGetLog()

/**
 * Mark the user-variable as set, check if it has been
 * set previously and issue a warning if set more than once ...
 */
void
check_and_mark_as_set ( LALUserVariable *varp )
{
  // output warning if this variable has been set before ...
  if ( varp->was_set ) {
    XLALPrintWarning ( "User-variable '%s' was set more than once!\n", varp->name ? varp->name : "(NULL)" );
  }

  varp->was_set = 1;

  return;
} // check_and_mark_as_set()

/* ========== DEPRECATED LAL INTERFACE FUNCTIONS, which have been replaced by XLAL functions,
 * These functions are just wrappers around the XLAL functions
 */
#define USERINPUTH_EXLAL        1
#define USERINPUTH_MSGEXLAL	"Failure in XLAL function"

/** \deprecated us XLALDestroyUserVars() instead */
void
LALDestroyUserVars (LALStatus *status)
{
  INITSTATUS(status);
  XLALDestroyUserVars();
  RETURN (status);
} // LALDestroyUserVars()


/** \deprecated use XLALUserVarReadCmdline() instead */
void
LALUserVarReadCmdline (LALStatus *status, int argc, char *argv[])
{
  INITSTATUS(status);
  if ( XLALUserVarReadCmdline(argc, argv) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALUserVarReadCmdline() failed with code %d\n", xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} // LALUserVarReadCmdline()

/** \deprecated use XLALUserVarCheckRequired() instead */
void
LALUserVarCheckRequired (LALStatus *status)
{
  INITSTATUS(status);
  if ( XLALUserVarCheckRequired() != XLAL_SUCCESS ) {
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} /* LALUserVarCheckRequired() */


/** \deprecated use XLALUserVarReadAllInput() instead */
void
LALUserVarReadAllInput (LALStatus *status, int argc, char *argv[])
{
  INITSTATUS(status);
  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    XLALPrintError ( "XLALUserVarReadAllInput() failed with code %d\n", xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} // LALReadUserInput()

/** \deprecated use XLALUserVarWasSet() instead */
INT4
LALUserVarWasSet (const void *cvar)
{
  return (XLALUserVarWasSet(cvar));
}

/** \deprecated use XLALUserVarGetLog() instead */
void
LALUserVarGetLog (LALStatus *status, CHAR **logstr,  UserVarLogFormat format)
{
  INITSTATUS(status);
  if ( ((*logstr) = XLALUserVarGetLog ( format )) == NULL ) {
    XLALPrintError ("UserVarLogFormat() failed.\n" );
    ABORT (status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL);
  }
  RETURN (status);
} /* LALUserVarGetLog() */

/**
 * \deprecated use XLALUserVarHelpString() instead
 */
void
LALUserVarHelpString (LALStatus *status,
		      CHAR **helpstring, /* output: allocated here! */
		      const CHAR *progname)
{
  INITSTATUS(status);
  if ( ((*helpstring) = XLALUserVarHelpString ( progname )) == NULL ) {
    XLALPrintError ("XLALUserVarHelpString() failed with code %d\n", xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
} /* LALUserVarHelpString() */

/** \deprecated use XLALRegisterREALUserVar() instead */
void
LALRegisterREALUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarCategory category,
			const CHAR *helpstr,
			REAL8 *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_REAL8, optchar, category, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALRegisterINTUserVar() instead */
void
LALRegisterINTUserVar (LALStatus *status,
		       const CHAR *name,
		       CHAR optchar,
		       UserVarCategory category,
		       const CHAR *helpstr,
		       INT4 *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_INT4, optchar, category, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALRegisterBOOLUserVar() instead */
void
LALRegisterBOOLUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarCategory category,
			const CHAR *helpstr,
			BOOLEAN *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_BOOLEAN, optchar, category, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALRegisterSTRINGUserVar() instead */
void
LALRegisterSTRINGUserVar (LALStatus *status,
			  const CHAR *name,
			  CHAR optchar,
			  UserVarCategory category,
			  const CHAR *helpstr,
			  CHAR **cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_STRING, optchar, category, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALRegisterSTRINGVectorUserVar() instead */
void
LALRegisterLISTUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarCategory category,
			const CHAR *helpstr,
			LALStringVector **cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_STRINGVector, optchar, category, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALUserVarReadCfgfile() instead */
void
LALUserVarReadCfgfile (LALStatus *status,
		       const CHAR *cfgfile) 	   /* name of config-file */
{

  INITSTATUS(status);
  if ( XLALUserVarReadCfgfile ( cfgfile ) != XLAL_SUCCESS ) {
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} // LALUserVarReadCfgfile()
