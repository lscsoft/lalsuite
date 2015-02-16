/*
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

#include "LALgetopt.h"

#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALString.h>
#include <lal/Date.h>
#include <lal/StringVector.h>

#define TRUE  (1==1)
#define FALSE (1==0)

/* Defines the type of a "user-variable": bool, int, real or string.
 * Should be used only internally !!
 */
typedef enum {
  UVAR_TYPE_BOOL,    /* boolean */
  UVAR_TYPE_INT4,    /* integer */
  UVAR_TYPE_REAL8,   /* float */
  UVAR_TYPE_STRING,  /* string */
  UVAR_TYPE_CSVLIST, /* list of comma-separated values */
  UVAR_TYPE_EPOCH,   /* time 'epoch', specified in either GPS or MJD(TT) format, translated into GPS */
  UVAR_TYPE_LONGITUDE,/* sky longitude (aka right-ascencion or RA), in either radians or hours:minutes:seconds format, translated into radians */
  UVAR_TYPE_LATITUDE, /* sky latitude (aka declination or DEC), in either radians or degrees:minutes:seconds format, translated into radians */
  UVAR_TYPE_LAST
} UserVarType;


/**
 * Array of descriptors for each UserVarType type
 */
static const struct {
  const char *const name;	/**< type name */
} UserVarTypeDescription[UVAR_TYPE_LAST] = {

  [UVAR_TYPE_BOOL]     = {"BOOLEAN"},
  [UVAR_TYPE_INT4]     = {"INT4"},
  [UVAR_TYPE_REAL8]    = {"REAL8"},
  [UVAR_TYPE_STRING]   = {"STRING"},
  [UVAR_TYPE_CSVLIST]  = {"CSVLIST"},
  [UVAR_TYPE_EPOCH]    = {"EPOCH"},
  [UVAR_TYPE_LONGITUDE]= {"LONGITUDE"},
  [UVAR_TYPE_LATITUDE] = {"LATITUDE"}
};

/**
 * Linked list to hold the complete information about the user-variables.
 */
typedef struct tagLALUserVariable {
  const CHAR *name;	/**< full name */
  UserVarType type;	/**< type: bool, int, float or string */
  CHAR optchar;		/**< cmd-line character */
  const CHAR *help;	/**< help-string */
  void *varp;		/**< pointer to the actual C-variable */
  UserVarFlag state;	/**< state (empty, default, set) */
  struct tagLALUserVariable *next; /* linked list */
} LALUserVariable;

/** The module-local linked list to hold the user-variables */
static LALUserVariable UVAR_vars;	/**< empty head */
static const CHAR *program_name;	/**< keep a pointer to the program name */

/* needed for command-line parsing */
extern char *LALoptarg;
extern int LALoptind, LALopterr, LALoptopt;

/* ---------- internal prototypes ---------- */

/* ----- XLAL interface ----- */
int XLALRegisterUserVar (const CHAR *name, UserVarType type, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, void *cvar);
CHAR *XLALUvarValue2String (LALUserVariable *uvar);
CHAR *XLAL_copy_string_unquoted ( const CHAR *in );
void check_and_mark_as_set ( LALUserVariable *varp );


// ==================== Function definitions ====================

/* these are type-specific wrappers to allow tighter type-checking! */
/** Register a user-variable of type REAL8, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterREALUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar )
{
  return XLALRegisterUserVar (name, UVAR_TYPE_REAL8, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type INT4, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterINTUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, INT4 *cvar )
{
  return XLALRegisterUserVar (name, UVAR_TYPE_INT4, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type BOOLEAN, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterBOOLUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, BOOLEAN *cvar )
{
  return XLALRegisterUserVar (name, UVAR_TYPE_BOOL, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type CHAR*, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterSTRINGUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, CHAR **cvar )
{
  return XLALRegisterUserVar (name, UVAR_TYPE_STRING, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of 'list' type LALStringVector, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterLISTUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, LALStringVector **cvar)
{
  return XLALRegisterUserVar ( name, UVAR_TYPE_CSVLIST, optchar, flag, helpstr, cvar );
}

/** Register a user-variable of 'EPOCH' type LIGOTimeGPS, allowing both GPS and MJD string inputs, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterEPOCHUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, LIGOTimeGPS *cvar)
{
  return XLALRegisterUserVar ( name, UVAR_TYPE_EPOCH, optchar, flag, helpstr, cvar );
}

/** Register a user-variable of 'LONGITUDE' type (REAL8), allowing both "hours:minutes:seconds" or radians as input */
int
XLALRegisterLONGITUDEUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar)
{
  return XLALRegisterUserVar ( name, UVAR_TYPE_LONGITUDE, optchar, flag, helpstr, cvar );
}

/** Register a user-variable of 'LATITUDE' type (REAL8), allowing both "degrees:minutes:seconds" or radians as input */
int
XLALRegisterLATITUDEUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar)
{
  return XLALRegisterUserVar ( name, UVAR_TYPE_LATITUDE, optchar, flag, helpstr, cvar );
}


/**
 * \ingroup UserInput_h
 * Internal function: Register a user-variable with the module.
 * Effectively put an appropriate entry into UVAR_vars
 *
 * Checks that long- and short-options are unique, an error is returned
 * if a previous option name collides.
 *
 * \note don't use this function directly, as it is not type-safe!!
 * ==> use one of the appropriate typed wrappers:  XLALRegisterREALUserVar(), XLALRegisterINTUserVar(), ...
 */
int
XLALRegisterUserVar ( const CHAR *name,		/**< name of user-variable to register */
                      UserVarType type,		/**< variable type (int,bool,string,real) */
                      CHAR optchar,		/**< optional short-option character */
                      UserVarFlag flag,	/**< sets state flag to this */
                      const CHAR *helpstr,	/**< help-string explaining this input-variable */
                      void *cvar		/**< pointer to the actual C-variabe to link to this user-variable */
                      )
{
  XLAL_CHECK ( (cvar != NULL) && (name != NULL), XLAL_EINVAL );

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
  ptr->state 	= flag;

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
      // is an allocated string here?
      if ( (ptr->type == UVAR_TYPE_STRING) && (*(CHAR**)(ptr->varp) != NULL) )
	{
	  XLALFree ( *(CHAR**)(ptr->varp) );
	  *(CHAR**)(ptr->varp) = NULL;
	}
      else if ( ptr->type == UVAR_TYPE_CSVLIST )
        {
          XLALDestroyStringVector ( *(LALStringVector**)ptr->varp );
          *(LALStringVector**)(ptr->varp) = NULL;
        }

      /* free list-entry behind us (except for the head) */
      if ( lastptr != NULL )
        {
          XLALFree ( lastptr );
          lastptr = NULL;
        }

      lastptr = ptr;

    } // while ptr->next

  if ( lastptr != NULL )
    {
      XLALFree ( lastptr );
      lastptr=NULL;
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
      if (ptr->type == UVAR_TYPE_BOOL) {	/* but for BOOL its optional */
	optstring[pos++] = ':';
      }
    } // while ptr->next
  optstring[pos] = '\0';

  // ---------- fill option-struct for long-options
  struct option *long_options = LALCalloc (1, (numvars+1) * sizeof(struct option));
  ptr = &UVAR_vars;	// start again from beginning: empty head
  pos = 0;
  while ( (ptr= ptr->next) != NULL)
    {
      if (ptr->name == NULL) {		/* if no long-option given, ignore */
	continue;
      }
      long_options[pos].name 	= ptr->name;
      long_options[pos].has_arg = (ptr->type == UVAR_TYPE_BOOL) ? optional_argument : required_argument;
      long_options[pos].flag 	= NULL;	// get val returned from getopt_long()
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
   * Bruce's notes: read getopt_long() source code, and in particular
   * _getopt_internal() to see what is initialized.
   */
  LALoptind = 0; 	// reset our local getopt(), getopt_long()

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

      switch (ptr->type)
	{
	case UVAR_TYPE_BOOL:
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
            XLAL_CHECK ( XLALParseStringValueToBOOLEAN ( (BOOLEAN*)(ptr->varp), LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
          }
	  break;

	case UVAR_TYPE_INT4:
          XLAL_CHECK ( XLALParseStringValueToINT4 ( (INT4*)(ptr->varp), LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
	  break;

	case UVAR_TYPE_REAL8:
          XLAL_CHECK ( XLALParseStringValueToREAL8 ( (REAL8*)(ptr->varp), LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
	  break;

	case UVAR_TYPE_STRING:
	  XLALFree ( *(CHAR**)(ptr->varp) ); // in case something allocated here before
	  XLAL_CHECK ( ( *(CHAR**)(ptr->varp) = XLAL_copy_string_unquoted ( LALoptarg )) != NULL, XLAL_EFUNC );
	  break;

	case UVAR_TYPE_CSVLIST:	// list of comma-separated string values
	  XLALDestroyStringVector ( *(LALStringVector**)(ptr->varp) );	// in case sth allocated here before
	  XLAL_CHECK ( (*(LALStringVector**)(ptr->varp) = XLALParseCSV2StringVector ( LALoptarg )) != NULL, XLAL_EFUNC );
	  break;

        case UVAR_TYPE_EPOCH:
          XLAL_CHECK ( XLALParseStringValueToEPOCH ( (LIGOTimeGPS *)(ptr->varp), LALoptarg ) != NULL, XLAL_EFUNC );
	  break;

        case UVAR_TYPE_LONGITUDE:
          XLAL_CHECK ( XLALParseStringValueToLONGITUDE ( (REAL8 *)(ptr->varp), LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
	  break;

        case UVAR_TYPE_LATITUDE:
          XLAL_CHECK ( XLALParseStringValueToLATITUDE ( (REAL8 *)(ptr->varp), LALoptarg ) == XLAL_SUCCESS, XLAL_EFUNC );
	  break;

	default:
	  XLALPrintError ( "%s: ERROR: unkown UserVariable-type encountered... points to a coding error!\n", __func__ );
	  XLAL_ERROR ( XLAL_EINVAL );
	  break;

	} // switch ptr->type

      check_and_mark_as_set ( ptr );

    } // while getopt_long()

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

  CHAR *stringbuf;
  LALParsedDataFile *cfg = NULL;
  XLAL_CHECK ( XLALParseDataFile ( &cfg, cfgfile ) == XLAL_SUCCESS, XLAL_EFUNC );

  // step through all user-variable: read those with names from config-file
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->name == NULL) {	// ignore name-less user-variable
	continue;
      }

      BOOLEAN wasRead = FALSE;

      switch (ptr->type)
	{
	case UVAR_TYPE_BOOL:
          XLAL_CHECK ( XLALReadConfigBOOLVariable ( ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

	case UVAR_TYPE_INT4:
	  XLAL_CHECK ( XLALReadConfigINT4Variable ( ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

	case UVAR_TYPE_REAL8:
	  XLAL_CHECK ( XLALReadConfigREAL8Variable(ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

	case UVAR_TYPE_STRING:
	  stringbuf = NULL;
	  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &stringbuf, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if ( wasRead )
	    {
	      XLALFree ( *(CHAR**)(ptr->varp) ); // if anything allocated here before
	      *(CHAR**)(ptr->varp) = stringbuf;
	      check_and_mark_as_set ( ptr );
	    }
	  break;

	case UVAR_TYPE_CSVLIST:
	  stringbuf = NULL;
	  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &stringbuf, cfg, NULL, ptr->name,&wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if ( wasRead )
	    {
	      XLALDestroyStringVector ( *(LALStringVector**)(ptr->varp) ); // if anything allocated here before
	      XLAL_CHECK ( (*(LALStringVector**)(ptr->varp) = XLALParseCSV2StringVector ( stringbuf )) != NULL, XLAL_EFUNC );
	      check_and_mark_as_set ( ptr );
	    }
	  break;

        case UVAR_TYPE_EPOCH:
	  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

        case UVAR_TYPE_LONGITUDE:
	  XLAL_CHECK ( XLALReadConfigLONGITUDEVariable ( ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

        case UVAR_TYPE_LATITUDE:
	  XLAL_CHECK ( XLALReadConfigLATITUDEVariable ( ptr->varp, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
	  if (wasRead) { check_and_mark_as_set ( ptr ); }
	  break;

	default:
          XLAL_ERROR ( XLAL_EFAILED, "ERROR: unkown UserVariable-type encountered...points to a coding error!\n" );
          break;

	} // switch ptr->type

    } // while ptr->next

  // ok, that should be it: check if there were more definitions we did not read
  XLAL_CHECK ( XLALCheckConfigReadComplete ( cfg, CONFIGFILE_WARN ) == XLAL_SUCCESS, XLAL_EFUNC );

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

  CHAR *helpstr_regular = NULL;
  CHAR *helpstr_developer = NULL;

  BOOLEAN showDeveloperOptions = (lalDebugLevel > 0);	// currently only output for lalDebugLevel > 0

  // ---------- ZEROTH PASS: find longest long-option and type names, for proper output formatting
  LALUserVariable *ptr = &UVAR_vars;
  UINT4 nameFieldLen = 0;
  UINT4 typeFieldLen = 0;
  BOOLEAN haveDevOpts = 0;
  while ( (ptr=ptr->next) != NULL )
    {
      if (ptr->state & UVAR_DEVELOPER) {
        haveDevOpts = 1;
      }
      if ( !showDeveloperOptions && (ptr->state & UVAR_DEVELOPER) ) {
	continue;	// skip developer-options if not requested
      }

      UINT4 len;
      // max variable name length
      if ( ptr->name != NULL )
        {
          len = strlen ( ptr->name );
          nameFieldLen = (len > nameFieldLen) ? len : nameFieldLen;
        }

      // max type name length
      len = strlen ( UserVarTypeDescription[ptr->type].name );
      typeFieldLen = (len > typeFieldLen) ? len : typeFieldLen;

    } // while ptr=ptr->next

  CHAR fmtStr[256];		// for building a dynamic format-string
  snprintf ( fmtStr, sizeof(fmtStr), "  %%s --%%-%ds   %%-%ds  %%s [%%s]\n", nameFieldLen, typeFieldLen );
  fmtStr[sizeof(fmtStr)-1]=0;

  CHAR defaultstr[256]; 	// for display of default-value
  CHAR strbuf[512];

  // ---------- provide header line: info about config-file reading
  snprintf (strbuf, sizeof(strbuf), "Usage: %s [@ConfigFile] [options], where options are:\n\n", progname);
  strbuf[sizeof(strbuf)-1] = 0;
  XLAL_CHECK_NULL ( (helpstr_regular = XLALStringDuplicate ( strbuf )) != NULL, XLAL_EFUNC );

  // ---------- MAIN LOOP: step through all user variables and add entry to appropriate help string
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL )
    {
      if ( ptr->state & UVAR_REQUIRED ) {
	strcpy (defaultstr, "REQUIRED");
      }
      else if ( ptr->state & UVAR_HELP ) {
        strcpy ( defaultstr, "");
      }
      else // write the current default-value into a string
	{
	  CHAR *valstr = NULL;
	  XLAL_CHECK_NULL ( (valstr = XLALUvarValue2String ( ptr )) != NULL, XLAL_EFUNC );
	  strncpy ( defaultstr, valstr, sizeof(defaultstr) );	// cut short for default-entry
	  defaultstr[sizeof(defaultstr)-1] = 0;
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
                 UserVarTypeDescription[ptr->type].name,
                 ptr->help ? ptr->help : "-NONE-",
                 defaultstr
                 );
      strbuf[sizeof(strbuf)-1] = 0;

      // now append new line to the appropriate helpstring
      if ( ptr->state & UVAR_DEVELOPER ) {
        if ( showDeveloperOptions ) {
          helpstr_developer = XLALStringAppend ( helpstr_developer, strbuf );
        }
      } else {
        helpstr_regular = XLALStringAppend ( helpstr_regular, strbuf );
      }

    } // while ptr=ptr->next

  CHAR *helpstr = NULL;
  XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, helpstr_regular )) != NULL, XLAL_EFUNC );

  // handle output of developer options, if requested
  if ( haveDevOpts )
    {
      if ( !showDeveloperOptions )
        {
          const char *str = "\n ---------- Use help with lalDebugLevel > 0 to also see all 'developer-options' ---------- \n\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
        }
      else
        {
          const char *str = "\n   ---------- The following are 'Developer'-options not useful for most users:----------\n\n";
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, str )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (helpstr = XLALStringAppend ( helpstr, helpstr_developer )) != NULL, XLAL_EFUNC );
        }
    } // if haveDevOpts

  XLALFree ( helpstr_regular );
  XLALFree ( helpstr_developer );

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

  BOOLEAN skipCheckRequired = FALSE;
  // ---------- check if help-string was requested
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL )
    {
      if ( (ptr->state & UVAR_HELP) && ( *((BOOLEAN*)ptr->varp) ) )
	{
	  CHAR *helpstring;
	  XLAL_CHECK ( ( helpstring = XLALUserVarHelpString(argv[0])) != NULL, XLAL_EFUNC );
          printf ("\n%s\n", helpstring);
	  XLALFree (helpstring);
	  return XLAL_SUCCESS;
	} // if help requested

      // check 'special' flag, which suppresses the CheckRequired test
      if ( (ptr->state & UVAR_SPECIAL) && (ptr->state & UVAR_WAS_SET) ) {
	skipCheckRequired = TRUE;
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
  if ( (ptr->state & UVAR_WAS_SET) != 0 ) {
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
      XLAL_CHECK ( ((ptr->state & UVAR_REQUIRED) == 0) || ((ptr->state & UVAR_WAS_SET) != 0), XLAL_EFAILED, "Required user-variable '%s' has not been specified!\n\n", ptr->name );
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
      if ( (ptr->state & UVAR_WAS_SET) == FALSE ) {	// skip unset variables
	continue;
      }

      CHAR *valstr;
      XLAL_CHECK_NULL ( (valstr = XLALUvarValue2String ( ptr )) != NULL, XLAL_EFUNC );

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
	  snprintf (append, sizeof(append), "--%s = %s :%s;", ptr->name, valstr, UserVarTypeDescription[ptr->type].name );
	  break;

	default:
          XLAL_ERROR_NULL ( XLAL_EINVAL, "Unknown format for recording user-input: '%i'\n", format );
	  break;
	} // switch (format)
      append[sizeof(append)-1] = 0;

      XLAL_CHECK_NULL ( (record = XLALStringAppend (record, append)) != NULL, XLAL_EFUNC );

      XLALFree (valstr);
    } // while ptr=ptr->next

  return record;

} // XLALUserVarGetLog()

/* Return the value of the given UserVariable as a string.
 * For INTERNAL use only!
 */
CHAR *
XLALUvarValue2String ( LALUserVariable *uvar )
{
  XLAL_CHECK_NULL ( uvar != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( uvar->varp != NULL, XLAL_EINVAL );

  char buf[512];
  char *retstr = NULL;

  switch ( uvar->type )
    {
    case UVAR_TYPE_BOOL:
      sprintf (buf, *(BOOLEAN*)(uvar->varp) ? "TRUE" : "FALSE");
      break;

    case UVAR_TYPE_INT4:
      sprintf (buf, "%" LAL_INT4_FORMAT, *(INT4*)(uvar->varp) );
      break;

    case UVAR_TYPE_REAL8:
    case UVAR_TYPE_LONGITUDE:
    case UVAR_TYPE_LATITUDE:
      if (*(REAL8*)(uvar->varp) == 0) {
	strcpy (buf, "0.0");	// makes it more explicit that's it a REAL
      } else {
	sprintf (buf, "%.16g", *(REAL8*)(uvar->varp) );
      }
      break;

    case UVAR_TYPE_STRING:
      if ( *(CHAR**)(uvar->varp) != NULL ) {
        snprintf (buf, sizeof(buf), "\"%s\"", *(CHAR**)(uvar->varp) );
        buf[sizeof(buf)-1] = 0;
      } else {
        strcpy (buf, "NULL");
      }
      break;

    case UVAR_TYPE_EPOCH:
      XLAL_CHECK_NULL ( XLALGPSToStr ( buf, (const LIGOTimeGPS *)(uvar->varp) ) != NULL, XLAL_EFUNC );
      strcat ( buf, "GPS" );	// postfix this with 'units' for explicitness (as opposed to 'MJD')
      break;

    case UVAR_TYPE_CSVLIST:
      if ( *(LALStringVector**)(uvar->varp) != NULL ) {
        XLAL_CHECK_NULL ( (retstr = XLALStringVector2CSV ( *(LALStringVector**)(uvar->varp) )) != NULL, XLAL_EFUNC );
      } else {
        strcpy (buf, "NULL");
      }
      break;

    default:
      XLAL_ERROR_NULL ( XLAL_EINVAL, "\nUnkown UserVariable-type encountered... this points to a coding error!\n" );
      break;

    } // switch uvar->type

  if ( retstr == NULL ) {
    XLAL_CHECK_NULL ( (retstr = XLALStringDuplicate ( buf )) != NULL, XLAL_EFUNC );
  }

  return retstr;

} // XLALUvarValue2String()


/**
 * Copy (and allocate) string 'in', possibly with quotes \" or \' removed.
 * If quotes are present at the beginning of 'in', they must have a matching
 * quote at the end of string, otherwise an error is printed and return=NULL
 */
CHAR *
XLAL_copy_string_unquoted ( const CHAR *in )
{
  XLAL_CHECK_NULL ( in != NULL, XLAL_EINVAL );


  CHAR opening_quote = 0;
  CHAR closing_quote = 0;
  UINT4 inlen = strlen ( in );

  if ( (in[0] == '\'') || (in[0] == '\"') ) {
    opening_quote = in[0];
  }
  if ( (inlen >= 2) && ( (in[inlen-1] == '\'') || (in[inlen-1] == '\"') ) ) {
    closing_quote = in[inlen-1];
  }

  // check matching quotes
  XLAL_CHECK_NULL ( opening_quote == closing_quote, XLAL_EINVAL, "Unmatched quotes in string [%s]\n", in );

  const CHAR *start = in;
  UINT4 outlen = inlen;
  if ( opening_quote )
    {
      start = in + 1;
      outlen = inlen - 2;
    }

  CHAR *ret;
  XLAL_CHECK_NULL ( (ret = LALCalloc (1, outlen + 1)) != NULL, XLAL_ENOMEM );
  strncpy ( ret, start, outlen );
  ret[outlen] = 0;

  return ret;

} // XLAL_copy_string_unquoted()

/**
 * Mark the user-variable as set, check if it has been
 * set previously and issue a warning if set more than once ...
 */
void
check_and_mark_as_set ( LALUserVariable *varp )
{
  // output warning if this variable has been set before ...
  if ( (varp->state & UVAR_WAS_SET) ) {
    XLALPrintWarning ( "User-variable '%s' was set more than once!\n", varp->name ? varp->name : "(NULL)" );
  }

  varp->state = (UserVarFlag)( varp->state |  UVAR_WAS_SET );

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
			UserVarFlag flag,
			const CHAR *helpstr,
			REAL8 *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_REAL8, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
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
		       UserVarFlag flag,
		       const CHAR *helpstr,
		       INT4 *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_INT4, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
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
			UserVarFlag flag,
			const CHAR *helpstr,
			BOOLEAN *cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_BOOL, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
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
			  UserVarFlag flag,
			  const CHAR *helpstr,
			  CHAR **cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_STRING, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALRegisterLISTUserVar() instead */
void
LALRegisterLISTUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarFlag flag,
			const CHAR *helpstr,
			LALStringVector **cvar)
{
  INITSTATUS(status);
  if ( XLALRegisterUserVar ( name, UVAR_TYPE_CSVLIST, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
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
