/*
 * Copyright (C) 2016 Karl Wette
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
#include <config.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_IOCTL_H
#include <sys/ioctl.h>
#endif

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
int XLALRegister ##UTYPE## UserVar ( CTYPE *cvar, const CHAR *name, CHAR optchar, UserVarCategory category, const CHAR *fmt, ... ) \
{                                                           \
  char helpstr[2048];                                       \
  va_list ap;                                               \
  va_start(ap, fmt);                                        \
  vsnprintf(helpstr, sizeof(helpstr), fmt, ap);             \
  va_end(ap);                                               \
  return XLALRegisterUserVar (cvar, name, UVAR_TYPE_ ## UTYPE, optchar, category, helpstr); \
}

// ---------- local type definitions ----------

// Define the type of a "user-variable": bool, int, real or string, ...
typedef enum {
  UVAR_TYPE_START=0, 		// internal start marker for range checking

  UVAR_TYPE_BOOLEAN, 		// boolean
  UVAR_TYPE_INT4,    		// 4-byte integer
  UVAR_TYPE_INT8,    		// 8-byte integer
  UVAR_TYPE_REAL8,   		// 8-byte float
  UVAR_TYPE_EPOCH,   		// time 'epoch', specified in either GPS or MJD(TT) format, translated into GPS
  UVAR_TYPE_RAJ,     		// sky equatorial longitude (aka right-ascencion or RA), in either radians or hours:minutes:seconds format, translated into radians
  UVAR_TYPE_DECJ,    		// sky equatorial latitude (aka declination or DEC), in either radians or degrees:minutes:seconds format, translated into radians
  UVAR_TYPE_STRING, 		// normal string

  UVAR_TYPE_REAL8Range,		// range of REAL8 values
  UVAR_TYPE_EPOCHRange,		// range of LIGOTimeGPS values
  UVAR_TYPE_RAJRange,		// range of RAJ values
  UVAR_TYPE_DECJRange,		// range of DECJ values

  UVAR_TYPE_STRINGVector,	// list of comma-separated strings
  UVAR_TYPE_REAL8Vector,	// list of comma-separated REAL8's
  UVAR_TYPE_INT4Vector,		// list of comma-separated INT4's

  UVAR_TYPE_END      	// internal end marker for range checking
} UserVarType;

//
// Linked list to hold the complete information about the user-variables.
//
typedef struct tagLALUserVariable {
  CHAR name[64];			// full name
  UserVarType type;			// variable type: BOOLEAN, INT4, REAL8, ...
  CHAR optchar;				// cmd-line character
  CHAR help[2048];			// help-string
  void *varp;				// pointer to the actual C-variable
  UserVarCategory category;		// category (optional, required, developer, ... )
  BOOLEAN was_set;			// was this set by the user: 0=no, 1=via cfg-file, 2=via cmdline
  struct tagLALUserVariable *next; 	// linked list
} LALUserVariable;

// ---------- local prototypes ----------
int XLALRegisterUserVar (void *cvar, const CHAR *name, UserVarType type, CHAR optchar, UserVarCategory category, const CHAR *help);
int XLALUserVarPrintUsage ( FILE *file );
int XLALUserVarPrintHelp ( FILE *file );
static void format_user_var_names( char *s );
static void fprint_wrapped( FILE *f, int line_width, const char *prefix, char *text );

// ----- define templated registration functions for all supported UVAR_TYPE_ 'UTYPES'
DEFN_REGISTER_UVAR(BOOLEAN,BOOLEAN);
DEFN_REGISTER_UVAR(INT4,INT4);
DEFN_REGISTER_UVAR(INT8,INT8);
DEFN_REGISTER_UVAR(REAL8,REAL8);
DEFN_REGISTER_UVAR(RAJ,REAL8);
DEFN_REGISTER_UVAR(DECJ,REAL8);
DEFN_REGISTER_UVAR(EPOCH,LIGOTimeGPS);
DEFN_REGISTER_UVAR(STRING,CHAR*);

DEFN_REGISTER_UVAR(REAL8Range,REAL8Range);
DEFN_REGISTER_UVAR(EPOCHRange,LIGOTimeGPSRange);
DEFN_REGISTER_UVAR(RAJRange,REAL8Range);
DEFN_REGISTER_UVAR(DECJRange,REAL8Range);

DEFN_REGISTER_UVAR(STRINGVector,LALStringVector*);
DEFN_REGISTER_UVAR(REAL8Vector,REAL8Vector*);
DEFN_REGISTER_UVAR(INT4Vector,INT4Vector*);

// ----- define helper types for casting
typedef int (*parserT)(void*, const char*);
typedef void (*destructorT)(void*);
typedef char *(*printerT)(const void*);

// ----- handy macro to simplify adding 'regular' entries for new UTYPES into UserVarTypeMap
#define REGULAR_MAP_ENTRY(UTYPE,DESTRUCTOR,FORMATHELP) \
  [UVAR_TYPE_##UTYPE] = { \
    .name = #UTYPE, \
    .format_help = FORMATHELP, \
    .parser = (parserT)XLALParseStringValueAs##UTYPE, \
    .printer = (printerT)XLALPrintStringValueOf##UTYPE, \
    .destructor = (destructorT)DESTRUCTOR \
  }

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
  const char *const format_help;		///< help string describing format of user variable
  int (*parser)(void*, const char*);		///< parser function to parse string as this type
  char *(*printer)(const void *);		///< 'printer' function returning string value for given type
  void (*destructor)(void*);			///< destructor for this variable type, NULL if none required
} UserVarTypeMap[UVAR_TYPE_END]
= {
  // either use 'manual' entries of the form
  // [UVAR_TYPE_\<UTYPE\>] = { "\<UTYPE\>",	(parserT)XLALParseStringValueAs\<UTYPE\>, (printerT)XLALPrintStringValueOf\<UTYPE\>, (destructorT)XLALDestroy\<UTYPE\> },
  // or the convenience macro for cases using 'standard' function names and API
  // REGULAR_MAP_ENTRY ( \<UTYPE\>, XLALDestroy\<UTYPE\> ),
  REGULAR_MAP_ENTRY ( BOOLEAN, NULL, "[TRUE|FALSE | YES|NO | 1|0]" ),
  REGULAR_MAP_ENTRY ( INT4, NULL, "<4-byte integer>" ),
  REGULAR_MAP_ENTRY ( INT8, NULL, "<8-byte integer>" ),
  REGULAR_MAP_ENTRY ( REAL8, NULL, "<8-byte real>" ),
  REGULAR_MAP_ENTRY ( STRING, XLALFree, "<string>" ),
  REGULAR_MAP_ENTRY ( EPOCH, NULL, "<seconds>[.<nano-seconds>][GPS|MJD]" ),
  REGULAR_MAP_ENTRY ( RAJ, NULL, "<radians>|<hours>:<minutes>:<seconds>" ),
  REGULAR_MAP_ENTRY ( DECJ, NULL, "<radians>|<degrees>:<minutes>:<seconds>" ),

  REGULAR_MAP_ENTRY ( REAL8Range, NULL, "<start>[,<end>|/<band>|~<plus-minus>] where <>=<8-byte real>" ),
  REGULAR_MAP_ENTRY ( EPOCHRange, NULL, "<start>[,<end>|/<band>|~<plus-minus>] where <>=<seconds>[.<nano-seconds>][GPS|MJD]" ),
  REGULAR_MAP_ENTRY ( RAJRange, NULL, "<start>[,<end>|/<band>|~<plus-minus>] where <>=<radians>|<hours>:<minutes>:<seconds>" ),
  REGULAR_MAP_ENTRY ( DECJRange, NULL, "<start>[,<end>|/<band>|~<plus-minus>] where <>=<radians>|<degrees>:<minutes>:<seconds>" ),

  REGULAR_MAP_ENTRY ( STRINGVector, XLALDestroyStringVector, "<string>,..." ),
  REGULAR_MAP_ENTRY ( REAL8Vector, XLALDestroyREAL8Vector, "<8-byte real>,..." ),
  REGULAR_MAP_ENTRY ( INT4Vector, XLALDestroyINT4Vector, "<4-byte integer>,..." )
};


// ---------- The module-local linked list to hold the user-variables
static LALUserVariable UVAR_vars;	// empty head
static CHAR *program_path = NULL;	// keep a pointer to the program path
static CHAR *program_name = NULL;	// keep a pointer to the program name


/**
 * An optional brief description of the program printed as part of the help page.
 */
const char *lalUserVarHelpBrief = NULL;


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
XLALRegisterUserVar ( void *cvar,		/**< pointer to the actual C-variabe to link to this user-variable */
                      const CHAR *name,		/**< name of user-variable to register */
                      UserVarType type,		/**< variable type (int,bool,string,real) */
                      CHAR optchar,		/**< optional short-option character */
                      UserVarCategory category,	/**< sets category to this */
                      const CHAR *help		/**< help-string explaining this input-variable */
                      )
{
  XLAL_CHECK ( cvar != NULL, XLAL_EINVAL );
  XLAL_CHECK ( name != NULL, XLAL_EINVAL );
  XLAL_CHECK ( strlen(name) < sizeof(UVAR_vars.name), XLAL_EINVAL, "User-variable name '%s' is too long", name );
  XLAL_CHECK ( (category > UVAR_CATEGORY_START) && (category < UVAR_CATEGORY_END), XLAL_EINVAL );
  XLAL_CHECK ( help != NULL, XLAL_EINVAL );
  XLAL_CHECK ( strlen(help) < sizeof(UVAR_vars.help), XLAL_EINVAL, "User-variable help '%s' is too long", help );

  // check that neither short- nor long-option are used by help
  XLAL_CHECK ( strcmp ( name, "help" ) != 0, XLAL_EINVAL, "Long-option name '--%s' is reserved for help!\n", name );
  XLAL_CHECK ( optchar != 'h', XLAL_EINVAL, "Short-option '-%c' is reserved for help!\n", optchar );

  // find end of uvar-list && check that neither short- nor long-option are taken already
  LALUserVariable *ptr = &UVAR_vars;
  while ( ptr->next != NULL )
    {
      ptr = ptr->next;

      // long-option name taken already?
      XLAL_CHECK ( strcmp ( name, ptr->name ) != 0, XLAL_EINVAL, "Long-option name '--%s' already taken!\n", name );
      // short-option character taken already?
      XLAL_CHECK ( (optchar == 0) || (ptr->optchar == 0) || (optchar != ptr->optchar), XLAL_EINVAL, "Short-option '-%c' already taken (by '--%s')!\n", optchar, ptr->name );

    } // while ptr->next

  // append new entry at the end
  XLAL_CHECK ( (ptr->next = XLALCalloc (1, sizeof(LALUserVariable))) != NULL, XLAL_ENOMEM );

  // set pointer to newly created entry
  ptr = ptr->next;

  // copy entry name, replacing '_' with '-' so that
  // e.g. uvar->a_long_option maps to --a-long-option
  XLALStringReplaceChar( strncpy( ptr->name, name, sizeof(ptr->name) - 1 ), '_', '-' );

  // copy entry help string
  strncpy( ptr->help, help, sizeof(ptr->help) - 1 );
  format_user_var_names( ptr->help );

  // fill in entry values
  ptr->type 	= type;
  ptr->optchar 	= optchar;
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
 *
 * If \p *should_exit is TRUE when this function returns, the
 * caller should exit immediately.
 */
int
XLALUserVarReadCmdline ( BOOLEAN *should_exit, int argc, char *argv[] )
{
  XLAL_CHECK ( should_exit != NULL, XLAL_EFAULT );
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL, "Input error, NULL argv[] pointer passed.\n" );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "Internal error, no UVAR memory allocated. Did you register any user-variables?" );

  *should_exit = 0;

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
  while ( 1 )
    {

      // call LALgetopt_long()
      char *old_argv0 = argv[0];
      argv[0] = program_name;   // use program_name for LALgetopt_long() error messages
      int longindex = -1;
      int c = LALgetopt_long(argc, argv, optstring, long_options, &longindex);
      argv[0] = old_argv0;
      if ( c == -1 ) {   // LALgetopt_long() is done
        break;
      }
      if ( c == '?' || c == ':' ) {   // LALgetopt_long() returned an error
        XLALUserVarPrintUsage( stderr );
        *should_exit = 1;
        return XLAL_SUCCESS;
      }

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
            if ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, LALoptarg ) != XLAL_SUCCESS )
              {
                XLALPrintError( "\n%s: could not parse value given to option " UVAR_FMT "\n\n", program_name, ptr->name );
                *should_exit = 1;
                return XLAL_SUCCESS;
              }
          }
	  break;

	default:
          // all other UVAR_TYPE_ types can be handled canonically: first destroy previous value, the parse new one
          if ( UserVarTypeMap [ ptr->type ].destructor != NULL )
            {
              UserVarTypeMap [ ptr->type ].destructor( *(char**)ptr->varp );
              *(char**)ptr->varp = NULL;
            } // if a destructor was registered
          if ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, LALoptarg ) != XLAL_SUCCESS )
            {
              XLALPrintError( "\n%s: could not parse value given to option " UVAR_FMT "\n\n", program_name, ptr->name );
              *should_exit = 1;
              return XLAL_SUCCESS;
            }
	  break;

	} // switch ptr->type

      switch ( ptr->was_set ) {
      case 0:    // this variable has not been set; mark as set on command line
        ptr->was_set = 2;
        break;
      case 1:    // this variable has been set in configuration file; print warning
        XLALPrintError ( "\n%s: option " UVAR_FMT " is overriding value set in configuration file!\n\n", program_name, ptr->name );
        break;
      default:   // this variable has been set before; error
        XLALUserVarCheck( should_exit, 0, "option " UVAR_FMT " was set more than once!", ptr->name );
        return XLAL_SUCCESS;
      }

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
      *should_exit = 1;
      return XLAL_SUCCESS;
    } // trailing non-option arguments found

  XLALFree (long_options);
  long_options=NULL;

  return XLAL_SUCCESS;

} // XLALUserVarReadCmdline()


/**
 * Read config-variables from cfgfile and parse into input-structure.
 * An error is reported if the config-file reading fails, but the
 * individual variable-reads are treated as optional
 *
 * If \p *should_exit is TRUE when this function returns, the
 * caller should exit immediately.
 */
int
XLALUserVarReadCfgfile ( BOOLEAN *should_exit, const CHAR *cfgfile )
{
  XLAL_CHECK ( should_exit != NULL, XLAL_EFAULT );
  XLAL_CHECK ( cfgfile != NULL, XLAL_EINVAL );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No memory allocated in UVAR_vars.next, did you register any user-variables?\n" );

  *should_exit = 0;

  LALParsedDataFile *cfg = NULL;
  XLAL_CHECK ( XLALParseDataFile ( &cfg, cfgfile ) == XLAL_SUCCESS, XLAL_EFUNC );

  // step through all user-variable: read those with names from config-file
  LALUserVariable *ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {

      XLAL_CHECK ( (ptr->type > UVAR_TYPE_START) && (ptr->type < UVAR_TYPE_END), XLAL_EFAILED, "Invalid UVAR_TYPE '%d' outside of [%d,%d]\n", ptr->type, UVAR_TYPE_START+1, UVAR_TYPE_END-1 );

      BOOLEAN wasRead;
      CHAR *valString = NULL;       // first read the value as a string
      XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &valString, cfg, NULL, ptr->name, &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( wasRead ) // if successful, parse this as the desired type
        {
          // destroy previous value, is applicable, then parse new one
          if ( UserVarTypeMap [ ptr->type ].destructor != NULL )
            {
              UserVarTypeMap [ ptr->type ].destructor( *(char**)ptr->varp );
              *(char**)ptr->varp = NULL;
            } // if a destructor was registered
          if ( UserVarTypeMap [ ptr->type ].parser( ptr->varp, valString ) != XLAL_SUCCESS )
            {
              XLALPrintError( "\n%s: could not parse value given to option " UVAR_FMT "\n\n", program_name, ptr->name );
              *should_exit = 1;
              return XLAL_SUCCESS;
            }
          XLALFree (valString);

          switch ( ptr->was_set ) {
          case 0:    // this variable has not been set; mark as set in configuration file
            ptr->was_set = 1;
            break;
          default:   // this variable has been set before; error
            XLALUserVarCheck( should_exit, 0, "configuration option `%s' was set more than once!", ptr->name );
            return XLAL_SUCCESS;
          }

        } // if wasRead

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
 * Print a one-line usage string
 */
int
XLALUserVarPrintUsage ( FILE *file )
{

  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  /* Print usage: only required and regular arguments are printed */
  fprintf( file, "\nUsage: %s [-h|--help] [@<config-file>]", program_name );
  for ( LALUserVariable *ptr = &UVAR_vars; (ptr=ptr->next) != NULL; )
    {
      switch ( ptr->category )
        {
        case UVAR_CATEGORY_DEVELOPER:
        case UVAR_CATEGORY_DEPRECATED:
        case UVAR_CATEGORY_DEFUNCT:
          continue;
        default:
          break;
        }
      fprintf( file, " " );
      if ( ptr->category != UVAR_CATEGORY_REQUIRED ) {
        fprintf( file, "[" );
      }
      if ( ptr->optchar != 0 ) {
        fprintf( file, "-%c|", ptr->optchar );
      }
      fprintf( file,"--%s", ptr->name);
      if ( ptr->category != UVAR_CATEGORY_REQUIRED ) {
        fprintf( file, "]" );
      }
    }
  fprintf( file,"\n\n");

  return XLAL_SUCCESS;

} // XLALUserVarPrintUsage()

/*
 * Format strings so that `--user_variables` have '_' replaced with '-'
 */
void
format_user_var_names( char *s )
{
  while ( ( s = strchr( s, '`' ) ) != NULL )
    {
      while ( *s != '\0' && *s != '\'' )
        {
          if ( *s == '_' ) *s = '-';
          ++s;
        }
    }
}

/*
 * Print text wrapped to a given line width
 */
void
fprint_wrapped( FILE *file, int line_width, const char *prefix, char *text )
{

  /* Adjust line width */
  line_width -= strlen(prefix) + 1;

  /* If line width is too short, just assume very long lines */
  if ( line_width < 1 ) {
    line_width = INT_MAX;
  }

  /* Iterate over text */
  char *pstart = text, *pbreak = NULL;
  for ( char *pend = text; *pend != '\0'; ++pend )
    {

      /* Record position of last space, in order to break line */
      if ( isspace( *pend ) ) {
        pbreak = pend;
      }

      /* If at end of line, or encountered a newline */
      if ( ( pend - pstart ) == line_width || *pend == '\n' ) {

        if ( pbreak != NULL ) {   /* If we have a space character */

          /* Print text up to before last space character */
          const char old_pbreak = *pbreak;
          *pbreak = '\0';
          fprintf( file, "%s%s\n", prefix, pstart );
          *pbreak = old_pbreak;

          /* Start from next non-printed character */
          pstart = pend = pbreak + 1;

          /* Reset space character */
          pbreak = NULL;

        } else {

          /* Print unbroken line, ending with hyphen */
          const char old_pend = *pend;
          *pend = '\0';
          fprintf( file, "%s%s-\n", prefix, pstart );
          *pend = old_pend;

          /* Start from next non-printed character */
          pstart = pend;

        }
      }

    }

  /* Print remaining text */
  if ( strlen( pstart ) > 0 ) {
    fprintf( file, "%s%s\n", prefix, pstart );
  }

}

/**
 * Print help page
 */
int
XLALUserVarPrintHelp ( FILE *file )
{

  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  /* Determine terminal line width, or set to zero otherwise */
  int line_width = 0;
#if defined(HAVE_ISATTY) && defined(HAVE_FILENO) && defined(HAVE_IOCTL)
#if HAVE_DECL_TIOCGWINSZ && HAVE_STRUCT_WINSIZE_WS_COL
  if ( isatty( fileno( file ) ) ) {
    struct winsize w;
    ioctl( fileno( file ), TIOCGWINSZ, &w );
    line_width = w.ws_col;
  }
#endif
#endif

  /* Pipe output through pager if possible */
  FILE *f = file;
#if defined(PAGER) && defined(HAVE_POPEN) && defined(HAVE_PCLOSE)
  f = popen(PAGER, "w");
  if ( f == NULL ) {
    f = file;
  }
#endif
  fflush( f );

  /* Print program name and synopsis of command line syntax */
  fprintf( f, "\nNAME\n" );
  fprintf( f, "       %s", program_name );
  if ( lalUserVarHelpBrief != NULL ) {
    fprintf( f, " - %s", lalUserVarHelpBrief );
  }
  fprintf( f, "\n" );
  fprintf( f, "\nSYNOPSIS\n" );
  fprintf( f, "       %s -h|--help\n", program_name );
  fprintf( f, "       %s [@<config-file>] [<options>...]\n", program_name );

  /* Print options in sections */
  const char* section_headers[] = { "OPTIONS", "DEVELOPER OPTIONS", "DEPRECATED OPTIONS" };
  for ( size_t section = 0; section < XLAL_NUM_ELEM(section_headers); ++section )
    {
      BOOLEAN print_section_header = 1;

      /* Go through all user variables */
      for ( LALUserVariable *ptr = &UVAR_vars; (ptr=ptr->next) != NULL; )
        {

          /* Decide what section to print option in, and with what formatting */
          size_t print_section = 0;
          BOOLEAN print_format_help = 1, print_default_value = 1;
          switch ( ptr->category )
            {
            case UVAR_CATEGORY_DEVELOPER:
              print_section = 1;
              break;
            case UVAR_CATEGORY_DEPRECATED:
              print_section = 2;
              print_format_help = print_default_value = 0;
              break;
            case UVAR_CATEGORY_DEFUNCT:
              continue;
            default:
              break;
            }
          XLAL_CHECK( print_section < XLAL_NUM_ELEM(section_headers), XLAL_EFAILED );

          if ( print_section == section )
            {

              /* Print section header */
              if ( print_section_header )
                {
                  fprintf( f, "\n%s\n", section_headers[section] );
                  print_section_header = 0;
                }

              /* Print option, format help, and default value */
              fprintf( f, "       " );
              if ( ptr->optchar != 0 )
                {
                  fprintf( f, "-%c, ", ptr->optchar );
                }
              fprintf( f, "--%s", ptr->name );
              if ( print_format_help )
                {
                  fprintf( f, "=%s", UserVarTypeMap [ ptr->type ].format_help );
                }
              if ( print_default_value )
                {
                  if ( ptr->category == UVAR_CATEGORY_REQUIRED )
                    {
                      fprintf( f, " [required]" );
                    }
                  else if ( ptr->category == UVAR_CATEGORY_NODEFAULT )
                    {
                      fprintf( f, " [optional]" );
                    }
                  else
                    {
                      char *valstr;
                      XLAL_CHECK( (valstr = UserVarTypeMap [ ptr->type ].printer( ptr->varp )) != NULL, XLAL_EFUNC );
                      fprintf( f, " [default: %s]",  valstr );
                      XLALFree( valstr );
                    }
                }
              fprintf( f, "\n" );

              /* Print option help string */
              fprint_wrapped( f, line_width, "           ", ptr->help );
              fprintf( f, "\n" );

            }

        }

    }

  /* Close pipe to pager if used */
  fflush( f );
#if defined(PAGER) && defined(HAVE_POPEN) && defined(HAVE_PCLOSE)
  if ( f != file ) {
    pclose( f );
  }
#endif

  return XLAL_SUCCESS;

} // XLALUserVarPrintHelp()


/**
 * Put all the pieces together, and basically does everything:
 * print help (if requested), get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 *
 * If \p *should_exit is TRUE when this function returns, the
 * program should exit immediately with a non-zero status.
 */
int
XLALUserVarReadAllInput ( BOOLEAN *should_exit, int argc, char *argv[] )
{
  XLAL_CHECK ( should_exit != NULL, XLAL_EFAULT );
  XLAL_CHECK ( argc > 0, XLAL_EINVAL );
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );
  XLAL_CHECK ( argv[0] != NULL, XLAL_EINVAL );
  XLAL_CHECK ( UVAR_vars.next != NULL, XLAL_EINVAL, "No UVAR memory allocated. Did you register any user-variables?" );

  *should_exit = 0;

  // keep a module-local pointer to the executable path/name
  program_path = argv[0];
  program_name = strrchr( program_path, '/' );
  if ( program_name == NULL ) {
    program_name = program_path;
  } else {
    ++program_name;
  }

  // ---------- manually parse command-line for help/usage arguments
  for ( INT4 i = 1; i < argc; i++ )
    {
      XLAL_CHECK( argv[i] != NULL, XLAL_EINVAL, "argc = %d, but argv[%d] == NULL!\n", argc, i );
      if ( strcmp( argv[i], "-h" ) == 0 )
        {
          XLALUserVarPrintUsage( stdout );
          *should_exit = 1;
          return XLAL_SUCCESS;
        }
      else if ( strcmp( argv[i], "--help" ) == 0 || strcmp( argv[i], "-help" ) == 0 )
        {
          XLALUserVarPrintHelp( stdout );
          *should_exit = 1;
          return XLAL_SUCCESS;
        }
    }

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
      XLAL_CHECK ( XLALUserVarReadCfgfile ( should_exit, cfgfile_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( *should_exit ) {
        return XLAL_SUCCESS;
      }
      XLALFree (cfgfile_name);
    }

  // ---------- now parse cmdline: overloads previous config-file settings
  XLAL_CHECK ( XLALUserVarReadCmdline ( should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( *should_exit ) {
    return XLAL_SUCCESS;
  }

  // ---------- handle special options that need some action ----------
  BOOLEAN skipCheckRequired = FALSE;
  for ( LALUserVariable *ptr = &UVAR_vars; (ptr=ptr->next) != NULL; )
    {

      // check 'special' category, which suppresses the CheckRequired test
      if ( (ptr->category == UVAR_CATEGORY_SPECIAL) && ptr->was_set ) {
	skipCheckRequired = TRUE;
      }

      // handle DEPRECATED options by outputting a warning (on error-level to make this very noticeable!)
      if ( ptr->category == UVAR_CATEGORY_DEPRECATED && ptr->was_set ) {
        XLALPrintError ("\n%s: option " UVAR_FMT " is DEPRECATED: %s\n\n", program_name, ptr->name, ptr->help );
      }

      // handle DEFUNCT options by throwing an error:
      XLALUserVarCheck( should_exit, ptr->category != UVAR_CATEGORY_DEFUNCT || !ptr->was_set, "option " UVAR_FMT " is DEFUNCT: %s", ptr->name, ptr->help );
      if ( *should_exit ) {
        return XLAL_SUCCESS;
      }

    } // while ptr = ptr->next

  // check that all required input-variables have been specified
  if ( !skipCheckRequired ) {

    // go through list of uvars
    for ( LALUserVariable *ptr = &UVAR_vars; (ptr=ptr->next) != NULL; )
      {
        XLALUserVarCheck( should_exit, ptr->category != UVAR_CATEGORY_REQUIRED || ptr->was_set, "required option " UVAR_FMT " has not been specified!", ptr->name );
        if ( *should_exit ) {
          return XLAL_SUCCESS;
        }
      }

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
 * If \p assertion is false, print the given error message, then the help usage;
 * \p should_exit is then set to true.
 */
void
XLALUserVarCheck( BOOLEAN *should_exit, const int assertion, const CHAR *fmt, ... )
{
  if ( !( *should_exit ) && !assertion ) {
    char buf[2048];
    va_list ap;
    va_start( ap, fmt );
    vsnprintf( buf, sizeof( buf ), fmt, ap );
    va_end( ap );
    format_user_var_names( buf );
    fprintf( stderr, "\n%s: %s\n", program_name, buf );
    XLALUserVarPrintUsage( stderr );
    fflush( stderr );
    *should_exit = 1;
  }
} // XLALUserVarCheck()

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
    XLAL_CHECK_NULL ( (record = XLALStringAppend ( record, program_path)) != NULL, XLAL_EFUNC );
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
LALUserVarReadCmdline (LALStatus *status, BOOLEAN *should_exit, int argc, char *argv[])
{
  INITSTATUS(status);
  if ( XLALUserVarReadCmdline(should_exit, argc, argv) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALUserVarReadCmdline() failed with code %d\n", xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} // LALUserVarReadCmdline()

/** \deprecated use XLALUserVarReadAllInput() instead */
void
LALUserVarReadAllInput (LALStatus *status, BOOLEAN *should_exit, int argc, char *argv[])
{
  INITSTATUS(status);
  if ( XLALUserVarReadAllInput ( should_exit, argc, argv ) != XLAL_SUCCESS ) {
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
  if ( XLALRegisterUserVar ( cvar, name, UVAR_TYPE_REAL8, optchar, category, helpstr ) != XLAL_SUCCESS ) {
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
  if ( XLALRegisterUserVar ( cvar, name, UVAR_TYPE_INT4, optchar, category, helpstr ) != XLAL_SUCCESS ) {
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
  if ( XLALRegisterUserVar ( cvar, name, UVAR_TYPE_BOOLEAN, optchar, category, helpstr ) != XLAL_SUCCESS ) {
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
  if ( XLALRegisterUserVar ( cvar, name, UVAR_TYPE_STRING, optchar, category, helpstr ) != XLAL_SUCCESS ) {
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
  if ( XLALRegisterUserVar ( cvar, name, UVAR_TYPE_STRINGVector, optchar, category, helpstr ) != XLAL_SUCCESS ) {
    XLALPrintError ("Call to XLALRegisterUserVar() failed: %d\n", xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN(status);
}

/** \deprecated use XLALUserVarReadCfgfile() instead */
void
LALUserVarReadCfgfile (LALStatus *status,
		       BOOLEAN *should_exit,
		       const CHAR *cfgfile) 	   /* name of config-file */
{

  INITSTATUS(status);
  if ( XLALUserVarReadCfgfile ( should_exit, cfgfile ) != XLAL_SUCCESS ) {
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }
  RETURN (status);
} // LALUserVarReadCfgfile()
