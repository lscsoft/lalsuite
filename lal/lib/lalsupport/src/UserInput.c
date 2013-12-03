/*
 * Copyright (C) 2010 Reinhard Prix (xlalified)
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

#include "getopt.h"

#include <lal/LALStdio.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>

#define TRUE  (1==1)
#define FALSE (1==0)

/* Defines the type of a "user-variable": bool, int, real or string.
 * Should be used only internally !!
 */
typedef enum {
  UVAR_BOOL,    /* boolean */
  UVAR_INT4,    /* integer */
  UVAR_REAL8,   /* float */
  UVAR_STRING,  /* string */
  UVAR_CSVLIST, /* list of comma-separated values */
  UVAR_LAST
} UserVarType;

/**
 * Linked list to hold the complete information about the user-variables.
 */
typedef struct tagLALUserVariable {
  const CHAR *name;	/**< full name */
  UserVarType type;	/**< type: bool, int, float or string */
  CHAR optchar;		/**< cmd-line character */
  const CHAR *help;	/**< help-string */
  void *varp;		/**< pointer to the actual C-variable */
  UserVarState state;	/**< state (empty, default, set) */
  struct tagLALUserVariable *next; /* linked list */
} LALUserVariable;

/** The module-local linked list to hold the user-variables */
static LALUserVariable UVAR_vars;	/**< empty head */
static const CHAR *program_name;	/**< keep a pointer to the program name */

/* needed for command-line parsing */
extern char *optarg;
extern int optind, opterr, optopt;

/* ---------- internal prototypes ---------- */

/* ----- XLAL interface ----- */
int XLALRegisterUserVar (const CHAR *name, UserVarType type, CHAR optchar, UserVarState flag, const CHAR *helpstr, void *cvar);
CHAR *XLALUvarValue2String (LALUserVariable *uvar);

CHAR *XLALUvarType2String (LALUserVariable *uvar);

CHAR *XLAL_copy_string_unquoted ( const CHAR *in );
void check_and_mark_as_set ( LALUserVariable *varp );


/* ----- LAL interface [deprecated] ----- */
static void RegisterUserVar (LALStatus *, const CHAR *name, UserVarType type, CHAR optchar,
			     UserVarState flag, const CHAR *helpstr, void *cvar);

/*---------- Function definitions ---------- */

/* these are type-specific wrappers to allow tighter type-checking! */
/** Register a user-variable of type REAL8, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterREALUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, REAL8 *cvar )
{
  return XLALRegisterUserVar (name, UVAR_REAL8, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type INT4, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterINTUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, INT4 *cvar )
{
  return XLALRegisterUserVar (name, UVAR_INT4, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type BOOLEAN, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterBOOLUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, BOOLEAN *cvar )
{
  return XLALRegisterUserVar (name, UVAR_BOOL, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of type CHAR*, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterSTRINGUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, CHAR **cvar )
{
  return XLALRegisterUserVar (name, UVAR_STRING, optchar, flag, helpstr, cvar);
}

/** Register a user-variable of 'list' type LALStringVector, see XLALRegisterUserVar() for API documentation */
int
XLALRegisterLISTUserVar ( const CHAR *name, CHAR optchar, UserVarState flag, const CHAR *helpstr, LALStringVector **cvar)
{
  return XLALRegisterUserVar ( name, UVAR_CSVLIST, optchar, flag, helpstr, cvar );
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
 * ==> use one of the 4 wrappers: XLALRegisterREALUserVar(),
 * XLALRegisterINTUserVar(), XLALRegisterBOOLUserVar(), XLALRegisterSTRINGUserVar().
 *
 */
int XLALRegisterUserVar ( const CHAR *name,	/**< name of user-variable to register */
                          UserVarType type,	/**< variable type (int,bool,string,real) */
                          CHAR optchar,		/**< optional short-option character */
                          UserVarState flag,	/**< sets state flag to this */
                          const CHAR *helpstr,	/**< help-string explaining this input-variable */
                          void *cvar		/**< pointer to the actual C-variabe to link to this user-variable */
                          )
{
  if (cvar == NULL || name == NULL ) {
    XLALPrintError ("%s: invalue NULL input 'cvar' or 'name'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  LALUserVariable *ptr = NULL;

  /* find end of uvar-list && check that neither short- nor long-option are taken already */
  ptr = &UVAR_vars;
  while ( ptr->next && (ptr = ptr->next) )
    {
      if ( name && ptr->name && !strcmp(name, ptr->name) ) {
        XLALPrintError ("%s: Long-option name '--%s' is already taken!\n", __func__, name );
        XLAL_ERROR ( XLAL_EINVAL );
      }
      if ( optchar && ptr->optchar && (optchar == ptr->optchar) ) {
        XLALPrintError ("%s: Short-option '-%c' is already taken (by '--%s')!\n", __func__, ptr->optchar, ptr->name );
        XLAL_ERROR ( XLAL_EINVAL );
      }
    } /* while user-var list */

  /* create new entry */
  if ( (ptr->next = XLALCalloc (1, sizeof(LALUserVariable))) == NULL ) {
    XLALPrintError ("%s: Failed to XLALCalloc (1, sizeof(LALUserVariable)\n", __func__ );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* set pointer to newly created entry (Note: the head remains empty!) */
  ptr = ptr->next;

  /* fill in values */
  ptr->name = name;
  ptr->type = type;
  ptr->optchar = optchar;
  ptr->help = helpstr;
  ptr->varp = cvar;
  ptr->state = flag;

  return XLAL_SUCCESS;

} /* XLALRegisterUserVar() */

/**
 * Free all memory associated with user-variable linked list
 */
void
XLALDestroyUserVars( void )
{
  LALUserVariable *ptr, *lastptr;

  /* step through user-variables: free list-entries and all allocated strings */
  ptr = &(UVAR_vars);
  lastptr = NULL;
  while ( (ptr=ptr->next) != NULL)
    {
      /* is an allocated string here? */
      if ( (ptr->type == UVAR_STRING) && (*(CHAR**)(ptr->varp) != NULL) )
	{
	  LALFree ( *(CHAR**)(ptr->varp) );
	  *(CHAR**)(ptr->varp) = NULL;		/* IMPORTANT: reset user-variable to NULL ! */
	}
      else if ( ptr->type == UVAR_CSVLIST ) {
	XLALDestroyStringVector ( *(LALStringVector**)ptr->varp );
	*(LALStringVector**)(ptr->varp) = NULL;
      }

      /* free list-entry behind us (except for the head) */
      if (lastptr) {
	LALFree (lastptr);
	lastptr=NULL;
      }

      lastptr = ptr;

    } /* while ptr->next */

  if (lastptr) {
    LALFree (lastptr);
    lastptr=NULL;
  }

  /* clean head */
  memset (&UVAR_vars, 0, sizeof(UVAR_vars));

  return;

} /* XLALDestroyUserVars() */


/**
 * Parse command-line into UserVariable array
 */
int
XLALUserVarReadCmdline ( int argc, char *argv[] )
{
  INT4 c;
  UINT4 pos;
  UINT4 numvars;
  LALUserVariable *ptr = NULL;
  char optstring[512] = "\0";	/* string of short-options, should be easily enough */
  struct option *long_options;
  int longindex = -1;
  CHAR *strp;
  LALStringVector *csv;

  if (!argv ) {
    XLALPrintError ("%s: Input error, NULL argv[] pointer passed.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !UVAR_vars.next ) {
    XLALPrintError ("%s: Internal error, no UVAR memory allocated. Did you register any user-variables?", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* build optstring of short-options */
  ptr = &UVAR_vars;	/* set to empty head */
  pos = 0;
  numvars = 0;

  /* special treatment of head, which could contain the debug-option:
   * NOTE: this one will *NOT* be treated by the remaining function,
   * (as the head is always skipped), but it has to be in the optstring
   * to avoid an error if specified on the command-line.
   * Treatment of debug-option reading has to be done separately using
   * UVARgetDebugLevel()
   */
  if ( (ptr->help != NULL) && (ptr->optchar != 0) )
    {
      numvars ++;
      optstring[pos++] = ptr->optchar;
      optstring[pos++] = ':';		/* requires an argument */
    } /* if debug-option */

  /* treat the remaining "normal" entries */
  while ( (ptr = ptr->next) != NULL )
    {
      numvars ++;			/* counter number of user-variables */
      if (ptr->optchar == 0)		/* if no short-option given, ignore */
	continue;
      optstring[pos++] = ptr->optchar;
      optstring[pos++] = ':';		/* everything but bool takes an argument */
      if (ptr->type == UVAR_BOOL)	/* but for BOOL its optional */
	optstring[pos++] = ':';
    } /* while ptr->next */
  optstring[pos] = '\0';

  /* fill option-struct for long-options */
  long_options = LALCalloc (1, (numvars+1) * sizeof(struct option));
  ptr = &UVAR_vars;	/* start again from beginning */
  pos = 0;
  while ( (ptr= ptr->next) != NULL)
    {
      if (ptr->name == NULL)		/* if no long-option given, ignore */
	continue;
      long_options[pos].name = ptr->name;
      long_options[pos].has_arg = (ptr->type == UVAR_BOOL) ? optional_argument : required_argument;
      long_options[pos].flag = NULL;	/* get val returned from getopt_long() */
      long_options[pos].val = 0;	/* we use longindex to find long-options */
      pos ++;
    } /* while ptr->next */
  /* null-terminate array */
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
  optind = 0; 	/* reset getopt(), getopt_long() */

  /* parse the command-line */
  while ( (c = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1 )
    {
      if (c == '?') {
	XLALPrintError ( "%s: ERROR: unkown command-line option encountered\n", __func__ );
	XLALPrintError ( "see '%s --help' for usage-help\n\n", argv[0]);
	XLAL_ERROR ( XLAL_EDOM );
      }
      if (c != 0) 	/* find short-option character */
	{
	  ptr = &UVAR_vars;
	  do {
	    if (c == ptr->optchar)
	      break;
	  } while ( (ptr=ptr->next) != NULL);
	} /* if short-option */
      else	/* find long-option: returned in longindex */
	{
	  ptr = &UVAR_vars;
	  while ( (ptr=ptr->next) != NULL)
	    if ( !strcmp (long_options[longindex].name, ptr->name) )
	      break;
	} /* if long-option */

      if (ptr == NULL) {	/* should not be possible: nothing found at all... */
	XLALPrintError ( "%s: ERROR: failed to find option.. this points to a coding-error!\n", __func__);
	XLAL_ERROR ( XLAL_EDOM );
      }

      /* if we found the debug-switch, ignore it (has been handled already */
      if (ptr == &UVAR_vars)
	continue;

      switch (ptr->type)
	{
	  INT2 ans;
	case UVAR_BOOL:
	  ans = -1;

	  /* subtlety with optional argument: it's not necessarily found in the *same* argv-entry
	   * eg, if no '=' was used, so we have to check for that case by hand: */

	  /* if the next entry is not an option, take it as an argument */
	  if (optarg == NULL && (optind<argc) && (argv[optind][0]!='-') && (argv[optind][0]!='@') )
	    optarg = argv[optind];

	  if ( optarg == NULL )	/* no argument found at all: defaults to TRUE */
	    {
	      ans = 1;
	    }
	  else	/* parse bool-argument: should be consistent with bool-parsing in ConfigFile!! */
	    {
	      /* get rid of case ambiguities */
	      if ( XLALStringToLowerCase (optarg) != XLAL_SUCCESS ) {
                XLAL_ERROR ( XLAL_EFUNC );
              }

	      if      ( !strcmp(optarg, "yes") || !strcmp(optarg, "true") || !strcmp(optarg,"1") )
		ans = 1;
	      else if ( !strcmp (optarg, "no") || !strcmp(optarg,"false") || !strcmp(optarg,"0") )
		ans = 0;
	      else {	/* failed to parse BOOL properly */
		XLALPrintError ( "%s: Illegal bool-value `%s`\n\n", __func__, optarg);
		XLAL_ERROR ( XLAL_EDOM );
	      }
	    } /* parse bool-argument */

	  /* only set if we properly parsed something */
	  if (ans != -1) {
	    *(BOOLEAN*)(ptr->varp)  = (BOOLEAN)ans;
	    check_and_mark_as_set ( ptr );
	  }

	  break;

	case UVAR_INT4:
	  if ( 1 != sscanf ( optarg, "%" LAL_INT4_FORMAT, (INT4*)(ptr->varp)) )
	    {
	      XLALPrintError ("%s: Illegal INT4 commandline argument to --%s: '%s'\n\n", __func__, ptr->name, optarg);
	      XLAL_ERROR ( XLAL_EDOM );
	    }

	  check_and_mark_as_set ( ptr );
	  break;

	case UVAR_REAL8:
	  if ( 1 != sscanf ( optarg, "%" LAL_REAL8_FORMAT, (REAL8*)(ptr->varp)) )
	    {
	      XLALPrintError ("%s: Illegal REAL8 commandline argument to --%s: '%s'\n\n", __func__, ptr->name, optarg);
	      XLAL_ERROR ( XLAL_EDOM );
	    }

	  check_and_mark_as_set ( ptr );
	  break;

	case UVAR_STRING:
	  if (!optarg) {	/* should not be possible, but let's be paranoid */
	    XLALPrintError ( "%s: optarg==NULL, something went badly wrong ...\n", __func__ );
            XLAL_ERROR ( XLAL_EFAULT );
	  }
	  strp = *(CHAR**)(ptr->varp);
	  if ( strp != NULL) 	 /* something allocated here before? */
	    XLALFree ( strp );
	  if ( (strp = XLAL_copy_string_unquoted ( optarg )) == NULL ) {
            XLALPrintError ("%s: XLAL_copy_string_unquoted() failed.\n", __func__ );
            XLAL_ERROR ( XLAL_EFUNC );
	  }
	  /* return value */
	  *(CHAR**)(ptr->varp) = strp;
	  check_and_mark_as_set ( ptr );
	  break;

	case UVAR_CSVLIST:	/* list of comma-separated values */
	  csv = *(LALStringVector**)(ptr->varp);
	  if ( csv != NULL) { 	/* something allocated here before? */
	    XLALDestroyStringVector ( csv );
          }
	  if ( (csv = XLALParseCSV2StringVector ( optarg )) == NULL ) {
            XLALPrintError ("%s: XLALParseCSV2StringVector() failed on '%s'\n", __func__, optarg );
            XLAL_ERROR ( XLAL_EFUNC );
	  }
	  /* return value */
	  *(LALStringVector**)(ptr->varp) = csv;
	  check_and_mark_as_set ( ptr );
	  break;

	default:
	  XLALPrintError ( "%s: ERROR: unkown UserVariable-type encountered... points to a coding error!\n", __func__ );
	  XLAL_ERROR ( XLAL_EINVAL );
	  break;

	} /* switch ptr->type */

    } /* while getopt_long() */

  // check if there's any non-option strings left (except for a config-file specification '@file')
  if ( (optind == argc - 1) && (argv[optind][0] == '@' ) ) {
    optind ++;	// advance counter in case of one config-file specification (only one allowed)
  }
  if ( optind < argc ) // still stuff left? ==> error
    {
      XLALPrintError ( "\nGot non-option ARGV-elements: [ ");
      while (optind < argc) {
        XLALPrintError ("%s ", argv[optind++]);
      }
      XLALPrintError(" ]\n");
      XLAL_ERROR ( XLAL_EDOM );
    } // non-option arguments found

  XLALFree (long_options);
  long_options=NULL;

  return XLAL_SUCCESS;

} /* XLALUserVarReadCmdline() */


/**
 * Read config-variables from cfgfile and parse into input-structure.
 * An error is reported if the config-file reading fails, but the
 * individual variable-reads are treated as optional
 */
int
XLALUserVarReadCfgfile ( const CHAR *cfgfile ) 	   /**< [in] name of config-file */
{
  LALParsedDataFile *cfg = NULL;
  CHAR *stringbuf, *strp;
  LALUserVariable *ptr;
  BOOLEAN wasRead;
  LALStringVector *csv;

  if ( !UVAR_vars.next ) {
    XLALPrintError ("%s: no memory allocated in UVAR_vars.next, did you register any user-variables?\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( XLALParseDataFile ( &cfg, cfgfile ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Call to XLALParseDataFile() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* step through all user-variable: read those with names from config-file */
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->name == NULL)	/* ignore name-less user-variable */
	continue;

      wasRead = FALSE;

      switch (ptr->type)
	{
	case UVAR_BOOL:
          if ( XLALReadConfigBOOLVariable(ptr->varp, cfg, NULL, ptr->name, &wasRead) != XLAL_SUCCESS ) {
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if (wasRead)
	    check_and_mark_as_set ( ptr );
	  break;
	case UVAR_INT4:
	  if ( XLALReadConfigINT4Variable(ptr->varp, cfg, NULL, ptr->name, &wasRead) != XLAL_SUCCESS ) {
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if (wasRead)
	    check_and_mark_as_set ( ptr );
	  break;
	case UVAR_REAL8:
	  if ( XLALReadConfigREAL8Variable(ptr->varp, cfg, NULL, ptr->name, &wasRead) != XLAL_SUCCESS ) {
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if (wasRead)
	    check_and_mark_as_set ( ptr );
	  break;
	case UVAR_STRING:
	  stringbuf = NULL;
	  if ( XLALReadConfigSTRINGVariable( &stringbuf, cfg, NULL, ptr->name, &wasRead) != XLAL_SUCCESS ) {
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if ( wasRead && stringbuf)	/* did we find something? */
	    {
	      strp = *(CHAR**)(ptr->varp);
	      if ( strp != NULL) /* something allocated here before? */
		XLALFree ( strp );
	      /* return value */
	      *(CHAR**)(ptr->varp) = stringbuf;
	      check_and_mark_as_set ( ptr );
	    } /* if stringbuf */
	  break;

	case UVAR_CSVLIST:
	  stringbuf = NULL;
	  if ( XLALReadConfigSTRINGVariable(&stringbuf, cfg, NULL, ptr->name,&wasRead) != XLAL_SUCCESS ) {
            XLAL_ERROR ( XLAL_EFUNC );
          }
	  if ( wasRead && stringbuf)	/* did we find something? */
	    {
	      csv = *(LALStringVector**)(ptr->varp);
	      if ( csv != NULL) { 	/* something allocated here before? */
		XLALDestroyStringVector ( csv );
              }
	      if ( (csv = XLALParseCSV2StringVector ( stringbuf )) == NULL ) {
                XLALPrintError ("%s: XLALParseCSV2StringVector() failed with code %d\n", __func__, xlalErrno );
                XLAL_ERROR ( XLAL_EFUNC );
	      }
	      XLALFree ( stringbuf );
	      *(LALStringVector**)(ptr->varp) = csv;
	      check_and_mark_as_set ( ptr );
	    } /* if stringbuf */
	  break;
	default:
	  XLALPrintError ("%s: ERROR: unkown UserVariable-type encountered...points to a coding error!\n", __func__);
          XLAL_ERROR ( XLAL_EFAILED );
          break;

	} /* switch ptr->type */

    } /* while ptr->next */

  /* ok, that should be it: check if there were more definitions we did not read */
  if ( XLALCheckConfigReadComplete (cfg, CONFIGFILE_WARN) != XLAL_SUCCESS ) {
    XLAL_ERROR ( XLAL_EFUNC );
  }

  XLALDestroyParsedDataFile(cfg);

  return XLAL_SUCCESS;
} /* XLALUserVarReadCfgfile() */


#define UVAR_MAXHELPLINE  512	/* max length of one help-line */
#define UVAR_MAXDEFSTR    100 	/* max length of default-string */
#define UVAR_MAXFMTLEN    128   /* max length of help-line format-string */

/**
 * Assemble all help-info from uvars into a help-string.
 */
CHAR *
XLALUserVarHelpString ( const CHAR *progname )
{
  CHAR strbuf[UVAR_MAXHELPLINE];	/* should be enough for one line...*/
  CHAR defaultstr[UVAR_MAXDEFSTR]; 	/* for display of default-value */
  CHAR fmtStr[UVAR_MAXFMTLEN];		/* for building a dynamic format-string */
  CHAR optstr[10];			/* display of opt-char */
  /* we need strings for UVAR_BOOL, UVAR_INT4, UVAR_REAL8, UVAR_STRING: */
  const CHAR *typestr[] = {"BOOL", "INT", "REAL", "STRING", "LIST"};
  LALUserVariable *ptr;
  LALUserVariable *helpptr = NULL;	/* pointer to help-option */
  CHAR *helpstr = NULL;
  size_t newlen = 0;
  BOOLEAN haveDevOpt = 0;	/* have we got any 'developer'-options */

  /* check input consistency */
  if ( !UVAR_vars.next ) {
    XLALPrintError ("%s: Internal error, no UVAR memory allocated. Did you register any user-variables?", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* prepare first lines of help-string: info about config-file reading */
  newlen = 0;
  sprintf (strbuf, "Usage: %s [@ConfigFile] [options], where options are:\n\n", progname);
  if ( (helpstr = XLALCalloc (1, strlen(strbuf) + 1)) == NULL) {
    XLALPrintError ("%s: failed to XLALCalloc(1,%d)\n", __func__, strlen(strbuf) + 1 );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  strcpy (helpstr, strbuf);
  newlen += strlen (strbuf) + 1;

  /* ZEROTH PASS: find longest long-option name, for proper output formatting */
  ptr = &UVAR_vars;
  UINT4 nameFieldLen = 0;
  while ( (ptr=ptr->next) != NULL)
    {
      if ( (lalDebugLevel == 0) && (ptr->state & UVAR_DEVELOPER) )
	continue;	/* skip developer-options if debugLevel = 0 */
      if ( ptr->name )
        {
          UINT4 len = strlen ( ptr->name );
          if ( len > nameFieldLen ) nameFieldLen = len;
        } /* if name */
    } /* while options */

  /* put together the help-string. Allocate memory on-the-fly... */
  ptr = &UVAR_vars;
  /* special treatment of debug-option in the head (if present) */
  if ( ptr->help && ptr->optchar )
    {
      snprintf ( fmtStr, UVAR_MAXFMTLEN, "  -%%c    %%-%ds   %%-6s   %%s [%%d]\n",  nameFieldLen );
      fmtStr[UVAR_MAXFMTLEN-1]=0;

      sprintf (strbuf, fmtStr, ptr->optchar, " ", typestr[ptr->type], ptr->help, *(INT4*)(ptr->varp) );
      newlen += strlen (strbuf);
      if ( (helpstr = XLALRealloc (helpstr, newlen)) == NULL ) {
        XLALPrintError ("%s: failed to XLALRealloc (helpstr, %d)\n", __func__, newlen );
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */
    }

  snprintf ( fmtStr, UVAR_MAXFMTLEN, "  %%s --%%-%ds   %%-6s   %%s [%%s]\n", nameFieldLen );
  fmtStr[UVAR_MAXFMTLEN-1]=0;

  /* FIRST PASS: treat all "normal" entries excluding DEVELOPER-options */
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->state & UVAR_DEVELOPER)
	continue;	/* skip developer-options to be treated a second pass */

      if (ptr->state & UVAR_REQUIRED)
	strcpy (defaultstr, "REQUIRED");
      else if (ptr->state & UVAR_HELP)
	{
	  helpptr = ptr;	/* keep a pointer to the help-option for later */
	  strcpy (defaultstr, "");
	}
      else /* write the current default-value into a string */
	{
	  CHAR *valstr = NULL;
	  if ( (valstr = XLALUvarValue2String(ptr)) == NULL ) {
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }
	  strncpy (defaultstr, valstr, UVAR_MAXDEFSTR);	/* cut short for default-entry */
	  defaultstr[UVAR_MAXDEFSTR-1] = 0;
	  XLALFree (valstr);
	  valstr=NULL;
	}

      if (ptr->optchar != 0)
	sprintf (optstr, "-%c,", ptr->optchar);
      else
	strcpy (optstr, "   ");

      snprintf (strbuf, UVAR_MAXHELPLINE,  fmtStr,
		   optstr,
		   ptr->name ? ptr->name : "-NONE-",
		   typestr[ptr->type],
		   ptr->help ? ptr->help : "-NONE-",
		   defaultstr);

      /* now increase allocated memory by the right amount */
      newlen += strlen (strbuf);
      if ( (helpstr = LALRealloc (helpstr, newlen)) == NULL ) {
        XLALPrintError ("%s: failed to LALRealloc (helpstr, %d)\n", __func__, newlen );
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */

    } /* while ptr->next */

  /* ---------- SECOND PASS through user-options:
   * show DEVELOPER-options only if lalDebugLevel >= 1
   */
  if ( lalDebugLevel == 0)	/* only give instructions as to how to see developer-options */
    {
      CHAR buf[256];
      if ( UVAR_vars.optchar && helpptr && helpptr->name )
	sprintf (buf, "(e.g. --%s -%c1)", helpptr->name, UVAR_vars.optchar);
      else
	sprintf (buf, " ");

      sprintf (strbuf, "\n ----- Hint: use help with lalDebugLevel > 0 %s to see all 'developer-options' ----- \n", buf);
      newlen += strlen (strbuf);
      if ( (helpstr = LALRealloc (helpstr, newlen)) == NULL ) {
        XLALPrintError ( "%f: LALRealloc (helpstr, %d) failed.\n", __func__, newlen );
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */
    }
  else	/* lalDebugLevel > 0 */
    {
      strcpy(strbuf,
	     "\n   ---------- The following are 'Developer'-options not useful "
	     "for most users:----------\n\n");
      newlen += strlen (strbuf);
      if ( (helpstr = LALRealloc (helpstr, newlen)) == NULL ) {
        XLALPrintError ( "%f: LALRealloc (helpstr, %d) failed.\n", __func__, newlen );
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */

      ptr = &UVAR_vars;
      while ( (ptr=ptr->next) != NULL )
	{
	  CHAR *valstr = NULL;

	  if ( ! (ptr->state & UVAR_DEVELOPER) )	/* only treat developer-options */
	    continue;

	  haveDevOpt = 1;

	  if ( (valstr = XLALUvarValue2String(ptr)) == NULL ) {
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }
	  strncpy (defaultstr, valstr, UVAR_MAXDEFSTR);	/* cut short for default-entry */
	  defaultstr[UVAR_MAXDEFSTR-1] = 0;
	  XLALFree (valstr);
	  valstr = NULL;

	  if (ptr->optchar != 0)
	    sprintf (optstr, "-%c,", ptr->optchar);
	  else
	    strcpy (optstr, "   ");

	  snprintf (strbuf, UVAR_MAXHELPLINE,  fmtStr,
		       optstr,
		       ptr->name ? ptr->name : "-NONE-",
		       typestr[ptr->type],
		       ptr->help ? ptr->help : "-NONE-",
		       defaultstr);

	  /* now increase allocated memory by the right amount */
	  newlen += strlen (strbuf);
          if ( (helpstr = LALRealloc (helpstr, newlen)) == NULL ) {
            XLALPrintError ( "%f: LALRealloc (helpstr, %d) failed.\n", __func__, newlen );
            XLAL_ERROR_NULL ( XLAL_ENOMEM );
          }

	  strcat (helpstr, strbuf);	/* add this line to the helpstring */

	} /* while ptr->next: 2nd PASS for DEVELOPER-options */

      if ( !haveDevOpt )	/* no developer-options found: say something */
	{
	  strcpy(strbuf, "   -- NONE --\n\n");
	  newlen += strlen (strbuf);
          if ( (helpstr = LALRealloc (helpstr, newlen)) == NULL ) {
            XLALPrintError ( "%f: LALRealloc (helpstr, %d) failed.\n", __func__, newlen );
            XLAL_ERROR_NULL ( XLAL_ENOMEM );
          }

	  strcat (helpstr, strbuf);	/* add this line to the helpstring */
	} /* if !haveDevOpt */

    } /* if lalDebugLevel: output developer-options */

  /* finished: return final helpstring */
  return helpstr;

} /* XLALUserVarHelpString() */


/**
 * Put all the pieces together, and basically does everything:
 * get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 */
int
XLALUserVarReadAllInput ( int argc, char *argv[] )
{
  INT4 i;
  CHAR* fname = NULL;
  CHAR *tmp;
  LALUserVariable *ptr;
  BOOLEAN skipCheckRequired = FALSE;

  if ( !UVAR_vars.next ) {
    XLALPrintError ("%s: Internal error, no UVAR memory allocated. Did you register any user-variables?", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  program_name = argv[0];	/* keep modul-local pointer to executable name */

  /* pre-process command-line: have we got a config-file ? */
  for (i=1; i < argc; i++)
    {
      tmp = argv[i];
      if ( *tmp == '@' )
	{
	  if (fname != NULL) {
            XLALPrintError ("%s: can handle only one config-file passed on commandline!\n", __func__ );
            XLAL_ERROR ( XLAL_EDOM );
	  }

	  tmp ++;
	  if ( (fname = XLALCalloc (1, strlen(tmp) + 5 )) == NULL) {
            XLALPrintError("%s: XLALCalloc (1, %s) failed.\n", __func__, strlen(tmp) + 5 );
            XLAL_ERROR ( XLAL_ENOMEM );
	  }
	  /* NOTE: if the filename given is not a relative or absolute path,
	   * we want to ensure it is interpreted relative to the CURRENT directory,
	   * NOT relative to LAL_DATA_PATH (if set), as we cannot rely on it containg
	   * the local directory.
	   *
	   * ==> therefore we ensure that the path is relative to "./" in that case.
	   */
	  if ( (tmp[0] != '.') && (tmp[0] != '/') )
	    sprintf (fname, "./%s", tmp);
	  else
	    strcpy (fname, tmp);

	} /* if argument starts with '@' */

    } /* for i < argc */


  /* if config-file specified, read from that first */
  if (fname) {
    if ( XLALUserVarReadCfgfile (fname) != XLAL_SUCCESS ) {
      XLALPrintError ("%s: XLALUserVarReadCfgfile (%s) failed with code %d\n", __func__, fname, xlalErrno );
      XLAL_ERROR ( XLAL_EFUNC );
    }
    XLALFree (fname);
    fname=NULL;
  }

  /* now do proper cmdline parsing: overloads config-file settings */
  if ( XLALUserVarReadCmdline (argc, argv) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALUserVarReadCmdline() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* now check if help-string was requested */
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      if ( (ptr->state & UVAR_HELP) && (ptr->state & UVAR_WAS_SET) )
	{
	  CHAR *helpstring = NULL;
	  if ( ( helpstring = XLALUserVarHelpString(argv[0])) == NULL ) {
            XLALPrintError ("%s: XLALUserVarHelpString() failed.\n", __func__);
            XLAL_ERROR ( XLAL_EFUNC );
          }
          printf ("\n%s\n", helpstring);
	  XLALFree (helpstring);
	  helpstring=NULL;
	  return XLAL_SUCCESS;
	} /* if help requested */

      /* check 'special' flag, which suppresses the CheckRequired test */
      if ( (ptr->state & UVAR_SPECIAL) && (ptr->state & UVAR_WAS_SET) )
	skipCheckRequired = TRUE;
    }
  /* check that all required input-variables have been specified */
  if ( !skipCheckRequired ) {
    if ( XLALUserVarCheckRequired() != XLAL_SUCCESS ) {
      XLAL_ERROR ( XLAL_EFUNC );
    }
  } /* !skipCheckRequired */

  return XLAL_SUCCESS;

} /* XLALUserVarReadAllInput() */



/**
 * Has this user-variable been set by the user?
 * returns TRUE/FALSE
 */
int
XLALUserVarWasSet (const void *cvar)
{
  LALUserVariable *ptr;

  if (!cvar) {
    XLALPrintError ("%s: invalid NULL input\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* find this varname in the list of user-variables */
  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL)
    if ( ptr->varp == cvar)
      break;

  if (ptr == NULL) {
    XLALPrintError ("%s: Variable passed to UVARwasSet is not a registered User-variable\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* we found it: has it been set by user? */
  return ( (ptr->state & UVAR_WAS_SET) != 0 );

} /* XLALUserVarWasSet() */


/**
 * Check that all required user-variables have been set successfully.
 * Print error if not
 */
int
XLALUserVarCheckRequired (void)
{
  LALUserVariable *ptr;

  /* go through list of uvars */
  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL)
    if ( (ptr->state & UVAR_REQUIRED) && !(ptr->state & UVAR_WAS_SET))
      {
	XLALPrintError ("%s: Required user-variable `%s` has not been specified!\n\n", __func__, ptr->name);
	XLAL_ERROR ( XLAL_EDOM );
      }

  return XLAL_SUCCESS;

} /* XLALUserVarCheckRequired() */


/**
 * Return a log-string representing the <em>complete</em> user-input.
 * <em>NOTE:</em> we only record user-variables that have been set
 * by the user.
 */
CHAR *
XLALUserVarGetLog ( UserVarLogFormat format 	/**< output format: return as config-file or command-line */
                    )
{
  LALUserVariable *ptr = NULL;
  CHAR *record = NULL;
  CHAR *valstr;		/* buffer to hold value-string */
  CHAR *append;
  CHAR *typestr=NULL;
  UINT4 len, appendlen;

  /* initialize return-string */
  record = XLALMalloc (1);
  record[0] = 0;
  len = 0;

  if ( format == UVAR_LOGFMT_CMDLINE )
    {
      len += strlen ( program_name );
      if ( (record = XLALRealloc (record, len+1)) == NULL ) {
        XLALPrintError ("%s: XLALRealloc (%d, %d) failed.\n", __func__, record, len+1);
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      strcat (record, program_name);
    }

  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) )   /* we skip the possible lalDebugLevel-entry for now (FIXME?) */
    {
      if ( (ptr->state & UVAR_WAS_SET) == FALSE )	/* skip unset variables */
	continue;

      valstr = NULL;
      typestr = NULL;
      if ( (valstr = XLALUvarValue2String ( ptr )) == NULL ) {
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      if ( (typestr = XLALUvarType2String ( ptr )) == NULL ) {
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

      appendlen = strlen(ptr->name) + strlen(valstr) + strlen(typestr) + 10;
      if ( (append = XLALMalloc(appendlen)) == NULL) {
        XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, appendlen );
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      switch (format)
	{
	case UVAR_LOGFMT_CFGFILE:
	  sprintf (append, "%s = %s ;\n", ptr->name, valstr);
	  break;
	case UVAR_LOGFMT_CMDLINE:
	  sprintf (append, " --%s=%s", ptr->name, valstr);
	  break;
	case UVAR_LOGFMT_PROCPARAMS:
	  sprintf (append, "--%s = %s :%s;", ptr->name, valstr, typestr);
	  break;
	default:
          XLALPrintError ("%s: Unknown format for recording user-input: '%s'\n", __func__, format );
          XLAL_ERROR_NULL ( XLAL_EINVAL );
	  break;
	} /* switch (format) */

      len += strlen(append);
      if ( (record = LALRealloc (record, len+1)) == NULL ) {
        XLALPrintError("%s: LALRealloc (%d, %d) failed.\n", __func__, record, len+1);
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcat (record, append);

      XLALFree (valstr);
      valstr=NULL;
      XLALFree(typestr);
      typestr = NULL;
      XLALFree (append);
      append=NULL;
    } /* while ptr->next */

  return record;

} /* XLALUserVarGetLog() */


/* Return the type of the given UserVariable as a string.
 * For INTERNAL use only!
 */
CHAR *
XLALUvarType2String ( LALUserVariable *uvar )
{
  CHAR *ret=NULL;
  CHAR buf[16];

  if ( !uvar ) {
    XLALPrintError ("%s: invalid NULL input\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  switch (uvar->type)
    {
    case UVAR_BOOL:
      sprintf(buf, "boolean");
      break;
    case UVAR_INT4:
      sprintf(buf, "int4");
      break;
    case UVAR_REAL8:
      sprintf(buf, "real8");
      break;
    case UVAR_STRING:
      sprintf(buf, "string");
      break;
    case UVAR_CSVLIST:
      sprintf(buf, "list");
      break;
    default:
      XLALPrintError ("%s: ERROR: unkown UserVariable-type encountered\n", __func__ );
      XLAL_ERROR_NULL ( XLAL_EINVAL );
      break;
    } /* switch */

  if ( (ret = XLALMalloc (strlen(buf) + 1)) == NULL) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, strlen(buf)+1);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  strcpy (ret, buf);

  return ret;

} /* XLALUvarType2String() */



/* Return the value of the given UserVariable as a string.
 * For INTERNAL use only!
 */
CHAR *
XLALUvarValue2String ( LALUserVariable *uvar )
{
  CHAR *str = NULL;
  CHAR *ptr;
  CHAR buf[512];	/* buffer for producing non-string values */

  if ( !uvar || !uvar->varp ) {
    XLALPrintError ("%s: illegal NULL in input uvar struct\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  switch (uvar->type)	      /* info on default-value for this variable */
    {
    case UVAR_BOOL:
      sprintf (buf, *(BOOLEAN*)(uvar->varp) ? "TRUE" : "FALSE");
      break;
    case UVAR_INT4:
      sprintf (buf, "%" LAL_INT4_FORMAT, *(INT4*)(uvar->varp) );
      break;
    case UVAR_REAL8:
      if (*(REAL8*)(uvar->varp) == 0)
	strcpy (buf, "0.0");	/* makes it more explicit that's it a REAL */
      else
	sprintf (buf, "%.16g", *(REAL8*)(uvar->varp) );
      break;
    case UVAR_STRING:
      ptr = *(CHAR**)(uvar->varp);
      if ( ptr != NULL )
	{
	  if ( (str = XLALMalloc ( strlen(ptr) + 3 )) == NULL) {
            XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, strlen(ptr)+3);
            XLAL_ERROR_NULL ( XLAL_ENOMEM );
	  }
	  sprintf (str, "\"%s\"", ptr);
	}
      else
	strcpy (buf, "NULL");
      break;

    case UVAR_CSVLIST:
      {
	UINT4 i, listlen, outlen = 0;
	LALStringVector *csv;
	str = NULL;
	csv = *(LALStringVector**)(uvar->varp);
	if ( csv != NULL )
	  {
	    listlen = csv->length;
	    for ( i=0; i < listlen; i++)
	      {
		outlen += strlen (csv->data[i]) + 3;
		str = XLALRealloc(str, outlen);
		if ( i==0 )
		  sprintf ( str, "\"" );
		else
		  strcat ( str, "," );
		strcat ( str, csv->data[i] );
	      } /* for i < listlen */
	    strcat(str,"\"" );
	  } /* if csv != NULL */
	else
	  strcpy (buf, "NULL");
      }
      break;
    default:
      XLALPrintError ("%s: ERROR: unkown UserVariable-type encountered... this points to a coding error!\n", __func__ );
      XLAL_ERROR_NULL ( XLAL_EINVAL );
      break;

    } /* switch uvar->type */

  if (str == NULL)
    {
      if ( (str = XLALMalloc (strlen(buf) + 1)) == NULL) {
        XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, strlen(buf)+1);
        XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }
      strcpy (str, buf);
    }

  return ( str );

} /* XLALUvarValue2String() */


/**
 * Copy (and allocate) string 'in', possibly with quotes \" or \' removed.
 * If quotes are present at the beginning of 'in', they must have a matching
 * quote at the end of string, otherwise an error is printed and return=NULL
 */
CHAR *
XLAL_copy_string_unquoted ( const CHAR *in )
{
  const CHAR *tmp;
  CHAR *out;
  CHAR opening_quote = 0;
  CHAR closing_quote = 0;
  UINT4 outlen;

  XLAL_CHECK_NULL ( in != NULL, XLAL_EINVAL );
  UINT4 inlen = strlen ( in );

  if ( (in[0] == '\'') || (in[0] == '\"') ) {
    opening_quote = in[0];
  }
  if ( (inlen >= 2) && ( (in[inlen-1] == '\'') || (in[inlen-1] == '\"') ) ) {
    closing_quote = in[inlen-1];
  }

  /* check matching quotes */
  XLAL_CHECK_NULL ( opening_quote == closing_quote, XLAL_EINVAL, "Unmatched quotes in string [%s]\n", in );

  if ( opening_quote )
    {
      tmp = in + 1;
      outlen = inlen - 2;
    }
  else
    {
      tmp = in;
      outlen = inlen;
    }

  XLAL_CHECK_NULL ( (out = LALCalloc (1, outlen + 1)) != NULL, XLAL_ENOMEM );

  strncpy ( out, tmp, outlen);
  out[outlen] = 0;

  return out;

} /* XLAL_copy_string_unquoted() */

/**
 * Mark the user-variable as set, check if it has been
 * set previously and issue a warning if set more than once ...
 */
void
check_and_mark_as_set ( LALUserVariable *varp )
{
  /* check if this variable had been set before ... */
  if ( (varp->state & UVAR_WAS_SET) )
    LogPrintf ( LOG_NORMAL, "Warning: user-variable '%s' was set more than once!\n", varp->name ? varp->name : "(NULL)" );

  varp->state = (UserVarState)( varp->state |  UVAR_WAS_SET );

  return;
} /* check_and_mark_as_set() */



/* ========== DEPRECATED LAL INTERFACE FUNCTIONS, which have been replaced by XLAL functions,
 * These functions are just wrappers around the XLAL functions
 */


/** \deprecated use XLALRegisterUserVar() instead */
static void
RegisterUserVar (LALStatus *status,
		 const CHAR *name,
		 UserVarType type,
		 CHAR optchar,
		 UserVarState flag,
		 const CHAR *helpstr,
		 void *cvar)
{
  const char *fn = __func__;

  INITSTATUS(status);

  ASSERT (cvar != NULL, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (name != NULL, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);

  if ( XLALRegisterUserVar ( name, type, optchar, flag, helpstr, cvar ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Call to XLALRegisterUserVar() failed: %d\n", fn, xlalErrno );
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }

  RETURN (status);

} /* LALRegisterUserVar() */


/** \deprecated us XLALDestroyUserVars() instead */
void
LALDestroyUserVars (LALStatus *status)
{

  INITSTATUS(status);

  XLALDestroyUserVars();

  RETURN(status);

} /* LALDestroyUserVars() */


/** \deprecated use XLALUserVarReadCmdline() instead */
void
LALUserVarReadCmdline (LALStatus *status, int argc, char *argv[])
{
  const char *fn = __func__;
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (argv, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (UVAR_vars.next, status, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  if ( XLALUserVarReadCmdline(argc, argv) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: Call to XLALUserVarReadCmdline() failed with code %d\n", fn, xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALUserVarReadCmdline() */

/** \deprecated use XLALUserVarCheckRequired() instead */
void
LALUserVarCheckRequired (LALStatus *status)
{
  INITSTATUS(status);

  if ( XLALUserVarCheckRequired() != XLAL_SUCCESS ) {
    ABORT (status, USERINPUTH_ENOTSET, USERINPUTH_MSGENOTSET);
  }

  RETURN (status);

} /* LALUserVarCheckRequired() */


/** \deprecated use XLALUserVarReadAllInput() instead */
void
LALUserVarReadAllInput (LALStatus *status, int argc, char *argv[])
{
  const char *fn = __func__;
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (UVAR_vars.next, status, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  if ( XLALUserVarReadAllInput ( argc, argv ) != XLAL_SUCCESS ) {
    XLALPrintError ( "%s: XLALUserVarReadAllInput() failed with code %d\n", fn, xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

} /* LALReadUserInput() */

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
  const char *fn = __func__;

  INITSTATUS(status);

  ASSERT (logstr, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (*logstr == NULL, status, USERINPUTH_ENONULL,USERINPUTH_MSGENONULL);

  if ( ((*logstr) = XLALUserVarGetLog ( format )) == NULL ) {
    XLALPrintError ("%s: UserVarLogFormat() failed.\n", fn );
    ABORT (status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL);
  }

  RETURN (status);

} /* LALUserVarGetLog() */


#if 0
/**
 * Return user log as a process-params table
 *
 * \param[out] **procPar the output ProcessParamsTable
 * \param[in] *progname  name of calling code
 */
void
LALUserVarGetProcParamsTable (LALStatus *status, ProcessParamsTable **out, CHAR *progname)
{
  LALUserVariable *ptr = NULL;
  CHAR *valstr=NULL;
  CHAR *typestr=NULL;
  ProcessParamsTable *this_proc_param=NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT (out, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (*out == NULL, status, USERINPUTH_ENONULL,USERINPUTH_MSGENONULL);
  ASSERT (progname, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);

  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) )   /* we skip the possible lalDebugLevel-entry */
    {
      if ( (ptr->state & UVAR_WAS_SET) == FALSE )	/* skip unset variables */
	continue;

      /* get value and type of the uservar */
      TRY ( UvarValue2String (status->statusPtr, &valstr, ptr), status);
      TRY ( UvarType2String (status->statusPtr, &typestr, ptr), status);

      /* *out is null in the first iteration of this loop in which case
	 we allocate memory for the header of the linked list, otherwise
	 allocate memory for the nodes */
      if (*out == NULL)
	this_proc_param = *out = (ProcessParamsTable *)LALCalloc( 1, sizeof(ProcessParamsTable) );
      else
	this_proc_param = this_proc_param->next =
	  (ProcessParamsTable *)LALCalloc( 1, sizeof(ProcessParamsTable) );

      /* copy the strings into the procparams table */
      snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", progname );
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", ptr->name );
      snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s", valstr );
      snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", typestr );

      LALFree (valstr);
      valstr=NULL;
      LALFree(typestr);
      typestr = NULL;

    } /* while ptr->next */


  DETATCHSTATUSPTR(status);
  RETURN (status);

} /* LALUserVarGetProcParamsTable() */
#endif


/**
 * \deprecated use XLALUserVarHelpString() instead
 */
void
LALUserVarHelpString (LALStatus *status,
		      CHAR **helpstring, /* output: allocated here! */
		      const CHAR *progname)
{
  const char *fn = __func__;

  INITSTATUS(status);

  ASSERT (UVAR_vars.next, status, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);
  ASSERT (helpstring != NULL, status, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT ( *helpstring == NULL, status, USERINPUTH_ENONULL, USERINPUTH_MSGENONULL);

  if ( ((*helpstring) = XLALUserVarHelpString ( progname )) == NULL ) {
    XLALPrintError ("%s: XLALUserVarHelpString() failed with code %d\n", fn, xlalErrno );
    ABORT ( status,  USERINPUTH_EXLAL,  USERINPUTH_MSGEXLAL );
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALUserVarHelpString() */

/** \deprecated use XLALRegisterREALUserVar() instead */
void
LALRegisterREALUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarState flag,
			const CHAR *helpstr,
			REAL8 *cvar)
{
  RegisterUserVar (status, name, UVAR_REAL8, optchar, flag, helpstr, cvar);
}

/** \deprecated use XLALRegisterINTUserVar() instead */
void
LALRegisterINTUserVar (LALStatus *status,
		       const CHAR *name,
		       CHAR optchar,
		       UserVarState flag,
		       const CHAR *helpstr,
		       INT4 *cvar)
{
  RegisterUserVar (status, name, UVAR_INT4, optchar, flag, helpstr, cvar);
}

/** \deprecated use XLALRegisterBOOLUserVar() instead */
void
LALRegisterBOOLUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarState flag,
			const CHAR *helpstr,
			BOOLEAN *cvar)
{
  RegisterUserVar (status, name, UVAR_BOOL, optchar, flag, helpstr, cvar);
}

/** \deprecated use XLALRegisterSTRINGUserVar() instead */
void
LALRegisterSTRINGUserVar (LALStatus *status,
			  const CHAR *name,
			  CHAR optchar,
			  UserVarState flag,
			  const CHAR *helpstr,
			  CHAR **cvar)
{
  RegisterUserVar (status, name, UVAR_STRING, optchar, flag, helpstr, cvar);
}

/** \deprecated use XLALRegisterLISTUserVar() instead */
void
LALRegisterLISTUserVar (LALStatus *status,
			const CHAR *name,
			CHAR optchar,
			UserVarState flag,
			const CHAR *helpstr,
			LALStringVector **cvar)
{
  RegisterUserVar ( status, name, UVAR_CSVLIST, optchar, flag, helpstr, cvar );
}

/** \deprecated use XLALUserVarReadCfgfile() instead */
void
LALUserVarReadCfgfile (LALStatus *status,
		       const CHAR *cfgfile) 	   /* name of config-file */
{

  INITSTATUS(status);

  ASSERT (UVAR_vars.next, status, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  if ( XLALUserVarReadCfgfile ( cfgfile ) != XLAL_SUCCESS ) {
    ABORT ( status, USERINPUTH_EXLAL, USERINPUTH_MSGEXLAL );
  }

  RETURN (status);

} /* LALUserVarReadCfgfile() */
