/************************************ <lalVerbatim file="UserInputCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UserInput.c}}
\label{ss:UserInput.c}

Convenient unified handling of user-input via config-file and/or command-line.

\subsubsection*{Prototypes}

\input{UserInputCP}

\subsubsection*{Description}

This module provides a very simple and convenient way to handle
user-input, wether it comes from the command-line or a config-file. 

The general procedure is the following:

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UserInputCV}}

******************************************************* </lalLaTeX> */
#include "getopt.h"

#include <lal/LALStdio.h>
#include <lal/UserInput.h>

NRCSID( USERINPUTC, "$Id$");

extern INT4 lalDebugLevel;

#define TRUE  (1==1)
#define FALSE (1==0)

/* --------------------
 * This structure defines a "user-variable", which can be read
 * automagically from command-line and/or config-file. 
 *     ** USED ONLY INTERNALLY !! ** 
 * --------------------     */
typedef enum {
  UVAR_BOOL,	/* boolean */
  UVAR_INT4,	/* integer */
  UVAR_REAL8,	/* float */
  UVAR_STRING	/* string */
} UserVarType;

typedef struct tagLALUserVariable {
  const CHAR *name;	/* full name */
  UserVarType type;	/* type: bool, int, float or string */
  CHAR optchar;		/* cmd-line character */
  const CHAR *help;	/* help-string */
  void *varp;		/* pointer to the actual C-variable */
  UserVarState state;	/* state (empty, default, set) */
  struct tagLALUserVariable *next; /* linked list */
} LALUserVariable;

/* this is the module-local linked list to store the user-variable */
static LALUserVariable UVAR_vars;	/* empty head */

/* needed for command-line parsing */
extern char *optarg;
extern int optind, opterr, optopt;

/* internal prototypes */
static void
RegisterUserVar (LALStatus *stat, 
		 const CHAR *name, 
		 UserVarType type, 
		 CHAR optchar, 
		 UserVarState flag,
		 const CHAR *helpstr, 
		 void *cvar);


/* these are type-specific wrappers to allow tighter type-checking! */
/* <lalVerbatim file="UserInputCP"> */
void
LALRegisterREALUserVar (LALStatus *stat, 
			const CHAR *name, 
			CHAR optchar, 
			UserVarState flag,
			const CHAR *helpstr, 
			REAL8 *cvar)
{ /* </lalVerbatim> */
  RegisterUserVar (stat, name, UVAR_REAL8, optchar, flag, helpstr, cvar);
}
/* <lalVerbatim file="UserInputCP"> */
void
LALRegisterINTUserVar (LALStatus *stat, 
		       const CHAR *name, 
		       CHAR optchar, 
		       UserVarState flag,
		       const CHAR *helpstr, 
		       INT4 *cvar)
{ /* </lalVerbatim> */
  RegisterUserVar (stat, name, UVAR_INT4, optchar, flag, helpstr, cvar);
} 

/* <lalVerbatim file="UserInputCP"> */
void
LALRegisterBOOLUserVar (LALStatus *stat, 
			const CHAR *name, 
			CHAR optchar, 
			UserVarState flag,
			const CHAR *helpstr, 
			BOOLEAN *cvar)
{ /* </lalVerbatim> */
  RegisterUserVar (stat, name, UVAR_BOOL, optchar, flag, helpstr, cvar);
} 

/* <lalVerbatim file="UserInputCP"> */
void
LALRegisterSTRINGUserVar (LALStatus *stat, 
			  const CHAR *name, 
			  CHAR optchar, 
			  UserVarState flag,
			  const CHAR *helpstr, 
			  CHAR **cvar)
{ /* </lalVerbatim> */
  RegisterUserVar (stat, name, UVAR_STRING, optchar, flag, helpstr, cvar);
} 




/*----------------------------------------------------------------------
 * "register" a user-variable with the module
 * effectively put an appropriate entry into UVAR_vars
 *
 * NOTE: don't use this directly, as it's not type-safe!!
 *      ==> use one of the 4 wrappers above!
 *
 *----------------------------------------------------------------------*/
static void
RegisterUserVar (LALStatus *stat, 
		 const CHAR *name, 
		 UserVarType type, 
		 CHAR optchar, 
		 UserVarState flag,
		 const CHAR *helpstr,
		 void *cvar)
{
  LALUserVariable *ptr;

  INITSTATUS( stat, "LALRegisterUserVar", USERINPUTC );

  ASSERT (cvar != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (name != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);

  /* find end of uvar-list */
  ptr = &UVAR_vars;
  while (ptr->next) 
    ptr = ptr->next;

  /* create new entry */
  ptr->next = LALCalloc (1, sizeof(LALUserVariable));
  if (ptr->next == NULL) {
    ABORT (stat, USERINPUTH_EMEM,  USERINPUTH_MSGEMEM);
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

  RETURN (stat);

} /* LALRegisterUserVar() */



/*--------------------------------------------------------------------------------
 * free memory associated with user-variable linked list
 *--------------------------------------------------------------------------------*/
/* <lalVerbatim file="UserInputCP"> */
void
LALDestroyUserVars (LALStatus *stat)
{ /* </lalVerbatim> */

  LALUserVariable *ptr, *lastptr;

  INITSTATUS( stat, "LALDestroyUserVars", USERINPUTC );

  /* step through user-variables: free list-entries and all allocated strings */
  ptr = &(UVAR_vars);
  lastptr = NULL;
  while ( (ptr=ptr->next) != NULL)
    {
      /* is an allocated string here? */
      if ( (ptr->type == UVAR_STRING) && (*(CHAR**)(ptr->varp) != NULL) ) 
	LALFree ( *(CHAR**)(ptr->varp) );

      /* free list-entry behind us (except for the head) */
      if (lastptr)
	LALFree (lastptr);

      lastptr = ptr;

    } /* while ptr->next */

  if (lastptr)
    LALFree (lastptr);

  /* clean head */
  memset (&UVAR_vars, 0, sizeof(UVAR_vars));
  
  RETURN(stat);
} /* LALDestroyUserVars() */


/*----------------------------------------------------------------------
 * parse command-line into UserVariable array
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="UserInputCP"> */
void
LALUserVarReadCmdline (LALStatus *stat, 
		     int argc, char *argv[]) 	  /* command-line contents */
{ /* </lalVerbatim> */
  INT4 c;
  UINT4 pos;
  UINT4 numvars;
  LALUserVariable *ptr = NULL;
  char optstring[512] = "\0";	/* string of short-options, should be easily enough */
  struct option *long_options;
  int longindex = -1;

  INITSTATUS( stat, "LALUserVarReadCmdline", USERINPUTC );

  ASSERT (argv, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  /* build optstring of short-options */
  ptr = &UVAR_vars;	/* set to empty head */
  pos = 0;
  numvars = 0;

  /* special treatment of head, which could contain the debug-option: 
   * NOTE: this one will *NOT* be treated by the remaining function, 
   * (as the head is always skipped), but it has to be in the optstring 
   * to avoid an error if specified on the command-line. 
   * Treatment of debug-option reading has to be done separately using
   * UVARgetDebugLevle() 
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

  /* parse the command-line */
  while ( (c = getopt_long(argc, argv, optstring, long_options, &longindex)) != -1 )
    {
      if (c == '?') {
	CHAR *helpstring = NULL;
	ATTATCHSTATUSPTR (stat);
	TRY (LALUserVarHelpString (stat->statusPtr, &helpstring, argv[0]), stat);
	printf ("\n%s\n", helpstring);
	LALFree (helpstring);
	ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
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
	LALPrintError ("WARNING: failed to find option.. this points to a coding-error!\n");
	ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
      }

      /* if we found the debug-switch, ignore it (has been handled already */
      if (ptr == &UVAR_vars)
	continue;

      switch (ptr->type)
	{
	  INT2 ret;
	case UVAR_BOOL:
	  ret = -1;

	  if (optarg == NULL)	/* no argument: counts a 'true' */
	    ret = 1;
	  else		/* parse bool-argument: should be consistent with bool-parsing in ConfigFile! */
	    {
	      if      ( !strcmp(optarg, "yes") || !strcmp(optarg, "true") || !strcmp(optarg,"1") )
		ret = 1;
	      else if ( !strcmp (optarg, "no") || !strcmp(optarg,"false") || !strcmp(optarg,"0") )
		ret = 0;
	      else {	/* failed to parse BOOL properly */
		LALPrintError ("illegal bool-value `%s`\n", optarg);
		ABORT (stat, USERINPUTH_EBOOL, USERINPUTH_MSGEBOOL);
	      }
	    } /* parse bool-argument */
	  
	  /* only set if we properly parsed something */
	  if (ret != -1) {
	    *(BOOLEAN*)(ptr->varp)  = (BOOLEAN)ret;
	    ptr->state |= UVAR_WAS_SET;
	  }

	  break;

	case UVAR_INT4:
	  *(INT4*)(ptr->varp) = (INT4) atoi (optarg);
	  ptr->state |= UVAR_WAS_SET;
	  break;

	case UVAR_REAL8:
	  *(REAL8*)(ptr->varp) = (REAL8) atof (optarg);
	  ptr->state |= UVAR_WAS_SET;
	  break;

	case UVAR_STRING:
	  if (!optarg) {	/* should not be possible, but let's be paranoid */
	    ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	  }
	  if ( *(CHAR**)(ptr->varp) != NULL) {	 /* something allocated here before? */
	    LALFree ( *(CHAR**)(ptr->varp) );
	  }

	  *(CHAR**)(ptr->varp) = LALCalloc (1, strlen(optarg) + 1);
	  if ( *(CHAR**)(ptr->varp) == NULL) {
	    ABORT (stat, USERINPUTH_EMEM, USERINPUTH_MSGEMEM);
	  }
	  strcpy ( *(CHAR**)(ptr->varp), optarg);
	  ptr->state |= UVAR_WAS_SET;
	  break; 

	default:
	  LALPrintError ("ERROR: unkown UserVariable-type encountered... this points to a coding error!\n");
	  ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	  break;

	} /* switch ptr->type */

    } /* while getopt_long() */

  LALFree (long_options);

  RETURN (stat);

} /* LALUserVarReadCmdline() */

/*----------------------------------------------------------------------
 * Read config-variables from cfgfile and parse into input-structure
 *
 * an error is reported if the config-file reading fails, but the 
 * individual variable-reads are treated as optional
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="UserInputCP"> */
void
LALUserVarReadCfgfile (LALStatus *stat, 
		       const CHAR *cfgfile) 	   /* name of config-file */
{/* </lalVerbatim> */
  LALParsedDataFile *cfg = NULL;
  CHAR *stringbuf;
  LALUserVariable *ptr;
  BOOLEAN wasRead;

  INITSTATUS( stat, "LALUserVarReadCfgfile", USERINPUTC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  TRY (LALParseDataFile (stat->statusPtr, &cfg, cfgfile), stat);

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
	  TRY( LALReadConfigBOOLVariable (stat->statusPtr, ptr->varp, cfg, ptr->name, &wasRead), stat);
	  if (wasRead)
	    ptr->state |= UVAR_WAS_SET;
	  break;
	case UVAR_INT4:
	  TRY( LALReadConfigINT4Variable (stat->statusPtr, ptr->varp, cfg, ptr->name, &wasRead), stat);
	  if (wasRead)
	    ptr->state |= UVAR_WAS_SET;
	  break;
	case UVAR_REAL8:
	  TRY( LALReadConfigREAL8Variable (stat->statusPtr, ptr->varp, cfg, ptr->name, &wasRead), stat);
	  if (wasRead)
	    ptr->state |= UVAR_WAS_SET;
	  break;
	case UVAR_STRING:
	  stringbuf = NULL;
	  TRY( LALReadConfigSTRINGVariable (stat->statusPtr, &stringbuf, cfg, ptr->name, &wasRead), stat);
	  if ( wasRead && stringbuf)	/* did we find something? */
	    {
	      if ( *(CHAR**)(ptr->varp) != NULL)	 /* something allocated here before? */
		LALFree ( *(CHAR**)(ptr->varp) );

	      *(CHAR**)(ptr->varp) = stringbuf;
	      ptr->state |= UVAR_WAS_SET;
	    } /* if stringbuf */
	  break;

	default:
	  LALPrintError ("ERROR: unkown UserVariable-type encountered... this points to a coding error!\n");
	  ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	  break;

	} /* switch ptr->type */

    } /* while ptr->next */

  /* ok, that should be it: check if there were more definitions we did not read */
  TRY (LALCheckConfigReadComplete (stat->statusPtr, cfg, CONFIGFILE_ERROR), stat);	/* be strict */

  TRY( LALDestroyParsedDataFile (stat->statusPtr, &cfg), stat);

  DETATCHSTATUSPTR(stat);
  RETURN (stat);

} /* LALUserVarReadCfgfile() */


/*----------------------------------------------------------------------
 * assemble all help-info from uvars into a help-string
 *----------------------------------------------------------------------*/
#define UVAR_MAXHELPLINE  512	/* max length of one help-line */
/* <lalVerbatim file="UserInputCP"> */
void
LALUserVarHelpString (LALStatus *stat, 
		      CHAR **helpstring, /* output: allocated here! */
		      const CHAR *progname)
{ /* </lalVerbatim> */

  CHAR strbuf[UVAR_MAXHELPLINE];	/* should be enough for one line...*/
  CHAR defaultstr[100]; /* for display of default-value */
  CHAR optstr[10];	/* display of opt-char */
  /* we need strings for UVAR_BOOL, UVAR_INT4, UVAR_REAL8, UVAR_STRING: */
  const CHAR *typestr[] = {"BOOL", "INT", "REAL", "STRING"}; 
  LALUserVariable *ptr;
  CHAR *helpstr = NULL;
  size_t newlen = 0;

  INITSTATUS (stat, "LALUserVarHelpString", USERINPUTC);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);
  ASSERT (helpstring != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT ( *helpstring == NULL, stat, USERINPUTH_ENONULL, USERINPUTH_MSGENONULL);

  /* prepare first lines of help-string: info about config-file reading */
  newlen = 0;
  sprintf (strbuf, "Usage: %s [@ConfigFile] [options], where options are:\n\n", progname);
  if ( (helpstr = LALCalloc (1, strlen(strbuf) + 1)) == NULL) {
    ABORT (stat, USERINPUTH_EMEM, USERINPUTH_MSGEMEM);
  }
  strcpy (helpstr, strbuf);
  newlen += strlen (strbuf) + 1;

  /* put together the help-string. Allocate memory on-the-fly... */
  ptr = &UVAR_vars;

  /* special treatment of debug-option in the head (if present) */
  if ( ptr->help && ptr->optchar )
    {
      sprintf (strbuf, "                   -%c     %-8s  %s    [%d] \n", 
	       ptr->optchar, typestr[ptr->type], ptr->help, *(INT4*)(ptr->varp) );
      newlen += strlen (strbuf);
      helpstr = LALRealloc (helpstr, newlen);
      if ( helpstr == NULL) {
	ABORT (stat, USERINPUTH_EMEM, USERINPUTH_MSGEMEM);
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */
    }

  /* treat the remaining "normal" entries */
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->state & UVAR_REQUIRED)
	strcpy (defaultstr, "REQUIRED");
      else if (ptr->state & UVAR_HELP)
	strcpy (defaultstr, "");
      else
	switch (ptr->type)	      /* info on default-value for this variable */
	  {
	  case UVAR_BOOL:
	    sprintf (defaultstr, *(BOOLEAN*)(ptr->varp) ? "True" : "False");
	    break;
	  case UVAR_INT4:
	    sprintf (defaultstr, "%d", *(INT4*)(ptr->varp) );
	    break;
	  case UVAR_REAL8:
	    if (*(REAL8*)(ptr->varp) == 0)
	      strcpy (defaultstr, "0.0");
	    else
	      sprintf (defaultstr, "%g", *(REAL8*)(ptr->varp) );
	    break;
	  case UVAR_STRING:
	    if ( *(CHAR**)(ptr->varp) )
	      sprintf (defaultstr, "\"%s\"", *(CHAR**)(ptr->varp) );
	    else
	      strcpy (defaultstr, "NULL");
	    break;
	    
	  default:
	    LALPrintError ("ERROR: unkown UserVariable-type encountered... this points to a coding error!\n");
	    ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	    break;

	  } /* switch ptr->type */
      
      if (ptr->optchar != 0)
	sprintf (optstr, "(-%c)", ptr->optchar);
      else
	strcpy (optstr, "");

      LALSnprintf (strbuf, UVAR_MAXHELPLINE,  " --%-14s %-5s   %-8s  %s   [%s] \n", 
		   ptr->name ? ptr->name : "-NONE-", 
		   optstr,
		   typestr[ptr->type], 
		   ptr->help ? ptr->help : "-NONE-",
		   defaultstr);

      /* now increase allocated memory by the right amount */
      newlen += strlen (strbuf);
      helpstr = LALRealloc (helpstr, newlen);
      if ( helpstr == NULL) {
	ABORT (stat, USERINPUTH_EMEM, USERINPUTH_MSGEMEM);
      }

      strcat (helpstr, strbuf);	/* add this line to the helpstring */

    } /* while ptr->next */

  *helpstring = helpstr;

  RETURN(stat);

} /* LALUserVarHelpString() */


/*----------------------------------------------------------------------
 * This function puts all the above pieces together, and basically does
 * everything: get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="UserInputCP"> */
void
LALUserVarReadAllInput (LALStatus *stat, int argc, char *argv[])
{ /* </lalVerbatim> */

  INT4 i;
  CHAR* fname = NULL;
  LALUserVariable *ptr;

  INITSTATUS( stat, "LALUserVarReadAllInput", USERINPUTC);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  ATTATCHSTATUSPTR (stat); 

  /* pre-process command-line: have we got a config-file ? */
  for (i=1; i < argc; i++)
    {
      if ( argv[i][0] == '@' )
	{
	  fname = LALCalloc (1, strlen(argv[i]+1) + 1 );
	  if (fname == NULL) {
	    ABORT (stat, USERINPUTH_EMEM,  USERINPUTH_MSGEMEM);
	  }
	  strcpy (fname, argv[i]+1);
	}
    } /* for i < argc */


  /* if config-file specified, read from that first */
  if (fname) {
    TRY (LALUserVarReadCfgfile (stat->statusPtr, fname), stat);
    LALFree (fname);
  }

  /* now do proper cmdline parsing: overloads config-file settings */
  TRY (LALUserVarReadCmdline (stat->statusPtr, argc, argv), stat);

  /* now check if help-string was requested */
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    if ( (ptr->state & UVAR_HELP) && (ptr->state & UVAR_WAS_SET) )
      {
	CHAR *helpstring = NULL;
	TRY (LALUserVarHelpString (stat->statusPtr, &helpstring, argv[0]), stat);
	printf ("\n%s\n", helpstring);
	LALFree (helpstring);
	DETATCHSTATUSPTR (stat);
	RETURN (stat);
      } /* if help requested */

  /* check that all required input-variables have been specified */
  TRY (LALUserVarCheckRequired (stat->statusPtr), stat);


  DETATCHSTATUSPTR (stat);
  RETURN (stat);
  
} /* LALReadUserInput() */


/*----------------------------------------------------------------------
 * Has this user-variable been set by the user?
 * return -1 on error, TRUE/FALSE otherwise
 *----------------------------------------------------------------------*/
/* <lalVerbatim file="UserInputCP"> */
INT4
LALUserVarWasSet (void *cvar)
{ /* </lalVerbatim> */

  LALUserVariable *ptr;

  if (!cvar)
    return (-1);


  /* find this varname in the list of user-variables */
  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL) 
    if ( ptr->varp == cvar)
      break;

  if (ptr == NULL) {
    LALPrintError ("Variable passed to UVARwasSet is not a registered Usser-variable\n");
    return (-1);
  }
  
  /* we found it: has it been set by user? */
  return (ptr->state & UVAR_WAS_SET);

} /* LALUserVarWasSet() */

/*----------------------------------------------------------------------
 * check that all required user-variables have been set successfully
 * print error if not 
 *----------------------------------------------------------------------*/
void
LALUserVarCheckRequired (LALStatus *stat)
{
  LALUserVariable *ptr;  

  INITSTATUS( stat, "LALUserVarCheckRequired", USERINPUTC);

  /* go through list of uvars */
  ptr = &UVAR_vars;
  while ( (ptr = ptr->next) != NULL) 
    if ( (ptr->state & UVAR_REQUIRED) && !(ptr->state & UVAR_WAS_SET))
      {
	LALPrintError ("\nRequired user-variable `%s` has not been specified!\n\n", ptr->name);
	ABORT (stat, USERINPUTH_ENOTSET, USERINPUTH_MSGENOTSET);
      }

  RETURN (stat);

} /* LALUserVarCheckRequired() */

/*----------------------------------------------------------------------
 * treat the delicate setting of lalDebuglevel
 *----------------------------------------------------------------------*/
void
LALGetDebugLevel (LALStatus *stat, int argc, char *argv[], CHAR optchar)
{
  static const char *help = "set lalDebugLevel";
  static INT4 defaultDebugLevel;
  INT4 i;
  CHAR *ptr;

  INITSTATUS( stat, "UVARgetDebugLevel", USERINPUTC);

  ASSERT (argv, stat,  USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (UVAR_vars.next == NULL, stat, USERINPUTH_EDEBUG,  USERINPUTH_MSGEDEBUG);
  ASSERT (UVAR_vars.varp == NULL, stat, USERINPUTH_EDEBUG,  USERINPUTH_MSGEDEBUG);

  /* "register" the debug-level variable in the head of the UVAR-list, 
   * to avoid any mallocs. We need this to show up in the help-string */
  UVAR_vars.name = NULL;
  UVAR_vars.type = UVAR_INT4;
  UVAR_vars.optchar = optchar;
  UVAR_vars.help = help;

  defaultDebugLevel = lalDebugLevel;
  UVAR_vars.varp = &defaultDebugLevel;	/* trick: use to store default-value (for help-string) */

  UVAR_vars.state = UVAR_OPTIONAL;
  UVAR_vars.next = NULL;
  
  /* the command-line has to be processed by hand for this... ! */
  for (i=1; i < argc; i++)
    {
      if ( (argv[i][0] == '-') && (argv[i][1] == optchar) )
	{
	  if (argv[i][2] != '\0')
	    ptr = argv[i]+2;
	  else
	    ptr = argv[i+1];

	  if ( (ptr==NULL) || (sscanf ( ptr, "%d", &lalDebugLevel) != 1) ) {
	    LALPrintError ("setting debug-level `-%c` requires an argument\n", optchar);
	    ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
	  }
	  break;
	} /* if debug-switch found */
      
    } /* for i < argc */

  RETURN (stat);

} /* LALGetDebugLevel() */
