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

/* --------------------
 * This structure defines a "user-variable", which can be read
 * automagically from command-line and/or config-file. 
 *     ** USED ONLY INTERNALLY !! ** 
 * --------------------     */
typedef struct tagLALUserVariable {
  const CHAR *name;	/* full name */
  LALUserVarType type;	/* type: bool, int, float or string */
  CHAR optchar;		/* cmd-line character */
  const CHAR *help;	/* help-string */
  void *varp;		/* pointer to the actual C-variable */
  struct tagLALUserVariable *next; /* linked list */
} LALUserVariable;

/* this is the module-local linked list to store the user-variable */
static LALUserVariable UVAR_vars;	/* empty head */

/* needed for command-line parsing */
extern char *optarg;
extern int optind, opterr, optopt;

/*----------------------------------------------------------------------
 * "register" a user-variable with the module
 * effectively put an appropriate entry into UVAR_vars
 *----------------------------------------------------------------------*/
void
LALRegisterUserVar (LALStatus *stat, 
		    const CHAR *name, 
		    LALUserVarType type, 
		    CHAR optchar, 
		    const CHAR *helpstr, 
		    void *cvar)
{
  LALUserVariable *ptr;

  INITSTATUS( stat, "LALRegisterUserVar", USERINPUTC );

  ASSERT (cvar != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);

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
	ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
      }
      if (c != 0) 	/* find short-option character */
	{
	  ptr = &UVAR_vars;
	  while ( (ptr=ptr->next) != NULL)
	    if (c == ptr->optchar)
	      break;
	} /* if short-option */
      else	/* find long-option: returned in longindex */
	{
	  ptr = &UVAR_vars;
	  while ( (ptr=ptr->next) != NULL)
	    if ( !strcmp (long_options[longindex].name, ptr->name) )
	      break;
	}
      if (ptr == NULL) {	/* should not be possible: nothing found at all... */
	LALPrintError ("WARNING: failed to find option.. this points to a coding-error!\n");
	ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
      }

      switch (ptr->type)
	{
	  BOOLEAN ret;
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
		ABORT (stat, USERINPUTH_EBOOL, USERINPUTH_MSGEBOOL);
	      }
	    } /* parse bool-argument */

	  *(BOOLEAN*)(ptr->varp)  = ret;
	  break;

	case UVAR_INT4:
	  *(INT4*)(ptr->varp) = (INT4) atoi (optarg);
	  break;

	case UVAR_REAL8:
	  *(REAL8*)(ptr->varp) = (REAL8) atof (optarg);
	  break;

	case UVAR_STRING:
	  if (!optarg) {	/* should not be possible, but let's be paranoid */
	    ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	  }
	  if ( *(CHAR**)(ptr->varp) != NULL) {	 /* something allocated here before? */
	    LALFree ( *(CHAR**)(ptr->varp) );
	  }

	  *(CHAR**)(ptr->varp) = LALMalloc (strlen(optarg) + 1);
	  if ( *(CHAR**)(ptr->varp) == NULL) {
	    ABORT (stat, USERINPUTH_EMEM, USERINPUTH_MSGEMEM);
	  }
	  strcpy ( *(CHAR**)(ptr->varp), optarg);
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
  LALConfigData *cfg = NULL;
  CHAR *stringbuf;
  LALUserVariable *ptr;

  INITSTATUS( stat, "LALUserVarReadCfgfile", USERINPUTC );
  ATTATCHSTATUSPTR (stat);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  TRY (LALLoadConfigFile (stat->statusPtr, &cfg, cfgfile), stat);

  /* step through all user-variable: read those with names from config-file */
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      if (ptr->name == NULL)	/* ignore name-less user-variable */
	continue;

      switch (ptr->type)
	{
	case UVAR_BOOL:
	  TRY( LALReadConfigBOOLVariable (stat->statusPtr, ptr->varp, cfg, ptr->name ), stat);
	  break;
	case UVAR_INT4:
	  TRY( LALReadConfigINT4Variable (stat->statusPtr, ptr->varp, cfg, ptr->name ), stat);
	  break;
	case UVAR_REAL8:
	  TRY( LALReadConfigREAL8Variable (stat->statusPtr, ptr->varp, cfg, ptr->name ), stat);
	  break;
	case UVAR_STRING:
	  stringbuf = NULL;
	  TRY( LALReadConfigSTRINGVariable (stat->statusPtr, &stringbuf, cfg, ptr->name ), stat);
	  if (stringbuf)	/* did we find something? */
	    {
	      if ( *(CHAR**)(ptr->varp) != NULL)	 /* something allocated here before? */
		LALFree ( *(CHAR**)(ptr->varp) );

	      *(CHAR**)(ptr->varp) = stringbuf;
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

  TRY( LALDestroyConfigData (stat->statusPtr, &cfg), stat);

  DETATCHSTATUSPTR(stat);
  RETURN (stat);

} /* LALUserVarReadCfgfile() */


/*----------------------------------------------------------------------
 * assemble all help-info from uvars into a help-string
 *----------------------------------------------------------------------*/
#define UVAR_MAXHELPLINE  512	/* max length of one help-line */
void
LALUserVarHelpString (LALStatus *stat, 
		      CHAR **helpstring) /* output: allocated here! */
{
  CHAR strbuf[UVAR_MAXHELPLINE];	/* should be enough for one line...*/
  CHAR defaultstr[100]; /* for display of default-value */
  /* we need strings for UVAR_BOOL, UVAR_INT4, UVAR_REAL8, UVAR_STRING: */
  const CHAR *typestr[] = {"BOOL", "INT", "FLOAT", "STRING"}; 
  LALUserVariable *ptr;
  CHAR *helpstr = NULL;
  size_t newlen = 0;

  INITSTATUS (stat, "LALUserVarHelpString", USERINPUTC);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);
  ASSERT (helpstring != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT ( *helpstring == NULL, stat, USERINPUTH_ENONULL, USERINPUTH_MSGENONULL);

  /* put together the help-string. Allocate memory on-the-fly... */
  ptr = &UVAR_vars;
  while ( (ptr=ptr->next) != NULL)
    {
      switch (ptr->type)
	{
	case UVAR_BOOL:
	  sprintf (defaultstr, *(BOOLEAN*)(ptr->varp) ? "True" : "False");
	  break;
	case UVAR_INT4:
	  sprintf (defaultstr, "%d", *(INT4*)(ptr->varp) );
	  break;
	case UVAR_REAL8:
	  sprintf (defaultstr, "%.3g", *(REAL8*)(ptr->varp) );
	  break;
	case UVAR_STRING:
	  if ( *(CHAR**)(ptr->varp) )
	    strcpy (defaultstr, *(CHAR**)(ptr->varp) );
	  else
	    strcpy (defaultstr, "\"\"");
	  break;

	default:
	  LALPrintError ("ERROR: unkown UserVariable-type encountered... this points to a coding error!\n");
	  ABORT (stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
	  break;

	} /* switch ptr->type */

      LALSnprintf (strbuf, UVAR_MAXHELPLINE,  "   --%-14s (-%c) \t%s\t%s  (Default: %s)\n", 
	       ptr->name ? ptr->name : "-NONE-", 
	       ptr->optchar, 
	       typestr[ptr->type], 
	       ptr->help ? ptr->help : "-NONE-",
	       defaultstr);

      /* now increase allocated memory by the right amount */
      newlen += strlen (strbuf) + 1;
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
void
LALUserVarReadAllInput (LALStatus *stat, int argc, char *argv[])
{
  INT4 i;
  CHAR* fname = NULL;

  INITSTATUS( stat, "LALUserVarReadAllInput", USERINPUTC);

  ASSERT (UVAR_vars.next, stat, USERINPUTH_ENOUVARS,  USERINPUTH_MSGENOUVARS);

  /* FIXME: debug-level should be parsed here, before ATTATCHSTATUSPTR(), because */
  /* !! this can only be done after lalDebugLevel has been set *definitively* !! */
  ATTATCHSTATUSPTR (stat); 

  /*----------------------------------------------------------------------
   * pre-process command-line: we want config-file (and debug-level, FIXME)
   */
  for (i=1; i < argc; i++)
    {
      /* was a config-file specified ? */
      if ( argv[i][0] == '@' )
	{
	  fname = LALMalloc ( strlen(argv[i]+1) + 1 );
	  if (fname == NULL) {
	    ABORT (stat, USERINPUTH_EMEM,  USERINPUTH_MSGEMEM);
	  }
	  strcpy (fname, argv[i]+1);
	}

    } /* for i < argc */
  /*----------------------------------------------------------------------*/

  /* if config-file specified, read from that first */
  if (fname) {
    TRY (LALUserVarReadCfgfile (stat->statusPtr, fname), stat);
    LALFree (fname);
  }

  /* now do proper cmdline parsing: overloads config-file settings */
  TRY (LALUserVarReadCmdline (stat->statusPtr, argc, argv), stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
  
} /* LALReadUserInput() */

