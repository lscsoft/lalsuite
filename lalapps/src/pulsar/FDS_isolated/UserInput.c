/************************************ <lalVerbatim file="UserInputCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{UserInput.c}}
\label{ss:UserInput.c}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UserInputCV}}

******************************************************* </lalLaTeX> */
#include <getopt.h>

#include "UserInput.h"

NRCSID( USERINPUTC, "$Id$");

extern INT4 lalDebugLevel;

/* needed for command-line parsing */
extern char *optarg;
extern int optind, opterr, optopt;


/*----------------------------------------------------------------------
 * parse command-line into UserVariable array
 *----------------------------------------------------------------------*/
void
LALReadCmdLineInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars)
{
  INT4 c;
  UINT4 index, numvars;
  UINT4 i, pos;
  char optstring[512] = "\0";	/* should be easily enough */
  struct option *long_options;

  INITSTATUS( stat, "ReadCmdlineInput", USERINPUTC );

  ASSERT (argv, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT (uvars, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);

  /* count number of user-variables */
  numvars = 0;
  while (uvars[numvars].name != NULL)
    numvars ++;

  /* build optstring */
  for (i=0; i< numvars; i++)
    {
      pos = strlen(optstring);
      optstring[pos++] = uvars[i].optchar;
      if (uvars[i].type != UVAR_BOOL)	/* everything but bool takes an argument */
	optstring[pos++] = ':';
      optstring[pos] = '\0';
    } /* for i < numvars */

  long_options = LALMalloc ( (numvars+1) * sizeof(struct option));
  
  for (i=0; i < numvars; i++)
    {
      long_options[i].name = uvars[i].name;
      long_options[i].has_arg = (uvars[i].type == UVAR_BOOL) ? no_argument : required_argument;
      long_options[i].flag = 0;
      long_options[i].val = uvars[i].optchar;
    }
  long_options[numvars].name = 0;
  long_options[numvars].has_arg = 0;
  long_options[numvars].flag = 0;
  long_options[numvars].val = 0;

  while ( (c = getopt_long(argc, argv, optstring, long_options, NULL)) != -1 )
    {
      index = 0;
      while ( long_options[index].name && ( c != long_options[index].val) )
	index ++;
      if (index >= numvars) {
	ABORT (stat, USERINPUTH_EOPT, USERINPUTH_MSGEOPT);
      }

      switch (uvars[index].type)
	{
	case UVAR_BOOL:
	  *(BOOLEAN*)(uvars[index].varp)  = 1;
	  break;
	case UVAR_INT4:
	  *(INT4*)(uvars[index].varp) = (INT4) atoi (optarg);
	  break;
	case UVAR_REAL8:
	  *(REAL8*)(uvars[index].varp) = (REAL8) atof (optarg);
	  break;
	case UVAR_CHAR:
	  if (!optarg) break;
	  if ( *(CHAR**)uvars[index].varp != NULL)	 /* something allocated here before? */
	    LALFree ( *(CHAR**)uvars[index].varp );

	  *(CHAR**)(uvars[index].varp) = LALMalloc (strlen(optarg) + 1);
	  strcpy ( *(CHAR**)uvars[index].varp, optarg);
	  break;
	} /* switch index */

    } /* while getopt() */

  LALFree (long_options);

  RETURN (stat);

} /* LALReadCmdLineInput() */

/*----------------------------------------------------------------------
 * Read config-variables from cfgfile and parse into input-structure
 *
 * an error is reported if the config-file reading fails, but the 
 * individual variable-reads are treated as optional
 *----------------------------------------------------------------------*/
void
LALReadCfgFileInput (LALStatus *stat, const CHAR *cfgfile, UserVariable *uvars)
{
  LALConfigData *cfg = NULL;
  UINT4 i, numvars;
  CHAR *stringbuf;
		  

  INITSTATUS( stat, "ReadCfgFile", USERINPUTC );
  ATTATCHSTATUSPTR (stat);

  TRY (LALLoadConfigFile (stat->statusPtr, &cfg, cfgfile), stat);


  /* count number of user-variables */
  numvars = 0;
  while (uvars[numvars].name != NULL)
    numvars ++;

  for (i=0; i< numvars; i++)
    {
      switch (uvars[i].type)
	{
	case UVAR_BOOL:
	  TRY( LALReadConfigBOOLVariable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_INT4:
	  TRY( LALReadConfigINT4Variable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_REAL8:
	  TRY( LALReadConfigREAL8Variable (stat->statusPtr, uvars[i].varp, cfg, uvars[i].name ), stat);
	  break;
	case UVAR_CHAR:
	  stringbuf = NULL;
	  TRY( LALReadConfigSTRINGVariable (stat->statusPtr, &stringbuf, cfg, uvars[i].name ), stat);
	  if (stringbuf)
	    {
	      if ( *(CHAR**)uvars[i].varp != NULL)	 /* something allocated here before? */
		LALFree ( *(CHAR**)uvars[i].varp );
	      *(CHAR**)uvars[i].varp = stringbuf;
	    }
	  break;
	} /* switch type */

    } /* for i < numvars */


  /* ok, that should be it: check if there were more definitions we did not read */
  TRY (LALCheckConfigReadComplete (stat->statusPtr, cfg, CONFIGFILE_ERROR), stat);	

  TRY( LALDestroyConfigData (stat->statusPtr, &cfg), stat);


  DETATCHSTATUSPTR(stat);
  RETURN (stat);

} /* LALReadCfgFileInput() */

/*--------------------------------------------------------------------------------
 * free memory associated with user-variable array (basically free all the strings 
 *--------------------------------------------------------------------------------*/
void 
LALFreeUserVars (LALStatus *stat, UserVariable *uvars)
{
  UINT4 i;

  INITSTATUS( stat, "FreeUserVars", USERINPUTC );

  ASSERT (uvars != NULL,  stat, USERINPUTH_ENULL,  USERINPUTH_MSGENULL);

  i = 0;
  while ( uvars[i].name != NULL )
    {
      if ( (uvars[i].type == UVAR_CHAR) && (uvars[i].varp != NULL) && (*(CHAR**)(uvars[i].varp) != NULL) )
	{
	  LALFree ( *((CHAR**)uvars[i].varp) );
	  *(CHAR**)uvars[i].varp = NULL; 	/* important */
	}
      
      i++;
    } /* while uvars != NULL */
  
  RETURN(stat);
} /* LALFreeUserVars() */

/*----------------------------------------------------------------------
 * assemble all help-info from uvars into a help-string
 *----------------------------------------------------------------------*/
void
LALGetUvarHelpString (LALStatus *stat, CHAR **helpstring, UserVariable *uvars)
{
  UINT4 numvars, i, mem;
  CHAR strbuf[512];	/* should be enough for one line...*/
  CHAR defaultstr[100];
  /* we need strings for UVAR_BOOL, UVAR_INT4, UVAR_REAL8, UVAR_CHAR: */
  const CHAR *typestr[] = {"BOOL", "INT", "FLOAT", "STRING"}; 

  INITSTATUS (stat, "GetUvarHelpString", USERINPUTC);

  ASSERT (helpstring != NULL, stat, USERINPUTH_ENULL, USERINPUTH_MSGENULL);
  ASSERT ( *helpstring == NULL, stat, USERINPUTH_ENONULL, USERINPUTH_MSGENONULL);

  /* count number of user-variables  and estimate lower-limit on memory requirements*/
  numvars = 0;
  mem = 0;
  while (uvars[numvars].name != NULL)
    {
      mem += strlen (uvars[numvars].name) + 10;
      if (uvars[numvars].help)
	mem += strlen (uvars[numvars].help) + 10;

      numvars ++;
    }/* while uvars[numvars] */
  
  if ( (*helpstring = LALCalloc (1, mem)) == NULL ) {
    ABORT (stat,  USERINPUTH_EMEM,  USERINPUTH_MSGEMEM);
  }
  
  for (i=0; i < numvars; i++)
    {
      switch (uvars[i].type)
	{
	case UVAR_BOOL:
	  sprintf (defaultstr, *(BOOLEAN*)uvars[i].varp ? "True" : "False");
	  break;
	case UVAR_INT4:
	  sprintf (defaultstr, "%d", *(INT4*)uvars[i].varp );
	  break;
	case UVAR_REAL8:
	  sprintf (defaultstr, "%.3g", *(REAL8*)uvars[i].varp);
	  break;
	case UVAR_CHAR:
	  if ( *(CHAR**)uvars[i].varp )
	    strcpy (defaultstr, *(CHAR**)uvars[i].varp );
	  else
	    defaultstr[0] = '\0';
	  break;
	} /* switch type */

      sprintf (strbuf, "   --%-14s (-%c) \t%s\t%s  (Default: %s)\n", 
	       uvars[i].name, uvars[i].optchar, typestr[uvars[i].type], uvars[i].help, defaultstr);

      strcat (*helpstring, strbuf);

    } /* for i < numvars */


  RETURN(stat);
} /* LALGetUvarHelpString() */


/*----------------------------------------------------------------------
 * This function puts all the above pieces together, and basically does
 * everything: get config-filename from cmd-line (if found),
 * then interpret config-file and then the command-line
 *----------------------------------------------------------------------*/
void
LALReadUserInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars)
{
  INT2 i;
#define MAXFILENAMELENGTH 512
  CHAR fname[MAXFILENAMELENGTH] = "";

  INITSTATUS( stat, "ReadUserInput", USERINPUTC);
  ATTATCHSTATUSPTR (stat);

  /*----------------------------------------------------------------------
   * pre-process command-line: we want config-file and debug-level 
   */
  for (i=1; i < argc; i++)
    {
      /* was a config-file specified ? */
      if ( argv[i][0] == '@' )
	{
	  strncpy (fname, argv[i]+1, MAXFILENAMELENGTH);
	  fname[MAXFILENAMELENGTH-1] = '\0'; /* terminate in case it was truncated */
	}

    } /* for cmd-line */
  /*----------------------------------------------------------------------*/

  /* if config-file specified, read from that first */
  if (fname[0] != '\0') {
    TRY (LALReadCfgFileInput (stat->statusPtr, fname, uvars ), stat);
  }

  /* now do proper cmdline parsing: overloads config-file settings */
  TRY (LALReadCmdLineInput (stat->statusPtr, argc, argv, uvars), stat);

  DETATCHSTATUSPTR (stat);
  RETURN (stat);
  
} /* LALReadUserInput() */

