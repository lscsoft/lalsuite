/*
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

/** \defgroup UserInput 
 * \ingroup support
 * \author Reinhard Prix
 * \date $Date$
 * \brief Module for general parsing of user-input from config-file and/or command-line.
 *
 * More documentation on how to use this module will appear here soon!!
 */
 
/** \file 
 * \ingroup UserInput
 * \author Reinhard Prix
 * \date $Date$
 * \brief Header file defining the API for the UserInput modules.
 */

/* <lalVerbatim file="UserInputHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/**************************************************** <lalLaTeX>
\section{Header \texttt{UserInput.h}}
\label{s:UserInput.h}

Module for general parsing of "user input" from config-file and/or command-line.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/UserInput.h>
\end{verbatim}

****************************************************** </lalLaTeX> */

#ifndef _USERINPUT_H  /* Double-include protection. */
#define _USERINPUT_H

#include <lal/ConfigFile.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( USERINPUTH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
/** \name Error codes */
/*@{*/
#define USERINPUTH_ENULL 	1
#define USERINPUTH_ENONULL	2
#define USERINPUTH_EMEM		3
#define USERINPUTH_EOPT		4
#define USERINPUTH_ENOUVARS	5
#define USERINPUTH_ECMDLARG     6
#define USERINPUTH_EUNKNOWN	7
#define USERINPUTH_ENOTSET	8
#define USERINPUTH_EDEBUG	9
#define USERINPUTH_EONECONFIG   10
#define USERINPUTH_ERECFORMAT   11


#define USERINPUTH_MSGENULL 	"Arguments contained an unexpected null pointer."
#define USERINPUTH_MSGENONULL	"Output pointer is not NULL"
#define USERINPUTH_MSGEMEM	"Out of memory"
#define USERINPUTH_MSGEOPT	"Unknown command-line option encountered"
#define USERINPUTH_MSGENOUVARS	"No user-variables have been registered!"
#define USERINPUTH_MSGECMDLARG   "Illegal command-line argument"
#define USERINPUTH_MSGEUNKNOWN	"Unknown user-variable"
#define USERINPUTH_MSGENOTSET	"Required user-variable was not set"
#define USERINPUTH_MSGEDEBUG	"lalDebugLevel can only be read before ANY mallocs(), even hidden.."
#define USERINPUTH_MSGEONECONFIG "Currently one ONE config-file can be specified using '@'"
#define USERINPUTH_MSGERECFORMAT   "Unknown format for recording user-input"
/*@}*/
/*************************************************** </lalErrTable> */

/* <lalLaTeX> 
\subsection*{Macros}

The following macros make it a bit easier to register the
user-variables in a constistent way. The convention followed is that
the C-variable corresponding to the user-variable \verb+Name+ is
called \verb+uvar_[Name]+. These macros register a user-variable "name"
of type \verb+REAL8+, \verb+INT4+, \verb+BOOLEAN+ or \verb+CHAR*+
respectively. 
</lalLaTeX> */
/* <lalVerbatim> 
regREALUserVar(stat,name,option,help)	
regINTUserVar(stat,name,option,help) 
regBOOLUserVar(stat,name,option,help)
regSTRINGUserVar(stat,name,option,help)
</lalVerbatim> */
/* <lalLaTeX> 
Some examples of use:
</lalLaTeX> */
/* <lalVerbatim> 
CHAR *uvar_inDataFilename = NULL;
REAL8 uvar_startTime;
BOOLEAN uvar_binaryoutput = FALSE;
INT4 uvar_nTsft;
regSTRINGUserVar(stat, inDataFilename,  'i', UVAR_OPTIONAL, "Name of input parameter file");
regREALUserVar(stat,   startTime,	'G', UVAR_REQUIRED, "Detector GPS time to start data");
regBOOLUserVar(stat,   binaryoutput,	'b', UVAR_OPTIONAL, "Output time-domain data in binary format");
regINTUserVar(stat,    nTsft,		'N', UVAR_REQUIRED, "Number of SFTs nTsft");
</lalVerbatim> */

#define LALregREALUserVar(status,name,option,flag,help) \
TRY(LALRegisterREALUserVar((status)->statusPtr, #name, option, flag, help,&(uvar_ ## name)), status)

#define LALregINTUserVar(status,name,option,flag,help) \
TRY(LALRegisterINTUserVar((status)->statusPtr, #name, option,flag, help,&(uvar_ ## name)), status)

#define LALregBOOLUserVar(status,name,option,flag,help) \
TRY(LALRegisterBOOLUserVar((status)->statusPtr, #name, option, flag, help, &(uvar_ ## name)),status)

#define LALregSTRINGUserVar(status,name,option,flag,help) \
TRY(LALRegisterSTRINGUserVar((status)->statusPtr, #name, option, flag, help, &(uvar_ ## name)),status)

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{UserInputHV}}
\newpage\input{UserInputC}
\newpage\input{UserInputTestC}
******************************************************* </lalLaTeX> */


/** State-flags: variable is optional, required, help, developer or was_set */
typedef enum {
  UVAR_OPTIONAL		= 0,	/**< not required, and hasn't been set */
  UVAR_REQUIRED 	= 1<<0,	/**< we require the user to set this variable */
  UVAR_HELP		= 1<<1,	/**< special variable: trigger output of help-string */
  UVAR_DEVELOPER	= 1<<2,	/**< OPTIONAL and hidden in help-output at lalDebugLevel==0 */
  UVAR_WAS_SET 		= 1<<7	/**< flag that this user-var has been set by user */
} UserVarState;

/** Format for logging User-input: configFile- or cmdLine-style.
 * This determines the format of the string returned from LALLogUserInput().
 */
typedef enum {
  UVAR_LOGFMT_CFGFILE,	/**< return UserVars as a config-file */
  UVAR_LOGFMT_CMDLINE,	/**< return UserVars as a command-line */
  UVAR_LOGFMT_LAST
} UserVarLogFormat;

/* Function prototypes */
void LALRegisterREALUserVar(LALStatus *, 
			    const CHAR *name, 
			    CHAR optchar, 
			    UserVarState flag,
			    const CHAR *helpstr, 
			    REAL8 *cvar);

void LALRegisterINTUserVar (LALStatus *,  
			    const CHAR *name, 
			    CHAR optchar, 
			    UserVarState flag, 
			    const CHAR *helpstr, 
			    INT4 *cvar);

void 
LALRegisterBOOLUserVar (LALStatus *, 
			const CHAR *name, 
			CHAR optchar, 
			UserVarState flag,
			const CHAR *helpstr, 
			BOOLEAN *cvar);

void
LALRegisterSTRINGUserVar (LALStatus *,
			  const CHAR *name,
			  CHAR optchar, 
			  UserVarState flag,
			  const CHAR *helpstr, 
			  CHAR **cvar);

void LALDestroyUserVars (LALStatus *);

void LALUserVarReadAllInput(LALStatus *, int argc, char *argv[]);
void LALUserVarReadCmdline (LALStatus *, int argc, char *argv[]);
void LALUserVarReadCfgfile (LALStatus *, const CHAR *cfgfile);

void LALUserVarHelpString (LALStatus *, CHAR **helpstring, const CHAR *progname);
void LALUserVarCheckRequired (LALStatus *);
INT4 LALUserVarWasSet (void *cvar);
void LALGetDebugLevel (LALStatus *, int argc, char *argv[], CHAR optchar);
void LALUserVarGetLog (LALStatus *, CHAR **logstr,  UserVarLogFormat format);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */



