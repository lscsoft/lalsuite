/************************************ <lalVerbatim file="UserInputHV">
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
#define USERINPUTH_ENULL 	1
#define USERINPUTH_ENONULL	2
#define USERINPUTH_EMEM		3
#define USERINPUTH_EOPT		4
#define USERINPUTH_ENOUVARS	5
#define USERINPUTH_EBOOL	6

#define USERINPUTH_MSGENULL 	"Arguments contained an unexpected null pointer."
#define USERINPUTH_MSGENONULL	"Output pointer is not NULL"
#define USERINPUTH_MSGEMEM	"Out of memory"
#define USERINPUTH_MSGEOPT	"Unknown command-line option encountered"
#define USERINPUTH_MSGENOUVARS	"No user-variables have been registered!"
#define USERINPUTH_MSGEBOOL	"Illegal BOOLEAN option-value"
/*************************************************** </lalErrTable> */

/****************************************************** <lalLaTeX>

\subsection*{Constants}
\idx[Type]{UserVarType}

The following constants are used to define the C-type of a user-variable.
</lalLaTeX> */
/* <lalVerbatim> */
typedef enum {
  UVAR_BOOL,	/* boolean */
  UVAR_INT4,	/* integer */
  UVAR_REAL8,	/* float */
  UVAR_STRING	/* string */
} LALUserVarType;
/* </lalVerbatim> */


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
CHAR *uvar_inDataFilename;
REAL8 uvar_startTime;
BOOLEAN uvar_binaryoutput;
INT4 uvar_nTsft;
regSTRINGUserVar(stat, inDataFilename,  'i', "Name of input parameter file");
regREALUserVar(stat,   startTime,	'G', "Detector GPS time to start data");
regBOOLUserVar(stat,   binaryoutput,	'b', "Output time-domain data in binary format");
regINTUserVar(stat,    nTsft,		'N', "Number of SFTs nTsft");
</lalVerbatim> */

#define regREALUserVar(stat,name,option,help) LALRegisterREALUserVar((stat), #name, option, help, &(uvar_ ## name))
#define regINTUserVar(stat,name,option,help) LALRegisterINTUserVar((stat), #name, option, help, &(uvar_ ## name))
#define regBOOLUserVar(stat,name,option,help)LALRegisterBOOLUserVar((stat), #name, option, help, &(uvar_ ## name))
#define regSTRINGUserVar(stat,name,option,help) LALRegisterSTRINGUserVar((stat), #name, option, help, &(uvar_ ## name))


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{UserInputHV}}
\newpage\input{UserInputC}
\newpage\input{UserInputTestC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALRegisterREALUserVar (LALStatus *stat, const CHAR *name, CHAR optchar, const CHAR *helpstr, REAL8 *cvar);
void LALRegisterINTUserVar (LALStatus *stat,  const CHAR *name, CHAR optchar, const CHAR *helpstr, INT4 *cvar);
void LALRegisterBOOLUserVar (LALStatus *stat, const CHAR *name, CHAR optchar, const CHAR *helpstr, BOOLEAN *cvar);
void LALRegisterSTRINGUserVar (LALStatus *stat,const CHAR *name,CHAR optchar, const CHAR *helpstr, CHAR **cvar);

void LALDestroyUserVars (LALStatus *stat);

void LALUserVarReadAllInput(LALStatus *stat, int argc, char *argv[]);
void LALUserVarReadCmdline (LALStatus *stat, int argc, char *argv[]);
void LALUserVarReadCfgfile (LALStatus *stat, const CHAR *cfgfile);

void LALUserVarHelpString (LALStatus *stat, CHAR **helpstring);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */



