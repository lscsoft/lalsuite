/************************************ <lalVerbatim file="UserInputHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/**************************************************** <lalLaTeX>
\section{Header \texttt{UserInput.h}}
\label{s:UserInput.h}

Header file for general parsing of "user input" via config-file and/or command-line.

\subsection*{Synopsis}
\begin{verbatim}
#include "UserInput.h"
\end{verbatim}

\noindent This module provides routines for 

******************************************************* </lalLaTeX> */

#ifndef _USERINPUT_H  /* Double-include protection. */
#define _USERINPUT_H

#include "ConfigFile.h"

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( USERINPUTH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define USERINPUTH_ENULL 	1
#define USERINPUTH_ENONULL	6
#define USERINPUTH_EMEM		9
#define USERINPUTH_EOPT		10

#define USERINPUTH_MSGENULL 	"Arguments contained an unexpected null pointer."
#define USERINPUTH_MSGENONULL	"Output pointer is not NULL"
#define USERINPUTH_MSGEMEM	"Out of memory"
#define USERINPUTH_MSGEOPT	"Unknown command-line option encountered"

/*************************************************** </lalErrTable> */

/****************************************************** <lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{LALConfigVar}}
\idx[Type]{LALConfigVar}

\end{description}

******************************************************* </lalLaTeX> */
typedef enum {
  UVAR_BOOL,
  UVAR_INT4,
  UVAR_REAL8,
  UVAR_CHAR
} UserVarType;

typedef struct {
  const CHAR *name;	/* full name */
  UserVarType type;	/* bool, int, float or char */
  CHAR optchar;		/* cmd-line character */
  const CHAR *help;
  void *varp;		/* pointer to the actual C-variable */
} UserVariable;


#define regUserVar(name,type,option,help) { #name, type, option, help, &(uvar_ ## name) }



/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{UserInputHV}}
\newpage\input{UserInputC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALReadUserInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars);
void LALReadCmdLineInput (LALStatus *stat, int argc, char *argv[], UserVariable *uvars);
void LALReadCfgFileInput (LALStatus *stat, const CHAR *cfgfile, UserVariable *uvars );
void LALFreeUserVars (LALStatus *stat, UserVariable *uvars);
void LALGetUvarHelpString (LALStatus *stat, CHAR **helpstring, UserVariable *uvars);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */



