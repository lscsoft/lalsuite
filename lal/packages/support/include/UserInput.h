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

#define USERINPUTH_MSGENULL 	"Arguments contained an unexpected null pointer."
#define USERINPUTH_MSGENONULL	"Output pointer is not NULL"
#define USERINPUTH_MSGEMEM	"Out of memory"
#define USERINPUTH_MSGEOPT	"Unknown command-line option encountered"

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
\subsection*{Types}

\subsubsection*{Structure \texttt{LALUserVariable}}
\idx[Type]{LALUserVariable}

This structure defines a "user-variable", which can be read
automagically from command-line and/or config-file.
</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
  const CHAR *name;	/* full name */
  LALUserVarType type;	/* type: bool, int, float or string */
  CHAR optchar;		/* cmd-line character */
  const CHAR *help;	/* help-string */
  void *varp;		/* pointer to the actual C-variable */
} LALUserVariable;
/* </lalVerbatim> */

/* this helps users to set up a consistent array of UserVariable's */
#define regUserVar(name,type,option,help) { #name, type, option, help, &(uvar_ ## name) }


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{UserInputHV}}
\newpage\input{UserInputC}
\newpage\input{UserInputTestC}
******************************************************* </lalLaTeX> */

/* Function prototypes */
void LALReadUserInput (LALStatus *stat, int argc, char *argv[], const LALUserVariable *uvars);
void LALReadCmdLineInput (LALStatus *stat, int argc, char *argv[], const LALUserVariable *uvars);
void LALReadCfgFileInput (LALStatus *stat, const CHAR *cfgfile, const LALUserVariable *uvars );
void LALFreeUserVars (LALStatus *stat, LALUserVariable *uvars);
void LALGetUvarHelpString (LALStatus *stat, CHAR **helpstring, const LALUserVariable *uvars);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */
#endif  /* Double-include protection. */



