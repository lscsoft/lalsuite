/************************************ <lalVerbatim file="ConfigFileTestCV">
Author: Reinhard Prix
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX> 

\subsection{Program \texttt{ConfigFileTest.c}}
\label{s:ConfigFileTest.c}

Tests the routines in \verb@ConfigFile.h@.

\subsubsection*{Usage}
\begin{verbatim}
ConfigFileTest 
\end{verbatim}

\subsubsection*{Description}

Do some standard-tests for the config-file reading routines. No
extensive error-condition checking is done here, we only check if the
basic functionality works.

\subsubsection*{Exit codes}

\input{ConfigFileErrors}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ConfigFileTestCV}}

</lalLaTeX> */

#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>

NRCSID (CONFIGFILETESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="ConfigFileErrors"> */
#define CONFIGFILETESTC_ENORM 		0
#define CONFIGFILETESTC_EFLOAT 		1
#define CONFIGFILETESTC_EINT 		2
#define CONFIGFILETESTC_EBOOL 		3
#define CONFIGFILETESTC_ESTRING 	4
#define CONFIGFILETESTC_ESUB	 	5

#define CONFIGFILETESTC_MSGENORM 	"Normal exit"
#define CONFIGFILETESTC_MSGEFLOAT 	"Read-in REAL8 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEINT 	"Read-in INT4 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEBOOL 	"Read-in BOOL variable is not what it should be..."
#define CONFIGFILETESTC_MSGESTRING 	"Read-in STRING-variable is not what it should be..."
#define CONFIGFILETESTC_MSGESUB	 	"Error occurred in sub-routine"

/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=3;


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, CONFIGFILETESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              CONFIGFILETESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( CONFIGFILETESTC_ESUB, CONFIGFILETESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return CONFIGFILETESTC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

#define TRUE (1==1)
#define FALSE (1==0)

int main(int argc, char *argv[]){ 
  static LALStatus       status;  
  static LALParsedDataFile *cfgdata;
  
  BOOLEAN testBool;
  CHAR *string1 = NULL;
  CHARVector *string2 = NULL;
  CHAR *string2b = NULL;
  CHAR *string3 = NULL;
  INT4 someint;
  REAL8 somefloat;
  BOOLEAN wasRead = FALSE;

  if ( argc > 1 )
    LALPrintError ("WARNING: commond-line arguments useless here \n");

  SUB (LALParseDataFile (&status, &cfgdata, "ConfigFileSample.cfg"), &status);

  SUB (LALReadConfigREAL8Variable  (&status, &somefloat, cfgdata, "float1", &wasRead), &status);
  SUB (LALReadConfigSTRINGVariable (&status, &string1,   cfgdata, "string1", &wasRead), &status);

  SUB (LALReadConfigINT4Variable   (&status, &someint,   cfgdata, "int1", &wasRead), &status);

  SUB (LALCHARCreateVector (&status, &string2, 35), &status);  
  SUB (LALReadConfigSTRINGNVariable(&status, string2,   cfgdata, "string2", &wasRead), &status);

  SUB (LALReadConfigSTRINGVariable(&status, &string2b,   cfgdata, "string2", &wasRead), &status);
  SUB (LALReadConfigSTRINGVariable(&status, &string3,   cfgdata, "string3", &wasRead), &status);

  SUB (LALReadConfigBOOLVariable   (&status, &testBool,  cfgdata, "testBool", &wasRead), &status);

  SUB (LALCheckConfigReadComplete (&status, cfgdata, CONFIGFILE_ERROR), &status);
  SUB (LALDestroyParsedDataFile (&status, &cfgdata), &status);

  /* now check the stuff got read-in correctly */
  if (somefloat != 1.0) {
    ERROR (CONFIGFILETESTC_EFLOAT, CONFIGFILETESTC_MSGEFLOAT, 0);
    return (CONFIGFILETESTC_EFLOAT);
  }
  if ( strcmp (string1, "some text. You can also use line-continuation") ) {
    LALPrintError ("read-in: '%s'\n", string1);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if (someint != 4) {
    ERROR (CONFIGFILETESTC_EINT, CONFIGFILETESTC_MSGEINT, 0);
    return (CONFIGFILETESTC_EINT);
  }
  if ( strcmp(string2->data, "this is also possible, and # here ") ) {
    LALPrintError ("read-in: '%s'\n", string2->data);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string2b, "this is also possible, and # here does nothing ")) {
    LALPrintError ("read-in: '%s'\n", string2b);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string3, "how about #quotes AND line-continuation?") ) {
    LALPrintError ("read-in: '%s'\n", string3);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }

  
  if ( testBool != 0 ) {
    ERROR (CONFIGFILETESTC_EBOOL, CONFIGFILETESTC_MSGEBOOL, 0);
    return (CONFIGFILETESTC_EBOOL);
  }

  LALFree (string1);
  LALCHARDestroyVector (&status, &string2);
  LALFree (string2b);
  LALFree (string3);

  LALCheckMemoryLeaks(); 

  return CONFIGFILETESTC_ENORM;
}


