/* <lalVerbatim file="LIGOLwXMLReadHV">
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LIGOLwXMLRead.h}}
\label{s:LIGOLwXMLRead.h}

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LIGOLwXMLRead.h>
\end{verbatim}

</lalLaTeX> */


#ifndef _LIGOLWXMLREAD_H
#define _LIGOLWXMLREAD_H

#include <stdio.h>
#include <stdlib.h>
#include <metaio.h>
#include <lal/LALDatatypes.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/LIGOMetadataTables.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID( LIGOLWXMLREADH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error conditions}
\input{LIGOLwXMLReadHE}

</lalLaTeX> */

/* <lalErrTable file="LIGOLwXMLReadHE"> */
#define LIGOLWXMLREADH_ENULL 1
#define LIGOLWXMLREADH_ENNUL 2
#define LIGOLWXMLREADH_EALOC 3
#define LIGOLWXMLREADH_EUTAB 4
#define LIGOLWXMLREADH_ENCOL 5
#define LIGOLWXMLREADH_ENTAB 6
#define LIGOLWXMLREADH_EPARS 7

#define LIGOLWXMLREADH_ENTAB 8
#define LIGOLWXMLREADH_EENDT 8
#define LIGOLWXMLREADH_ETMSM 9
#define LIGOLWXMLREADH_ETNOP 10

#define LIGOLWXMLREADH_MSGENULL "Null pointer"
#define LIGOLWXMLREADH_MSGENNUL "Non-null pointer"
#define LIGOLWXMLREADH_MSGEALOC "Memory allocation error"
#define LIGOLWXMLREADH_MSGEUTAB "Unknown metadata table type"
#define LIGOLWXMLREADH_MSGENCOL "Unable to find table column"
#define LIGOLWXMLREADH_MSGENTAB "Requested table not found in file"
#define LIGOLWXMLREADH_MSGEPARS "Error parsing table"

#define LIGOLWXMLREADH_MSGENTAB "No table type specified"
#define LIGOLWXMLREADH_MSGEENDT "Ending a table without an beginning a table"
#define LIGOLWXMLREADH_MSGETMSM "Table type mismatch"
#define LIGOLWXMLREADH_MSGETNOP "Table not begun for writing"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Structures}

\subsubsection*{Type \texttt{tagMetaTableDirectory}}
\idx[Type]{tagMetaTableDirectory}
</lalLaTeX>
#endif

/* <lalVerbatim> */
typedef struct
tagMetaTableDirectory
{
  CHAR *name;
  INT4   pos;
  INT4   idx;
}
MetaTableDirectory;
/* </lalVerbatim> */

/* <lalLaTeX>
\newpage\input{LIGOLwXMLReadC}
</lalLaTeX> */

void
LALCreateMetaTableDir(
    LALStatus              *status,
    MetaTableDirectory    **tableDir,
    const MetaioParseEnv    env,
    MetadataTableType       table
    );

void
LALSnglBurstTableFromLIGOLw (
    LALStatus          *status,
    SnglBurstTable    **eventHead,
    CHAR               *fileName
    );

int
LALSnglInspiralTableFromLIGOLw (
    SnglInspiralTable **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    );

/* these functions need to be lalified, but they are in support... */

int
SnglInspiralTableFromLIGOLw (
    SnglInspiralTable **eventHead,
    CHAR               *fileName,
    INT4                startEvent,
    INT4                stopEvent
    );

int
InspiralTmpltBankFromLIGOLw (
    InspiralTemplate   **bankHead,
    CHAR                *fileName,
    INT4                startTmplt,
    INT4                stopTmplt
    );

int
SimInspiralTableFromLIGOLw (
    SimInspiralTable   **simHead,
    CHAR                *fileName,
    INT4                 startTime,
    INT4                 endTime
    );

int
SearchSummaryTableFromLIGOLw (
    SearchSummaryTable **sumHead,
    CHAR                *fileName
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LIGOLWXMLREAD_H */
