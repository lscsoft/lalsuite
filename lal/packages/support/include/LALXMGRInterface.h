/*----------------------------------------------------------------------- 
 * 
 * File Name: LALXMGRInterface.h
 *
 * Author: Brady, P.R, and Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="LALXMGRInterfaceHV">
Author: Brady P., R., and Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{LALXMGRInterface.h}}
\label{s:LALXMGRInterface.h}

Provides protypes, structures and functions to allow visualisation of
the events generated \texttt{findchirp} and the \texttt{inspiral} shared
object.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/LALXMGRInterface.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _LALXMGRINTERFACEH_H
#define _LALXMGRINTERFACEH_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (LALXMGRINTERFACEH, "$Id$");

#if 0
<lalLaTeX> 
\newpage\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define LALXMGRINTERFACEH_ENULL 1
#define LALXMGRINTERFACEH_ENNUL 2
#define LALXMGRINTERFACEH_EALOC 3
#define LALXMGRINTERFACEH_EOPEN 4
#define LALXMGRINTERFACEH_EFCLO 5
#define LALXMGRINTERFACEH_ENGRA 6
#define LALXMGRINTERFACEH_MSGENULL "Null pointer"
#define LALXMGRINTERFACEH_MSGENNUL "Non-null pointer"
#define LALXMGRINTERFACEH_MSGEALOC "Memory allocation error"
#define LALXMGRINTERFACEH_MSGEOPEN "Error opening file"
#define LALXMGRINTERFACEH_MSGEFCLO "Error closing file"
#define LALXMGRINTERFACEH_MSGENGRA "Already have max number of graphs in array"
/* </lalErrTable> */


/*
 *
 * typedefs of structures used by findchip view functions
 *
 */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

typedef enum
{
  xmgrSymbolNone  = 0,
  xmgrSymbolDot   = 1,
  xmgrSymbolPlus  = 9,
  xmgrSymbolCross = 10
}
XMGRSymbol;

typedef enum
{
  xmgrLineNone,
  xmgrLineSolid,
  xmgrLineDotted,
  xmgrLineDashed
}
XMGRLine;

typedef enum
{
  xmgrColorWhite,
  xmgrColorBlack,
  xmgrColorRed,
  xmgrColorGreen,
  xmgrColorBlue
}
XMGRColor;

typedef struct
tagXMGRDataSet
{
  XMGRSymbol    symbol;
  XMGRColor     symbolColor;
  REAL4         symbolSize;
  XMGRLine      line;
  XMGRColor     lineColor;
  REAL4         lineWidth;
  CHARVector   *name;
  REAL8Vector  *x;
  REAL8Vector  *y;
}
XMGRDataSet;

typedef struct
tagXMGRDataSetVector
{
  UINT4         length;
  XMGRDataSet  *data;
}
XMGRDataSetVector;

typedef struct
tagXMGRAxisParams
{
  CHARVector   *label;
  CHARVector   *format;
  REAL4         min;
  REAL4         max;
  REAL4         tickMajor;
  REAL4         tickMinor;
}
XMGRAxisParams;

typedef struct
tagXMGRGraph
{
  CHARVector                   *type;
  CHARVector                   *title;
  REAL4                         viewx[2];
  REAL4                         viewy[2];
  XMGRAxisParams               *xaxis;
  XMGRAxisParams               *yaxis;
  XMGRDataSetVector            *setVector;
}
XMGRGraph;

typedef struct
tagXMGRGraphVector
{
  UINT4         length;
  XMGRGraph    *data;
}
XMGRGraphVector;

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{LALXMGRInterfaceHV}}
</lalLaTeX> 
#endif


/*
 *
 * function prototypes
 *
 */

void 
LALXMGROpenFile ( 
    LALStatus          *status,
    FILE              **fp,
    CHAR               *title,
    CHAR               *fileName
    );

void
LALXMGRCloseFile ( 
    LALStatus          *status,
    FILE               *fp
    );

void
LALXMGRTimeTitle (
    LALStatus          *status,
    CHARVector         *title,
    LIGOTimeGPS        *startGPS,
    LIGOTimeGPS        *stopGPS,
    CHAR               *comment
    );

void
LALXMGRCreateGraph (
    LALStatus          *status,
    XMGRGraphVector    *graphVec
    );

void
LALXMGRGPSTimeToTitle(
    LALStatus          *status,
    CHARVector         *title,
    LIGOTimeGPS        *startGPS,
    LIGOTimeGPS        *stopGPS,
    CHAR               *comment
    );

#if 0
<lalLaTeX>
\newpage\input{LALXMGRInterfaceC}
</lalLaTeX> 
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _LALXMGRINTERFACEH_H */
