/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Benjamin Owen
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

/*-----------------------------------------------------------------------
 *
 * File Name: LALXMGRInterface.h
 *
 * Author: Brady, P.R, and Brown, D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brady P., R., and Brown, D. A.
\file
\ingroup pulsarTODO

\heading{Header \ref LALXMGRInterface.h}
\latexonly\tag{s_LALXMGRInterface_h}\endlatexonly

Provides protypes, structures and functions to allow visualisation of
the events generated \c findchirp and the \c inspiral shared
object.

\heading{Synopsis}

\code
#include <lal/LALXMGRInterface.h>
\endcode


*/

#ifndef _LALXMGRINTERFACEH_H
#define _LALXMGRINTERFACEH_H

#include <lal/LALDatatypes.h>
#include <lal/Date.h>
#include <lal/TwoDMesh.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**\name Error Codes */ /*@{*/
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
/*@}*/


/*
 *
 * typedefs of structures used by findchip view functions
 *
 */


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

void
LALXMGRPlotMesh(
    LALStatus          *status,
    TwoDMeshNode       *head,
    FILE               *fp,
    TwoDMeshParamStruc *mesh
    );





#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALXMGRINTERFACEH_H */
