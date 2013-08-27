/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Benjamin Owen, Reinhard Prix
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
 * File Name: LALXMGRInterface.c
 *
 * Author: Brady, P. R., Brown, D. A., and Owen, B. J.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
\author Brady, P. R., Brown, D. A., and Owen, B. J.
\file
\ingroup pulsarTODO

\heading{Module \ref LALXMGRInterface.c}
\latexonly\tag{ss_LALXMGRInterface_c}\endlatexonly

Functions for creating XMGR graphs from LAL structures and functions.


*/

#include <math.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALXMGRInterface.h>
#include <lal/TwoDMesh.h>

/* ---------------------------------------------------------------------- */
void
LALXMGROpenFile (
    LALStatus          *status,
    FILE              **fp,
    CHAR               *title,
    CHAR               *fileName
    )
{
  const char    xmgrHeader[] =
    "# ACE/gr parameter file\n#\n"
    "@version 40102\n"
    "@page layout free\n"
    "@ps linewidth begin 1\n"
    "@ps linewidth increment 2\n"
    "@hardcopy device 1\n"
    "@page 5\n"
    "@page inout 5\n"
    "@link page off\n"
    "@default linestyle 1\n"
    "@default linewidth 1\n"
    "@default color 1\n"
    "@default char size 1.000000\n"
    "@default font 4\n"
    "@default font source 0\n"
    "@default symbol size 1.000000\n"
    "@timestamp off\n"
    "@timestamp 0.03, 0.03\n"
    "@timestamp linewidth 1\n"
    "@timestamp color 1\n"
    "@timestamp rot 0\n"
    "@timestamp font 4\n"
    "@timestamp char size 1.000000\n"
    "@timestamp def \"Wed Jan  2 19:55:05 2002\"\n";

  const CHAR    titleHeader[] =
    "@with string\n"
    "@    string on\n"
    "@    string loctype view\n"
    "@    string 0.10, 0.95\n"
    "@    string linewidth 2\n"
    "@    string color 1\n"
    "@    string rot 0\n"
    "@    string font 4\n"
    "@    string just 0\n"
    "@    string char size 0.750000\n";

  INITSTATUS(status);

  ASSERT( fp, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );
  ASSERT( !(*fp), status,
      LALXMGRINTERFACEH_ENNUL, LALXMGRINTERFACEH_MSGENNUL );
  ASSERT( fileName, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );

  if ( ! (*fp = LALFopen( fileName, "w" )) )
  {
    ABORT( status, LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGEOPEN );
  }

  fprintf( *fp, "%s", xmgrHeader );

  if ( title )
  {
    fprintf( *fp, "%s", titleHeader );
    fprintf( *fp,  "@    string def \"%s\"\n", title );
  }

  RETURN( status );
}


/* ---------------------------------------------------------------------- */
void
LALXMGRCloseFile (
    LALStatus          *status,
    FILE               *fp
    )
{
  INITSTATUS(status);

  ASSERT( fp, status, LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );

  if ( fclose( fp ) )
  {
    ABORT( status, LALXMGRINTERFACEH_EFCLO, LALXMGRINTERFACEH_MSGEFCLO );
  }

  RETURN( status );
}

/* ---------------------------------------------------------------------- */
void
LALXMGRCreateGraph (
    LALStatus          *status,
    XMGRGraphVector    *graphVec
    )
{
  XMGRGraph            *graph;
  XMGRGraph            *newGraph;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( graphVec, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );

  graph = graphVec->data;


  /*
   *
   * allocate enough memory in the graph array for the graphs
   *
   */


  if ( graphVec->length > 5 )
  {
    ABORT( status, LALXMGRINTERFACEH_ENGRA, LALXMGRINTERFACEH_MSGENGRA );
  }

  graph = (XMGRGraph *)
    LALRealloc( graph, ++(graphVec->length) * sizeof(XMGRGraph) );
  if ( ! graphVec->data )
  {
    ABORT( status, LALXMGRINTERFACEH_EALOC, LALXMGRINTERFACEH_MSGEALOC );
  }

  /* zero the array element we have just created */
  newGraph = graph + graphVec->length - 1;
  memset( newGraph, 0, sizeof(XMGRGraph) );


  /*
   *
   * allocate memory for the parameter structures in the graph
   *
   */


  /* always create xy graphs at the moment: may want to extend this */
  LALCHARCreateVector( status->statusPtr, &(newGraph->type), (UINT4) 3 );
  CHECKSTATUSPTR( status );

  sprintf( newGraph->type->data, "xy" );

  newGraph->xaxis = (XMGRAxisParams *) LALCalloc( 1, sizeof(XMGRAxisParams) );
  newGraph->yaxis = (XMGRAxisParams *) LALCalloc( 1, sizeof(XMGRAxisParams) );

  if ( ! graph->xaxis || ! graph->yaxis )
  {
    ABORT( status, LALXMGRINTERFACEH_EALOC, LALXMGRINTERFACEH_MSGEALOC );
  }


  /*
   *
   * set and reset the viewport coords depending on number of graphs
   *
   */


  {
    UINT4       i;
    REAL4       w, h;   /* width and height of graphs */
    REAL4       xshift = 0.02;

    /* we use the same algorithm from dataviewer (with a few modifications) */
    for ( i = 0; i < graphVec->length; ++i )
    {
      graph[i].viewx[0] = 0.0;
      graph[i].viewy[0] = 0.0;
    }
    switch ( graphVec->length )
    {
      case 1:
        graph[0].viewx[0] = 0.06; graph[0].viewy[0] = 0.028;
        w = 0.90;
        h = 0.88;
        break;
      case 2:
        graph[0].viewx[0] = 0.06; graph[0].viewy[0] = 0.026;
        graph[1].viewx[0] = 0.06; graph[1].viewy[0] = 0.51;
        w = 0.90;
        h = 0.42;
        break;
      case 3:
        graph[3].viewx[0] = 0.54; graph[3].viewy[0] = 0.57;
      case 4:
        graph[0].viewx[0] = 0.04; graph[0].viewy[0] = 0.10;
        graph[1].viewx[0] = 0.04; graph[1].viewy[0] = 0.57;
        graph[2].viewx[0] = 0.54; graph[2].viewy[0] = 0.10;
        w = 0.45;
        h = 0.35;
        break;
      case 5:
        graph[5].viewx[0] = 0.56; graph[5].viewy[0] = 0.62;
      case 6:
        graph[0].viewx[0] = 0.10; graph[0].viewy[0] = 0.02;
        graph[1].viewx[0] = 0.10; graph[1].viewy[0] = 0.32;
        graph[2].viewx[0] = 0.10; graph[2].viewy[0] = 0.62;
        graph[3].viewx[0] = 0.56; graph[3].viewy[0] = 0.02;
        graph[4].viewx[0] = 0.56; graph[4].viewy[0] = 0.32;
        w = 0.37;
        h = 0.26;
        break;
      default:
        ABORT( status, LALXMGRINTERFACEH_ENGRA, LALXMGRINTERFACEH_MSGENGRA );
    }

    for ( i = 0; i < graphVec->length; ++i )
    {
      graph[i].viewy[1] = graph[i].viewy[1] + h;
      graph[i].viewx[1] = graph[i].viewx[0] + w;
      graph[i].viewx[0] += xshift;
    }
  }


  /* XXX */


  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* ---------------------------------------------------------------------- */
void
LALXMGRGPSTimeToTitle(
    LALStatus          *status,
    CHARVector         *title,
    LIGOTimeGPS        *startGPS,
    LIGOTimeGPS        *stopGPS,
    CHAR               *comment
    )
{
  CHARVector           *startString = NULL;
  CHARVector           *stopString  = NULL;
  struct tm             thisDate;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( title, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );
  ASSERT( startGPS, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );
  ASSERT( stopGPS, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );
  ASSERT( comment, status,
      LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );

  LALCHARCreateVector( status->statusPtr, &startString, (UINT4) 64 );
  CHECKSTATUSPTR( status );
  LALCHARCreateVector( status->statusPtr, &stopString, (UINT4) 64 );
  CHECKSTATUSPTR( status );

  XLALGPSToUTC(&thisDate, startGPS->gpsSeconds);
  strftime(startString->data, startString->length, "%Y-%m-%d %H:%M:%S UTC %a", &thisDate);

  XLALGPSToUTC(&thisDate, stopGPS->gpsSeconds);
  strftime(stopString->data, stopString->length, "%Y-%m-%d %H:%M:%S UTC %a", &thisDate);

  snprintf( title->data, title->length * sizeof(CHAR),
	    "%s from %s to %s", comment, startString->data, stopString->data );

  LALCHARDestroyVector( status->statusPtr, &startString );
  CHECKSTATUSPTR( status );
  LALCHARDestroyVector( status->statusPtr, &stopString);
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void LALXMGRPlotMesh( LALStatus          *status,
                      TwoDMeshNode       *head,
                      FILE               *fp,
                      TwoDMeshParamStruc *mesh )
{
  INT4          i;
  INT4          set = 0;
  TwoDMeshNode *node;
  REAL4 xlast = 0, ylast = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( fp, status, LALXMGRINTERFACEH_ENULL, LALXMGRINTERFACEH_MSGENULL );
  ASSERT( head, status, LALXMGRINTERFACEH_ENULL,
          LALXMGRINTERFACEH_MSGENULL );

  /* Plot boundary. */
  fprintf( fp, "@target s%d\n@type xy\n", set );
  for( i = 0; i < 1000; i++ )
  {
    REAL4 x = mesh->domain[0] + i*(mesh->domain[1]-mesh->domain[0])/1000;
    REAL4 y[2];
    mesh->getRange( status->statusPtr, y, x, mesh->rangeParams );
    fprintf( fp, "%e %e\n", x, y[1] );
    if (i == 0 ) {
      xlast = x;
      ylast = y[1];
    }
  }
  for( i = 1000; i >= 0; i-- )
  {
    REAL4 x = mesh->domain[0] + i*(mesh->domain[1]-mesh->domain[0])/1000;
    REAL4 y[2];
    mesh->getRange( status->statusPtr, y, x, mesh->rangeParams );
    fprintf( fp, "%e %e\n", x, y[0] );
  }
  fprintf (fp, "%e %e\n", xlast, ylast );

  set += 2;

  /* Plot mesh points. */
  fprintf( fp, "@s%d symbol 9\n@s%d symbol size 0.33\n", set, set );
  fprintf( fp, "@s%d line type 0\n", set );
  fprintf( fp, "@target s%d\n@type xy\n", set );
  for( node = head; node; node = node->next )
  {
    fprintf( fp, "%e %e\n", node->x, node->y );
  }



  set++;

#if 0
  /* Plot ellipses if requested. */
  for( node = head; node; node = node->next )
  {
  }
#endif

  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALXMGRPlotMesh() */
