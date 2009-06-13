/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/********************************* <lalVerbatim file="TwoDMeshPlotCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TwoDMeshPlot.c}}
\label{ss:TwoDMeshPlot.c}

Plots a hierarchical mesh of templates on an 2-dimensional parameter
space.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TwoDMeshPlotCP}
\idx{LALPlotTwoDMesh()}

\subsubsection*{Description}

This routine creates a PostScript plot of the parameter mesh list
pointed to by \verb@mesh@, using the plotting parameters given in
\verb@*params@.  The PostScript output is printed to the writable
output stream \verb@*stream@ using \verb@fprintf()@.

\subsubsection*{Algorithm}

The algorithm is set up so that it requires only one pass through the
list.  After defining PostScript macros to plot mesh points, mesh
tiles, and mismatch ellipses, the routine then defines a macro to plot
the boundary.  Since some PostScript interpreters will fail if a macro
contains too many objects, the boundary-plotting macro may be split
into several macros.

\verb@LALPlotTwoDMesh()@ then calls a static (but LAL-compliant)
subroutine \verb@LALMakeMeshMacro()@ to create one or more macros to
plot the mesh points, tiles, or ellipses, as required by
\verb@*params@.  This subroutine takes a pointer to the head of a list
of mesh points as input, and traverses down the list, calling itself
recursively on any submeshes it encounters (if
\verb@params->maxLevels@ permits).

While plotting the boundary and other mesh objects,
\verb@LALPlotTwoDMesh()@ and \verb@LALMakeMeshMacro()@ keep track of
the bounding box surrounding all plotted objects.  This is used either
to set a bounding box for the overall plot, or to adjust the scale of
the plot, depending on \verb@params->autoscale@.  If the resulting
bounding box is larger than a single $8.5''\times11''$ page,
\verb@LALPlotTwoDMesh()@ will divide the plot area up into pages of
this side, calling the plotting macros on each page.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()             LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoDMeshPlotCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>
#include "TwoDMeshPlot.h"

NRCSID( TWODEMESHPLOTC, "$Id$" );

/* Local constants. */
#define TWODMESHPLOTC_MAXOBJ 797 /* Maximum number of objects in a
                                    PostScript macro */

/* Local datatypes. */
typedef struct tagMeshMacroParamStruc {
  UINT4 nMacro;                  /* running count of macro number */
  UINT4 nObj;                    /* running count number of objects */
  UINT4 level;                   /* current recursion level */
  TwoDMeshPlotStruc *plotParams; /* top-level plotting parameters */
} MeshMacroParamStruc;

/* Local prototypes. */
static void
LALMakeMeshMacro( LALStatus           *stat,
		  FILE                *stream,
		  TwoDMeshNode        *mesh,
		  MeshMacroParamStruc *params );

static void
AdjustBBox( REAL4 x, REAL4 y, TwoDMeshPlotStruc *params );


/* <lalVerbatim file="TwoDMeshPlotCP"> */
void
LALPlotTwoDMesh( LALStatus         *stat,
		 FILE              *stream,
		 TwoDMeshNode      *mesh,
		 TwoDMeshPlotStruc *params )
{ /* </lalVerbatim> */
  UINT4 i;          /* an index */
  UINT4 nObj = 0;   /* counter of number of objects boundary macro */
  UINT4 nMacro = 0; /* counter of number of boundary macros */
  UINT4 nPage = 0;  /* number of pages plotted */
  REAL4 bBox[4];    /* bounding box in plot coordinates */
  REAL4 xOff, yOff; /* horizontal and vertical offsets */
  MeshMacroParamStruc macroParams; /* parameters for
				      LALMakeMeshMacro() */

  INITSTATUS( stat, "LALPlotTwoDMesh", TWODEMESHPLOTC );
  ATTATCHSTATUSPTR( stat );

  /* Check that arguments and their fields exist. */
  ASSERT( stream, stat, TWODMESHPLOTH_ENUL, TWODMESHPLOTH_MSGENUL );
  ASSERT( mesh, stat, TWODMESHPLOTH_ENUL, TWODMESHPLOTH_MSGENUL );
  ASSERT( params, stat, TWODMESHPLOTH_ENUL, TWODMESHPLOTH_MSGENUL );
  if ( params->nLevels ) {
    ASSERT( params->plotPoints, stat, TWODMESHPLOTH_ENUL,
	    TWODMESHPLOTH_MSGENUL );
    ASSERT( params->plotTiles, stat, TWODMESHPLOTH_ENUL,
	    TWODMESHPLOTH_MSGENUL );
    ASSERT( params->plotEllipses, stat, TWODMESHPLOTH_ENUL,
	    TWODMESHPLOTH_MSGENUL );
    ASSERT( params->params, stat, TWODMESHPLOTH_ENUL,
	    TWODMESHPLOTH_MSGENUL );
  }

  /* Perform some setup. */
  if ( params->autoscale )
    memcpy( bBox, params->bBox, 4*sizeof(REAL4) );
  params->bBox[0] = params->bBox[1] = LAL_REAL4_MAX;
  params->bBox[2] = params->bBox[3] = -LAL_REAL4_MAX;
  params->cosTheta = cos( LAL_PI_180*params->theta );
  params->sinTheta = sin( LAL_PI_180*params->theta );
  params->clip = ( ( params->clipBox[2] > params->clipBox[0] ) &&
		   ( params->clipBox[3] > params->clipBox[1] ) );

  /* Write the PostScript header. */
  fprintf( stream,
	   "%%!PS-Adobe-1.0\n"
	   "%%%%Creator: LALPlotTwoDMesh()\n"
	   "%%%%Title: mesh.ps\n"
	   "%%%%BoundingBox: %i %i %i %i\n"
	   "%%%%EndComments\n\n",
	   TWODMESHPLOTH_XMARG, TWODMESHPLOTH_YMARG,
	   TWODMESHPLOTH_XMARG + TWODMESHPLOTH_XSIZE,
	   TWODMESHPLOTH_YMARG + TWODMESHPLOTH_YSIZE );

  /* Write PostScript macros for plotting mesh points.  The macros are
     called simply as "point[N]", where [N] is the recursive submesh
     level. */
  for ( i = 0; i < params->nLevels; i++ )
    if ( params->plotPoints[i] > 0 )
      fprintf( stream,
	       "/point%u { gsave currentpoint translate %f %f scale\n"
	       "  auto auto scale %f rotate\n"
	       "  newpath 0 0 %u 0 360 arc closepath fill grestore }"
	       " def\n",
	       i, 1.0/params->xScale, 1.0/params->yScale,
	       -params->theta, params->plotPoints[i] );
    else if ( params->plotPoints[i] < 0 )
      fprintf( stream,
	       "/point%u { gsave currentpoint translate %f %f scale\n"
	       "  auto auto scale %f rotate\n"
	       "  newpath 0 0 r%u 0 360 arc closepath stroke grestore }"
	       " def\n",
	       i, 1.0/params->xScale, 1.0/params->yScale,
	       -params->theta, i );

  /* Write PostScript macro for plotting ellipses.  The macro is
     called as "[axis1] [axis2] [angle] ellipse", where [axis1] and
     [axis2] are the two principal axis lengths, and [angle] is the
     angle in degrees counterclockwise from the x-axis to the first
     principal axis. */
  fprintf( stream,
	   "/ellipse { gsave currentpoint translate rotate scale\n"
	   "  newpath 0 0 1 0 360 arc closepath stroke grestore } def\n" );

  /* Write PostScript macro for plotting tiles.  The macro is called
     as "[dx] [dy1] [dy2] tile", where [dx] is the half-width of the
     tile, and [dy1], [dy2] are the heights of the corners of the tile
     relative to the centre. */
  fprintf( stream,
	   "/tile { gsave currentpoint translate 3 copy\n"
	   "  dup 4 1 roll exch 4 1 roll moveto sub 0 exch rlineto\n"
	   "  2 copy add neg 4 3 roll -2 mul exch rlineto\n"
	   "  sub neg 0 exch rlineto closepath stroke grestore } def\n" );

  /* Write PostScript macro to clip the x-y area, if necessary. */
  if ( params->clip )
    fprintf( stream,
	     "/xyclip { %f %f moveto %f %f lineto %f %f lineto\n"
	     "  %f %f lineto closepath clip } def\n",
	     params->clipBox[0], params->clipBox[1],
	     params->clipBox[0], params->clipBox[3],
	     params->clipBox[2], params->clipBox[3],
	     params->clipBox[2], params->clipBox[1] );

  /* Write PostScript macros for plotting boundary, if necessary. */
  if ( params->nBoundary >= 2 ) {
    REAL4 x, x0, dx;   /* x coordinate, initial value, and increment */
    REAL4 *yBound;     /* array of y-values of boundary points */
    yBound = (REAL4 *)LALMalloc( 2*params->nBoundary*sizeof(REAL4) );
    if ( yBound == NULL ) {
      ABORT( stat, TWODMESHPLOTH_EMEM, TWODMESHPLOTH_MSGEMEM );
    }

    /* Fill array of boundary points. */
    x0 = params->params->domain[0];
    dx = ( params->params->domain[1] - x0 )/( params->nBoundary - 1 );
    for ( i = 0; i < params->nBoundary - 1; i++ ) {
      x = x0 + i*dx;
      (params->params->getRange)( stat->statusPtr, yBound + 2*i, x,
				  params->params->rangeParams );
      BEGINFAIL( stat )
	LALFree( yBound );
      ENDFAIL( stat );
    }
    x = params->params->domain[1];
    (params->params->getRange)( stat->statusPtr, yBound + 2*i, x,
				params->params->rangeParams );
    BEGINFAIL( stat )
      LALFree( yBound );
    ENDFAIL( stat );

    /* Write macro. */
    fprintf( stream,
	     "/boundary%u {\n"
	     "%f %f moveto\n", nMacro, x0, yBound[1] );
    nObj = 3;
    for ( i = 1; i < params->nBoundary - 1; i++ ) {
      x = x0 + i*dx;
      fprintf( stream, "%f %f lineto\n", x, yBound[2*i+1] );
      AdjustBBox( x, yBound[2*i+1], params );
      nObj += 3;
      if ( nObj > TWODMESHPLOTC_MAXOBJ ) {
	fprintf( stream,
		 "stroke } def\n"
		 "/boundary%u {\n"
		 "%f %f moveto\n", ++nMacro, x, yBound[2*i+1] );
	nObj = 3;
      }
    }
    x = params->params->domain[1];
    fprintf( stream, "%f %f lineto\n", x, yBound[2*i+1] );
    fprintf( stream, "%f %f lineto\n", x, yBound[2*i] );
    AdjustBBox( x, yBound[2*i+1], params );
    AdjustBBox( x, yBound[2*i], params );
    nObj += 6;
    for ( i = params->nBoundary - 2; i < (UINT4)( -1 ); i-- ) {
      x = x0 + i*dx;
      fprintf( stream, "%f %f lineto\n", x, yBound[2*i] );
      AdjustBBox( x, yBound[2*i], params );
      nObj += 3;
      if ( nObj > TWODMESHPLOTC_MAXOBJ ) {
	fprintf( stream,
		 "stroke } def\n"
		 "/boundary%u {\n"
		 "%f %f moveto\n", ++nMacro, x, yBound[2*i] );
	nObj = 3;
      }
    }
    fprintf( stream,
	     "%f %f lineto\n"
	     "stroke } def\n", x0, yBound[1] );
    nMacro++;
    LALFree( yBound );
  }

  /* Set up parameters for LALMakeMeshMacro(). */
  macroParams.nMacro = 0;
  macroParams.nObj = 0;
  macroParams.level = 0;
  macroParams.plotParams = params;

  /* Write PostScript macro for plotting the mesh.  The routine
     MakeMeshMacro() traverses the linked list, adding a line to the
     macro for each node, and calling itself recursively on any
     submeshes. */
  fprintf( stream, "\n/mesh%u {\n", macroParams.nMacro );
  TRY( LALMakeMeshMacro( stat->statusPtr, stream, mesh,
			 &macroParams ), stat );
  fprintf( stream, "} def\n" );

  /* Increment macro counter only if the last macro list is not
     empty. */
  if ( macroParams.nObj > 0 )
    macroParams.nMacro++;

  /* Autoscale the axes, if necessary. */
  if ( params->autoscale ) {
    REAL4 xScaleFac = params->bBox[2] - params->bBox[0];
    REAL4 yScaleFac = params->bBox[3] - params->bBox[1];
    if ( ( xScaleFac == 0.0 ) && ( yScaleFac == 0.0 ) ) {
      ABORT( stat, TWODMESHPLOTH_ENOPLOT, TWODMESHPLOTH_MSGENOPLOT );
    }
    xScaleFac = ( bBox[2] - bBox[0] )/xScaleFac;
    yScaleFac = ( bBox[3] - bBox[1] )/yScaleFac;
    if ( yScaleFac < xScaleFac )
      xScaleFac = yScaleFac;
    params->xScale *= xScaleFac;
    params->yScale *= xScaleFac;
    for ( i = 0; i < 4; i++ )
      params->bBox[i] *= xScaleFac;
    fprintf( stream, "/auto %f def\n", 1.0/xScaleFac );
  } else
    fprintf( stream, "/auto 1 def\n" );

  /* Set up coordinate system. */
  if ( params->bBox[2] > params->bBox[0] ) {
    bBox[0] = params->bBox[0];
    bBox[2] = params->bBox[2];
  } else {
    bBox[0] = params->bBox[2];
    bBox[2] = params->bBox[0];
  }
  if ( params->bBox[3] > params->bBox[1] ) {
    bBox[3] = params->bBox[3];
    bBox[1] = params->bBox[1];
  } else {
    bBox[3] = params->bBox[1];
    bBox[1] = params->bBox[3];
  }
  nPage = 0;

  /* Set the global graphics state. */
  fprintf( stream, "\n0 setlinewidth 0 setgray\n" );

  /* Define an overall clipping region for all pages. */
  fprintf( stream, "%i %i moveto %i %i lineto %i %i lineto\n"
	   "%i %i lineto closepath clip newpath\n",
	   TWODMESHPLOTH_XMARG, TWODMESHPLOTH_YMARG,
	   TWODMESHPLOTH_XMARG,
	   TWODMESHPLOTH_YMARG + TWODMESHPLOTH_YSIZE,
	   TWODMESHPLOTH_XMARG + TWODMESHPLOTH_XSIZE,
	   TWODMESHPLOTH_YMARG + TWODMESHPLOTH_YSIZE,
	   TWODMESHPLOTH_XMARG + TWODMESHPLOTH_XSIZE,
	   TWODMESHPLOTH_YMARG );

  /* Plot macros on each page. */
  for ( yOff = params->bBox[1] - TWODMESHPLOTH_YMARG;
	yOff < params->bBox[3] - TWODMESHPLOTH_YMARG - 1;
	yOff += TWODMESHPLOTH_YSIZE )
    for ( xOff = params->bBox[0] - TWODMESHPLOTH_XMARG;
	  xOff < params->bBox[2] - TWODMESHPLOTH_XMARG - 1;
	  xOff += TWODMESHPLOTH_XSIZE ) {
      fprintf( stream, "\n"
	       "%%%%Page: %u\n"
	       "gsave %f %f translate %f rotate %f %f scale",
	       ++nPage, -xOff, -yOff, params->theta, params->xScale,
	       params->yScale );
      if ( params->clip )
	fprintf( stream, " xyclip\n" );
      else
	fprintf( stream, "\n" );
      for ( i = 0; i < nMacro; i++ )
	fprintf( stream, "boundary%u\n", i );
      for ( i = 0; i < macroParams.nMacro; i++ )
	fprintf( stream, "mesh%u\n", i );
      fprintf( stream, "showpage grestore\n" );
    }

  /* Finished plotting.  Restore params->bBox to its original setting,
     if necessary. */
  fprintf( stream, "\n%%%%EOF\n" );
  if ( params->autoscale )
    memcpy( params->bBox, bBox, 4*sizeof(REAL4) );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


static void
LALMakeMeshMacro( LALStatus           *stat,
		  FILE                *stream,
		  TwoDMeshNode        *mesh,
		  MeshMacroParamStruc *params )
     /* This routine writes lines to the macro to draw the mesh
        points, tiles, and ellipses (as required), for the list
        pointed to by mesh.  It calls itself recursively on any
        submeshes, if the current recursion level is less than the
        maximum requested recursion level.  Along the way it adjusts
        the bounding box for the figure. */
{
  UINT4 rLevel;         /* current recursion depth */
  BOOLEAN plotPoints;   /* whether to plot mesh points at this level */
  BOOLEAN plotEllipses; /* whether to plot ellipses at this level */
  BOOLEAN plotTiles;    /* whether to plot tiles at this level */
  BOOLEAN plotAny;      /* whether to plot anything at this level */

  /* Pointer to the metric function, its optional arguments, and the
     mismatch value, for the current recursion level. */
  void (*getMetric)( LALStatus *, REAL4 [3], REAL4 [2], void *) = NULL;
  void *metricParams = NULL;
  REAL4 mThresh = 0.0;

  INITSTATUS( stat, "LALMakeMeshMacro", TWODEMESHPLOTC );
  ATTATCHSTATUSPTR( stat );

  ASSERT( params, stat, TWODMESHPLOTH_ENUL, TWODMESHPLOTH_MSGENUL );
  ASSERT( params->plotParams, stat, TWODMESHPLOTH_ENUL,
	  TWODMESHPLOTH_MSGENUL );

  /* Do nothing if we're beyond the maximum recursion level. */
  if ( ( rLevel = params->level ) >= params->plotParams->nLevels ) {
    DETATCHSTATUSPTR( stat );
    RETURN( stat );
  }

  /* Otherwise, determine what needs to be plotted at this level. */
  plotPoints = ( params->plotParams->plotPoints[rLevel] != 0 );
  plotTiles = params->plotParams->plotTiles[rLevel];
  plotEllipses = params->plotParams->plotEllipses[rLevel];
  plotAny = ( plotPoints || plotTiles || plotEllipses );
  if ( plotEllipses ) {
    getMetric = params->plotParams->params[rLevel].getMetric;
    metricParams = params->plotParams->params[rLevel].metricParams;
    mThresh = params->plotParams->params[rLevel].mThresh;
  }

  /* Move down the list, plotting what needs to be plotted. */
  while ( mesh != NULL ) {

    if ( plotAny ) {
      fprintf( stream, "%f %f moveto", mesh->x, mesh->y );
      if ( plotPoints ) {
	fprintf( stream, " point%u", rLevel );
	AdjustBBox( mesh->x, mesh->y, params->plotParams );
	params->nObj += 1;
      }
      if ( plotTiles ) {
	fprintf( stream, " %f %f %f tile", mesh->dx, mesh->dy[0],
		 mesh->dy[1] );
	AdjustBBox( mesh->x - mesh->dx, mesh->y - mesh->dy[0],
		    params->plotParams );
	AdjustBBox( mesh->x - mesh->dx, mesh->y - mesh->dy[1],
		    params->plotParams );
	AdjustBBox( mesh->x + mesh->dx, mesh->y + mesh->dy[0],
		    params->plotParams );
	AdjustBBox( mesh->x + mesh->dx, mesh->y + mesh->dy[1],
		    params->plotParams );
	params->nObj += 4;
      }
      if ( plotEllipses ) {
	REAL4 position[2];  /* location of current mesh point */
	REAL4 metric[3];    /* value of the metric at position */
	REAL4 axes[2];      /* lengths of principal axes of ellipse */
	REAL4 theta;        /* angle from x-axis to axis[0] */
	REAL4 cost, sint;   /* cosine and sine of theta */
	REAL4 lambda;       /* eigenvalue of metric */
	REAL4 term1, term2; /* temporary variables */

	/* Get metric and angle. */
	position[0] = mesh->x;
	position[1] = mesh->y;
	TRY( (getMetric)( stat->statusPtr, metric, position,
			  metricParams ), stat );
	theta = 0.5*atan2( -2.0*metric[2], metric[1] - metric[0] );
	cost = cos( theta );
	sint = sin( theta );

	/* Get first principal axis. */
	term1 = fabs( metric[0]*cost + metric[2]*sint );
	term2 = fabs( metric[2]*cost + metric[1]*sint );
	if ( term1 > term2 ) {
	  term2 /= term1;
	  lambda = term1*sqrt( 1.0 + term2*term2 );
	} else {
	  term1 /= term2;
	  lambda = term2*sqrt( 1.0 + term1*term1 );
	}
	if ( lambda <= 0.0 ) {
	  ABORT( stat, TWODMESHPLOTH_EMETRIC,
		 TWODMESHPLOTH_MSGEMETRIC );
	}
	axes[0] = sqrt( mThresh / lambda );

	/* Get second principal axis. */
	term1 = fabs( metric[0]*sint - metric[2]*cost );
	term2 = fabs( metric[2]*sint - metric[1]*cost );
	if ( term1 > term2 ) {
	  term2 /= term1;
	  lambda = term1*sqrt( 1.0 + term2*term2 );
	} else {
	  term1 /= term2;
	  lambda = term2*sqrt( 1.0 + term1*term1 );
	}
	if ( lambda <= 0.0 ) {
	  ABORT( stat, TWODMESHPLOTH_EMETRIC,
		 TWODMESHPLOTH_MSGEMETRIC );
	}
	axes[1] = sqrt( mThresh / lambda );

	/* Plot ellipse. */
	fprintf( stream, " %f %f %f ellipse", axes[0], axes[1],
		 theta*(REAL4)( LAL_180_PI ) );
	AdjustBBox( mesh->x + axes[0]*cost - axes[1]*sint,
		    mesh->y + axes[0]*sint + axes[1]*cost,
		    params->plotParams );
	AdjustBBox( mesh->x - axes[0]*cost - axes[1]*sint,
		    mesh->y - axes[0]*sint + axes[1]*cost,
		    params->plotParams );
	AdjustBBox( mesh->x - axes[0]*cost + axes[1]*sint,
		    mesh->y - axes[0]*sint - axes[1]*cost,
		    params->plotParams );
	AdjustBBox( mesh->x + axes[0]*cost + axes[1]*sint,
		    mesh->y + axes[0]*sint - axes[1]*cost,
		    params->plotParams );
	params->nObj += 4;
      }
      fprintf( stream, "\n" );

      /* Start a new macro if necessary. */
      if ( params->nObj > TWODMESHPLOTC_MAXOBJ ) {
	fprintf( stream,
		 "} def\n"
		 "\n/mesh%u {\n", ++( params->nMacro ) );
	params->nObj = 0;
      }
    }

    /* Plot any submeshes, if necessary. */
    if ( ( mesh->subMesh != NULL ) &&
	 ( rLevel < params->plotParams->nLevels - 1 ) ) {
      params->level += 1;
      TRY( LALMakeMeshMacro( stat->statusPtr, stream, mesh->subMesh,
			     params ), stat );
      params->level -= 1;
    }

    /* Move on to next node. */
    mesh = mesh->next;
  }
  /* End of loop over list. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


static void
AdjustBBox( REAL4 x, REAL4 y, TwoDMeshPlotStruc *params )
     /* This routine expands the params->bBox field to include the
        point with x-y coordinates given by position. */
{
  if ( !params->clip ||
       ( ( x > params->clipBox[0] ) && ( x < params->clipBox[2] ) &&
	 ( y > params->clipBox[1] ) && ( y < params->clipBox[3] ) ) ) {
    REAL4 xp = x*params->xScale*params->cosTheta -
      y*params->yScale*params->sinTheta;
    REAL4 yp = x*params->xScale*params->sinTheta +
      y*params->yScale*params->cosTheta;
    if ( params->bBox[0] > xp )
      params->bBox[0] = xp;
    if ( params->bBox[2] < xp )
      params->bBox[2] = xp;
    if ( params->bBox[1] > yp )
      params->bBox[1] = yp;
    if ( params->bBox[3] < yp )
      params->bBox[3] = yp;
  }
  return;
}
