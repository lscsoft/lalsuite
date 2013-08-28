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

/**
 * \author Creighton, T. D.
 * \file
 * \ingroup TwoDMesh_h
 * \brief Creates or destroys a hierarchical mesh of templates on an 2-dimensional parameter space.
 *
 * \heading{Description}
 *
 * The routine <tt>LALCreateTwoDMesh()</tt> lays out an unevenly-spaced
 * mesh on a 2-dimensional parameter space, according to the method
 * presented in \ref TwoDMesh_h and detailed in
 * TwoDMeshInternal.c.  The parameter \c mesh is a handle to
 * the head of the newly-created linked list of mesh points, while
 * \c params points to the parameter structure used to create the
 * list.  On completion, <tt>params->nOut</tt> is set to the number of mesh
 * points created.
 *
 * The routine <tt>LALDestroyTwoDMesh()</tt> destroys the list pointed to
 * by <tt>*mesh</tt>, including all sub-meshes, and sets
 * <tt>*mesh</tt>=\c NULL.  If <tt>*mesh</tt> is already \c NULL,
 * nothing is done (this is \e not an erroneous usage).  If
 * \c nFree\f$\neq\f$\c NULL, then <tt>*nFree</tt> is set to the number
 * of nodes freed.
 *
 * The routine <tt>LALRefineTwoDMesh()</tt> creates a heirarchical search
 * mesh by inserting copies of the nodes in the list pointed to by
 * \c fineMesh into the \c subMesh fields of appropriate nodes in
 * the list pointed to by \c coarseMesh.  The contents of the
 * \c fineMesh list are untouched.  If a \c fineMesh tile does
 * not overlap with any \c cosarseMesh tile, a warning is generated,
 * but this is not treated as an error.  If an internal error does occur,
 * the refinement will be left in a state of partial completion; there is
 * just too much overhead involved in maintaining an uncorrupted copy of
 * the \c coarseMesh list for it to be worthwhile.
 *
 * \heading{Algorithm}
 *
 * \c LALCreateTwoDMesh() simply creates a dummy node to serve as the
 * head of the linked list, and calls <tt>LALTwoDMesh()</tt> in
 * TwoDMeshInternal.c to attach a mesh to it.  The details of the
 * algorithm are given in TwoDMeshInternal.c.
 *
 * <tt>LALDestroyTwoDMesh()</tt> navigates down the linked list of mesh
 * points, destroying them as it goes.  It calls itself recursively on
 * any non-empty sub-meshes to destroy them too.
 *
 * <tt>LALRefineTwoDMesh()</tt> moves along the \c fineMesh list; for
 * each node in the list, it searches the \c coarseMesh list for the
 * any tile that overlaps with the fine mesh tile.  It then \e copies
 * the fine mesh node (and its submesh, if any) into the coarse mesh
 * node's \c subMesh list, using <tt>LALTwoDNodeCopy()</tt> in
 * TwoDMeshInternal.c.  Although it uses more memory, this
 * recursive copy routine is preferred over simple relinking, so as to
 * avoid any possible memory leaks: destroying the coarse mesh list will
 * leave the fine mesh list intact, and vice-versa.
 *
 * To create a \f$>2\f$~level hierarchical search mesh, build it from the
 * bottom up: call <tt>LALRefineTwoDMesh()</tt> to add the finest mesh to
 * the next finest, add that to the next finest, and so on up to the
 * coarsest mesh.
 *
 * \heading{Uses}
 * \code
 * lalDebugLevel               XLALPrintError()
 * LALWarning()                LALInfo()
 * LALTwoDMesh()               LALTwoDNodeCopy()
 * LALFree()
 * \endcode
 *
 * \heading{Notes}
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>

void
LALCreateTwoDMesh( LALStatus          *stat,
		   TwoDMeshNode       **mesh,
		   TwoDMeshParamStruc *params )
{
  TwoDMeshNode head;     /* dummy head node */
  TwoDMeshNode *headPtr; /* pointer to above */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check that input parameters exist, but that the mesh handle
     points to a null pointer.  Parameter fields will be checked
     within the subroutine. */
  ASSERT( mesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( !(*mesh), stat, TWODMESHH_EOUT, TWODMESHH_MSGEOUT );

  /* Ben wants a warning if the widthMaxFac or widthRetryFac are
     larger than is reasonable. */
#ifndef NDEBUG
  if ( lalDebugLevel&LALWARNING ) {
    REAL4 retry = params->widthRetryFac;
    if ( params->widthMaxFac > LAL_SQRT2 )
      LALWarning( stat, "widthMaxFac > sqrt(2)" );
    if ( retry > 1.0 && retry*retry > params->widthMaxFac )
      LALWarning( stat, "widthRetryFac > sqrt(widthMaxFac)" );
  }
#endif

  /* Create the list using LALTwoDMesh(). */
  params->nOut = 0;
  head.next = NULL;
  headPtr = &head;
#ifndef NDEBUG
  if ( lalDebugLevel&LALINFO )
    LALInfo( stat, "Generating mesh\n" );
#endif
  TRY( LALTwoDMesh( stat->statusPtr, &headPtr, params ), stat );
#ifndef NDEBUG
  if ( lalDebugLevel&LALINFO )
    if ( ( params->nIn == 0 ) || ( params->nOut < params->nIn ) ) {
      XLALPrintError( "\n" );
      LALInfo( stat, "Mesh complete" );
      XLALPrintError( "\tnumber of mesh points: %u\n", params->nOut );
    }
#endif

  /* Update the output, and exit. */
  *mesh = head.next;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



void
LALDestroyTwoDMesh( LALStatus    *stat,
		    TwoDMeshNode **mesh,
		    UINT4        *nFree )
{
  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check that all parameters exist. */
  ASSERT( mesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  if ( nFree )
    *nFree = 0;

  /* Free everything, recursively freeing sub-meshes if necessary. */
  while ( *mesh ) {
    UINT4 nSub = 0;             /* nodes freed from sub-meshes */
    TwoDMeshNode *last = *mesh; /* pointer to previous node */
    if ( last->subMesh ) {
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &(last->subMesh),
			       &nSub ), stat );
    }
    *mesh = last->next;
    LALFree( last );
    if ( nFree )
      *nFree += nSub + 1;
  }

  /* If we got here without sigsegving, we're done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



void
LALRefineTwoDMesh( LALStatus    *stat,
		   TwoDMeshNode *coarseMesh,
		   TwoDMeshNode *fineMesh )
{
  BOOLEAN found;      /* whether a fine point is in any coarse tile */
  UINT4 lost = 0;     /* number of fine points not found */
  TwoDMeshNode *here; /* pointer to coarse mesh list */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  ASSERT( coarseMesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( fineMesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Scan through fine mesh.  We never return to the head of the fine
     mesh, so we can directly increment the parameter *finemesh. */
  for ( ; fineMesh != NULL; fineMesh = fineMesh->next ) {
    /* Centre, width, height, and slope of fine tile */
    REAL4 xFine = fineMesh->x;
    REAL4 yFine = fineMesh->y;
    REAL4 dxFine = fineMesh->dx;
    REAL4 dyFine = 0.5*( fineMesh->dy[1] - fineMesh->dy[0] );
    REAL4 mFine = 0.5*( fineMesh->dy[1] + fineMesh->dy[0] )/dxFine;

    /* For each fine mesh tile, scan through coarse mesh, and look for
       ones that overlap in their x domains. */
    found = 0;
    for ( here = coarseMesh; here != NULL; here = here->next ) {
      REAL4 x = here->x - xFine;
      if ( fabs( x ) <= dxFine + here->dx ) {

	/* Transform the coarse mesh tile into sheared coordinates
           centred on the fine mesh tile. */
	REAL4 y = here->y + mFine*x - yFine;
	REAL4 dy0 = here->dy[0] + mFine*here->dx;
	REAL4 dy1 = here->dy[1] + mFine*here->dx;

	/* See if there is any possibility of overlap. */
	if ( ( fabs( y ) <= dyFine + fabs( dy0 ) ) ||
	     ( fabs( y ) <= dyFine + fabs( dy1 ) ) ) {

	  /* We check for overlap on the left and right sides of the
             common domain of the two tiles.  On either side, the
             coarse tile can be either completely below the fine tile
             (-1), completely above the fine tile (+1), or overlapping
             it (0).  We store this information in two INT2's,
             below. */
	  INT2 overlap[2] = { 0, 0 };

	  /* Compute height and slope of coarse tile in the sheared
             coordinates of the fine mesh tile. */
	  REAL4 dy = 0.5*( dy1 - dy0 );
	  REAL4 m = 0.5*( dy1 + dy0 );

	  /* Find leftmost point of overlap of the two tiles relative
             to the coarse mesh, and test the range of the coarse-mesh
             tile at that point. */
	  REAL4 xOver = -here->dx;
	  if ( xOver < -x - dxFine )
	    xOver = -x - dxFine;
	  if ( -dy + m*xOver <= dyFine ) {
	    if ( dy + m*xOver >= -dyFine )
	      overlap[0] = 0;
	    else
	      overlap[0] = -1;
	  } else
	    overlap[0] = 1;

	  /* Find rightmost point of overlap of the two tiles relative
             to the coarse mesh, and test the range of the coarse-mesh
             tile at that point. */
	  if ( overlap[0] ) {
	    xOver = here->dx;
	    if ( xOver > -x + dxFine )
	      xOver = -x + dxFine;
	    if ( -dy + m*xOver <= dyFine ) {
	      if ( dy + m*xOver >= -dyFine )
		overlap[1] = 0;
	      else
		overlap[1] = -1;
	    } else
	      overlap[1] = 1;
	  }

	  /* The two tiles overlap if either side has a value of 0 or
             if the two sides have opposite sign. */
	  overlap[0] *= overlap[1];
	  if ( !overlap[0] ) {

	    /* This is it!  Copy the fine mesh node and add it into
               the coarse submesh. */
	    TwoDMeshNode *copy = NULL;
	    TRY( LALTwoDNodeCopy( stat->statusPtr, &copy, fineMesh ),
		 stat );
	    copy->next = coarseMesh->subMesh;
	    coarseMesh->subMesh = copy;
	    found = 1;
	  }
	}
      }
    }
    /* If no coarse tile overlapped, make a note of this. */
#ifndef NDEBUG
    if ( !found ) {
      lost++;
      if ( lalDebugLevel&LALINFO ) {
	LALInfo( stat, "Fine mesh tile has no overlapping coarse tile" );
	XLALPrintError( "\tlocation: (%f,%f)\n", xFine, yFine );
      }
    }
#endif
  }
  /* If any fine mesh tiles were lost, warn the user. */
#ifndef NDEBUG
  if ( lalDebugLevel&LALWARNING )
    if ( lost > 0 ) {
      LALWarning( stat, "Some fine mesh tiles were lost" );
      XLALPrintError( "\tnumber lost = %u\n", lost );
    }
#endif

  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
