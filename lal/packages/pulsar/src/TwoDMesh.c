/************************************* <lalVerbatim file="TwoDMeshCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TwoDMesh.c}}
\label{ss:TwoDMesh.c}

Creates or destroys a hierarchical mesh of templates on an
2-dimensional parameter space.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TwoDMeshCP}
\index{\texttt{LALCreateTwoDMesh()}}
\index{\texttt{LALDestroyTwoDMesh()}}
\index{\texttt{LALRefineTwoDMesh()}}

\subsubsection*{Description}

The routine \verb@LALCreateTwoDMesh()@ lays out an unevenly-spaced
mesh on a 2-dimensional parameter space, according to the method
presented in \verb@TwoDMesh.h@ and detailed in
\verb@TwoDMeshInternal.c@.  The parameter \verb@mesh@ is a handle to
the head of the newly-created linked list of mesh points, while
\verb@params@ points to the parameter structure used to create the
list.  On completion, \verb@params->nOut@ is set to the number of mesh
points created.

The routine \verb@LALDestroyTwoDMesh()@ destroys the list pointed to
by \verb@*mesh@, including all sub-meshes, and sets
\verb@*mesh@=\verb@NULL@.  If \verb@*mesh@ is already \verb@NULL@,
nothing is done (this is \emph{not} an erroneous usage).  If
\verb@nFree@$\neq$\verb@NULL@, then \verb@*nFree@ is set to the number
of nodes freed.

The routine \verb@LALRefineTwoDMesh()@ creates a heirarchical search
mesh by inserting copies of the nodes in the list pointed to by
\verb@fineMesh@ into the \verb@subMesh@ fields of appropriate nodes in
the list pointed to by \verb@coarseMesh@.  The contents of the
\verb@fineMesh@ list are untouched.  If a \verb@fineMesh@ tile does
not overlap with any \verb@cosarseMesh@ tile, a warning is generated,
but this is not treated as an error.  If an internal error does occur,
the refinement will be left in a state of partial completion; there is
just too much overhead involved in maintaining an uncorrupted copy of
the \verb@coarseMesh@ list for it to be worthwhile.

\subsubsection*{Algorithm}

\verb@LALCreateTwoDMesh@ simply creates a dummy node to serve as the
head of the linked list, and calls \verb@LALTwoDMesh()@ in
\verb@TwoDMeshInternal.c@ to attach a mesh to it.  The details of the
algorithm are given in \verb@TwoDMeshInternal.c@.

\verb@LALDestroyTwoDMesh()@ navigates down the linked list of mesh
points, destroying them as it goes.  It calls itself recursively on
any non-empty sub-meshes to destroy them too.

\verb@LALRefineTwoDMesh()@ moves along the \verb@fineMesh@ list; for
each node in the list, it searches the \verb@coarseMesh@ list for the
any tile that overlaps with the fine mesh tile.  It then \emph{copies}
the fine mesh node (and its submesh, if any) into the coarse mesh
node's \verb@subMesh@ list, using \verb@LALTwoDNodeCopy()@ in
\verb@TwoDMeshInternal.c@.  Although it uses more memory, this
recursive copy routine is preferred over simple relinking, so as to
avoid any possible memory leaks: destroying the coarse mesh list will
leave the fine mesh list intact, and vice-versa.

To create a $>2$~level hierarchical search mesh, build it from the
bottom up: call \verb@LALRefineTwoDMesh()@ to add the finest mesh to
the next finest, add that to the next finest, and so on up to the
coarsest mesh.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel               LALPrintError()
LALWarning()                LALInfo()
LALTwoDMesh()               LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoDMeshCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>

NRCSID( TWODMESHC, "$Id$" );

/* <lalVerbatim file="TwoDMeshCP"> */
void
LALCreateTwoDMesh( LALStatus          *stat,
		   TwoDMeshNode       **mesh,
		   TwoDMeshParamStruc *params )
{ /* </lalVerbatim> */
  TwoDMeshNode head;     /* dummy head node */
  TwoDMeshNode *headPtr; /* pointer to above */

  INITSTATUS( stat, "LALCreateTwoDMesh", TWODMESHC );
  ATTATCHSTATUSPTR( stat );

  /* Check that input parameters exist, but that the mesh handle
     points to a null pointer.  Parameter fields will be checked
     within the subroutine. */
  ASSERT( mesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( !(*mesh), stat, TWODMESHH_EOUT, TWODMESHH_MSGEOUT );

  /* Create the list using LALTwoDMesh(). */
  params->nOut = 0;
  head.next = NULL;
  headPtr = &head;
  TRY( LALTwoDMesh( stat->statusPtr, &headPtr, params ), stat );
#ifndef NDEBUG
  if ( lalDebugLevel&LALINFO )
    if ( ( params->nIn == 0 ) || ( params->nOut < params->nIn ) ) {
      LALInfo( stat, "Mesh complete" );
      LALPrintError( "\tnumber of mesh points: %u\n", params->nOut );
    }
#endif

  /* Update the output, and exit. */
  *mesh = head.next;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TwoDMeshCP"> */
void
LALDestroyTwoDMesh( LALStatus    *stat,
		    TwoDMeshNode **mesh,
		    UINT4        *nFree )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALDestroyTwoDMesh", TWODMESHC );
  ATTATCHSTATUSPTR( stat );

  /* Check that all parameters exist. */
  ASSERT( mesh, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  if ( nFree )
    *nFree = 0;

  /* Free everything, recursively freeing sub-meshes if necessary. */
  while ( *mesh ) {
    UINT4 nSub = 0;             /* nodes freed from sub-meshes */
    TwoDMeshNode *last = *mesh; /* pointer to previous node */
    if ( last->subMesh )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &(last->subMesh),
			       &nSub ), stat );
    *mesh = last->next;
    LALFree( last );
    if ( nFree )
      *nFree -= nSub + 1;
  }

  /* If we got here without sigsegving, we're done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TwoDMeshCP"> */
void
LALRefineTwoDMesh( LALStatus    *stat,
		   TwoDMeshNode *coarseMesh,
		   TwoDMeshNode *fineMesh )
{ /* </lalVerbatim> */
  BOOLEAN found;      /* whether a fine point is in any coarse tile */
  UINT4 lost = 0;     /* number of fine points not found */
  TwoDMeshNode *here; /* pointer to coarse mesh list */

  INITSTATUS( stat, "LALDestroyTwoDMesh", TWODMESHC );
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
	  INT2 overlap[0] = { 0, 0 };

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
	LALPrintError( "\tlocation: (%f,%f)\n", xFine, yFine );
      }
    }
#endif
  }
  /* If any fine mesh tiles were lost, warn the user. */
#ifndef NDEBUG
  if ( lalDebugLevel&LALWARNING )
    if ( lost > 0 ) {
      LALWarning( stat, "Some fine mesh tiles were lost" );
      LALPrintError( "\tnumber lost = %u\n", lost );
    }
#endif

  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
