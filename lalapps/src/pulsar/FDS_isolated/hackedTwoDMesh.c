/***************************** <lalVerbatim file="TwoDMeshInternalCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>
#include "hackedTwoDMesh.h"

NRCSID( HTWODMESHC, "$Id$" );

/* Whether or not to track progress internally. */
static UINT4 columnNo;

/* Local constants. */
#define TWODMESHINTERNALC_WMAXFAC   (1.189207115)
#define TWODMESHINTERNALC_WRETRYFAC LAL_SQRT2

/* Local macros.  These macros replace repeated code blocks in the
routines, and operate within the scope of variables declared within
those routines.  Required variables are documented below. */


/* This macro terminates the current routine if the column is found to
be too wide.  It frees the linked list, decrements the node count,
sets the output flag, prints as an informational message the node
count where the error occured, and returns from the current routine.
In addition to the passed parameters, the folowing external variables
are required:

TwoDMeshNode **tail: The list (*tail)->next is destroyed.

LALStatus *stat: Used in reporting the error.  stat->statusPtr is also
                 used in destroying the list.

INT4 lalDebugLevel: Used in reporting the error.

TwoDColumnParamStruc *column: The flag column->tooWide is set.

TwoDMeshParamStruc *params: The field params->nOut is decremented. */

#define TOOWIDERETURN                                                \
do {                                                                 \
  UINT4 nFree;                                                       \
  if ( lalDebugLevel&LALINFO ) {                                     \
    LALInfo( stat, "Column too wide" );                              \
    LALPrintError( "\tnode count %u\n", params->nOut );              \
  }                                                                  \
  TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),        \
			   &nFree ), stat );                         \
  params->nOut -= nFree;                                             \
  column->tooWide = 1;                                               \
  DETATCHSTATUSPTR( stat );                                          \
  RETURN( stat );                                                    \
} while (0)


/* This macro computes the maximum half-width of a mismatch ellipse
given the metric and a mismatch value.  If the metric is not
positive-definite, an error returned.  In addition to the passed
parameters, the folowing external variables are required:

TwoDMeshNode **tail: The list (*tail)->next is destroyed on an error.

LALStatus *stat: Used in reporting an error.  stat->statusPtr is also
                 used in destroying the list.  */

#define GETWIDTH( dx, metric, mismatch )                             \
do {                                                                 \
  REAL8 det = (REAL8)(metric)[0]*(REAL8)(metric)[1] - (REAL8)(metric)[2]*(REAL8)(metric)[2];     \
  if ( ( (metric)[0] <= 0.0 ) || ( (metric)[1] <= 0.0 ) ||           \
       ( det <= 0.0 ) ) {                                            \
    TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*(tail))->next),    \
			     NULL ), stat );                         \
    ABORT( stat, TWODMESHH_EMETRIC, TWODMESHH_MSGEMETRIC );          \
  }                                                                  \
  (dx) = (REAL8) sqrt( (REAL8)(metric)[1]*(mismatch)/det );                         \
} while (0)


/* This macro computes the positions of the right-hand corners of a
tile given a tile half width, the metric, and a mismatch value.  If
the metric is not positive-definite, then an error is returned.  If
the ellipse is not sufficiently wider than the requested width, then a
flag is set and the current subroutine will return.  In addition to
the passed parameters, the folowing external variables are required:

TwoDMeshNode **tail: The list (*tail)->next is destroyed on an error.

LALStatus *stat: Used in reporting an error.  stat->statusPtr is
                     also used in destroying the list.

REAL4 widthMaxFac: The factor by which the maximum ellipse half-width
                   must exceed the given column half-width. */

#define GETSIZE( dy, dx, metric, mismatch )                          \
do {                                                                 \
  REAL8 det = (REAL8)(metric)[0]*(metric)[1] - (REAL8)(metric)[2]*(metric)[2];     \
  REAL8 disc;                                                        \
  if ( ( metric[0] <= 0.0 ) || ( metric[1] <= 0.0 ) ||               \
       ( det <= 0.0 ) ) {                                            \
    TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*(tail))->next),    \
			     NULL ), stat );                         \
    ABORT( stat, TWODMESHH_EMETRIC, TWODMESHH_MSGEMETRIC );          \
  }                                                                  \
  if ( widthMaxFac*(dx) > (REAL8)sqrt( (REAL8)(metric)[1]*(mismatch)/det ) )       \
    TOOWIDERETURN;                                                   \
  disc = sqrt( (REAL8)(metric)[1]*(mismatch) - det*(dx)*(dx) );             \
  (dy)[0] = (REAL8)( - (REAL8)metric[2]*dx - disc ) / metric[1];                    \
  (dy)[1] = (REAL8)( - (REAL8)metric[2]*dx + disc ) / metric[1];                    \
} while (0)


/* <lalVerbatim file="TwoDMeshCP"> */
void
hackedLALCreateTwoDMesh( LALStatus          *stat,
			 TwoDMeshNode       **mesh,
			 TwoDMeshParamStruc *params )
{ /* </lalVerbatim> */
  TwoDMeshNode head;     /* dummy head node */
  TwoDMeshNode *headPtr; /* pointer to above */

  INITSTATUS( stat, "LALCreateTwoDMesh", HTWODMESHC );
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
    REAL8 retry = params->widthRetryFac;
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
  TRY( hackedLALTwoDMesh( stat->statusPtr, &headPtr, params ), stat );
#ifndef NDEBUG
  if ( lalDebugLevel&LALINFO )
    if ( ( params->nIn == 0 ) || ( params->nOut < params->nIn ) ) {
      LALPrintError( "\n" );
      LALInfo( stat, "Mesh complete" );
      LALPrintError( "\tnumber of mesh points: %u\n", params->nOut );
    }
#endif

  /* Update the output, and exit. */
  *mesh = head.next;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );

} /* hackedLALCreateTwoDMesh() */



/* <lalVerbatim file="TwoDMeshInternalCP"> */
void
hackedLALTwoDMesh( LALStatus          *stat,
		   TwoDMeshNode       **tail,
		   TwoDMeshParamStruc *params )
{ /* </lalVerbatim> */
  TwoDColumnParamStruc column; /* parameters for current column */
  TwoDMeshNode *here;          /* current tail of linked list */

  /* Default parameter values: */
  REAL8 widthRetryFac = TWODMESHINTERNALC_WRETRYFAC;
  REAL8 maxColumnFac = 0.0;
  UINT4 nIn = (UINT4)( -1 );

  INITSTATUS( stat, "LALTwoDMesh", HTWODMESHC );
  ATTATCHSTATUSPTR( stat );

  /* Check that all parameters exist. */
  ASSERT( tail, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( *tail, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params->getRange, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params->getMetric, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Check that **tail really is the tail. */
  ASSERT( !( (*tail)->next ), stat, TWODMESHH_EOUT, TWODMESHH_MSGEOUT );

  /* Reassign default parameter values if necessary. */
  if ( params->widthRetryFac > 1.0 )
    widthRetryFac = params->widthRetryFac;
  if ( params->maxColumns > 0 )
    maxColumnFac = 1.0 / params->maxColumns;
  if ( params->nIn > 0 )
    nIn = params->nIn;

  /* Set the clipping area to something irrelevant. */
  column.leftClip[1] = column.rightClip[1] = 0.9*LAL_REAL4_MAX;
  column.leftClip[0] = column.rightClip[0] = -column.leftClip[1];

  /* Locate the first column's right-hand edge. */
  column.domain[0] = params->domain[0];
  TRY( (params->getRange)( stat->statusPtr, column.leftRange, column.domain[0], params->rangeParams ),
       stat );

  if (lalDebugLevel >= 3)
    {
      columnNo = 0;
      LALPrintError( "      Node count    Column count\n" );
    }

  /* Main loop: add columns until we're past the end of the space. */
  here = *tail;
  while ( column.domain[0] < params->domain[1] ) {
    REAL4 position[2]; /* position in parameter space */
    REAL4 metric[3];   /* components of metric at position */
    REAL8 w1, w2;      /* bottom and top widths of column */

    /* Estimate column width. */
    position[0] = column.domain[0];
    position[1] = column.leftRange[0];
    (params->getMetric)( stat->statusPtr, metric, position, params->metricParams );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
    ENDFAIL( stat );
    GETWIDTH( w1, metric, params->mThresh );
    position[1] = column.leftRange[1];
    (params->getMetric)( stat->statusPtr, metric, position, params->metricParams );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
    ENDFAIL( stat );
    GETWIDTH( w2, metric, params->mThresh );
    if ( w2 < w1 )
      w1 = w2;
    w1 *= LAL_SQRT2;

    /* Loop to try successively smaller column widths. */
    do {
      /* Make sure width is not too small or too big. */
      if ( maxColumnFac*( params->domain[1] - params->domain[0] ) > w1 ) {
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
	ABORT( stat, TWODMESHH_EWIDTH, TWODMESHH_MSGEWIDTH );
      }
      column.domain[1] = (REAL8)( (REAL8)column.domain[0] + (REAL8)w1);
      if ( column.domain[1] > params->domain[1] ) {
	column.domain[1] = params->domain[1];
	w1 = (REAL8)((REAL8)column.domain[1] - (REAL8)column.domain[0]);
      }

      /* Set remaining column parameters. */
      (params->getRange)( stat->statusPtr, column.rightRange,
			  column.domain[1], params->rangeParams );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
      ENDFAIL( stat );
      column.tooWide = 0;

      /* Call LALTwoDColumn() to place the column. */
      hackedLALTwoDColumn( stat->statusPtr, &here, &column, params );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
      ENDFAIL( stat );

      /* See if we've reached the maximum number of mesh points. */
      if ( params->nOut >= nIn ) {
	*tail = here;
	DETATCHSTATUSPTR( stat );
	RETURN( stat );
      }

      /* If necessary, repeat with a narrower column. */
      w1 = (REAL8)( (REAL8) w1 / widthRetryFac );
    } while ( column.tooWide );

    /* Otherwise, go on to the next column. */
    column.domain[0] = column.domain[1];
    column.leftRange[0] = column.rightRange[0];
    column.leftRange[1] = column.rightRange[1];
    if (lalDebugLevel >= 3) 
      {
	LALPrintError( "\r%16u%16u", params->nOut, columnNo++ );
      }
  }

  /* We're done.  Update the *tail pointer and exit. */
  if (lalDebugLevel >= 3) 
    {
      LALPrintError( "\n" );
    }

  *tail = here;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TwoDMeshInternalCP"> */
void
hackedLALTwoDColumn( LALStatus            *stat,
	       TwoDMeshNode         **tail,
	       TwoDColumnParamStruc *column,
	       TwoDMeshParamStruc   *params )
{ /* </lalVerbatim> */
  BOOLEAN tiled = 0;    /* whether tiles were placed on the centreline */
  REAL4 position[2];    /* current top of column */
  REAL8 dx;             /* half-width of column */
  REAL8 myy0, myy1;         /* temporary variables storing y-coordinates */
  REAL4 centreRange[2]; /* centreline of column parameter space */
  REAL8 centreClip[2];  /* centre of clip boundary */
  REAL8 leftTiled[2];   /* left side of region tiled */
  REAL8 rightTiled[2];  /* right side of region tiled */
  REAL4 metric[3];      /* current metric components */
  TwoDMeshNode *here;   /* current node in list */

  /* Default parameter values: */
  REAL8 widthMaxFac = TWODMESHINTERNALC_WMAXFAC;
  UINT4 nIn = (UINT4)( -1 );

  INITSTATUS( stat, "LALTwoDColumn", HTWODMESHC );
  ATTATCHSTATUSPTR( stat );

  /* Check that all parameters exist. */
  ASSERT( tail, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( *tail, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( column, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params->getRange, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params->getMetric, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Check that **tail really is the tail. */
  ASSERT( !( (*tail)->next ), stat, TWODMESHH_EOUT, TWODMESHH_MSGEOUT );
  here = *tail;

  /* Reassign default parameter values if necessary. */
  if ( params->widthMaxFac > 1.0 )
    widthMaxFac = params->widthMaxFac;
  if ( params->nIn > 0 )
    nIn = params->nIn;

  /* Set the boundaries of the regions that no longer need tiling. */
  centreClip[0] = (REAL8)0.5*column->leftClip[0] + (REAL8)0.5*column->rightClip[0];
  centreClip[1] = (REAL8)0.5*column->leftClip[1] + (REAL8)0.5*column->rightClip[1];
  leftTiled[0] = column->leftClip[1];
  leftTiled[1] = column->leftClip[0];
  rightTiled[0] = column->rightClip[1];
  rightTiled[1] = column->rightClip[0];

  /* Get the width and heights of this column. */
  position[0] =  0.5*( (REAL8)column->domain[1] + column->domain[0] );
  dx = 0.5*( (REAL8)column->domain[1] - column->domain[0] );
  TRY( (params->getRange)( stat->statusPtr, centreRange, position[0], params->rangeParams ), stat );

  /* Add the column of tiles along the centreline, if the parameter
     space intersects the clipping area along the centreline. */
  position[1] = centreClip[0];
  if ( position[1] < centreRange[0] )
    position[1] = centreRange[0];
  if ( position[1] <= centreRange[1] ) {

    /* Add base tile of column. */
    tiled = 1;
    TRY( (params->getMetric)( stat->statusPtr, metric, position, params->metricParams ), stat );
    here->next = (TwoDMeshNode *)LALMalloc( sizeof(TwoDMeshNode) );
    if ( here == NULL ) {
      ABORT( stat, TWODMESHH_EMEM, TWODMESHH_MSGEMEM );
    }
    memset( here->next, 0, sizeof(TwoDMeshNode) );
    params->nOut++;
    if (lalDebugLevel >= 3) 
      {
	LALPrintError( "\r%16u", params->nOut );
      }
    GETSIZE( here->next->dy, dx, metric, params->mThresh );
    here->next->y = position[1];
    here = here->next;
    here->x = position[0];
    here->dx = dx;
    here->next = here->subMesh = NULL;
    if ( params->nOut >= nIn ) {
      *tail = here;
      DETATCHSTATUSPTR( stat );
      RETURN( stat );
    }

    /* Determine the region that we've covered. */
    myy0 = (REAL8)here->y + here->dy[0];
    myy1 = (REAL8)here->y - here->dy[1];
    if ( leftTiled[0] > myy1 )
      leftTiled[0] = myy1;
    if ( rightTiled[0] > myy0 )
      rightTiled[0] = myy0;
    leftTiled[1] = (REAL8)here->y - here->dy[0];
    rightTiled[1] = (REAL8)here->y + here->dy[1];
    position[1] = 0.5*leftTiled[1] + 0.5*rightTiled[1];

    /* Continue stacking tiles until we reach the top. */
    while ( ( position[1] < centreRange[1] ) && ( position[1] < centreClip[1] ) ) {
      (params->getMetric)( stat->statusPtr, metric, position, params->metricParams );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
      ENDFAIL( stat );
      here->next = (TwoDMeshNode *)LALMalloc( sizeof(TwoDMeshNode) );
      if ( here == NULL ) {
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
	ABORT( stat, TWODMESHH_EMEM, TWODMESHH_MSGEMEM );
      }
      memset( here->next, 0, sizeof(TwoDMeshNode) );
      params->nOut++;
      if (lalDebugLevel >= 3) 
	{
	  LALPrintError( "\r%16u", params->nOut );
	}
      GETSIZE( here->next->dy, dx, metric, params->mThresh );
      myy0 = (REAL8)here->dy[1] - here->next->dy[0];
      myy1 = (REAL8)here->next->dy[1] - here->dy[0];
      if ( myy0 > myy1 )
	myy0 = myy1;
      if ( myy0 <= 0.0 )
	TOOWIDERETURN;
      here->next->y = (REAL8)here->y + myy0;
      here = here->next;
      if ( here->y > centreRange[1] )
	here->y = centreRange[1];
      here->x = position[0];
      here->dx = dx;
      here->next = here->subMesh = NULL;
      if ( params->nOut >= nIn ) {
	*tail = here;
	DETATCHSTATUSPTR( stat );
	RETURN( stat );
      }

      /* Extend the covered region upwards. */
      leftTiled[1] = (REAL8)here->y - here->dy[0];
      rightTiled[1] = (REAL8)here->y + here->dy[1];
      position[1] = 0.5*leftTiled[1] + 0.5*rightTiled[1];
    }
  }

  /* Centreline stacking is complete.  Now check for exposed corners
     of the parameter space, and call LALTwoDColumn() recursively. */

  /* Check bottom corners. */
  myy0 = 0.5*leftTiled[0] + 0.5*rightTiled[0];

  /* Bottom-left: */
  if ( ( ( column->leftClip[0] < leftTiled[0] ) || ( centreClip[0] < myy0 ) ) &&
       ( column->leftRange[0] < leftTiled[0] ) &&( ( column->leftRange[1] > column->leftClip[0] ) ||
	 ( centreRange[1] > centreClip[0] ) ) ) {
    TwoDColumnParamStruc column2;
    column2.domain[0] = column->domain[0];
    column2.domain[1] = position[0];
    memcpy( column2.leftRange, column->leftRange, 2*sizeof(REAL4) );
    memcpy( column2.leftClip, column->leftClip, 2*sizeof(REAL4) );
    memcpy( column2.rightRange, centreRange, 2*sizeof(REAL4) );
    memcpy( column2.rightClip, centreClip, 2*sizeof(REAL4) );
    if ( ( leftTiled[0] < column2.leftClip[1] ) && ( myy0 < column2.rightClip[1] ) ) {
      column2.leftClip[1] = leftTiled[0];
      column2.rightClip[1] = myy0;
    }
    LALTwoDColumn( stat->statusPtr, &here, &column2, params );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
    ENDFAIL( stat );
    if ( params->nOut >= nIn ) {
      *tail = here;
      DETATCHSTATUSPTR( stat );
      RETURN( stat );
    }
    if ( column2.tooWide )
      TOOWIDERETURN;
  }

  /* Bottom-right: */
  if ( ( ( column->rightClip[0] < rightTiled[0] ) || ( centreClip[0] < myy0 ) ) &&
       ( column->rightRange[0] < rightTiled[0] ) && ( ( column->rightRange[1] > column->rightClip[0] ) ||
	 ( centreRange[1] > centreClip[0] ) ) ) {
    TwoDColumnParamStruc column2;
    column2.domain[1] = column->domain[1];
    column2.domain[0] = position[0];
    memcpy( column2.rightRange, column->rightRange, 2*sizeof(REAL4) );
    memcpy( column2.rightClip, column->rightClip, 2*sizeof(REAL4) );
    memcpy( column2.leftRange, centreRange, 2*sizeof(REAL4) );
    memcpy( column2.leftClip, centreClip, 2*sizeof(REAL4) );
    if ( ( rightTiled[0] < column2.rightClip[1] ) && ( myy0 < column2.leftClip[1] ) ) {
      column2.rightClip[1] = rightTiled[0];
      column2.leftClip[1] = myy0;
    }
    LALTwoDColumn( stat->statusPtr, &here, &column2, params );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
    ENDFAIL( stat );
    if ( params->nOut >= nIn ) {
      *tail = here;
      DETATCHSTATUSPTR( stat );
      RETURN( stat );
    }
    if ( column2.tooWide )
      TOOWIDERETURN;
  }

  /* Check top corners. */
  if ( tiled ) {
    myy0 = 0.5*leftTiled[1] + 0.5*rightTiled[1];

    /* Top-left: */
    if ( ( ( column->leftClip[1] > leftTiled[1] ) ||
	   ( centreClip[1] > myy0 ) ) &&
	 ( column->leftRange[1] > leftTiled[1] ) &&
	 ( ( column->leftRange[0] < column->leftClip[1] ) ||
	   ( centreRange[0] < centreClip[1] ) ) ) {
      TwoDColumnParamStruc column2;
      column2.domain[0] = column->domain[0];
      column2.domain[1] = position[0];
      memcpy( column2.leftRange, column->leftRange, 2*sizeof(REAL4) );
      memcpy( column2.leftClip, column->leftClip, 2*sizeof(REAL4) );
      memcpy( column2.rightRange, centreRange, 2*sizeof(REAL4) );
      memcpy( column2.rightClip, centreClip, 2*sizeof(REAL4) );
      if ( ( leftTiled[1] > column2.leftClip[0] ) &&
	   ( myy0 > column2.rightClip[0] ) ) {
	column2.leftClip[0] = leftTiled[1];
	column2.rightClip[0] = myy0;
      }
      LALTwoDColumn( stat->statusPtr, &here, &column2, params );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
      ENDFAIL( stat );
      if ( params->nOut >= nIn ) {
	*tail = here;
	DETATCHSTATUSPTR( stat );
	RETURN( stat );
      }
      if ( column2.tooWide )
	TOOWIDERETURN;
    }

    /* Top-right: */
    if ( ( ( column->rightClip[1] > rightTiled[1] ) ||
	   ( centreClip[1] > myy0 ) ) &&
	 ( column->rightRange[1] > rightTiled[1] ) &&
	 ( ( column->rightRange[0] < column->rightClip[1] ) ||
	   ( centreRange[0] < centreClip[1] ) ) ) {
      TwoDColumnParamStruc column2;
      column2.domain[1] = column->domain[1];
      column2.domain[0] = position[0];
      memcpy( column2.rightRange, column->rightRange, 2*sizeof(REAL4) );
      memcpy( column2.rightClip, column->rightClip, 2*sizeof(REAL4) );
      memcpy( column2.leftRange, centreRange, 2*sizeof(REAL4) );
      memcpy( column2.leftClip, centreClip, 2*sizeof(REAL4) );
      if ( ( rightTiled[1] > column2.rightClip[0] ) && ( myy0 > column2.leftClip[0] ) ) {
	column2.rightClip[0] = rightTiled[1];
	column2.leftClip[0] = myy0;
      }
      LALTwoDColumn( stat->statusPtr, &here, &column2, params );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next), NULL ), stat );
      ENDFAIL( stat );
      if ( params->nOut >= nIn ) {
	*tail = here;
	DETATCHSTATUSPTR( stat );
	RETURN( stat );
      }
      if ( column2.tooWide )
	TOOWIDERETURN;
    }
  }

  /* Everything worked fine, so update *tail and exit. */
  *tail = here;
  column->tooWide = 0;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

