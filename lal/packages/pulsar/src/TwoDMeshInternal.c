/***************************** <lalVerbatim file="TwoDMeshInternalCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TwoDMeshInternal.c}}
\label{ss:TwoDMeshInternal.c}

Low-level routines to place a mesh of templates on an 2-dimensional
parameter space.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TwoDMeshInternalCP}
\idx{LALTwoDMesh()}
\idx{LALTwoDColumn()}

\subsubsection*{Description}

These are low-level ``internal'' routines called by
\verb@LALCreateTwoDMesh()@ to lay out an unevenly-spaced mesh on a
2-dimensional parameter space, according to the method discussed in
\verb@TwoDMesh.h@.  They are made globally available to allow greater
control to users attempting to tile complicated parameter spaces.

\verb@LALTwoDMesh()@ places a mesh on the parameter space specified by
\verb@*params@.  On successful completion, the linked list of mesh
points is attached to \verb@(*tail)->next@ (which must previously have
been \verb@NULL@), updates \verb@*tail@ to point to the new tail of
the list, and increases \verb@params->nOut@ by the number of mesh
points added.  (This is useful for tiling several parameter regions
independently with successive calls to \verb@LALTwoDMesh()@.)  On an
error, \verb@**tail@ is left unchanged, and \verb@params->nOut@
indicates where the error occurred.

\verb@LALTwoDColumn()@ places a single column of such a mesh,
according to the additional column restrictions in \verb@*column@.
Again, on success, the mesh points are attached to
\verb@(*tail)->next@, \verb@*tail@ is updated, and \verb@params->nOut@
increased.  If the column specified by \verb@*column@ is deemed to be
too wide for a single column of templates, then \verb@column->tooWide@
is set to 1, \verb@*tail@ and \verb@params->nOut@ are \emph{not}
updated, and the function returns normally (i.e.\ not with an error
code set).  Other more fatal errors are treated as for
\verb@LALTwoDMesh()@, above.

\verb@LALTwoDNodeCopy()@ creates a copy of the node \verb@*old@ and
points \verb@*new@ to the copy.  If \verb@old->subMesh@ exists, each
node in the submesh is copied into its corresponding place by
recursive calls to \verb@LALTwoDNodeCopy()@.  On an error, the copy is
destroyed and \verb@*new@ left unchanged.

\subsubsection*{Algorithm}

\paragraph{\texttt{LALTwoDMesh()}:} This routine starts placing mesh
points at the left side of the parameter region, $x=x_\mathrm{min}$.
It calls \verb@params->getRange()@ to get the bottom and top of the
left edge of the first column.  It also calls \verb@params->getMetric@
at these two corners, estimates the optimum width for the first
column, and uses \verb@params->getRange()@ again to get the two
corners of the right edge of the column.  It then calls the subroutine
\verb@LALTwoDColumn()@ (below) to place the mesh points in this
column.

If \verb@LALTwoDColumn()@ reports that the estimated column width was
too large, \verb@LALCreateTwoDMesh()@ tries again with the width
reduced by the factor \verb@params->widthRetryFac@.  This continues
until the estimated number of columns exceeds
\verb@params->maxColumns@; i.e.\ until the current column width is
less than $(x_\mathrm{max}-x_\mathrm{min})/$\verb@params->maxColumns@.
If this occurs, the linked list is destroyed using
\verb@LALDestroyTwoDMesh()@, and an error is returned.

Otherwise, if \verb@LALTwoDColumn()@ returns success (and does not
complain about the column width), \verb@LALCreateTwoDMesh()@ gets the
width and heights of the next column, and calls \verb@LALTwoDColumn()@
again.  This continues until eventually the right edge of a column
lies beyond $x_\mathrm{max}$.  This last column is squeezed so that
its right edge lies exactly at $x_\mathrm{max}$; once it is filled,
the mesh is deemed complete, and no further columns are generated.

\paragraph{\texttt{LALTwoDColumn()}:} This routine first locates the
centreline of the column, and uses \verb@params->getRange()@ to see
how much of this centreline is taken up by the requested parameter
region, restricted by any clipping area specified in \verb@*column@.
If any region of this centreline remains uncovered,
\verb@LALTwoDColumn()@ places a tile at the bottom of the uncovered
portion, and stacks more tiles upward until it reaches the top of the
uncovered line.  The sizes and orientations of each tile are
determined from calls to \verb@params->getMetric@.

While it is doing this, \verb@LALTwoDColumn()@ keeps track of the
bottom and top corners of the newly-covered region.  If it finds that
the top corners are not increasing monotonically, this is usually an
indication that the metric is changing too rapidly, or that the tiles
are getting too thin and tilted.  Often this can be corrected by using
narrower (and taller) tiles, so \verb@LALTwoDColumn()@ reports this as
a ``column too wide'' result: it sets \verb@column->tooWide@, frees
everythin attached to \verb@**tail@ and reduced \verb@params->nOut@
accordingly, then returns.  This is also done if
\verb@LALTwoDColumn()@ ever determines that the maximum width of a
mismatch ellipse is less than \verb@params->widthMaxFac@ times the
current column width.

Having successfully stacked up the centreline of the current column,
\verb@LALTwoDColumn()@ then checks to see whether corners of the
parameter region extend above or below the top and bottom of the
newly-tiled region on either side of the centreline, by looking at the
values in \verb@column->leftRange@ and \verb@column->rightRange@.  If
a corner remains uncovered, \verb@LALTwoDColumn()@ calls itself
recursively on a half-width column on the appropriate side, setting
the clipping area of the subroutine call to exclude the region already
covered.  In principle it could call itself up to four times (once for
each column), and each recursive call could in turn call itself
recursively in order to cover a particularly steep or complicated
boundary.  However, in most cases at most two additional tiles need to
be placed (one on a top corner, one on a bottom corner).  If you're
concerned about a runaway process, set \verb@params->nIn@ to stop
generation after a given number of tiles.  If a recursive call reports
the column is too wide, this information is passed up the calling
chain.

Once the centreline and any corners have been successfully covered,
\verb@LALTwoDColumn()@ updates \verb@*tail@ to the new tail of the
list, and returns.

\paragraph{\texttt{LALTwoDNodeCopy()}:} This routine works by a simple
recursive algorithm.  First, a copy of the node is allocated and the
output handle is pointed to it.  Next, all non-pointer fields are
copied over.  Then, if \verb@new->subMesh@ exists,
\verb@LALTwoDNodeCopy()@ navigates its way along the list, calling
itself recursively on each node, and attaching copies of each node to
a corresponding list in \verb@(*new)->subMesh@.  If any errors are
detected, \verb@*new@ is destroyed via \verb@LALDestroyTwoDMesh()@,
restoring it to \verb@NULL@.

\paragraph{Computing tile sizes:} Given a positive-definite
2-dimensional metric $\mathsf{g}_{ab}$, the elliptical contour
corresponding to a maximum mismatch $m_\mathrm{thresh}$ is given by
the set of points offset from the centre point by amounts $(\Delta
x,\Delta y)$, where:
$$
m_\mathrm{thresh} = g_{xx}(\Delta x)^2 + g_{yy}(\Delta y)^2
	+ 2g_{xy}\Delta x\Delta y \; .
$$
Thus for a tile of some half-width $\Delta x$, the heights of the two
right-hand corners of the tile relative to its centre are:
$$
\Delta y = \frac{-g_{xy}\Delta x \pm\sqrt{g_{yy}m_\mathrm{thresh} -
	( g_{xx}g_{yy} - g_{xy}^2 )(\Delta x)^2}}{g_{yy}} \; ,
$$
and the maximum half-width of a tile is:
$$
\Delta x_\mathrm{max} = \sqrt{\frac{g_{yy}m_\mathrm{thresh}}
	{g_{xx}g_{yy} - g_{xy}^2}} \; .
$$
The positive-definiteness of the metric ensures that $g_{xx}>0$,
$g_{yy}>0$, and $g_{xx}g_{yy}>g_{xy}^2$.  Furthermore, if one
maximizes the proper area of a tile with respect to $\Delta x$, one
finds that the \emph{optimal} tile half-width is:
$$
\Delta x_\mathrm{opt} = \frac{\Delta x_\mathrm{max}}{\sqrt{2}} \; .
$$

When estimating the width for the next column, \verb@LALTwoDMesh()@
computes $\Delta x_\mathrm{opt}$ at both the bottom and the top of the
column and uses the smaller value (it is almost always better to
underestimate $\Delta x$ than to overestimate it).  In
\verb@LALTwoDColumn()@, tile heights are computed using the column
half-width $\Delta x$ and the value of the metric components at its
particular location.

We also note that the width of a column is computed using the metric
evaluated at the edge joining it to the previous column, and the
height of a tile is computed using the metric evaluated at the edge
joining it to the previous tile.  In principle it might be more
accurate to refine our estimate of the column width or tile height by
re-evaluating the metric at their centres, but this may be a
significant excess computational burden.  Furthermore, if the metric
varies enough that the estimated width or height changes significantly
over that distance, then the quadratic approximation to the match
function is breaking down, and we shouldn't be treating the
constant-mismatch contour as an ellipse.  The routines here do not do
any sophisticated sanity-checking, though.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALInfo()                   LALPrintError()
LALMalloc()                 LALDestroyTwoDMesh()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoDMeshInternalCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>

NRCSID( TWODMESHINTERNALC, "$Id$" );

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
  REAL4 det = (metric)[0]*(metric)[1] - (metric)[2]*(metric)[2];     \
  if ( ( (metric)[0] <= 0.0 ) || ( (metric)[1] <= 0.0 ) ||           \
       ( det <= 0.0 ) ) {                                            \
    TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*(tail))->next),    \
			     NULL ), stat );                         \
    ABORT( stat, TWODMESHH_EMETRIC, TWODMESHH_MSGEMETRIC );          \
  }                                                                  \
  (dx) = sqrt( (metric)[1]*(mismatch)/det );                         \
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
  REAL4 det = (metric)[0]*(metric)[1] - (metric)[2]*(metric)[2];     \
  REAL4 disc;                                                        \
  if ( ( metric[0] <= 0.0 ) || ( metric[1] <= 0.0 ) ||               \
       ( det <= 0.0 ) ) {                                            \
    TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*(tail))->next),    \
			     NULL ), stat );                         \
    ABORT( stat, TWODMESHH_EMETRIC, TWODMESHH_MSGEMETRIC );          \
  }                                                                  \
  if ( widthMaxFac*(dx) > sqrt( (metric)[1]*(mismatch)/det ) )       \
    TOOWIDERETURN;                                                   \
  disc = sqrt( (metric)[1]*(mismatch) - det*(dx)*(dx) );             \
  (dy)[0] = ( -metric[2]*dx - disc ) / metric[1];                    \
  (dy)[1] = ( -metric[2]*dx + disc ) / metric[1];                    \
} while (0)


/* <lalVerbatim file="TwoDMeshInternalCP"> */
void
LALTwoDMesh( LALStatus          *stat,
	     TwoDMeshNode       **tail,
	     TwoDMeshParamStruc *params )
{ /* </lalVerbatim> */
  TwoDColumnParamStruc column; /* parameters for current column */
  TwoDMeshNode *here;          /* current tail of linked list */

  /* Default parameter values: */
  REAL4 widthRetryFac = TWODMESHINTERNALC_WRETRYFAC;
  REAL4 maxColumnFac = 0.0;
  UINT4 nIn = (UINT4)( -1 );

  INITSTATUS( stat, "LALTwoDMesh", TWODMESHINTERNALC );
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
  TRY( (params->getRange)( stat->statusPtr, column.leftRange,
			   column.domain[0], params->rangeParams ),
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
    REAL4 w1, w2;      /* bottom and top widths of column */

    /* Estimate column width. */
    position[0] = column.domain[0];
    position[1] = column.leftRange[0];
    (params->getMetric)( stat->statusPtr, metric, position,
			 params->metricParams );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
			       NULL ), stat );
    ENDFAIL( stat );
    GETWIDTH( w1, metric, params->mThresh );
    position[1] = column.leftRange[1];
    (params->getMetric)( stat->statusPtr, metric, position,
			 params->metricParams );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
			       NULL ), stat );
    ENDFAIL( stat );
    GETWIDTH( w2, metric, params->mThresh );
    if ( w2 < w1 )
      w1 = w2;
    w1 *= LAL_SQRT2;

    /* Loop to try successively smaller column widths. */
    do {
      /* Make sure width is not too small or too big. */
      if ( maxColumnFac*( params->domain[1] - params->domain[0] )
	   > w1 ) {
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
	ABORT( stat, TWODMESHH_EWIDTH, TWODMESHH_MSGEWIDTH );
      }
      column.domain[1] = column.domain[0] + w1;
      if ( column.domain[1] > params->domain[1] ) {
	column.domain[1] = params->domain[1];
	w1 = column.domain[1] - column.domain[0];
      }

      /* Set remaining column parameters. */
      (params->getRange)( stat->statusPtr, column.rightRange,
			  column.domain[1], params->rangeParams );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
      ENDFAIL( stat );
      column.tooWide = 0;

      /* Call LALTwoDColumn() to place the column. */
      LALTwoDColumn( stat->statusPtr, &here, &column, params );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
      ENDFAIL( stat );

      /* See if we've reached the maximum number of mesh points. */
      if ( params->nOut >= nIn ) {
	*tail = here;
	DETATCHSTATUSPTR( stat );
	RETURN( stat );
      }

      /* If necessary, repeat with a narrower column. */
      w1 /= widthRetryFac;
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
LALTwoDColumn( LALStatus            *stat,
	       TwoDMeshNode         **tail,
	       TwoDColumnParamStruc *column,
	       TwoDMeshParamStruc   *params )
{ /* </lalVerbatim> */
  BOOLEAN tiled = 0;    /* whether tiles were placed on the centreline */
  REAL4 position[2];    /* current top of column */
  REAL4 dx;             /* half-width of column */
  REAL4 myy0, myy1;         /* temporary variables storing y-coordinates */
  REAL4 centreRange[2]; /* centreline of column parameter space */
  REAL4 centreClip[2];  /* centre of clip boundary */
  REAL4 leftTiled[2];   /* left side of region tiled */
  REAL4 rightTiled[2];  /* right side of region tiled */
  REAL4 metric[3];      /* current metric components */
  TwoDMeshNode *here;   /* current node in list */

  /* Default parameter values: */
  REAL4 widthMaxFac = TWODMESHINTERNALC_WMAXFAC;
  UINT4 nIn = (UINT4)( -1 );

  INITSTATUS( stat, "LALTwoDColumn", TWODMESHINTERNALC );
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
  centreClip[0] = 0.5*column->leftClip[0] + 0.5*column->rightClip[0];
  centreClip[1] = 0.5*column->leftClip[1] + 0.5*column->rightClip[1];
  leftTiled[0] = column->leftClip[1];
  leftTiled[1] = column->leftClip[0];
  rightTiled[0] = column->rightClip[1];
  rightTiled[1] = column->rightClip[0];

  /* Get the width and heights of this column. */
  position[0] = 0.5*( column->domain[1] + column->domain[0] );
  dx = 0.5*( column->domain[1] - column->domain[0] );
  TRY( (params->getRange)( stat->statusPtr, centreRange, position[0],
			   params->rangeParams ), stat );

  /* Add the column of tiles along the centreline, if the parameter
     space intersects the clipping area along the centreline. */
  position[1] = centreClip[0];
  if ( position[1] < centreRange[0] )
    position[1] = centreRange[0];
  if ( position[1] < centreRange[1] ) {

    /* Add base tile of column. */
    tiled = 1;
    TRY( (params->getMetric)( stat->statusPtr, metric, position,
			      params->metricParams ), stat );
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
    myy0 = here->y + here->dy[0];
    myy1 = here->y - here->dy[1];
    if ( leftTiled[0] > myy1 )
      leftTiled[0] = myy1;
    if ( rightTiled[0] > myy0 )
      rightTiled[0] = myy0;
    leftTiled[1] = here->y - here->dy[0];
    rightTiled[1] = here->y + here->dy[1];
    position[1] = 0.5*leftTiled[1] + 0.5*rightTiled[1];

    /* Continue stacking tiles until we reach the top. */
    while ( ( position[1] < centreRange[1] ) &&
	    ( position[1] < centreClip[1] ) ) {
      (params->getMetric)( stat->statusPtr, metric, position,
			   params->metricParams );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
      ENDFAIL( stat );
      here->next = (TwoDMeshNode *)LALMalloc( sizeof(TwoDMeshNode) );
      if ( here == NULL ) {
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
	ABORT( stat, TWODMESHH_EMEM, TWODMESHH_MSGEMEM );
      }
      memset( here->next, 0, sizeof(TwoDMeshNode) );
      params->nOut++;
      if (lalDebugLevel >= 3) 
	{
	  LALPrintError( "\r%16u", params->nOut );
	}
      GETSIZE( here->next->dy, dx, metric, params->mThresh );
      myy0 = here->dy[1] - here->next->dy[0];
      myy1 = here->next->dy[1] - here->dy[0];
      if ( myy0 > myy1 )
	myy0 = myy1;
      if ( myy0 <= 0.0 )
	TOOWIDERETURN;
      here->next->y = here->y + myy0;
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
      leftTiled[1] = here->y - here->dy[0];
      rightTiled[1] = here->y + here->dy[1];
      position[1] = 0.5*leftTiled[1] + 0.5*rightTiled[1];
    }
  }

  /* Centreline stacking is complete.  Now check for exposed corners
     of the parameter space, and call LALTwoDColumn() recursively. */

  /* Check bottom corners. */
  myy0 = 0.5*leftTiled[0] + 0.5*rightTiled[0];

  /* Bottom-left: */
  if ( ( ( column->leftClip[0] < leftTiled[0] ) ||
	 ( centreClip[0] < myy0 ) ) &&
       ( column->leftRange[0] < leftTiled[0] ) &&
       ( ( column->leftRange[1] > column->leftClip[0] ) ||
	 ( centreRange[1] > centreClip[0] ) ) ) {
    TwoDColumnParamStruc column2;
    column2.domain[0] = column->domain[0];
    column2.domain[1] = position[0];
    memcpy( column2.leftRange, column->leftRange, 2*sizeof(REAL4) );
    memcpy( column2.leftClip, column->leftClip, 2*sizeof(REAL4) );
    memcpy( column2.rightRange, centreRange, 2*sizeof(REAL4) );
    memcpy( column2.rightClip, centreClip, 2*sizeof(REAL4) );
    if ( ( leftTiled[0] < column2.leftClip[1] ) &&
	 ( myy0 < column2.rightClip[1] ) ) {
      column2.leftClip[1] = leftTiled[0];
      column2.rightClip[1] = myy0;
    }
    LALTwoDColumn( stat->statusPtr, &here, &column2, params );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
			       NULL ), stat );
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
  if ( ( ( column->rightClip[0] < rightTiled[0] ) ||
	 ( centreClip[0] < myy0 ) ) &&
       ( column->rightRange[0] < rightTiled[0] ) &&
       ( ( column->rightRange[1] > column->rightClip[0] ) ||
	 ( centreRange[1] > centreClip[0] ) ) ) {
    TwoDColumnParamStruc column2;
    column2.domain[1] = column->domain[1];
    column2.domain[0] = position[0];
    memcpy( column2.rightRange, column->rightRange, 2*sizeof(REAL4) );
    memcpy( column2.rightClip, column->rightClip, 2*sizeof(REAL4) );
    memcpy( column2.leftRange, centreRange, 2*sizeof(REAL4) );
    memcpy( column2.leftClip, centreClip, 2*sizeof(REAL4) );
    if ( ( rightTiled[0] < column2.rightClip[1] ) &&
	 ( myy0 < column2.leftClip[1] ) ) {
      column2.rightClip[1] = rightTiled[0];
      column2.leftClip[1] = myy0;
    }
    LALTwoDColumn( stat->statusPtr, &here, &column2, params );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
			       NULL ), stat );
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
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
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
      if ( ( rightTiled[1] > column2.rightClip[0] ) &&
	   ( myy0 > column2.leftClip[0] ) ) {
	column2.rightClip[0] = rightTiled[1];
	column2.leftClip[0] = myy0;
      }
      LALTwoDColumn( stat->statusPtr, &here, &column2, params );
      BEGINFAIL( stat )
	TRY( LALDestroyTwoDMesh( stat->statusPtr, &((*tail)->next),
				 NULL ), stat );
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


/* <lalVerbatim file="TwoDMeshInternalCP"> */
void
LALTwoDNodeCopy( LALStatus    *stat,
		 TwoDMeshNode **new,
		 TwoDMeshNode *old )
{ /* </lalVerbatim> */
  TwoDMeshNode *tail;      /* current tail of old->subMesh */
  TwoDMeshNode **tailCopy; /* pointer to copy of *tail */

  INITSTATUS( stat, "LALTwoDNodeCopy", TWODMESHINTERNALC );
  ATTATCHSTATUSPTR( stat );

  /* Check parameters. */
  ASSERT( old, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( new, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( !(*new), stat, TWODMESHH_EOUT, TWODMESHH_MSGEOUT );

  /* Copy top-level node. */
  *new = (TwoDMeshNode *)LALMalloc( sizeof(TwoDMeshNode) );
  if ( *new == NULL ) {
    ABORT( stat, TWODMESHH_EMEM, TWODMESHH_MSGEMEM );
  }
  **new = *old;
  (*new)->next = (*new)->subMesh = NULL;

  /* Recursively copy submesh. */
  tail = old->next;
  tailCopy = &((*new)->next);
  while ( tail != NULL ) {
    LALTwoDNodeCopy( stat->statusPtr, tailCopy, tail );
    BEGINFAIL( stat )
      TRY( LALDestroyTwoDMesh( stat->statusPtr, new, NULL ), stat );
    ENDFAIL( stat );
    tail = tail->next;
    tailCopy = &((*tailCopy)->next);
  }

  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
