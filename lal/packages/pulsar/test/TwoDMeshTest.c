/********************************* <lalVerbatim file="TwoDMeshTestCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{TwoDMeshTest.c}}
\label{ss:TwoDMeshTest.c}

Creates a 2-dimensional template mesh for linearly-changing mismatch
ellipses.

\subsubsection*{Usage}
\begin{verbatim}
TwoDMeshTest [-o outfile] [-p psfile flags] [-d debug] [-b dx1 dy1 dx2 dy2 ]
             [-e a b c] [-x dadx dbdx dcdx] [-y dady dbdy dcdy]
\end{verbatim}

\subsubsection*{Description}

This test program creates a template mesh for a parameter space with a
constant mismatch metric.  The following option flags are accepted:
\begin{itemize}
\item[\texttt{-o}] Writes the output mesh list to the file
\verb@outfile@.  If absent, no output is written.
\item[\texttt{-p}] Plots the output mesh in a PostScript file
\verb@psfile@, using plot flags \verb@flags@ (see below).  If absent,
no plot is made.
\item[\texttt{-d}] Sets the debug level to \verb@debug@.  If
absent, a debug level of zero is used.
\item[\texttt{-b}] Sets the parameter space boundary to be a
parallelogram defined by the vectors (\verb@dx1@,\verb@dy1@) and
(\verb@dx2@,\verb@dy2@) from the origin.  If absent, the region is
taken to be a unit square.
\item[\texttt{-e}] Sets the parameters of the mismatch ellipse at the
origin: its principal axis lengths are \verb@a@ and \verb@b@ units,
and the angle from the $x$-axis to the first principal axis is
\verb@c@ radians.  If absent, the values \verb@a@=0.1, \verb@b@=0.05,
and \verb@c@=1 are assumed.
\item[\texttt{-x}] Sets the rates of change in the $x$-direction of
\verb@a@, \verb@b@, and \verb@c@ (above) to \verb@dadx@, \verb@dbdx@,
and \verb@dcdx@, respectively.  If absent, the rates are taken to be
zero.
\item[\texttt{-y}] Sets the rates of change in the $y$-direction of
\verb@a@, \verb@b@, and \verb@c@ (above) to \verb@dady@, \verb@dbdy@,
and \verb@dcdy@, respectively.  If absent, the rates are taken to be
zero.
\end{itemize}

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
#define TWODMESHTESTC_ENORM   0
#define TWODMESHTESTC_ESUB    1
#define TWODMESHTESTC_EARG    2
#define TWODMESHTESTC_EBAD    3
#define TWODMESHTESTC_EMEM    4
#define TWODMESHTESTC_EFILE   5
#define TWODMESHTESTC_EMETRIC 6

#define TWODMESHTESTC_MSGENORM   "Normal exit"
#define TWODMESHTESTC_MSGESUB    "Subroutine failed"
#define TWODMESHTESTC_MSGEARG    "Error parsing arguments"
#define TWODMESHTESTC_MSGEBAD    "Bad argument value"
#define TWODMESHTESTC_MSGEMEM    "Memory allocation error"
#define TWODMESHTESTC_MSGEFILE   "Could not open file"
#define TWODMESHTESTC_MSGEMETRIC "Ellipse axis length is zero or negative within the specified region"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

The test program reads the input arguments and creates a parameter
structure \verb@*params@ to be passed to \verb@LALCreateTwoDMesh()@.
In particular, it computes the domain of the parameter space, and
defines functions and parameter lists to compute the range in $y$ at
any $x$, and the metric at any point $(x,y)$.  If PostScript output is
requested, it is generated using \verb@LALPlotTwoDMesh()@, using the
value of the command-line number \verb@flags@ to set the plotting
parameters.  Each of these functions is discussed below.

\paragraph{Parameter ranges:} The parameter space is defined to be a
parallelogram with one corner on the origin, and two sides defined by
vectors $(x_1,y_1)$ and $(x_2,y_2)$.  Without loss of generality we
assume that $x_1<x_2$.  The functions defining the boundaries are
denoted $y_{a,b}(x)$, and we make no assumption about their signs or
relative order.  The algorithm used then depends on the signs of $x_1$
and $x_2$.

\medskip\noindent
If $x_1=x_2=0$, then the parameter space is singular, and no mesh need
be generated.

\medskip\noindent
If $x_1=0$ and $x_2\neq0$, then the domain is $[0,x_2]$, and the
boundary functions are:
\begin{eqnarray}
y_a(x) & = & y_2x/x_2 \nonumber\\
y_b(x) & = & y_1 + y_2x/x_2 \nonumber
\end{eqnarray}

\noindent
If $x_2=0$ and $x_1\neq0$, then the domain is $[x_1,0]$, and the above
equations for $y_{a,b}(x)$ simply have 1 and 2 reversed.

\medskip\noindent
If $x_1$ and $x_2$ have the same sign, then the domain is
$[0,x_1+x_2]$ if $x_1$ and $x_2$ are positive, and $[x_1+x_2,0]$
otherwise.  The boundary functions are:
\begin{eqnarray}
y_a(x) & = & \left\{\begin{array}{c@{\qquad}c}
	y_1x/x_1             & x\mathrm{~between~}0\mathrm{~and~}x_1 \\
	y_1 + y_2(x-x_1)/x_2 & x\mathrm{~between~}x_1\mathrm{~and~}x_1+x_2
	\end{array}\right.\nonumber\\
y_b(x) & = & \left\{\begin{array}{c@{\qquad}c}
	y_2x/x_2             & x\mathrm{~between~}0\mathrm{~and~}x_2 \\
	y_2 + y_1(x-x_2)/x_1 & x\mathrm{~between~}x_2\mathrm{~and~}x_1+x_2
	\end{array}\right.\nonumber
\end{eqnarray}

\noindent
If $x_1$ and $x_2$ have opposite sign, the domain is $[x_1,x_2]$ if
$x_1<0$, and $[x_2,x_1]$ otherwise.  The boundary functions are:
\begin{eqnarray}
y_a(x) & = & \left\{\begin{array}{c@{\qquad}c}
	y_1x/x_1 & x\mathrm{~between~}0\mathrm{~and~}x_1 \\
	y_2x/x_2 & x\mathrm{~between~}0\mathrm{~and~}x_2
	\end{array}\right.\nonumber\\
y_b(x) & = & \left\{\begin{array}{c@{\qquad}c}
	y_1 + y_2(x-x_1)/x_2 & x\mathrm{~between~}x_1\mathrm{~and~}x_1+x_2 \\
	y_2 + y1(x-x_2)/x_1  & x\mathrm{~between~}x_2\mathrm{~and~}x_1+x_2
	\end{array}\right.\nonumber
\end{eqnarray}

The main program sorts the input parameters so that $x_1\leq x_2$,
stores them in a 4-dimensional array, and assigns a \verb@void@
pointer to that array.  It also computes the domain.  The routine
\verb@LALTwoDRangeTest()@ takes a value of $x$ and the \verb@void@
pointer, computes the values of $y_a(x)$ and $y_b(x)$ according to the
algorithm above, sorts them, and returns them ordered from lower to
higher.

\paragraph{Metric values:} The main program takes the input parameters
\verb@a@, \verb@b@, \verb@c@, \verb@dadx@, \verb@dbdx@, and
\verb@dcdx@, stores them in a 9-dimensional array, and assigns a
\verb@void@ pointer to it.  The routine \verb@LALTwoDMetricTest()@
takes a position $(x,y)$ and the \verb@void@ pointer, and computes the
``local'' value of the principal axis
$a=$\verb@a@$+x\times$\verb@dadx@$+y\times$\verb@dady@, and similarly
for $b$ and $c$.  If that ellipse corresponds to the
$m_\mathrm{thresh}$ mismatch level contour, then the eigenvalues of
the corresponding metric are $\lambda_1=m_\mathrm{thresh}/a^2$ and
$\lambda_2=m_\mathrm{thresh}/b^2$.  The metric components are thus:
\begin{eqnarray}
g_{xx} & = & \lambda_1\cos^2(c) + \lambda_2\sin^2(c) \;,\nonumber\\
g_{yy} & = & \lambda_1\sin^2(c) + \lambda_2\cos^2(c) \;,\nonumber\\
g_{xy} \quad = \quad g_{yx} & = & (\lambda_1-\lambda_2)\cos(c)\sin(c)
	\;.\nonumber
\end{eqnarray}
The routine assumes that the values of $a$, $b$, and $c$ refer to an
$m_\mathrm{thresh}=1$ mismatch ellipse.  It computes and returns
$g_{xx}$, $g_{yy}$, and $g_{xy}$ in a 3-dimensional array.

\paragraph{PostScript flags:} The parameter \verb@flags@ is an
unsigned integer whose lowest-order bits contain parameters to be
passed to \verb@LALPlotTwoDMesh()@.  The bits and their meanings are:
\begin{description}
\item[bit 0:] 1 if mesh points will be plotted, 0 otherwise.
\item[bit 1:] 1 if mesh tiles will be plotted, 0 otherwise.
\item[bit 2:] 1 if mismatch ellipses will be plotted, 0 otherwise.
\item[bit 3:] 1 if the boundary will be plotted, 0 otherwise.
\end{description}
Thus a value of 15 will plot everything, while a value of 9 will just
plot the mesh points and the boundary.  A value of zero suppresses the
plot.

If mesh points are to be plotted, they will be filled circles $1/72''$
(1~point) in diameter.  The parameter space will be rotated so that
the longer of the diagonals of the parallelogram will be vertical, and
scaled to fit on one $8.5''\times11''$ page.  That is, if
$||(x_1+x_2,y_1+y_2)||\geq||(x_1-x_2,y_1-y_2)||$, the rotation angle
of the coordinate axes will be
$\theta=\pi/2-\arctan\!2(y_1+y_2,x_1+x_2)$, or
$\theta=\pi/2-\arctan\!2(y_1-y_2,x_1-x_2)$ otherwise.  We note that
the function $\arctan\!2(y,x)$ returns the argument of the complex
number $x+iy$ in the range $[-\pi,\pi]$.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALCreateTwoDMesh()             LALDestroyTwoDMesh()
LALPlotTwoDMesh()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoDMeshTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TwoDMesh.h>
#include "TwoDMeshPlot.h"

NRCSID( TWODMESHTESTC, "$Id$" );

/* Default parameter settings. */
int lalDebugLevel = 0;
#define DX1 (1.0)
#define DY1 (0.0)
#define DX2 (0.0)
#define DY2 (1.0)
#define A_DEFAULT (0.1)
#define B_DEFAULT (0.05)
#define C_DEFAULT (1.0)
#define DADX (0.0)
#define DBDX (0.0)
#define DCDX (0.0)
#define DADY (0.0)
#define DBDY (0.0)
#define DCDY (0.0)

/* Other numerical constants. */
#define MISMATCH (1.0)  /* arbitrary mismatch threshold used */

/* Usage format string. */
#define USAGE "Usage: %s TwoDMeshTest [-o outfile] [-p psfile flags] [-d debug]\n\t[-b dx1 dy1 dx2 dy2 ] [-e a b c]\n\t[-x dadx dbdx dcdx] [-y dady dbdy dcdy]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, TWODMESHTESTC, statement ? statement :    \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 TWODMESHTESTC, (statement) );                       \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( TWODMESHTESTC_ESUB, TWODMESHTESTC_MSGESUB,                  \
         "Function call \"" #func "\" failed:" );                    \
  return TWODMESHTESTC_ESUB;                                         \
}                                                                    \
else (void)(0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

/* Local prototypes. */
void
LALRangeTest( LALStatus *stat, REAL4 range[2], REAL4 x, void *params );

void
LALMetricTest( LALStatus *stat,
	       REAL4 metric[3],
	       REAL4 position[2],
	       void *params );

int
main(int argc, char **argv)
{
  INT4 arg;                  /* argument counter */
  static LALStatus stat;     /* top-level status structure */
  CHAR *outfile = NULL;      /* name of output file */
  CHAR *psfile = NULL;       /* name of PostScript output file */
  UINT2 flags;               /* plotting flags */
  REAL4 rangeParams[4];      /* LALRangeTest() parameters */
  REAL4 metricParams[9];     /* LALMetricTest() parameters */
  TwoDMeshParamStruc params; /* LALCreateTwoDMesh() parameters */
  TwoDMeshNode *mesh = NULL; /* head of mesh list */
  FILE *fp;                  /* output file pointer */

  /* Boundary parameters: */
  REAL4 dx1 = DX1, dy1 = DY1, dx2 = DX2, dy2 = DY2;
  /* Ellipse parameters: */
  REAL4 a = A_DEFAULT, b = B_DEFAULT, c = C_DEFAULT;
  REAL4 dadx = DADX, dbdx = DBDX, dcdx = DCDX;
  REAL4 dady = DADY, dbdy = DBDY, dcdy = DCDY;


  /* DEBUG: Do nothing at present. */
  return 0;


  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse PostScript file option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	psfile = argv[arg++];
	flags = atoi( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	lalDebugLevel = atoi( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse boundary parameters option. */
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 4 ) {
	arg++;
	dx1 = atof( argv[arg++] );
	dy1 = atof( argv[arg++] );
	dx2 = atof( argv[arg++] );
	dy2 = atof( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse ellipse parameters option. */
    else if ( !strcmp( argv[arg], "-e" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	a = atof( argv[arg++] );
	b = atof( argv[arg++] );
	c = atof( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse ellipse variation in x option. */
    else if ( !strcmp( argv[arg], "-x" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	dadx = atof( argv[arg++] );
	dbdx = atof( argv[arg++] );
	dcdx = atof( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse ellipse variation in y option. */
    else if ( !strcmp( argv[arg], "-y" ) ) {
      if ( argc > arg + 3 ) {
	arg++;
	dady = atof( argv[arg++] );
	dbdy = atof( argv[arg++] );
	dcdy = atof( argv[arg++] );
      }else{
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return TWODMESHTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /* Set up range function and parameters. */
  if ( ( dx1 == 0.0 ) && ( dx2 == 0.0 ) ) {
    ERROR( TWODMESHTESTC_EBAD, TWODMESHTESTC_MSGEBAD,
	   "dx1 = dx2 = 0:" );
    return TWODMESHTESTC_EBAD;
  }
  if ( dx1 < dx2 ) {
    rangeParams[0] = dx1;
    rangeParams[1] = dy1;
    rangeParams[2] = dx2;
    rangeParams[3] = dy2;
  } else {
    rangeParams[0] = dx2;
    rangeParams[1] = dy2;
    rangeParams[2] = dx1;
    rangeParams[3] = dy1;
  }
  params.getRange = LALRangeTest;
  params.rangeParams = (void *)( rangeParams );
  if ( dx1*dx2 < 0.0 ) {
    if ( dx1 < dx2 ) {
      params.domain[0] = dx1;
      params.domain[1] = dx2;
    } else {
      params.domain[0] = dx2;
      params.domain[1] = dx1;
    }
  } else {
    if ( dx1 < 0.0 ) {
      params.domain[0] = dx1 + dx2;
      params.domain[1] = 0.0;
    } else {
      params.domain[0] = 0.0;
      params.domain[1] = dx1 + dx2;
    }
  }

  /* Set up metric function and parameters. */
  if ( ( a <= 0.0 ) ||
       ( a + dadx*dx1 + dady*dy1 <= 0.0 ) ||
       ( a + dadx*dx2 + dady*dy2 <= 0.0 ) ||
       ( a + dadx*( dx1 + dx2 ) + dady*( dy1 + dy2 ) <= 0.0 ) ) {
    ERROR( TWODMESHTESTC_EMETRIC, TWODMESHTESTC_MSGEMETRIC,
	   "axis a:" );
    return TWODMESHTESTC_EBAD;
  }
  if ( ( b <= 0.0 ) ||
       ( b + dbdx*dx1 + dbdy*dy1 <= 0.0 ) ||
       ( b + dbdx*dx2 + dbdy*dy2 <= 0.0 ) ||
       ( b + dbdx*( dx1 + dx2 ) + dbdy*( dy1 + dy2 ) <= 0.0 ) ) {
    ERROR( TWODMESHTESTC_EMETRIC, TWODMESHTESTC_MSGEMETRIC,
	   "axis b:" );
    return TWODMESHTESTC_EBAD;
  }
  metricParams[0] = a;
  metricParams[1] = b;
  metricParams[2] = c;
  metricParams[3] = dadx;
  metricParams[4] = dbdx;
  metricParams[5] = dcdx;
  metricParams[6] = dady;
  metricParams[7] = dbdy;
  metricParams[8] = dcdy;
  params.getMetric = LALMetricTest;
  params.metricParams = (void *)( metricParams );

  /* Set up remaining mesh creation parameters. */
  params.mThresh = MISMATCH;
  params.widthMaxFac = 0.0;
  params.widthRetryFac = 0.0;
  params.maxColumns = 0;
  params.nIn = 0;

  /* Create mesh. */
  SUB( LALCreateTwoDMesh( &stat, &mesh, &params ), &stat );

  /* Print mesh list to a file, if requested. */
  if ( outfile ) {
    TwoDMeshNode *here; /* current node in mesh list */

    if ( !( fp = fopen( outfile, "w" ) ) ) {
      ERROR( TWODMESHTESTC_EFILE, "- " TWODMESHTESTC_MSGEFILE,
	     outfile );
      return TWODMESHTESTC_EFILE;
    }
    for ( here = mesh; here; here = here->next )
      fprintf( fp, "%f %f %f %f %f\n", here->x, here->y, here->dx,
	       here->dy[0], here->dy[1] );
    fclose( fp );
  }

  /* Make a PostScript plot of the mesh, if requested. */
  if ( psfile && flags ) {
    REAL4 xSum = dx1 + dx2, xDiff = dx1 - dx2;
    REAL4 ySum = dy1 + dy2, yDiff = dy1 - dy2;
    INT2 plotPoints = flags & 1;
    BOOLEAN plotTiles = flags & 2;
    BOOLEAN plotEllipses = flags & 4;
    TwoDMeshPlotStruc plotParams;

    if ( !( fp = fopen( psfile, "w" ) ) ) {
      ERROR( TWODMESHTESTC_EFILE, "- " TWODMESHTESTC_MSGEFILE,
	     psfile );
      return TWODMESHTESTC_EFILE;
    }
    if ( xSum*xSum + ySum*ySum > xDiff*xDiff + yDiff*yDiff )
      plotParams.theta = LAL_PI_2 - atan2( ySum, xSum );
    else
      plotParams.theta = LAL_PI_2 - atan2( yDiff, xDiff );
    plotParams.xScale = plotParams.yScale = 10.0;
    plotParams.bBox[0] = 36.0;
    plotParams.bBox[1] = 36.0;
    plotParams.bBox[2] = 576.0;
    plotParams.bBox[3] = 756.0;
    plotParams.autoscale = 1;
    memset( plotParams.clipBox, 0, 4*sizeof(REAL4) );
    plotParams.nLevels = 1;
    if ( flags & 8 )
      plotParams.nBoundary = (UINT4)( fabs( 0.1*( dx1 + dx2 )/a ) );
    else
      plotParams.nBoundary = 0;
    plotParams.plotPoints = &plotPoints;
    plotParams.plotTiles = &plotTiles;
    plotParams.plotEllipses = &plotEllipses;
    plotParams.params = &params;
    SUB( LALPlotTwoDMesh( &stat, fp, mesh, &plotParams ), &stat );
    fclose( fp );
  }

  /* Free the mesh, and exit. */
  SUB( LALDestroyTwoDMesh( &stat, &mesh, NULL ), &stat );
  LALCheckMemoryLeaks();
  INFO( TWODMESHTESTC_MSGENORM );
  return TWODMESHTESTC_ENORM;
}


void
LALRangeTest( LALStatus *stat, REAL4 range[2], REAL4 x, void *params )
{
  REAL4 *xy = (REAL4 *)( params ); /* params recast to its proper type */
  REAL4 ya, yb;                    /* unsorted range values */

  INITSTATUS( stat, "LALRangeTest", TWODMESHTESTC );
  ASSERT( range, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Case 1: one of the side vectors is vertical. */
  if ( xy[0] == 0.0 ) {
    ya = xy[3]*( x/xy[2] );
    yb = ya + xy[1];
  } else if ( xy[2] == 0.0 ) {
    ya = xy[1]*( x/xy[0] );
    yb = ya + xy[3];
  }

  /* Case 2: Both side vectors point in the same direction (either
     left or right). */
  else if ( xy[0]*xy[2] > 0.0 ) {
    REAL4 dx1 = x - xy[0];
    REAL4 dx2 = x - xy[2];
    if ( x*dx1 < 0.0 )
      ya = xy[1]*( x/xy[0] );
    else
      ya = xy[1] + xy[3]*( dx1/xy[2] );
    if ( x*dx2 < 0.0 )
      yb = xy[3]*( x/xy[2] );
    else
      yb = xy[3] + xy[1]*( dx2/xy[0] );
  }

  /* Case 3: One side vector points left and the other points
     right. */
  else {
    REAL4 dx1 = x - xy[0];
    REAL4 dx2 = x - xy[2];
    REAL4 dx12 = x - xy[0] - xy[2];
    if ( x*dx1 < 0.0 )
      ya = xy[1]*( x/xy[0] );
    else
      ya = xy[3]*( x/xy[2] );
    if ( dx1*dx12 < 0.0 )
      yb = xy[1] + xy[3]*( dx1/xy[2] );
    else
      yb = xy[3] + xy[1]*( dx2/xy[0] );
  }

  /* Sort and return the range values. */
  if ( ya < yb ) {
    range[0] = ya;
    range[1] = yb;
  } else {
    range[0] = yb;
    range[1] = ya;
  }
  RETURN( stat );
}


void
LALMetricTest( LALStatus *stat,
	       REAL4 metric[3],
	       REAL4 position[2],
	       void *params )
{
  REAL4 *abc = (REAL4 *)( params ); /* params recast to proper type */
  REAL4 a, b, c;                    /* axis lengths and angle */
  REAL4 lambda1, lambda2;           /* metric eigenvalues */
  REAL4 cosc, sinc;                 /* cosine and sine of c */

  INITSTATUS( stat, "LALMetricTest", TWODMESHTESTC );
  ASSERT( metric, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( position, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Compute axis lengths and angle at current position. */
  a = abc[0] + position[0]*abc[3] + position[1]*abc[6];
  b = abc[1] + position[0]*abc[4] + position[1]*abc[7];
  c = abc[2] + position[0]*abc[5] + position[1]*abc[8];
  if ( a*b == 0.0 ) {
    ABORT( stat, TWODMESHTESTC_EMETRIC, TWODMESHTESTC_MSGEMETRIC );
  }

  /* Compute eigenvalues and trigonometric functions. */
  lambda1 = MISMATCH/( a*a );
  lambda2 = MISMATCH/( b*b );
  cosc = cos( c );
  sinc = sin( c );

  /* Return metric components. */
  metric[0] = lambda1*cosc*cosc + lambda2*sinc*sinc;
  metric[1] = lambda1*sinc*sinc + lambda2*cosc*cosc;
  metric[2] = ( lambda1 - lambda2 )*cosc*sinc;
  RETURN( stat );
}
