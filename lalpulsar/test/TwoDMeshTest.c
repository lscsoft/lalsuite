/*
*  Copyright (C) 2007 Teviet Creighton
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
 * \ingroup TwoDMeshPlot_h
 * \brief Creates a 2-dimensional template mesh for linearly-changing mismatch ellipses.
 *
 * ### Usage ###
 *
 * \code
 * TwoDMeshTest [-o outfile] [-p psfile flags] [-d debug] [-m mismatch nmax cmax]
 * [-i metricfile rangefile] [-b x1 y1 x2 y2 ] [-e a b c]
 * [-x dadx dbdx dcdx] [-y dady dbdy dcdy]
 * \endcode
 *
 * ### Description ###
 *
 * This test program creates a template mesh for a parameter space with
 * an arbitrary mismatch metric.  The following option flags are
 * accepted:
 * <ul>
 * <li><b>-o</b> Writes the output mesh list to the file
 * \c outfile.  If absent, no output is written.</li>
 * <li><b>-p</b> Plots the output mesh in a PostScript file
 * \c psfile, using plot flags \c flags (see below).  If absent,
 * no plot is made.</li>
 * <li><b>-d</b> Sets the debug level to \c debug.  If
 * absent, a debug level of zero is used.</li>
 * <li><b>-m</b> Sets the maximum mismatch to \c mismatch,
 * maximum number of mesh points to \c nmax and the maximum estimated
 * number of columns to \c cmax.  If \c mismatch is not in the
 * range (0,1], it is taken to be 1.  If \c nmax or \c cmax is
 * non-positive, it is ignored (no maximum).  If this option is not
 * given, <b>-m 1 0 0</b> is assumed.</li>
 * <li><b>-i</b> Determines the metric and the parameter space
 * boundary from \c REAL4Grid structures stored in the files
 * \c metricfile and \c rangefile, read using the generic parser
 * LALSReadGrid().  The formats for these grids are discussed
 * below.  If present, this option \e overrides the <b>-b</b>,
 * <b>-e</b>, <b>-x</b>, and <b>-y</b> options, below.  If absent, these
 * options, or their defaults, will be used.</li>
 * <li><b>-b</b> Sets the parameter space boundary to be a
 * parallelogram defined by the vectors (\c x1,\c y1) and
 * (\c x2,\c y2) from the origin.  If absent, the region is
 * taken to be a unit square.</li>
 * <li><b>-e</b> Sets the parameters of the mismatch ellipse at the
 * origin: its principal axis lengths are \c a and \c b units,
 * and the angle from the \f$x\f$-axis to the first principal axis is
 * \c c radians.  If absent, <b>-e 0.1 0.05 1</b> is assumed.</li>
 * <li><b>-x</b> Sets the rates of change in the \f$x\f$-direction of
 * \c a, \c b, and \c c (above) to \c dadx, \c dbdx,
 * and \c dcdx, respectively.  If absent, the rates are taken to be
 * zero.</li>
 * <li><b>-y</b> Sets the rates of change in the \f$y\f$-direction of
 * \c a, \c b, and \c c (above) to \c dady, \c dbdy,
 * and \c dcdy, respectively.  If absent, the rates are taken to be
 * zero.</li>
 * </ul>
 *
 * ### Algorithm ###
 *
 * The test program reads the input arguments and creates a parameter
 * structure <tt>*params</tt> to be passed to LALCreateTwoDMesh().
 * In particular, it computes the domain of the parameter space, and
 * defines functions and parameter lists to compute the range in \f$y\f$ at
 * any \f$x\f$, and the metric at any point \f$(x,y)\f$.  If PostScript output is
 * requested, it is generated using LALPlotTwoDMesh(), using the
 * value of the command-line number \c flags to set the plotting
 * parameters.  Each of these functions is discussed below.
 *
 * \par Metric and range grid files:
 * If the <b>-i</b> option was
 * given, the metric and parameter ranges are read from the files named
 * by the \c metricfile and \c rangefile arguments.  These two
 * files must be in a format parseable by LALSReadGrid(); see the
 * documentation of that routine for more details.
 *
 * The \c REAL4Grid extracted from \c metricfile must
 * have grid dimension 2 and data dimension 3: the grid dimensions refer
 * to the \f$(x,y)\f$ coordinates of the points where the metric is
 * evaluated, while the third dimension must have length 3, storing the
 * metric components \f$g_{xx}\f$, \f$g_{yy}\f$, and \f$g_{xy}\f$ (in that order) at
 * each point.  Within the \ref TwoDMesh_h routines, this metric grid
 * is interpolated using LALInterpolateMetricGrid().
 *
 * The \c REAL4Grid extracted from \c rangefile must have grid
 * dimension 1 and data dimension 2: the grid dimension refers to an \f$x\f$
 * coordinate, and the second dimension must have length 3, storing the
 * lower and upper boundaries \f$y_1(x)\f$ and \f$y_2(x)\f$ at each sampled value
 * of \f$x\f$.  Within the \ref TwoDMesh_h routines, this range grid is
 * interpolated using LALInterpolateRangeGrid().
 *
 * If the <b>-i</b> option is \e not given, then the parameter
 * boundary and metric are determined by internal routines, with default
 * settings that can be overridden using command-line options <b>-p</b>,
 * <b>-e</b>, <b>-x</b>, and <b>-y</b>.
 *
 * \par Parameter ranges:
 * The parameter space boundary can be
 * specified by input parameters \c x1\f$=x_1\f$, \c x2\f$=x_2\f$,
 * \c y1\f$=y_1\f$, and \c y2\f$=y_2\f$.  The parameter space is then
 * defined to be a parallelogram with one corner on the origin, and two
 * sides defined by vectors \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$.  Without loss of
 * generality we assume that \f$x_1<x_2\f$.  The functions defining the
 * boundaries are denoted \f$y_{a,b}(x)\f$, and we make no assumption about
 * their signs or relative order.  The algorithm used then depends on the
 * signs of \f$x_1\f$ and \f$x_2\f$.
 *
 * If \f$x_1=x_2=0\f$, then the parameter space is singular, and no mesh need
 * be generated.
 *
 * If \f$x_1=0\f$ and \f$x_2\neq0\f$, then the domain is \f$[0,x_2]\f$, and the
 * boundary functions are:
 * \f{eqnarray}{
 * y_a(x) & = & y_2x/x_2 \nonumber\\
 * y_b(x) & = & y_1 + y_2x/x_2 \nonumber
 * \f}
 *
 * If \f$x_2=0\f$ and \f$x_1\neq0\f$, then the domain is \f$[x_1,0]\f$, and the above
 * equations for \f$y_{a,b}(x)\f$ simply have 1 and 2 reversed.
 *
 * If \f$x_1\f$ and \f$x_2\f$ have the same sign, then the domain is
 * \f$[0,x_1+x_2]\f$ if \f$x_1\f$ and \f$x_2\f$ are positive, and \f$[x_1+x_2,0]\f$
 * otherwise.  The boundary functions are:
 * \f{eqnarray}{
 * y_a(x) & = & \left\{\begin{array}{c@{\qquad}c}
 * y_1x/x_1             & x\mathrm{~between~}0\mathrm{~and~}x_1 \\
 * y_1 + y_2(x-x_1)/x_2 & x\mathrm{~between~}x_1\mathrm{~and~}x_1+x_2
 * \end{array}\right.\nonumber\\
 * y_b(x) & = & \left\{\begin{array}{c@{\qquad}c}
 * y_2x/x_2             & x\mathrm{~between~}0\mathrm{~and~}x_2 \\
 * y_2 + y_1(x-x_2)/x_1 & x\mathrm{~between~}x_2\mathrm{~and~}x_1+x_2
 * \end{array}\right.\nonumber
 * \f}
 *
 * If \f$x_1\f$ and \f$x_2\f$ have opposite sign, the domain is \f$[x_1,x_2]\f$ if
 * \f$x_1<0\f$, and \f$[x_2,x_1]\f$ otherwise.  The boundary functions are:
 * \f{eqnarray}{
 * y_a(x) & = & \left\{\begin{array}{c@{\qquad}c}
 * y_1x/x_1 & x\mathrm{~between~}0\mathrm{~and~}x_1 \\
 * y_2x/x_2 & x\mathrm{~between~}0\mathrm{~and~}x_2
 * \end{array}\right.\nonumber\\
 * y_b(x) & = & \left\{\begin{array}{c@{\qquad}c}
 * y_1 + y_2(x-x_1)/x_2 & x\mathrm{~between~}x_1\mathrm{~and~}x_1+x_2 \\
 * y_2 + y1(x-x_2)/x_1  & x\mathrm{~between~}x_2\mathrm{~and~}x_1+x_2
 * \end{array}\right.\nonumber
 * \f}
 *
 * The main program sorts the input parameters so that \f$x_1\leq x_2\f$,
 * stores them in a 4-dimensional array, and assigns a \c void
 * pointer to that array.  It also computes the domain.  The routine
 * LALTwoDRangeTest() takes a value of \f$x\f$ and the \c void
 * pointer, computes the values of \f$y_a(x)\f$ and \f$y_b(x)\f$ according to the
 * algorithm above, sorts them, and returns them ordered from lower to
 * higher.
 *
 * \par Metric values:
 * The main program takes the input parameters
 * \c a, \c b, \c c, \c dadx, \c dbdx, and
 * \c dcdx, stores them in a 9-dimensional array, and assigns a
 * \c void pointer to it.  The routine LALTwoDMetricTest()
 * takes a position \f$(x,y)\f$ and the \c void pointer, and computes the
 * ``local'' value of the principal axis
 * \f$a=\f$\c a\f$+x\times\f$\c dadx\f$+y\times\f$\c dady, and similarly
 * for \f$b\f$ and \f$c\f$.  If that ellipse corresponds to the
 * \f$m_\mathrm{thresh}\f$ mismatch level contour, then the eigenvalues of
 * the corresponding metric are \f$\lambda_1=m_\mathrm{thresh}/a^2\f$ and
 * \f$\lambda_2=m_\mathrm{thresh}/b^2\f$.  The metric components are thus:
 * \f{eqnarray}{
 * g_{xx} & = & \lambda_1\cos^2(c) + \lambda_2\sin^2(c) \;,\nonumber\\
 * g_{yy} & = & \lambda_1\sin^2(c) + \lambda_2\cos^2(c) \;,\nonumber\\
 * g_{xy} \quad = \quad g_{yx} & = & (\lambda_1-\lambda_2)\cos(c)\sin(c)
 * \;.\nonumber
 * \f}
 * The routine assumes that the values of \f$a\f$, \f$b\f$, and \f$c\f$ refer to an
 * \f$m_\mathrm{thresh}=1\f$ mismatch ellipse.  It computes and returns
 * \f$g_{xx}\f$, \f$g_{yy}\f$, and \f$g_{xy}\f$ in a 3-dimensional array.
 *
 * \par PostScript flags:
 * The parameter \c flags is an
 * unsigned integer whose lowest-order bits contain parameters to be
 * passed to LALPlotTwoDMesh().  The bits and their meanings are:
 * <dl>
 * <dt>bit 0:</dt><dd> 1 if mesh points will be plotted, 0 otherwise.</dd>
 * <dt>bit 1:</dt><dd> 1 if mesh tiles will be plotted, 0 otherwise.</dd>
 * <dt>bit 2:</dt><dd> 1 if mismatch ellipses will be plotted, 0 otherwise.</dd>
 * <dt>bit 3:</dt><dd> 1 if the boundary will be plotted, 0 otherwise.</dd>
 * </dl>
 * Thus a value of 15 will plot everything, while a value of 9 will just
 * plot the mesh points and the boundary.  A value of zero suppresses the
 * plot.
 *
 * If mesh points are to be plotted, they will be filled circles \f$1/72''\f$
 * (1~point) in diameter.  The parameter space will be rotated so that
 * the longer of the diagonals of the parallelogram will be vertical, and
 * scaled to fit on one \f$8.5''\times11''\f$ page.  That is, if
 * \f$||(x_1+x_2,y_1+y_2)||\geq||(x_1-x_2,y_1-y_2)||\f$, the rotation angle
 * of the coordinate axes will be
 * \f$\theta=\pi/2-\arctan\!2(y_1+y_2,x_1+x_2)\f$, or
 * \f$\theta=\pi/2-\arctan\!2(y_2-y_1,x_2-x_1)\f$ otherwise.  We note that
 * the function \f$\arctan\!2(y,x)\f$ returns the argument of the complex
 * number \f$x+iy\f$ in the range \f$[-\pi,\pi]\f$.
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * LALPrintError()                 LALCheckMemoryLeaks()
 * LALCreateTwoDMesh()             LALDestroyTwoDMesh()
 * LALSReadGrid()                  LALSDestroyGrid()
 * LALPlotTwoDMesh()
 * \endcode
 *
 * ### Notes ###
 *
 */

/** \name Error Codes */ /*@{*/
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
#define TWODMESHTESTC_MSGEMETRIC "Axis length is zero or negative within specified region"
/*@}*/

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Grid.h>
#include <lal/StreamInput.h>
#include <lal/TwoDMesh.h>
#include "TwoDMeshPlot.h"

/** \cond DONT_DOXYGEN */

/* Default parameter settings. */
#define X1 (1.0)
#define Y1 (0.0)
#define X2 (0.0)
#define Y2 (1.0)
#define A_DEFAULT (0.1)
#define B_DEFAULT (0.05)
#define C_DEFAULT (1.0)
#define DADX (0.0)
#define DBDX (0.0)
#define DCDX (0.0)
#define DADY (0.0)
#define DBDY (0.0)
#define DCDY (0.0)
#define MISMATCH (1.0)

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-p psfile flags] [-d debug]\n" \
"\t[-m mismatch nmax cmax] [-b x1 y1 x2 y2 ] [-e a b c]\n"        \
"\t[-x dadx dbdx dcdx] [-y dady dbdy dcdy] [-i metricfile rangefile]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$", statement ? statement :    \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                       \
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
  INT4 arg;                     /* argument counter */
  static LALStatus stat;        /* top-level status structure */
  REAL4 rangeParams[4];         /* LALRangeTest() params. */
  REAL4 metricParams[9];        /* LALMetricTest() params. */
  REAL4Grid *metricGrid = NULL; /* LALInterpolateMetric() params. */
  REAL4Grid *rangeGrid = NULL;  /* LALInterpolateRangeGrid() params. */
  TwoDMeshParamStruc params;    /* LALCreateTwoDMesh() params. */
  TwoDMeshNode *mesh = NULL;    /* head of mesh list */

  /* Command-line arguments. */
  CHAR *outfile = NULL;                       /* output filename */
  CHAR *psfile = NULL;                        /* PostScript filename */
  UINT2 flags = 0;                            /* PostScript flags */
  CHAR *metricfile = NULL, *rangefile = NULL; /* input filenames */
  REAL4 mismatch = MISMATCH;                  /* maximum mismatch */
  UINT4 nmax = 0, cmax = 0;                   /* maximum nodes/columns */
  REAL4 x_1 = X1, y_1 = Y1, x_2 = X2, y_2 = Y2;   /* boundary params. */
  REAL4 a = A_DEFAULT, b = B_DEFAULT, c = C_DEFAULT; /* ellipse params. */
  REAL4 dadx = DADX, dbdx = DBDX, dcdx = DCDX;       /* ellipse x gradient */
  REAL4 dady = DADY, dbdy = DBDY, dcdy = DCDY;       /* ellipse y gradient */


  /******************************************************************
   * ARGUMENT PARSING                                               *
   ******************************************************************/

  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse PostScript file option. */
    else if ( !strcmp( argv[arg], "-p" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	psfile = argv[arg++];
	flags = atoi( argv[arg++] );
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse maximum numbers option. */
    else if ( !strcmp( argv[arg], "-m" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	nmax = atoi( argv[arg++] );
	cmax = atoi( argv[arg++] );
	mismatch = atof( argv[arg++] );
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse boundary parameters option. */
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 4 ) {
	arg++;
	x_1 = atof( argv[arg++] );
	y_1 = atof( argv[arg++] );
	x_2 = atof( argv[arg++] );
	y_2 = atof( argv[arg++] );
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
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
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
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
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
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
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Parse metric and range grid input option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 2 ) {
	arg++;
	metricfile = argv[arg++];
	rangefile = argv[arg++];
      } else {
	ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return TWODMESHTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( TWODMESHTESTC_EARG, TWODMESHTESTC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return TWODMESHTESTC_EARG;
    }
  } /* End of argument parsing loop. */

  /******************************************************************
   * SETUP (INTERNAL METRIC/RANGE FUNCTIONS)                        *
   ******************************************************************/

  if ( !metricfile ) {
    REAL4 axisMin, axisTemp; /* Min. ellipse size, and a temp value */

    /* Set up range function and parameters. */
    if ( ( x_1 == 0.0 ) && ( x_2 == 0.0 ) ) {
      ERROR( TWODMESHTESTC_EBAD, TWODMESHTESTC_MSGEBAD,
	     "x_1 = x_2 = 0:" );
      return TWODMESHTESTC_EBAD;
    }
    if ( x_1 < x_2 ) {
      rangeParams[0] = x_1;
      rangeParams[1] = y_1;
      rangeParams[2] = x_2;
      rangeParams[3] = y_2;
    } else {
      rangeParams[0] = x_2;
      rangeParams[1] = y_2;
      rangeParams[2] = x_1;
      rangeParams[3] = y_1;
    }
    params.getRange = LALRangeTest;
    params.rangeParams = (void *)( rangeParams );
    if ( x_1*x_2 <= 0.0 ) {
      if ( x_1 < x_2 ) {
	params.domain[0] = x_1;
	params.domain[1] = x_2;
      } else {
	params.domain[0] = x_2;
	params.domain[1] = x_1;
      }
    } else {
      if ( x_1 < 0.0 ) {
	params.domain[0] = x_1 + x_2;
	params.domain[1] = 0.0;
      } else {
	params.domain[0] = 0.0;
	params.domain[1] = x_1 + x_2;
      }
    }

    /* Check that metric will be positive everywhere in the region. */
    axisMin = a;
    axisTemp = a + dadx*x_1 + dady*y_1;
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    axisTemp = a + dadx*x_2 + dady*y_2;
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    axisTemp = a + dadx*( x_1 + x_2 ) + dady*( y_1 + y_2 );
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    if ( axisMin <= 0.0 ) {
      ERROR( TWODMESHTESTC_EMETRIC, TWODMESHTESTC_MSGEMETRIC,
	     "axis a:" );
      return TWODMESHTESTC_EBAD;
    }
    axisTemp = b;
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    axisTemp = b + dbdx*x_1 + dbdy*y_1;
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    axisTemp = b + dbdx*x_2 + dbdy*y_2;
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    axisTemp = b + dbdx*( x_1 + x_2 ) + dbdy*( y_1 + y_2 );
    if ( axisTemp < axisMin ) axisMin = axisTemp;
    if ( axisMin <= 0.0 ) {
      ERROR( TWODMESHTESTC_EMETRIC, TWODMESHTESTC_MSGEMETRIC,
	     "axis b:" );
      return TWODMESHTESTC_EBAD;
    }

    /* Set up metric function and parameters. */
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
  }

  /******************************************************************
   * SETUP (METRIC/RANGE GRID)                                      *
   ******************************************************************/

  else {
    FILE *fp = fopen( metricfile, "r" ); /* input file pointer */
    if ( !fp ) {
      ERROR( TWODMESHTESTC_EFILE, "- " TWODMESHTESTC_MSGEFILE,
	     metricfile );
      return TWODMESHTESTC_EFILE;
    }
    SUB( LALSReadGrid( &stat, &metricGrid, fp ), &stat );
    fclose( fp );
    fp = fopen( rangefile, "r" );
    if ( !fp ) {
      ERROR( TWODMESHTESTC_EFILE, "- " TWODMESHTESTC_MSGEFILE,
	     rangefile );
      return TWODMESHTESTC_EFILE;
    }
    SUB( LALSReadGrid( &stat, &rangeGrid, fp ), &stat );
    fclose( fp );
    params.getMetric = LALInterpolateMetricGrid;
    params.metricParams = (void *)( metricGrid );
    params.getRange = LALInterpolateRangeGrid;
    params.rangeParams = (void *)( rangeGrid );
    params.domain[0] = rangeGrid->offset->data[0];
    params.domain[1] = rangeGrid->offset->data[0]
      + rangeGrid->interval->data[0]
      *( rangeGrid->data->dimLength->data[0] - 1 );
  }

  /******************************************************************
   * MESH CREATION AND OUTPUT                                       *
   ******************************************************************/

  /* Set up remaining mesh creation parameters. */
  params.mThresh = mismatch;
  params.widthMaxFac = 0.0;
  params.widthRetryFac = 0.0;
  params.maxColumns = cmax;
  params.nIn = nmax;

  /* Create mesh. */
  SUB( LALCreateTwoDMesh( &stat, &mesh, &params ), &stat );

  /* Print mesh list to a file, if requested. */
  if ( outfile ) {
    FILE *fp = fopen( outfile, "w" ); /* output file pointer */
    TwoDMeshNode *here;               /* current node in mesh list */

    if ( !fp ) {
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
    FILE *fp = fopen( psfile, "w" );  /* PostScript file pointer */
    INT2 plotPoints = flags & 1;      /* whether to plot mesh points */
    BOOLEAN plotTiles = flags & 2;    /* whether to plot mesh tiles */
    BOOLEAN plotEllipses = flags & 4; /* whether to plot mesh ellipses */
    TwoDMeshPlotStruc plotParams;     /* plotting parameter structure */

    if ( !fp ) {
      ERROR( TWODMESHTESTC_EFILE, "- " TWODMESHTESTC_MSGEFILE,
	     psfile );
      return TWODMESHTESTC_EFILE;
    }

    /* For a parallelogram region specified on the command line, find
       the rotation angle for best fit.  Otherwise, just plot it
       straight. */
    if ( !metricfile ) {
      REAL4 theta1, theta2;
      REAL4 xSum = x_1 + x_2, xDiff = x_2 - x_1;
      REAL4 ySum = y_1 + y_2, yDiff = y_2 - y_1;
      theta1 = LAL_180_PI*atan2( (REAL4)( TWODMESHPLOTH_YSIZE ),
				 (REAL4)( TWODMESHPLOTH_XSIZE ) );
      if ( xSum*xSum + ySum*ySum >= xDiff*xDiff + yDiff*yDiff )
	theta1 -= LAL_180_PI*atan2( ySum, xSum );
      else
	theta1 -= LAL_180_PI*atan2( yDiff, xDiff );
      theta2 = theta1 + LAL_180_PI*atan2( y_1, x_1 );
      while ( theta2 < -180.0 ) theta2 += 360.0;
      while ( theta2 > 180.0 ) theta2 -= 360.0;
      if ( theta2 > 90.0 )
	theta1 -= theta2 - 90.0;
      else if ( ( -90.0 < theta2 ) && ( theta2 < 0.0 ) )
	theta1 -= theta2 + 90.0;
      theta2 = theta1 + LAL_180_PI*atan2( y_2, x_2 );
      while ( theta2 < -180.0 ) theta2 += 360.0;
      while ( theta2 > 180.0 ) theta2 -= 360.0;
      if ( theta2 > 90.0 )
	theta1 -= theta2 - 90.0;
      else if ( ( -90.0 < theta2 ) && ( theta2 < 0.0 ) )
	theta1 -= theta2 + 90.0;
      plotParams.theta = theta1;
    } else
      plotParams.theta = 0.0;

    /* Set remaining parameters, and make the plot. */
    plotParams.xScale = plotParams.yScale = 100.0;
    plotParams.bBox[0] = 36.0;
    plotParams.bBox[1] = 36.0;
    plotParams.bBox[2] = 576.0;
    plotParams.bBox[3] = 756.0;
    plotParams.autoscale = 1;
    memset( plotParams.clipBox, 0, 4*sizeof(REAL4) );
    plotParams.nLevels = 1;
    if ( flags & 8 )
      plotParams.nBoundary = 2 +
	(UINT4)( ( params.domain[1] - params.domain[0] )/a );
    else
      plotParams.nBoundary = 0;
    plotParams.plotPoints = &plotPoints;
    plotParams.plotTiles = &plotTiles;
    plotParams.plotEllipses = &plotEllipses;
    plotParams.params = &params;
    SUB( LALPlotTwoDMesh( &stat, fp, mesh, &plotParams ), &stat );
    fclose( fp );
  }

  /* Free memory, and exit. */
  if ( metricfile ) {
    SUB( LALSDestroyGrid( &stat, &metricGrid ), &stat );
    SUB( LALSDestroyGrid( &stat, &rangeGrid ), &stat );
  }
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

  /* NOTE: It is assumed and required that xy[0] <= xy[2]. */

  INITSTATUS(stat);
  ASSERT( range, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  ASSERT( params, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* Case 1: one of the side vectors is vertical. */
  if ( xy[0] == 0.0 ) {
    ya = xy[3]*( x/xy[2] );
    yb = ya + xy[1];
  }
  else if ( xy[2] == 0.0 ) {
    ya = xy[1]*( x/xy[0] );
    yb = ya + xy[3];
  }

  /* Case 2: Both side vectors point in the same direction (either
     left or right). */
  else if ( ( xy[0] > 0.0 ) || ( xy[2] < 0.0 ) ) {
    if ( fabs( x ) < fabs( xy[0] ) )
      ya = xy[1]*( x/xy[0] );
    else
      ya = xy[1] + xy[3]*( ( x - xy[0] )/xy[2] );
    if ( fabs( x ) < fabs( xy[2] ) )
      yb = xy[3]*( x/xy[2] );
    else
      yb = xy[3] + xy[1]*( ( x - xy[2] )/xy[0] );
  }

  /* Case 3: The first side vector points left and the second points
     right.  (The reverse is impossible, given our assumed ordering.)  */
  else {
    if ( x < 0.0 )
      ya = xy[1]*( x/xy[0] );
    else
      ya = xy[3]*( x/xy[2] );
    if ( x - xy[0] - xy[2] < 0.0 )
      yb = xy[1] + xy[3]*( ( x - xy[0] )/xy[2] );
    else
      yb = xy[3] + xy[1]*( ( x - xy[2] )/xy[0] );
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

  INITSTATUS(stat);
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

/** \endcond */
