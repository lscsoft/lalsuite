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
\author Creighton, T. D.
\file
\ingroup FlatMesh_h
\brief Creates a template mesh for an arbitrary but constant \f$n\f$-dimensional mismatch metric.

\heading{Program <tt>FlatMeshTest.c</tt>}

\heading{Usage}
\code
FlatMeshTest [-o outfile] [-d debuglevel] [-m mismatch]
             [eigenvectorfile inversefile rangefile]
\endcode

\heading{Description}

This test program creates a template mesh for a parameter space with a
constant mismatch metric.  The following option flags are accepted:
<ul>
<li><b>-o</b> Writes the output mesh to the file \c outfile.</li>
<li><b>-d</b> Sets the debug level to \c debuglevel.</li>
<li><b>-m</b> Sets the maximum allowed mismatch to
\c mismatch, a positive number less than 1.</li>
</ul>
Once the above options are processed, any remaining command-line
arguments must be the names of three files containing information
about the eigenvectors of the metric and the desired search range;
these files are described below.  They are read using the function
LALSReadVectorSequence().  If the <b>-o</b> option is not
specified, results are written to \c stdout; if other options or
arguments are not specified, the information is taken from
<tt>\#define</tt>d constants.

\heading{\c eigenvectorfile:} This file contains the
eigenvectors of the \f$n\f$-dimensional mismatch metric \f$\mathsf{g}_{ab}\f$
described in \ref FlatMesh_h.  The file format is simply \f$n\f$ lines
each containing \f$n\f$ whitespace-separated numbers in any standard
floating-point format.  Each line lists the components of a particular
eigenvector; the eigenvector must be normalized so that its squared
magnitude is 1 over the corresponding eigenvalue.

\heading{\c inversefile:} This file also consists of \f$n\f$ lines
each with \f$n\f$ floating-point numbers.  It is simply the matrix inverse
of the contents of \c eigenvectorfile taken as an \f$n\times n\f$
matrix.

\heading{\c rangefile:} This file consists of two lines of \f$n\f$
floating-point numbers; these specify two opposite corners of a
rectilinear region in parameter space to be covered by the mesh.
Additional lines will be ignored.


\heading{Algorithm}

For the most part this test program simply reads the input arguments
and files, passes them to the function LALCreateFlatMesh()
using LALRectIntersect() to define the parameter-space
boundary, and prints the resulting mesh.  However, there are two
additional bits of processing that deserve comment.

The rows of the matrix in \c eigenvectorfile are already of the
form \f$\mathsf{e}^i_{(j)}/\sqrt{\lambda_{(j)}}\f$, as discussed in
\ref FlatMesh_h.  To get the proper orthonormalized transformation
matrix, one must simply multiply each element by
\f$2m_\mathrm{thresh}/\sqrt{n}\f$.  Similarly, the inverse transformation
matrix elements should be \e divided by this number.

In order to ensure \e complete coverage of the desired parameter
space, \c FlatMeshTest extends the boundaries of the rectilinear
region specified in \c rangefile to include any mesh point whose
patch volume touches on the desired search region.  If
\f$\mathsf{M}^a{}_b\f$ is the renormalized transformation matrix described
above, then the sum of the magnitudes of the components along a
column, \f$\Delta x_j=\sum_i|M^i{}_j|\f$ represents the maximum extent of
a mesh point's patch in the \f$j^\mathrm{th}\f$ dimension.  The algorithm
in \c FlatMeshTest extends the rectangular search region by half
this amount in each direction to ensure that any patch touching on the
desired search volume is included.  This assumes that the boundary of
the search region is ``soft''; i.e.\ that no harm will come of
stepping slightly outside it.

\heading{Uses}
\code
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALCalloc()                     LALFree()
LALCreateFlatMesh()             LALSReadVectorSequence()
LALSCreateVectorSequence()      LALSDestroyVectorSequence()
LALSCreateVector()              LALSDestroyVector()
\endcode

\heading{Notes}

*/

/** \name Error Codes */ /*@{*/
#define FLATMESHTESTC_ENORM 0
#define FLATMESHTESTC_ESUB  1
#define FLATMESHTESTC_EARG  2
#define FLATMESHTESTC_EMEM  3
#define FLATMESHTESTC_EDIM  4
#define FLATMESHTESTC_ELEN  5
#define FLATMESHTESTC_EFILE 6

#define FLATMESHTESTC_MSGENORM "Normal exit"
#define FLATMESHTESTC_MSGESUB  "Subroutine failed"
#define FLATMESHTESTC_MSGEARG  "Error parsing arguments"
#define FLATMESHTESTC_MSGEMEM  "Memory allocation error"
#define FLATMESHTESTC_MSGEDIM  "Inconsistent parameter space dimension"
#define FLATMESHTESTC_MSGELEN  "Too few points specified"
#define FLATMESHTESTC_MSGEFILE "Could not open file"
/*@}*/

#include <math.h>
#include <stdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FlatMesh.h>
#include <lal/StreamInput.h>

/** \cond DONT_DOXYGEN */

/* Default parameter settings. */
#define MISMATCH 0.1
#define DIM 2
REAL4 defaultMatrix[] = { 1.0, 0.0, 0.0, 1.0 };
REAL4 defaultMatrixInv[] = { 1.0, 0.0, 0.0, 1.0 };
REAL4 defaultCorners[] = { 0.0, 0.0, 0.5, 1.0 };

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel] [-m mismatch]\n\
                    [eigenvectorfile inversefile rangefile]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, "$Id$", statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 "$Id$", (statement) );                    \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( FLATMESHTESTC_ESUB, FLATMESHTESTC_MSGESUB,            \
         "Function call \"" #func "\" failed:" );                    \
  return FLATMESHTESTC_ESUB;                                      \
}                                                                    \
else (void)(0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

int
main(int argc, char **argv)
{
  INT4 arg;
  UINT4 dim;                 /* dimension of parameter space */
  static LALStatus stat;     /* top-level status structure */
  CHAR *outfile = NULL;      /* name of output file */
  CHAR *eigenfile = NULL;    /* name of eigenvector matrix file */
  CHAR *inversefile = NULL;  /* name of inverse matrix file */
  CHAR *rangefile = NULL;    /* name of parameter range file */
  REAL4 mismatch = MISMATCH; /* maximum mismatch level */
  REAL4VectorSequence *matrix = NULL;    /* tranformation matrix */
  REAL4VectorSequence *matrixInv = NULL; /* inverse tranformation */
  REAL4VectorSequence *corners = NULL;   /* corners of serach region */
  REAL4VectorSequence *mesh = NULL;      /* mesh of parameter values */
  FILE *fp;                  /* input/output file pointer */


  /* Parse argument list.  arg stores the current position. */
  arg = 1;
  while ( arg < argc ) {
    /* Parse output file option. */
    if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	outfile = argv[arg++];
      }else{
	ERROR( FLATMESHTESTC_EARG, FLATMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return FLATMESHTESTC_EARG;
      }
    }
    /* Parse debug level option. */
    else if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
      }else{
	ERROR( FLATMESHTESTC_EARG, FLATMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return FLATMESHTESTC_EARG;
      }
    }
    /* Parse mismatch level option. */
    else if ( !strcmp( argv[arg], "-m" ) ) {
      if ( argc > arg + 1 ) {
	arg++;
	mismatch = fabs( atof( argv[arg++] ) );
      }else{
	ERROR( FLATMESHTESTC_EARG, FLATMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return FLATMESHTESTC_EARG;
      }
    }
    /* Check for unrecognized options. */
    else if ( argv[arg][0] == '-' ) {
      ERROR( FLATMESHTESTC_EARG, FLATMESHTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return FLATMESHTESTC_EARG;
    }
    /* Parse remaining parameters. */
    else {
      if ( argc > arg + 2 ) {
	eigenfile = argv[arg++];
	inversefile = argv[arg++];
	rangefile = argv[arg++];
      }else{
	ERROR( FLATMESHTESTC_EARG, FLATMESHTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return FLATMESHTESTC_EARG;
      }
    }
  } /* End of argument parsing loop. */

  /* If input files have been specified... */
  if ( eigenfile ) {

    /* Read input files into vector sequences. */
    if ( !( fp = LALOpenDataFile( eigenfile ) ) ) {
      ERROR( FLATMESHTESTC_EFILE, "- " FLATMESHTESTC_MSGEFILE,
	     eigenfile );
      return FLATMESHTESTC_EFILE;
    }
    SUB( LALSReadVectorSequence( &stat, &matrix, fp ), &stat );
    fclose( fp );
    if ( !( fp = LALOpenDataFile( inversefile ) ) ) {
      ERROR( FLATMESHTESTC_EFILE, "- " FLATMESHTESTC_MSGEFILE,
	     inversefile );
      return FLATMESHTESTC_EFILE;
    }
    SUB( LALSReadVectorSequence( &stat, &matrixInv, fp ), &stat );
    fclose( fp );
    if ( !( fp = LALOpenDataFile( rangefile ) ) ) {
      ERROR( FLATMESHTESTC_EFILE, "- " FLATMESHTESTC_MSGEFILE,
	     rangefile );
      return FLATMESHTESTC_EFILE;
    }
    SUB( LALSReadVectorSequence( &stat, &corners, fp ), &stat );
    fclose( fp );

    /* Determine dimension, and check consistency. */
    dim = matrix->length;
    if ( matrix->vectorLength != dim ) {
      ERROR( FLATMESHTESTC_EDIM, FLATMESHTESTC_MSGEDIM,
	     eigenfile );
      return FLATMESHTESTC_EDIM;
    }
    if ( ( matrixInv->length != dim ) ||
	 ( matrixInv->vectorLength != dim ) ) {
      ERROR( FLATMESHTESTC_EDIM, FLATMESHTESTC_MSGEDIM, inversefile );
      return FLATMESHTESTC_EDIM;
    }
    if ( corners->vectorLength != dim ) {
      ERROR( FLATMESHTESTC_EDIM, FLATMESHTESTC_MSGEDIM, rangefile );
      return FLATMESHTESTC_EDIM;
    }
    if ( corners->length < 2 ) {
      ERROR( FLATMESHTESTC_ELEN, FLATMESHTESTC_MSGELEN, rangefile );
      return FLATMESHTESTC_ELEN;
    }

  } else {
    /* If no input files are specified, get data from #defined
       constants.  No consistency checks are required. */
    CreateVectorSequenceIn in;
    dim = DIM;
    in.length = dim;
    in.vectorLength = dim;
    SUB( LALSCreateVectorSequence( &stat, &matrix, &in ), &stat );
    memcpy( matrix->data, defaultMatrix, dim*dim*sizeof(REAL4) );
    SUB( LALSCreateVectorSequence( &stat, &matrixInv, &in ), &stat );
    memcpy( matrixInv->data, defaultMatrixInv, dim*dim*sizeof(REAL4) );
    in.length = 2;
    SUB( LALSCreateVectorSequence( &stat, &corners, &in ), &stat );
    memcpy( corners->data, defaultCorners, 2*dim*sizeof(REAL4) );
  }

  /* Apply mismatch threshold to the transformation matrices. */
  {
    UINT4 i;
    REAL4 adjust = 2.0*mismatch/sqrt( (REAL4)(dim) );
    REAL4 *data;  /* pointer to matrix data */

    i = matrix->length*matrix->vectorLength;
    data = matrix->data;
    while ( i-- )
      *(data++) *= adjust;

    adjust = 1.0/adjust;
    i = matrixInv->length*matrixInv->vectorLength;
    data = matrixInv->data;
    while ( i-- )
      *(data++) *= adjust;
  }

  /* Extend the range boundary to ensure edge coverage. */
  {
    UINT4 i, j;    /* indecies */
    INT2 direct;  /* sign of direction from first corner to second */
    REAL4 *data;  /* pointer to matrix data */
    REAL4 *width; /* maximum width of a patch in each dimension */

    /* Allocate local memory. */
    width = (REAL4 *)LALCalloc( dim, sizeof(REAL4) );
    if ( !width ) {
      ERROR( FLATMESHTESTC_EMEM, FLATMESHTESTC_MSGEMEM, 0 );
      return FLATMESHTESTC_EMEM;
    }

    /* Determine patch width. */
    for ( data = matrix->data, i = 0; i < dim; i++ )
      for ( j = 0; j < dim; j++, data++ )
	width[j] += fabs( *data );

    /* Extend each corner by 0.5*width in the appropriate
       direction. */
    for ( data = corners->data, i = 0; i < dim; i++, data++ ) {
      direct = ( data[0] < data[dim] ) ? -1 : 1;
      data[0] += 0.5*direct*width[i];
      data[dim] -= 0.5*direct*width[i];
    }

    /* Free local memory. */
    LALFree( width );
  }

  /* Generate mesh using LALFlatMesh() and LALRectIntersect(). */
  {
    /* Set up parameter structure for LALFlatMesh. */
    static FlatMeshParamStruc params;
    params.matrix = matrix;
    params.matrixInv = matrixInv;
    params.controlPoints = corners;
    params.intersection = LALRectIntersect;
    SUB( LALSCreateVector( &stat, &(params.xMin), dim ), &stat );
    memcpy( params.xMin->data, corners->data, dim*sizeof(REAL4) );
    SUB( LALSCreateVector( &stat, &(params.xMax), dim ), &stat );
    memcpy( params.xMax->data, corners->data+dim, dim*sizeof(REAL4) );

    /* Compute the mesh, and clean up local memory. */
    SUB( LALCreateFlatMesh( &stat, &mesh, &params ), &stat );
    SUB( LALSDestroyVector( &stat, &(params.xMin) ), &stat );
    SUB( LALSDestroyVector( &stat, &(params.xMax) ), &stat );
    SUB( LALSDestroyVectorSequence( &stat, &matrix ), &stat );
    SUB( LALSDestroyVectorSequence( &stat, &matrixInv ), &stat );
    SUB( LALSDestroyVectorSequence( &stat, &corners ), &stat );
  }

  /* Print result. */
  {
    UINT4 i;
    REAL4 *data;
    if ( outfile ) {
      if ( !( fp = fopen( outfile, "w" ) ) ) {
        ERROR( FLATMESHTESTC_EFILE, "- " FLATMESHTESTC_MSGEFILE,
            outfile );
        return FLATMESHTESTC_EFILE;
      }
    } else
      fp = stdout;
    i = mesh->length;
    data = mesh->data;
    while ( i-- ) {
      UINT4 j = mesh->vectorLength;
      while ( j-- )
        fprintf( fp, "%10.3e ", *(data++) );
      fprintf( fp, "\n" );
    }
    if ( outfile )
      fclose( fp );
  }

  /* Free the mesh, and exit. */
  SUB( LALSDestroyVectorSequence( &stat, &mesh ), &stat );
  LALCheckMemoryLeaks();
  INFO( FLATMESHTESTC_MSGENORM );
  return FLATMESHTESTC_ENORM;
}

/** \endcond */
