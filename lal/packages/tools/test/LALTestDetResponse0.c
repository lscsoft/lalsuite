/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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

/* <lalVerbatim file="LALTestDetResponse0CV">
   Author: David Chin <dwchin@umich.edu> +1-734-709-9119
   $Id$
   </lalVerbatim>
*/

/*
<lalLaTeX>

\subsection{Program {\texttt{LALTestDetResponse0.c}}
\label{ss:LALTestDetResponse0.c}

\subsubsection*{Usage}

\subsubsection*{Description}

Performs zeroth-order test of \texttt{LALComputeDetAMResponse()} and
\texttt{LALComputeDetAMResponseSeries()}.

\subsubsection*{Exit codes}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALComputeDetAMResponse()
\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX>
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <lal/LALConfig.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/Units.h>

#include <lal/PrintFTSeries.h>
#include <lal/StreamOutput.h>

NRCSID( LALTESTDETRESPONSE0C, "$Id$" );

#define FALSE 0
#define TRUE  1
#define LALDR_MATRIXSIZE 3

#define TESTDR_MIN(a, b) (((a) < (b)) ? (a) : (b))

/* these two constants are for the sky grid */
#define NUM_DEC 21
#define NUM_RA  24

const INT4 lim = NUM_RA * NUM_DEC;
const INT4 declim = (NUM_DEC-1)/2;

typedef REAL8 LALDR_3Vector[3];
typedef REAL8 LALDR_33Matrix[3][3];

typedef REAL4 skygrid_t[NUM_RA * NUM_DEC];



/*
 * globals
 */
int        lalDebugLevel = 0;
BOOLEAN    verbose_p     = FALSE;
int        verbose_level = 0;
const INT4 oneBillion    = 1000000000;

const REAL8 zero_tolerance = 0.;
REAL8 real8_tolerance = (REAL8)2.*LAL_REAL8_EPS;
REAL4 real4_tolerance = (REAL4)2.*LAL_REAL4_EPS;

LALDetector det_north_pole;
LALDetector det_south_pole;
LALDetector det_green_equator;
LALDetector det_green_tropic_of_cancer;
LALDetector det_foo_tropic_of_cancer;

/* source specified by RA, Dec; the numbers represent degrees */
LALSource   src_0_0_p;
LALSource   src_0_0_c;
LALSource   src_0_90_p;
LALSource   src_0_90_c;
LALSource   src_0_45_p;
LALSource   src_0_45_c;

/* this #if is so I can use Emacs's hide-ifdef-mode to make this
   block invisible, making this file a little easier to scroll through */
#if 1
/* static void REAL4VectorSubtraction(REAL4Vector * pA,
                                   REAL4Vector * pB,
                                   REAL4Vector *pAminusB); */
/*static REAL4 REAL4VectorRMS(REAL4Vector *pVector); */
static void PrintLALDetector(LALDetector * const detector);
static void PrintDetResponse(const LALDetAMResponse * const response,
                             const char * const title);


/* static void make_me_an_Sarray_sequence(LALStatus *status,
                                       REAL4ArraySequence **sequence,
                                       UINT4 rows, UINT4 cols, UINT4 length); */

/* static void print_diagnostics(LALStatus * const status,
                              REAL4ArraySequence *sequence); */


static BOOLEAN almost_equal_real4_p(REAL4 a, REAL4 b, REAL4 tolerance);

static BOOLEAN almost_equal_real8_p(REAL8 a, REAL8 b, REAL8 tolerance);

static BOOLEAN almost_equal_real4_relative_p(REAL4 computed, REAL4 expected,
                                             REAL4 tolerance);

static BOOLEAN matrix_ok_p(LALDR_33Matrix * computed,
                           LALDR_33Matrix * expected,
                           REAL8 tolerance);
static BOOLEAN vector_ok_p(LALDR_3Vector * computed,
                           LALDR_3Vector * expected,
                           REAL8 tolerance);

static BOOLEAN vector_relative_ok_p(LALDR_3Vector * computed,
                                    const LALDR_3Vector * expected,
                                    REAL8 tolerance);

static void print_m_results_maybe(const char * title,
                                  LALDR_33Matrix * computed,
                                  LALDR_33Matrix * expected);

static void print_v_results_maybe(const char * title,
                                  LALDR_3Vector * computed,
                                  LALDR_3Vector * expected);

static void print_s_results_maybe(const char * title,
                                  REAL8 computed, REAL8 expected);

static int print_separator_maybe(void);

static int print_small_separator_maybe(void);

static int print_passed_maybe(void);


static BOOLEAN detresponse_ok_p(LALStatus * status,
                                LALDetAndSource * det_and_src,
                                LALGPSandAcc * gps_and_acc,
                                LALDetAMResponse * expected_resp,
                                REAL4 tolerance);

static void handle_detresponse_test(BOOLEAN passed_p, int line);

static BOOLEAN frdetector_ok_p(LALFrDetector * computed,
                               const LALFrDetector * expected);

static BOOLEAN detector_ok_p(LALDetector * computed,
                             const LALDetector * expected);
#endif

static REAL8 deg_to_rad(REAL8 degrees)
{
  return degrees * (REAL8)LAL_PI / (REAL8)180.;
}

static REAL8 rad_to_deg(REAL8 radians)
{
  return radians * (REAL8)180. / (REAL8)LAL_PI;
}


/* axis for LALDR_EulerRotation() */
typedef enum { xAxis = 1, yAxis = 2, zAxis = 3 } LALDR_Axis_t;

/* This #if is so I can hide this block in Emacs's hide-ifdef-mode */
#if 1
static void
LALDR_Set3Vector(LALDR_3Vector * v,
                 REAL8 v1, REAL8 v2, REAL8 v3)
{
  (*v)[0] = v1;
  (*v)[1] = v2;
  (*v)[2] = v3;
}

/*
 * Cross product of two 3-vectors:
 *    result = a x b
 */
static void
LALDR_CrossProd3Vector(LALDR_3Vector * result,
                       LALDR_3Vector * a,
                       LALDR_3Vector * b)
{
  (*result)[0] =  (*a)[1]*(*b)[2] - (*a)[2]*(*b)[1];
  (*result)[1] = -(*a)[0]*(*b)[2] + (*a)[2]*(*b)[0];
  (*result)[2] =  (*a)[0]*(*b)[1] - (*a)[1]*(*b)[0];

  return;
} /* END: LALDR_CrossProd3Vector() */



/*
 * Dot product of two 3-vectors
 */
static REAL8
LALDR_DotProd3Vector(LALDR_3Vector * a,
                     LALDR_3Vector * b)
{
  INT4 i;
  REAL8 result = 0.;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    result += (*a)[i] * (*b)[i];

  return result;
} /* END: LALDR_DotProd3Vector() */



static void
LALDR_OuterProd3Vector(LALDR_33Matrix * a,
                       LALDR_3Vector * u,
                       LALDR_3Vector * v)
{
  INT4 i;
  INT4 j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*a)[i][j] = (*u)[i] * (*v)[j];

  return;
}



/*
 * Scalar product of two 3x3 matrices
 */
static REAL8
LALDR_DotProd33Matrix(LALDR_33Matrix * a, LALDR_33Matrix * b)
{
  INT4 i, j;
  REAL8 result = 0.;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      result += (*a)[i][j] * (*b)[i][j];

  return result;
} /* END: LALDR_DotProd33Matrix() */




/*
 * Sets all elements of a 3x3 matrix
 */
static void
LALDR_Set33Matrix(LALDR_33Matrix * matrix,
                  REAL8 a11, REAL8 a12, REAL8 a13,
                  REAL8 a21, REAL8 a22, REAL8 a23,
                  REAL8 a31, REAL8 a32, REAL8 a33)
{
  (*matrix)[0][0] = a11;  (*matrix)[0][1] = a12;  (*matrix)[0][2] = a13;

  (*matrix)[1][0] = a21;  (*matrix)[1][1] = a22;  (*matrix)[1][2] = a23;

  (*matrix)[2][0] = a31;  (*matrix)[2][1] = a32;  (*matrix)[2][2] = a33;

  return;
} /* END: LALDR_Set33Matrix() */



#if 0 /* NOT USED */
/*
 * Copy matrix source to matrix target
 */
static void
LALDR_Copy33Matrix(LALDR_33Matrix * target, LALDR_33Matrix * source)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*target)[i][j] = (*source)[i][j];

  return;
} /* END: LALDR_Copy33Matrix() */
#endif



/*
 * Zero matrix
 */
static void
LALDR_Zero33Matrix(LALDR_33Matrix * matrix)
{
  LALDR_Set33Matrix(matrix,
                    0., 0., 0.,
                    0., 0., 0.,
                    0., 0., 0.);
  return;
}



/*
 * Matrix multiply
 */
static void
LALDR_Multiply33Matrix(LALDR_33Matrix * product,
                       LALDR_33Matrix * matrixL,
                       LALDR_33Matrix * matrixR)
{
  /* loop counters */
  INT4 i, j, k;

  /*
   * Zero out output matrix
   */
  LALDR_Set33Matrix(product,
                    0., 0., 0.,
                    0., 0., 0.,
                    0., 0., 0.);

  /*
   * Multiply the matrices: matrixL * matrixR
   */
  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (k = 0; k < LALDR_MATRIXSIZE; ++k)
      for (j = 0; j < LALDR_MATRIXSIZE; ++j)
        (*product)[i][k] += (*matrixL)[i][j] * (*matrixR)[j][k];

  return;
}



/*
 * Scalar multiply
 */
static void
LALDR_ScalarMult33Matrix(LALDR_33Matrix * result,
                         REAL8 coefficient,
                         LALDR_33Matrix * matrix)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*result)[i][j] = coefficient * (*matrix)[i][j];

  return;
}



/*
 * Add matrix
 */
static void
LALDR_Add33Matrix(LALDR_33Matrix * result,
                  LALDR_33Matrix * matrix1,
                  LALDR_33Matrix * matrix2)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*result)[i][j] = (*matrix1)[i][j] + (*matrix2)[i][j];

  return;
}



#if 0 /* NOT USED */
/*
 * Subtract matrices (M1 - M2)
 */
static void
LALDR_Subtract33Matrix(LALDR_33Matrix * result,
                       LALDR_33Matrix * matrix1,
                       LALDR_33Matrix * matrix2)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*result)[i][j] = (*matrix1)[i][j] - (*matrix2)[i][j];

  return;
}
#endif



/*
 * Transpose matrix
 */
static void
LALDR_Transpose33Matrix(LALDR_33Matrix * transpose,
                        LALDR_33Matrix * matrix)
{
  INT4 i, j;

  /*
   * Zero out output matrix
   */
  LALDR_Set33Matrix(transpose,
                    0., 0., 0.,
                    0., 0., 0.,
                    0., 0., 0.);

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*transpose)[i][j] = (*matrix)[j][i];

  return;
}



#if 0 /* NOT USED */
/*
 * The L2 norm of a matrix
 */
static REAL8
LALDR_L2Norm33Matrix(LALDR_33Matrix * matrix)
{
    INT4   i, j;
    REAL8  l2norm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            l2norm += (*matrix)[i][j] * (*matrix)[i][j];

    l2norm = sqrt((double)l2norm);

    return l2norm;
}
#endif




#if 1
/*
 * The RMS norm of a matrix: RMS sum of all elements.
 */
static REAL8
LALDR_RMSNorm33Matrix(LALDR_33Matrix * matrix)
{
    INT4   i, j;
    REAL8  rmsnorm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            rmsnorm += (*matrix)[i][j] * (*matrix)[i][j];

    rmsnorm = sqrt((double)rmsnorm / 9.);

    return rmsnorm;
}
#endif



#if 0 /* NOT USED */
/*
 * The "infinity" norm of a matrix: max over all elems
 */
static REAL8
LALDR_InfNorm33Matrix(LALDR_33Matrix * matrix)
{
    INT4 i, j;
    REAL8 infnorm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            if (fabs((double)((*matrix)[i][j])) > infnorm)
                infnorm = fabs((*matrix)[i][j]);

    return infnorm;
}
#endif



/*
 * Print out matrix
 */
static void
LALDR_Print33Matrix(LALDR_33Matrix * matrix,
                    const CHAR  *varname,
                    UINT4        format,
                    FILE        *file,
                    const CHAR  *graph_title)
{
  /* counters */
  INT4 i, j;

  REAL8 max;


  /*
   * Human-readable format
   */
  if (format == 0)
    {
      fprintf(file, "%s:\n", varname);
      for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        {
          for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            fprintf(file, "% 20.14e\t", (*matrix)[i][j]);

          fprintf(file, "\n");
        }
    }


  /*
     * Maple V.3 format: makes a patchcontour plot
     */
  if (format == 1)
    {
      max = 0.;
      for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
          if (fabs((*matrix)[i][j]) > max)
            max = fabs((*matrix)[i][j]);


      /* load linalg and plots packages for Maple; plot options */
      fprintf(file, "with(linalg): with(plots):\n");
      fprintf(file, "setoptions3d(shading=ZHUE,style=PATCHCONTOUR,");
      fprintf(file, "projection=0.5,axes=BOXED,");
      fprintf(file, "title=`%s (abs. max = %10.5e)`):\n",
              graph_title, fabs(max));

      /* input "varname" is name of variable */
      fprintf(file, "%s := array(1 .. 3, 1 .. 3,[",
              varname);

      for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
          {
            /* last entry doesn't have trailing comma */
            if (i == 2 && j == 2)
              fprintf(file, "(3, 3)=%20.14e)",
                      fabs((*matrix)[i][j]));
            else
              fprintf(file, "(%d, %d)=%20.14e,",
                      j, i, fabs((*matrix)[i][j]));
          }

      fprintf(file, "]):\n");
      fprintf(file, "matrixplot(%s);\n", varname);
    }

  fflush(file);

  return;
}



static void
LALDR_Print3Vector(LALDR_3Vector * vector,
                   const CHAR    * varname,
                   FILE          * file)
{
  /* counter */
  INT4 i;

  fprintf(file, "%s:\n", varname);
  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
      fprintf(file, "% 20.14e\t", (*vector)[i]);

  fprintf(file, "\n");
  fflush(file);
  return;
}






#if 0 /* NOT USED */
/*
 * Returns a rotation matrix representing rotation by angle "theta" about
 * axis "axis"
 *
 * theta in RADIANS
 */
static void
LALDR_EulerRotation(LALDR_33Matrix * rotationMatrix,
                    REAL8        theta,
                    LALDR_Axis_t axis)
{
  REAL8 cosTheta, sinTheta;

  cosTheta = cos((double)theta);
  sinTheta = sin((double)theta);

  switch (axis)
    {
    case zAxis:
      LALDR_Set33Matrix(rotationMatrix,
                        cosTheta,  sinTheta, 0.,
                        -sinTheta, cosTheta, 0.,
                        0.,        0.,       1.);
      break;

    case yAxis:
      LALDR_Set33Matrix(rotationMatrix,
                        cosTheta, 0., -sinTheta,
                        0.,       1.,  0.,
                        sinTheta, 0., cosTheta);
      break;

    case xAxis:
      LALDR_Set33Matrix(rotationMatrix,
                        1., 0.,        0.,
                        0.,  cosTheta, sinTheta,
                        0., -sinTheta, cosTheta);
      break;
    }

  return;
}
#endif
#endif


REAL4 skygrid_avg(const skygrid_t response);
void  skygrid_square(skygrid_t square, const skygrid_t input);
REAL4 skygrid_rms(const skygrid_t input);
void  skygrid_sqrt(skygrid_t result, const skygrid_t input);
INT4  skygrid_copy(skygrid_t dest, const skygrid_t src);
void  skygrid_print(const char * comments, const skygrid_t input,
                    const char * filename);
void  skygrid_fabs(skygrid_t absgrid, const skygrid_t input);
void  skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_subtract(skygrid_t sum, const skygrid_t a, const skygrid_t b);
void  skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b);

void  setup_global_detectors(LALStatus *status);
void  set_source_params(LALSource * source, const char *name, REAL8 ra_rad,
                        REAL8 dec_rad, REAL8 orien_rad);
void  print_source_maybe(const LALSource * source);
void  setup_global_sources(void);
void  crab_pulsar_test(LALStatus * status);

void find_zero_gmst(LALStatus * status);

typedef enum
  {
    gwpol_scalar = 0,
    gwpol_plus  = 1,
    gwpol_cross = 2
  }
GWPolarization;

/* response function in local horizon coordinates (from GRASP)
 * psi = source orientation
 * theta = source polar angle = Pi/2 - altitude
 * phi = source azimuth
 * pol = selects which response, F+ or Fx, to compute */
REAL4 resp_local(REAL8 psi, REAL8 theta, REAL8 phi, GWPolarization pol);

/* error-handled fopen */
FILE *xfopen(const char *path, const char *mode);
int  xfclose(FILE *stream);

/* wrapped laldr_strlcpy() to guarantee NUL termination */
char *laldr_strlcpy(char *dst, const char *src, size_t len);

/* static int local_strncasecmp(const char *, const char *, size_t); */

/*
 * Test modules
 */
void fudge_factor_test(LALStatus *status);
BOOLEAN passed_special_locations_tests_p(LALStatus *status);
BOOLEAN passed_almost_equal_tests_p(void);

/* FIXME */
BOOLEAN passed_matrix_test_p(void);

/* Yes, I do mean for the following block to be #if'ed out */
#if 0
/*
 * This computes the effective longitude and latitude of a detector, given
 * an LALDetector.  It stuffs its results in a SkyPosition structure in
 * RADIANS.
 */
static void
LALDR_GetEffectiveLoc(SkyPosition *eff_loc, const LALDetector *detector)
{
    REAL8 theta_x, phi_x;
    REAL8 theta_y, phi_y;
    REAL8 sinTheta;

    LALDR_3Vector e_xarm;
    LALDR_3Vector e_yarm;

    LALDR_3Vector normal;

    SkyPosition delta_loc;


    /*
     * First, we need the unit vectors representing the arm directions
     */

    /* polar angle, x-arm */
    theta_x = LAL_PI_2 - detector->frDetector.xArmAltitudeRadians;

    /* azimuthal angle. I want it measure anti-clockwise from S, the LAL
     * convention measures anti-clockwise of E */
    phi_x   = LAL_PI_2 + detector->frDetector.xArmAzimuthRadians;

    /* and the unit vector is: */
    sinTheta  = sin(theta_x);
    e_xarm[0] = sinTheta * cos(phi_x);
    e_xarm[1] = sinTheta * sin(phi_x);
    e_xarm[2] = cos(theta_x);

    /* polar angle, y-arm */
    theta_y = LAL_PI_2 - detector->frDetector.yArmAltitudeRadians;

    /* azimuthal angle. I want it measure anti-clockwise from S, the LAL
     * convention measures anti-clockwise of E */
    phi_y   = LAL_PI_2 + detector->frDetector.yArmAzimuthRadians;

    /* and the unit vector is: */
    sinTheta  = sin(theta_y);
    e_yarm[0] = sinTheta * cos(phi_y);
    e_yarm[1] = sinTheta * sin(phi_y);
    e_yarm[2] = cos(theta_y);

    /* so, the normal to the plane of the detector is */
    LALDR_CrossProd3Vector(normal, e_xarm, e_yarm);

    /* then, polar angles of this normal vector will tell me the residuals
     * of the latitude and longitude, though not directly:
     *
     * say   a = decrease in longitude
     *       b = decrease in latitude
     *
     * an active rotation formed by
     *
     *       R_x = [ 1    0        0
     *               0  cos(a)  -sin(a)
     *               0  sin(a)   cos(a) ]
     *
     *       R_y = [ cos(b)  0   sin(b)
     *                 0     1     0
     *              -sin(b)  0   cos(b) ]
     *
     * and acting on e_z = [0 0 1] should give me the normal vector as computed
     * by doing the crossproduct of the arms above.  We then solve for the
     * angles a and b:
     *
     *  normal = R_y * R_x * e_z
     *         = [ sin(b)*cos(a)    -sin(a)    cos(b)*cos(a) ]
     *
     *  so,
     *         a = arcsin(-normal[2])
     *         b = arctan(normal[1] / normal[3])
     */
    delta_loc.longitude = asin(-normal[1]);
    delta_loc.latitude  = atan2(normal[0], normal[2]);

    /* And thus, the effective location is */
    eff_loc->longitude = deg_to_rad(detector->frDetector.vertexLongitudeDegrees)
        - delta_loc.longitude;
    eff_loc->latitude  = deg_to_rad(detector->frDetector.vertexLatitudeDegrees)
        - delta_loc.latitude;

    return;

} /* END: GetEffectiveLoc() */
#endif



int main(int argc, char *argv[])
{
  static LALStatus  status;
  LALSource         pulsar;
  LALFrDetector     frdet;    /* Framelib detector info */
  LALDetector       detector;
  LIGOTimeGPS       gps;
  LALDate           utcDate;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
  LALGPSandAcc      gps_and_acc;
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;
  LALDetAMResponse  expected_resp;

  LALDetAMResponseSeries    am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series,
    circ_series, sum_series;
  REAL4                     mean;
  /* REAL4Vector               diffVector; */
  LALTimeIntervalAndNSample time_info;

  LALTimeInterval           interval;

  INT4  k;
  INT4  i, j, count;
  INT4  cnt;
  REAL4 tolerance;

  REAL8 tmpgmst;
  REAL8 gmst1;
  LALMSTUnitsAndAcc tmp_uandacc;

  skygrid_t plus;
  skygrid_t cross;
  skygrid_t sqsum;
  skygrid_t plus_sq_time_avg;
  skygrid_t cross_sq_time_avg;
  skygrid_t sum_of_sq_time_avg;
  skygrid_t tmpskygrid;
  skygrid_t tmpskygrid2;
  FILE     *file_plus_sq_avg = NULL;
  FILE     *file_cross_sq_avg = NULL;
  FILE     *file_plus_at_0_0 = NULL;
  FILE     *file_cross_at_0_0 = NULL;
  FILE     *file_plus_at_2_10 = NULL;
  FILE     *file_cross_at_2_10 = NULL;
  FILE     *file_cross_at_4_15 = NULL;
  FILE     *file_plus_at_4_15 = NULL;
  FILE     *file_cross_at_m4_15 = NULL;
  FILE     *file_plus_at_m4_15 = NULL;
  FILE     *file_sum_sq_avg = NULL;
  FILE     *file_sum_sq = NULL;
  FILE     *file_theta = NULL;
  FILE     *file_phi   = NULL;

  if (argc >= 3)
    {
      lalDebugLevel = atoi(argv[1]);
      verbose_level = atoi(argv[2]);
      verbose_p     = verbose_level;
    }

  setup_global_detectors(&status);
  setup_global_sources();


  /* this section is just for finding out a time when GMST is equal to 0 */
  /* and the answer is:  GPS = 13675020:943750000 */
  if (lalDebugLevel == 69)
    {
      find_zero_gmst(&status);
      goto conclusion;
    }


  /* "crab pulsar"  test, to compare with numbers produced by Greg Mendell
   * and Malik Rakhmanov */
  if (lalDebugLevel == 13)
    {
      lalDebugLevel = 1;  /* avoid verbose debug output */
      crab_pulsar_test(&status);
      goto conclusion;
    }

  if (!passed_almost_equal_tests_p())
    {
      exit(2);
    }

  /*
   * TEST 0: Test of matrix/vector manipulations
   */
  if (!passed_matrix_test_p())
    {
      fprintf(stderr, "Matrix test failed");
      exit(3);
    }

  /*********************************************/

  if (verbose_p)
    {
      printf("TEST OF LALCreateDetector()\n");
      printf("---------------------------\n\n");
      printf("Manual inspection required.\n");
    }

  (void)laldr_strlcpy(frdet.name, "Reference", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(0.);
  frdet.vertexLatitudeRadians  = deg_to_rad(0.);
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(90.);
  frdet.yArmAzimuthRadians     = deg_to_rad(0.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (verbose_p)
    {
      PrintLALDetector(&detector);
    }

  print_separator_maybe();

  /*********************************************/

  /* parameters come from LIGO-T980044-10 (Althouse, et al.)
   * and we note that the suffix number 10 is greater than A....
   * I know consistency is the hobgoblin etc. etc. but this is a bit
   * ridiculous.
   *
   * The referenced paper gives:
   *   LHO -- x-arm azimuth  = N 35.9994deg W
   *                altitude = -6.195e-04 rad
   *       -- y-arm azimuth  = S 54.0006deg W
   *                altitude = +1.25e-05 rad
   *
   *   LLO -- x-arm azimuth  = S 72.2835deg W
   *                altitude = -3.121e-04 rad
   *       -- y-arm azimuth  = S 17.7165deg E
   *                altitude = -6.107e-04 rad
   *
   * Which gives us the following in conventional bearing notation,
   * i.e. degrees clockwise from North
   *
   * 1. LHO -- x-arm azimuth = 324.0006deg =
   *           y-arm azimuth = 234.0006deg =
   *    (check: x-arm azi - y-arm azi = 90.0000deg)
   *
   * 2. LLO -- x-arm azimuth = 252.2835deg =
   *           y-arm azimuth = 162.2835deg =
   *    (check: x-arm azi - y-arm azi = 90.0000deg)
   *
   */

  /* First, pass a LALFrDetector using Frame spec for azimuths (East of
     North) */
  (void)laldr_strlcpy(frdet.name, "LHO, from FrDetector struct (Frame spec)",
                LALNameLength);
  frdet.vertexLongitudeRadians = (REAL8)deg_to_rad(-119. - 25./60. - 27.5657/3600.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(46. + 27./60. + 18.528/3600.);
  frdet.vertexElevation        = 142.554;
  frdet.xArmAltitudeRadians    = -6.195e-04;
  frdet.yArmAltitudeRadians    =  1.250e-05;
  frdet.xArmAzimuthRadians     = deg_to_rad(324.0006);
  frdet.yArmAzimuthRadians     = deg_to_rad(234.0006);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (!detector_ok_p(&detector, &(lalCachedDetectors[LALDetectorIndexLHODIFF])))
    {
      if (verbose_p)
        {
          fprintf(stderr, "WARNING: LHO computed w/ Frame spec != LHO cached\n");
        }
    }


  if (verbose_p)
    {
      printf("LHO tensor converted from LALFrDetector (Frame spec):\n");
      PrintLALDetector(&detector);
      printf("\n- - - - -\n");
    }

  /* Second, pass a LALFrDetector using LAL (package-tools) spec for azi
     (counterclockwise from East) */
  (void)laldr_strlcpy(frdet.name, "LHO, from FrDetector struct (package-tools spec)",
          LALNameLength);
  frdet.vertexLongitudeRadians  = (REAL8)deg_to_rad(-119. - 25./60. - 27.5657/3600.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(46. + 27./60. + 18.528/3600.);
  frdet.vertexElevation        = 142.554;
  frdet.xArmAltitudeRadians    = -6.195e-04;
  frdet.yArmAltitudeRadians    =  1.250e-05;
  frdet.xArmAzimuthRadians     = deg_to_rad(125.9994);
  frdet.yArmAzimuthRadians     = deg_to_rad(215.9994);
  /* check: y-arm azi - x-arm azi = 90deg */

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (!detector_ok_p(&detector, &(lalCachedDetectors[LALDetectorIndexLHODIFF])))
    {
      if (verbose_p)
        fprintf(stderr, "WARNING: LHO computed w/ package-tools spec != LHO cached\n");
    }


  if (verbose_p)
    {
      printf("LHO tensor converted from LALFrDetector (package-tools spec):\n");
      PrintLALDetector(&detector);
      printf("\n- - -\n");
    }

  /* Third, use data from the CreateDetector (DetectorSite.h) doco */
  (void)laldr_strlcpy(frdet.name, "LHO, from FrDetector struct (numbers from doco)",
          LALNameLength);
  frdet.vertexLongitudeRadians  = (REAL8)deg_to_rad(-119. - 25./60. - 27.5657/3600.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(46. + 27./60. + 18.528/3600.);
  frdet.vertexElevation        = 142.554;
  frdet.xArmAltitudeRadians    = -6.195e-04;
  frdet.yArmAltitudeRadians    =  1.250e-05;
  frdet.xArmAzimuthRadians     = deg_to_rad(324.0006);
  frdet.yArmAzimuthRadians     = deg_to_rad(234.0006);
  /* check: y-arm azi - x-arm azi = -90deg ; ok, i think. */

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (!detector_ok_p(&detector, &(lalCachedDetectors[LALDetectorIndexLHODIFF])))
    {
      if (verbose_p)
        fprintf(stderr,
                "WARNING: LHO computed w/ numbers from LAL documentation != LHO cached\n");
    }

  if (verbose_p)
    {
      printf("LHO tensor converted from LALFrDetector (numbers from doco):\n");
      PrintLALDetector(&detector);
      printf("\n- - -\n");
    }



  /* Fourth, look at the cached detector for comparison */
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];

  if (verbose_p)
    {
      printf("LHO tensor, cached by DetectorSite:\n");
      PrintLALDetector(&detector);
    }

  print_separator_maybe();


  /*********************************************/

  /*
   * Let's try to make a trivial detector...
   * My "greenwich_equator" from the old AM code is at location (0,0),
   * with arms pointing south and west.
   */
  (void)laldr_strlcpy(frdet.name, "TRIVIAL 1", LALNameLength);
  frdet.vertexLongitudeRadians = 0.;
  frdet.vertexLatitudeRadians  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(180.);
  frdet.yArmAzimuthRadians     = deg_to_rad( 90.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr,
              "LALTestDetResponse0: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }


  if (verbose_p)
    {
      printf("TRIVIAL 1 (converted from FrDetector):\n");
      PrintLALDetector(&detector);
    }

  print_separator_maybe();


  /*********************************************/

  tolerance = 0.05;

  if (verbose_p)
    {
      printf("TEST OF LALDetAMResponse()\n");
      printf("--------------------------\n\n");
      printf("Tolerance =  % 15.8e\n\n", tolerance);
    }

  /*
   * TEST 1: test LALDetAMResponse() at particular points
   */

  /* Set a GPS time that's close to 0h GMST1. (Found this by trial and
   * error.)  */
  gps.gpsSeconds     =     61094;
  gps.gpsNanoSeconds = 640000000;
  gps_and_acc.gps = gps;
  gps_and_acc.accuracy = accuracy;

  /* Set up a source at (RA=0, Dec=0, orientation=0, at time GMST1=0) */
  (void)laldr_strlcpy(pulsar.name, "TEST PULSAR 1", LALNameLength);
  pulsar.equatorialCoords.longitude = 0.;  /* RA */
  pulsar.equatorialCoords.latitude  = 0.;  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = LAL_PI_2;  /* orientation */


  /* Stuff Detector and Source into structure */
  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*** expect (1, 0) */
  expected_resp.plus   = 1.;
  expected_resp.cross  = 0.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (0, -1) */
  pulsar.orientation = -LAL_PI_4;
  expected_resp.plus   =  0.;
  expected_resp.cross  = -1.;
  expected_resp.scalar =  0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (-1, 0) */
  pulsar.orientation = 0;
  expected_resp.plus   = -1.;
  expected_resp.cross  =  0.;
  expected_resp.scalar =  0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (0.5, -sqrt(3)/2.) */
  pulsar.orientation = -LAL_PI/3.;
  expected_resp.plus   =  0.5;
  expected_resp.cross  = -sqrt(3.)/2.;
  expected_resp.scalar =  0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /* switch detector to something less trivial */
  (void)laldr_strlcpy(frdet.name, "TRIVIAL 2", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(0.);
  frdet.vertexLatitudeRadians  = deg_to_rad(15.);
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(180.);
  frdet.yArmAzimuthRadians     = deg_to_rad( 90.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (verbose_p)
    PrintLALDetector(&detector);

  /* switch source to be overhead the detector */
  (void)laldr_strlcpy(pulsar.name, "TEST PULSAR 2", LALNameLength);
  pulsar.equatorialCoords.longitude = deg_to_rad(0.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(15.);
  pulsar.orientation                = -LAL_PI_2;

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*** expect (1, 0 ) */
  expected_resp.plus = 1.;
  expected_resp.cross = 0.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (0, -1) */
  pulsar.orientation = -LAL_PI_2/2.;
  expected_resp.plus = 0.;
  expected_resp.cross = -1.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (-1, 0) */
  pulsar.orientation = 0.;
  expected_resp.plus = -1.;
  expected_resp.cross = 0.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (0.5, -sqrt(3)/2) */
  pulsar.orientation = -LAL_PI/3.;
  expected_resp.plus = 0.5;
  expected_resp.cross = -sqrt(3.)/2.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();


  /*** expect () */
  pulsar.orientation = -LAL_PI_2;
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude = 0.;

  (void)laldr_strlcpy(frdet.name, "TRIVIAL 1", LALNameLength);
  frdet.vertexLongitudeRadians = 0.;
  frdet.vertexLatitudeRadians  = LAL_PI_2;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(180.);
  frdet.yArmAzimuthRadians     = deg_to_rad( 90.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr,
              "LALTestDetResponse0: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  utcDate.unixDate.tm_sec = 46;
  utcDate.unixDate.tm_min = 20;
  utcDate.unixDate.tm_hour = 8;
  utcDate.unixDate.tm_mday = 17;
  utcDate.unixDate.tm_mon  = LALMONTH_MAY;
  utcDate.unixDate.tm_year = 1994 - 1900;

  /*  accuracy = LALLEAPSEC_LOOSE; */
  LALUTCtoGPS(&status, &gps, &utcDate, &accuracy);

  tmp_uandacc.units = MST_RAD;
  tmp_uandacc.accuracy = accuracy;

  LALGPStoGMST1(&status, &tmpgmst, &gps, &tmp_uandacc);

  if (verbose_p)
    printf("GMST1 = % 14.9e rad.\n", tmpgmst);

  gps_and_acc.gps = gps;
  gps_and_acc.accuracy = accuracy;

  expected_resp.plus = 0.5;
  expected_resp.cross = 0.;
  expected_resp.scalar = 0.;

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  if (verbose_p)
  {
    PrintLALDetector(&detector);
    fflush(stdout);
  }

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp, tolerance),
                          __LINE__);

  print_small_separator_maybe();



  /* switch detector to LHO */
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];

  /* switch source */
  (void)laldr_strlcpy(pulsar.name, "TEST PULSAR 3", LALNameLength);
  pulsar.equatorialCoords.longitude = deg_to_rad(16.037547 * 15.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(46.475430);
  pulsar.orientation                = -LAL_PI_2;

  utcDate.unixDate.tm_sec = 46;
  utcDate.unixDate.tm_min = 20;
  utcDate.unixDate.tm_hour = 8;
  utcDate.unixDate.tm_mday = 17;
  utcDate.unixDate.tm_mon  = LALMONTH_MAY;
  utcDate.unixDate.tm_year = 1994 - 1900;

  accuracy = LALLEAPSEC_LOOSE;
  LALUTCtoGPS(&status, &gps, &utcDate, &accuracy);

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /* expect (3.08260644404358e-01, -9.51301793267616e-01 */
  expected_resp.plus  =  3.08260644404358e-01;
  expected_resp.cross = -9.51301793267616e-01;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps_and_acc, &expected_resp,
                                           tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /* change to even less trivial detector */
  (void)laldr_strlcpy(frdet.name, "TRIVIAL 3", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(15.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (verbose_p)
    PrintLALDetector(&detector);


  /* now have to choose a source whose RA is 15 degrees towards the East;
     whaddaya know? that's 1 hour */
  pulsar.equatorialCoords.longitude = deg_to_rad(15.);

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*
   * Compute a time series AM response
   */
  if (verbose_p)
    printf("Starting vector test\n");

  /* fake detector */
  detector.location[0] = 0.;
  detector.location[1] = 0.;
  detector.location[2] = LAL_AWGS84_SI;
  detector.response[0][0] = 0.;
  detector.response[1][1] = 0.5;
  detector.response[2][2] = -0.5;
  detector.response[0][1] = detector.response[1][0] = 0.;
  detector.response[0][2] = detector.response[2][0] = 0.;
  detector.response[1][2] = detector.response[2][1] = 0.;
  detector.type = LALDETECTORTYPE_ABSENT;
  (void)laldr_strlcpy(detector.frDetector.name, "FAKE", LALNameLength);


  /* Make a fake detector by specifying frame format detector */
  (void)laldr_strlcpy(frdet.name, "FAKE FAKE, NOT THE REAL MCCOY", LALNameLength);
  frdet.vertexLongitudeRadians = (REAL8)deg_to_rad(0.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(90.); /* @ N pole */
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(90.);
  frdet.yArmAzimuthRadians     = deg_to_rad( 0.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  /* override - try the cached LHO */
  /* detector = lalCachedDetectors[LALDetectorIndexLHODIFF]; */

  if (verbose_p)
    {
      printf("Series test...\n");
      PrintLALDetector(&detector);
    }

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;
  circ_series.data = NULL;
  sum_series.data = NULL;

  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(&status, &(am_response_series.pCross->data), 1);
  LALSCreateVector(&status, &(am_response_series.pScalar->data), 1);
  LALSCreateVector(&status, &(circ_series.data), 1);
  LALSCreateVector(&status, &(sum_series.data), 1);

  if (lalDebugLevel > 0)
    {
      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCros->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
      printf("circ_series.data->length = %d\n", circ_series.data->length);
      printf("sum_series.data->length = %d\n", sum_series.data->length);
    }

  time_info.epoch.gpsSeconds     = 61094;
  time_info.epoch.gpsNanoSeconds = 640000000;
  time_info.deltaT               = 60;
  /* time_info.nSample              = 17*24*60; */
  time_info.nSample              = 24*60;
  time_info.accuracy             = LALLEAPSEC_STRICT;

  LALComputeDetAMResponseSeries(&status,
                                &am_response_series,
                                &det_and_pulsar,
                                &time_info);

  if (status.statusCode && verbose_p)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponseSeries, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      return status.statusCode;
    }

  if (lalDebugLevel > 0)
    {
      printf("Done computing AM response vectors\n");

      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCross->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);

      printf("TimeSeries data written to files plus_series.txt, ");
      printf("cross_series.txt, and scalar_series.txt\n");

      LALSPrintTimeSeries(am_response_series.pPlus, "plus_series.txt");
      LALSPrintTimeSeries(am_response_series.pCross, "cross_series.txt");
      LALSPrintTimeSeries(am_response_series.pScalar, "scalar_series.txt");

      /* compute circular response */
      /* need to resize vector */
      LALSDestroyVector(&status, &(circ_series.data));
      LALSCreateVector(&status, &(circ_series.data),
                       am_response_series.pPlus->data->length);
      printf("circ_series.data->length = %d\n", circ_series.data->length);

      LALSDestroyVector(&status, &(sum_series.data));
      LALSCreateVector(&status, &(sum_series.data),
                       am_response_series.pPlus->data->length);
      printf("sum_series.data->length = %d\n", sum_series.data->length);

      circ_series.epoch = am_response_series.pPlus->epoch;
      circ_series.deltaT = am_response_series.pPlus->deltaT;
      circ_series.f0 = am_response_series.pPlus->f0;
      circ_series.sampleUnits = lalDimensionlessUnit;

      sum_series.epoch = am_response_series.pPlus->epoch;
      sum_series.deltaT = am_response_series.pPlus->deltaT;
      sum_series.f0 = am_response_series.pPlus->f0;
      sum_series.sampleUnits = lalDimensionlessUnit;

      mean = 0.;
      for (k = 0; k < (INT4)(circ_series.data->length); ++k)
        {
          /* sqrt(Fplus^2 + Fcross^2) */
          circ_series.data->data[k] =
            sqrt(am_response_series.pPlus->data->data[k] *
                 am_response_series.pPlus->data->data[k] +
                 am_response_series.pCross->data->data[k] *
                 am_response_series.pCross->data->data[k]);

          /* (Fplus + Fcross)^2 */
          sum_series.data->data[k] =
            (am_response_series.pPlus->data->data[k]
             + am_response_series.pCross->data->data[k])
            * (am_response_series.pPlus->data->data[k]
               + am_response_series.pCross->data->data[k]);

          mean += sum_series.data->data[k];
        }

      mean /= (REAL8)(sum_series.data->length);

      printf("mean = % 14.7e\n", mean);

      LALSPrintTimeSeries(&circ_series, "circ_series.txt");
      LALSPrintTimeSeries(&sum_series, "sum_series.txt");
    }

  if (lalDebugLevel & 8)
    {
      printf("plus: (");
      for (k = 0; k < (INT4)(time_info.nSample); ++k)
        {
          printf("% 14.6e, ", am_response_series.pPlus->data->data[k]);
        }
      printf(")\n");

      printf("cross: (");
      for (k = 0; k < (INT4)(time_info.nSample); ++k)
        {
          printf("% 14.6e, ", am_response_series.pCross->data->data[k]);
        }
      printf(")\n");

      printf("scalar: (");
      for (k = 0; k < (INT4)(time_info.nSample); ++k)
        {
          printf("% 14.6e, ", am_response_series.pScalar->data->data[k]);
        }
      printf(")\n");


      /* print out quadrature sum of plus- and cross-response */
      printf("sqrt(PLUS^2 + CROSS^2): (");
      for (k = 0; k < (INT4)(time_info.nSample); ++k)
        {
          printf("% 1.6e, ", sqrt(am_response_series.pPlus->data->data[k] *
                                  am_response_series.pPlus->data->data[k] +
                                  am_response_series.pCross->data->data[k] *
                                  am_response_series.pCross->data->data[k]));
        }
      printf(")\n");
    }


  print_separator_maybe();


  /*
   * Loop over whole sky
   */
  if (verbose_p && lalDebugLevel)
    {
      printf("ALOHA\n");
      printf("time_info.nSample = %d\n", time_info.nSample);
    }

  /* use Livingston */
  /* detector = lalCachedDetectors[LALDetectorIndexLLODIFF]; */

  if (lalDebugLevel >= 1)
    {
      gmst1 = 0.;

      printf("\nStarting whole-sky test...\n");
      det_and_pulsar.pDetector = &det_north_pole;
      PrintLALDetector(det_and_pulsar.pDetector);
      count = 0;
      printf("NUM_RA = %d; NUM_DEC = %d\n", NUM_RA, NUM_DEC);

      file_plus_sq_avg    = xfopen("plus_sq_avg.txt", "w");
      file_cross_sq_avg   = xfopen("cross_sq_avg.txt", "w");
      file_plus_at_0_0    = xfopen("plus_at_0_0.txt", "w");
      file_cross_at_0_0   = xfopen("cross_at_0_0.txt", "w");
      file_plus_at_2_10   = xfopen("plus_at_2_10.txt", "w");
      file_cross_at_2_10  = xfopen("cross_at_2_10.txt", "w");
      file_plus_at_4_15   = xfopen("plus_at_4_15.txt", "w");
      file_cross_at_4_15  = xfopen("cross_at_4_15.txt", "w");
      file_plus_at_m4_15  = xfopen("plus_at_m4_15.txt", "w");
      file_cross_at_m4_15 = xfopen("cross_at_m4_15.txt", "w");
      file_sum_sq_avg     = xfopen("sum_sq_avg.txt", "w");
      file_sum_sq         = xfopen("sum_sq.txt", "w");
      file_theta          = xfopen("theta.txt", "w");
      file_phi            = xfopen("phi.txt", "w");

      printf("Done opening files.\n");

      gps.gpsSeconds     = time_info.epoch.gpsSeconds;
      gps.gpsNanoSeconds = time_info.epoch.gpsNanoSeconds;
      LALFloatToInterval(&status, &interval, &(time_info.deltaT));

      printf("Time info set.\n");

      /* Set a GPS time that's close to 0h GMST1. Found this by trial and
       * error:
       * GPS = 13675020:943728537; gmst1 =   2.68743469376486e-10
       * Later, need to use a Science Run time period. */
      gps_and_acc.gps.gpsSeconds     =  13675020;
      gps_and_acc.gps.gpsNanoSeconds = 943728537;
      interval.seconds               =       600;
      interval.nanoSeconds           =         0;

      printf("N sample = %d\n", time_info.nSample);

      pulsar.orientation = deg_to_rad(45.);

      for (k = 0; k < (int)time_info.nSample; ++k)
        {
          LALMSTUnitsAndAcc uandacc;
          uandacc.units    = MST_RAD;
          uandacc.accuracy = gps_and_acc.accuracy;
          LALGPStoGMST1(&status, &gmst1, &(gps_and_acc.gps), &uandacc);

          if (verbose_level & 16)
            printf("GRAR: k = %6d; gmst1 = % 20.14e\n", k, gmst1);

          for (j = 0; j < NUM_RA; ++j)
            {
              pulsar.equatorialCoords.longitude =
                (REAL8)j/(REAL8)NUM_RA * ((REAL8)LAL_TWOPI); /* RA */

              for (i = -declim; i <= declim; ++i)
                {
                  cnt = j*NUM_DEC + i + declim;

                  if (verbose_level & 16)
                    printf("OY: k = %6d; j = %6d; i = %6d\n", k, j, i);

                  pulsar.equatorialCoords.latitude =
                    asin((REAL8)i/(REAL8)declim);

                  if (k == 0 && j == 0 && i == -declim)
                    printf("FOO: gmst1 = % 20.14e\n", gmst1);
                  LALComputeDetAMResponse(&status, &am_response,
                                          &det_and_pulsar, &gps_and_acc);

                  plus[cnt]  = am_response.plus;
                  cross[cnt] = am_response.cross;
                  sqsum[cnt] = (plus[cnt] * plus[cnt])
                    + (cross[cnt] * cross[cnt]);

                  if (i == 0 && j == 0)
                    {
                      fprintf(file_plus_at_0_0, "%4d % 14.9e %9d % 14.9e\n",
                              k, gmst1, gps_and_acc.gps.gpsSeconds, plus[cnt]);
                      fprintf(file_cross_at_0_0, "%4d % 14.9e %9d % 14.9e\n",
                              k, gmst1, gps_and_acc.gps.gpsSeconds,
                              cross[cnt]);
                      fprintf(file_sum_sq, "%4d % 14.9e %9d % 14.9e\n",
                              k, gmst1, gps_and_acc.gps.gpsSeconds,
                              sqsum[cnt]);
                    }

                  if (i == 2 && j == 10)
                    {
                      fprintf(file_plus_at_2_10, "% 14.9e % 14.9e\n",
                              gmst1, plus[cnt]);
                      fprintf(file_cross_at_2_10, "% 14.9e % 14.9e\n",
                              gmst1, cross[cnt]);
                    }

                  if (i == 4 && j == 15)
                    {
                      fprintf(file_plus_at_4_15, "% 14.9e % 14.9e\n", gmst1,
                              plus[cnt]);
                      fprintf(file_cross_at_4_15, "% 14.9e % 14.9e\n",
                              gmst1, cross[cnt]);
                    }

                  if (i == -4 && j == 15)
                    {
                      fprintf(file_plus_at_m4_15, "% 14.9e\n", plus[cnt]);
                      fprintf(file_cross_at_m4_15, "% 14.9e\n", cross[cnt]);
                    }
                }
            }

          /* print out avg. square over sky for each time step */
          skygrid_square(tmpskygrid, plus);
          fprintf(file_plus_sq_avg, "% 14.9e\n", skygrid_avg(tmpskygrid));

          skygrid_square(tmpskygrid2, cross);
          fprintf(file_cross_sq_avg, "% 14.9e\n", skygrid_avg(tmpskygrid2));

          skygrid_scalar_mult(tmpskygrid, tmpskygrid, 1./time_info.nSample);
          skygrid_add(plus_sq_time_avg, plus_sq_time_avg, tmpskygrid);

          skygrid_scalar_mult(tmpskygrid2, tmpskygrid2, 1./time_info.nSample);
          skygrid_add(cross_sq_time_avg, cross_sq_time_avg, tmpskygrid2);

          skygrid_add(sum_of_sq_time_avg, plus_sq_time_avg, cross_sq_time_avg);

          if (k == 0)
            {
              skygrid_print("LAL-computed response to plus", plus,
                            "plus.txt");
              skygrid_print("LAL-computed response to cross", cross,
                            "cross.txt");

              skygrid_square(tmpskygrid, plus);
              skygrid_square(tmpskygrid2, cross);
              skygrid_add(tmpskygrid, tmpskygrid, tmpskygrid2);
              printf("BAR: gmst1 = % 20.15e\n", gmst1);
              printf("avg(F+^2 + Fx^2) = % 14.8e\n", skygrid_avg(tmpskygrid));
            }

          /* increment observation time */
          LALIncrementGPS(&status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                          &interval);
        }
      printf("\n");

      skygrid_print("LAL-computed time avg. of square of response to plus",
                    plus_sq_time_avg, "plus_sq_time_avg.txt");
      skygrid_print("LAL-computed time avg. of square of response to cross",
                    cross_sq_time_avg, "cross_sq_time_avg.txt");
      skygrid_print("LAL-computed time avg. of sum of squares of response to plus and cross",
                    sum_of_sq_time_avg, "sum_of_sq_time_avg.txt");

      fprintf(file_plus_sq_avg, "\n");
      fprintf(file_cross_sq_avg, "\n");
      fprintf(file_plus_at_0_0, "\n");
      fprintf(file_cross_at_0_0, "\n");
      fprintf(file_plus_at_2_10, "\n");
      fprintf(file_cross_at_2_10, "\n");
      fprintf(file_plus_at_4_15, "\n");
      fprintf(file_cross_at_4_15, "\n");
      fprintf(file_plus_at_m4_15, "\n");
      fprintf(file_cross_at_m4_15, "\n");
      fprintf(file_sum_sq_avg, "\n");
      fprintf(file_sum_sq, "\n");

      if (verbose_p)
        printf("avg(F+^2 + Fx^2) = % 14.8e\n", skygrid_avg(sqsum));

      xfclose(file_plus_sq_avg);
      xfclose(file_cross_sq_avg);
      xfclose(file_plus_at_0_0);
      xfclose(file_cross_at_0_0);
      xfclose(file_plus_at_2_10);
      xfclose(file_cross_at_2_10);
      xfclose(file_plus_at_4_15);
      xfclose(file_cross_at_4_15);
      xfclose(file_plus_at_m4_15);
      xfclose(file_cross_at_m4_15);
      xfclose(file_sum_sq_avg);
      xfclose(file_sum_sq);
      xfclose(file_theta);
      xfclose(file_phi);
    }

  if (verbose_p)
  {
    printf("GRASP plus = % 14.9e\n", resp_local(0., 0., 0., gwpol_plus));
    printf("     cross = % 14.9e\n", resp_local(0., 0., 0., gwpol_cross));

    printf("GRASP plus = % 14.9e\n", resp_local(LAL_PI_4, 0., 0., gwpol_plus));
    printf("     cross = % 14.9e\n", resp_local(LAL_PI_4, 0., 0., gwpol_cross));

    fflush(stdout);
  }

  fudge_factor_test(&status);

  if (verbose_p)
    printf("\n\nGOODBYE.\n");

  /*
   * Housekeeping
   */
  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));
  LALSDestroyVector(&status, &(circ_series.data));
  LALSDestroyVector(&status, &(sum_series.data));

 conclusion:
  LALCheckMemoryLeaks();

  return 0;
} /* END: main() */


#if 0 /* NOT USED */
/*
 * subtracts two REAL4Vectors; user must do all allocation beforehand
 */
static void REAL4VectorSubtraction(REAL4Vector *pA,
                                   REAL4Vector *pB,
                                   REAL4Vector *pAminusB)
{
  UINT4 i;

  /* Check for compatible dimensions */
  if ((pA->length != pB->length) || (pAminusB->length != pA->length) ||
      (pAminusB->length != pB->length))
    {
      fprintf(stderr, "VectorSubtraction: ERROR: incompatible dimensions\n");
      exit(13);
    }

  for (i = 0; i < pA->length; ++i)
    {
      pAminusB->data[i] = pA->data[i] - pB->data[i];
    }

    return;
}
#endif

#if 0 /* NOT USED */
static REAL4 REAL4VectorRMS(REAL4Vector *pVector)
{
  UINT4 i;
  REAL4 result = 0.;

  for (i = 0; i < pVector->length; ++i)
    {
      result += pVector->data[i] * pVector->data[i];
    }

  result /= pVector->length;

  return sqrt(result);
}
#endif


static void
PrintLALDetector(LALDetector * const detector)
{
  printf( "Detector  = \n");
  printf("{\n");
  printf("  { % 22.15e, % 22.15e, % 22.15e },\n",
         detector->location[0],
         detector->location[1],
         detector->location[2]);
  printf("  {\n");
  printf("    { % 22.15e, % 22.15e, % 22.15e },\n",
         detector->response[0][0],
         detector->response[0][1],
         detector->response[0][2]);
  printf( "    { % 22.15e, % 22.15e, % 22.15e },\n",
          detector->response[1][0],
          detector->response[1][1],
          detector->response[1][2]);
  printf( "    { % 22.15e, % 22.15e, % 22.15e }\n  },\n",
          detector->response[2][0],
          detector->response[2][1],
          detector->response[2][2]);
  printf("  {\n");
  printf("    \"%s\",\n", detector->frDetector.name);
  printf("    vertex longitude = % 22.15e,\n",
         detector->frDetector.vertexLongitudeRadians);
  printf("    vertex latitude  = % 22.15e,\n",
         detector->frDetector.vertexLatitudeRadians);
  printf("    vertex elevation = % 22.15e,\n",
         detector->frDetector.vertexElevation);
  printf("    X-arm altitude   = % 22.15e,\n",
         detector->frDetector.xArmAltitudeRadians);
  printf("    X-arm azimuth    = % 22.15e,\n",
         detector->frDetector.xArmAzimuthRadians);
  printf("    Y-arm altitude   = % 22.15e,\n",
         detector->frDetector.yArmAltitudeRadians);
  printf("    Y-arm azimuth    = % 22.15e\n  }\n}\n",
         detector->frDetector.yArmAzimuthRadians);
  return;
}




static BOOLEAN matrix_ok_p(LALDR_33Matrix * computed,
                           LALDR_33Matrix * expected,
                           REAL8 tolerance)
{
  INT4 i, j;
  BOOLEAN retval;

  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
      {
        if (!almost_equal_real8_p((REAL8)((*computed)[i][j]),
                                  (REAL8)((*expected)[i][j]),
                                  tolerance))
          {
            if (verbose_p)
              {
                LALDR_Print33Matrix(expected,
                                    "INFO: matrix_ok_p(): expected",
                                    0, stdout, "");
                LALDR_Print33Matrix(computed,
                                    "INFO: matrix_ok_p(): computed",
                                    0, stdout, "");
              }

            retval = FALSE;
          }
        else
          {
            retval = TRUE;
          }
      }

  return retval;
}



static BOOLEAN vector_ok_p(LALDR_3Vector * computed,
                           LALDR_3Vector * expected,
                           REAL8 tolerance)
{
  INT4 i;
  BOOLEAN retval;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    {
      if (!almost_equal_real8_p((*computed)[i], (*expected)[i], tolerance))
        {
          if (verbose_p)
            {
              LALDR_Print3Vector(computed, "INFO: vector_ok_p(): computed",
                                 stdout);
              LALDR_Print3Vector(expected, "INFO: vector_ok_p(): expected",
                                 stdout);
            }

          retval = FALSE;
        }
      else
        {
          retval = TRUE;
        }
    }

  return retval;
}




static BOOLEAN vector_relative_ok_p(LALDR_3Vector * computed,
                                    const LALDR_3Vector * expected,
                                    REAL8 tolerance)
{
  /* do this term-by-term rather than using a metric. bah. */
  INT4 i;
  REAL8 relative_diff;
  BOOLEAN retval;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    {
      relative_diff = fabs((*computed)[i]/(*expected)[i] - 1.);

      if (relative_diff <= tolerance)
        {
          retval = TRUE;
        }
      else
        {
          retval = FALSE;
        }
    }

  return retval;
}




static void print_m_results_maybe(const char * title,
                                  LALDR_33Matrix * computed,
                                  LALDR_33Matrix * expected)
{
  size_t i;
  size_t title_length = strlen(title);

  if (verbose_p)
    {
      printf("%s:\n", title);
      for (i = 0; i < title_length + 1; ++i)
        printf("-");
      printf("\n");
      LALDR_Print33Matrix(computed, "Computed", 0, stdout, "");
      LALDR_Print33Matrix(expected, "Expected", 0, stdout, "");
    }
}



static void print_v_results_maybe(const char * title,
                                  LALDR_3Vector * computed,
                                  LALDR_3Vector * expected)
{
  size_t i;
  size_t title_length = strlen(title);

  if (verbose_p)
    {
      printf("%s:\n", title);
      for (i = 0; i < title_length + 1; ++i)
        printf("-");
      printf("\n");
      LALDR_Print3Vector(computed, "Computed", stdout);
      LALDR_Print3Vector(expected, "Expected", stdout);
    }
}


static void print_s_results_maybe(const char * title,
                                  REAL8 computed,
                                  REAL8 expected)
{
  size_t i;
  size_t title_length = strlen(title);

  if (verbose_p)
    {
      printf("%s:\n", title);
      for (i = 0; i < title_length + 1; ++i)
        printf("-");
      printf("\n");
      printf("Computed = %22.14e\n", computed);
      printf("Expected = %22.14e\n", expected);
    }
}



static int print_separator_maybe(void)
{
  if (verbose_p)
    return printf("\n* * * * * * * * * *\n\n");
  else
    return 0;
}



static int print_small_separator_maybe(void)
{
  if (verbose_p)
    return printf("\n- - - - -\n\n");
  else
    return 0;
}



static int print_passed_maybe(void)
{
  if (verbose_p)
    return printf("\nPASS\n");
  else
    return 0;
}



static BOOLEAN almost_equal_real4_p(REAL4 a, REAL4 b, REAL4 tolerance)
{
  if (lalDebugLevel >= 8)
    {
      printf("almost_equal_real4_p() info:\n");
      printf("  a = % 16.10e\n  b = % 16.10e\n", a, b);
      printf("  (REAL4)fabs(a - b) = % 16.10e\n", (REAL4)fabs(a - b));
      printf("  tolerance = % 16.10e\n\n", tolerance);
    }
  if (tolerance == 0.)
    return (a == b);
  else
    return ((REAL4)fabs(a - b) <= tolerance);
}



static BOOLEAN almost_equal_real8_p(REAL8 a, REAL8 b, REAL8 tolerance)
{
  if (tolerance == 0.)
    return (a == b);
  else
    return ((REAL8)fabs(a - b) <= tolerance);
}



static BOOLEAN almost_equal_real4_relative_p(REAL4 computed, REAL4 expected,
                                             REAL4 tolerance)
{
  REAL4 relative_err;

  if (tolerance == 0.)
    {
      return (computed == expected);
    }
  else
    {
      relative_err = fabs(computed - expected)/fabs(expected);
      if (relative_err <= tolerance)
        return TRUE;
      else
        return FALSE;
    }
}



static BOOLEAN detresponse_ok_p(LALStatus * status,
                                LALDetAndSource * det_and_src,
                                LALGPSandAcc * gps_and_acc,
                                LALDetAMResponse * expected_resp,
                                REAL4 tolerance)

{
  LALDetAMResponse computed_resp;
  BOOLEAN resp_plus_ok_p;
  BOOLEAN resp_cross_ok_p;
  BOOLEAN resp_scalar_ok_p;
  BOOLEAN result;
  REAL4   computed_circ_resp, expected_circ_resp;

  LALComputeDetAMResponse(status, &computed_resp, det_and_src, gps_and_acc);

  expected_circ_resp = sqrt((expected_resp->plus) * (expected_resp->plus)
                            + (expected_resp->cross) * (expected_resp->cross));

  computed_circ_resp = sqrt(computed_resp.plus * computed_resp.plus
                            + computed_resp.cross * computed_resp.cross);


  if (verbose_p)
    {
      PrintDetResponse(&computed_resp, "Computed");
      PrintDetResponse(expected_resp, "Expected");
    }

  if (expected_resp->plus == 0.)
    {
      resp_plus_ok_p = almost_equal_real4_p(computed_resp.plus,
                                            expected_resp->plus,
                                            tolerance);
    }
  else
    {
      resp_plus_ok_p = almost_equal_real4_relative_p(computed_resp.plus,
                                                     expected_resp->plus,
                                                     tolerance);
    }

  if (lalDebugLevel >= 2)
    {
      printf("INFO: detresponse_ok_p(): resp_plus_ok_p = %d\n",
             resp_plus_ok_p);
      printf("                          computed = % 14.10e\n",
             computed_resp.plus);
      printf("                          expected = % 14.10e\n",
             expected_resp->plus);
    }

  if (expected_resp->cross == 0.)
    resp_cross_ok_p = almost_equal_real4_p(computed_resp.cross,
                                           expected_resp->cross,
                                           tolerance);
  else
    resp_cross_ok_p = almost_equal_real4_relative_p(computed_resp.cross,
                                                     expected_resp->cross,
                                                     tolerance);

  if (lalDebugLevel >= 2)
    {
      printf("INFO: detresponse_ok_p(): resp_cross_ok_p = %d\n",
             resp_cross_ok_p);
      printf("                          computed = % 14.10e\n",
             computed_resp.cross);
      printf("                          expected = % 14.10e\n",
             expected_resp->cross);
    }

  resp_scalar_ok_p = almost_equal_real4_p(computed_resp.scalar,
                                          expected_resp->scalar,
                                          tolerance);

  result = (resp_plus_ok_p && resp_cross_ok_p && resp_scalar_ok_p);

  if (result == FALSE && verbose_p)
    {
      printf("INFO: detresponse_ok_p(): circ_resp -- computed = % 10.6e\n",
             computed_circ_resp);
      printf("                                       expected = % 10.6e\n",
             expected_circ_resp);
    }

  return result;
}



void handle_detresponse_test(BOOLEAN passed_p, int line)
{
  if (!passed_p)
    {
      fprintf(stderr, "ERROR in LALComputeDetAMResponse(): computed result != expected result; line %d\n", line);
      exit (1);
    }
  else
    {
      print_passed_maybe();
    }
}




static void PrintDetResponse(const LALDetAMResponse * const response,
                             const char * const title)
{
  REAL4 circ_response = sqrt((response->plus) * (response->plus)
                             + (response->cross) * (response->cross));

  printf("%s:\n", title);
  printf("      plus = % 15.8e  }  circular = % 15.8e\n",
         response->plus, circ_response);
  printf("     cross = % 15.8e  }  \n", response->cross);
  printf("    scalar = % 15.8e\n", response->scalar);
}



static BOOLEAN frdetector_ok_p(LALFrDetector * computed,
                               const LALFrDetector * expected)
{
  /* let's not bother with comparing the name */

  return (almost_equal_real8_p(computed->vertexLongitudeRadians,
                               expected->vertexLongitudeRadians,
                               2.*real8_tolerance)
          && almost_equal_real8_p(computed->vertexLatitudeRadians,
                                  expected->vertexLatitudeRadians,
                                  2.*real8_tolerance)
          && almost_equal_real8_p(computed->vertexElevation,
                                  expected->vertexElevation,
                                  2.*real8_tolerance)
          && almost_equal_real8_p(computed->xArmAltitudeRadians,
                                  expected->xArmAltitudeRadians,
                                  2.*real8_tolerance)
          && almost_equal_real8_p(computed->xArmAzimuthRadians,
                                  expected->xArmAzimuthRadians,
                                  2.*real8_tolerance)
          && almost_equal_real8_p(computed->yArmAltitudeRadians,
                                  expected->yArmAltitudeRadians,
                                  2.*real8_tolerance)
          && almost_equal_real8_p(computed->yArmAzimuthRadians,
                                  expected->yArmAzimuthRadians,
                                  2.*real8_tolerance));
}




static BOOLEAN detector_ok_p(LALDetector * computed,
                             const LALDetector * expected)
{
  LALDR_33Matrix tmp_computed;
  LALDR_33Matrix tmp_expected;

  LALDR_Set33Matrix(&tmp_computed,
                    computed->response[0][0], computed->response[0][1],
                    computed->response[0][2],
                    computed->response[1][0], computed->response[1][1],
                    computed->response[1][2],
                    computed->response[2][0], computed->response[2][1],
                    computed->response[2][2]);

  LALDR_Set33Matrix(&tmp_expected,
                    expected->response[0][0], expected->response[0][1],
                    expected->response[0][2],
                    expected->response[1][0], expected->response[1][1],
                    expected->response[1][2],
                    expected->response[2][0], expected->response[2][1],
                    expected->response[2][2]);


  return (vector_relative_ok_p(&(computed->location), &(expected->location),
                                1.e-4)
          && matrix_ok_p(&tmp_computed, &tmp_expected, 1.e-4)
          && computed->type == expected->type
          && frdetector_ok_p(&(computed->frDetector), &(expected->frDetector)));
}



REAL4 skygrid_avg(const skygrid_t response)
{
  INT4 i, j, cnt;
  REAL4 retval = 0.;

  for (j = 0; j < NUM_RA; ++j)
    {
      for (i = -declim+1; i <= declim-1; ++i)
        {
          cnt = j*NUM_DEC + i + declim;
          retval += response[cnt];
        }
    }

  /*
  for (i = 0; i < lim; ++i)
    retval += response[i];
  */

  retval /= lim;

  return retval;
}




void skygrid_square(skygrid_t square, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    square[i] = (input[i]) * (input[i]);

}




void skygrid_sqrt(skygrid_t result, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    result[i] = (REAL4)sqrt((double)(input[i]));
}



REAL4 skygrid_rms(const skygrid_t input)
{
  skygrid_t tmpgrid;

  skygrid_square(tmpgrid, input);
  return (REAL4)(sqrt(skygrid_avg(tmpgrid)));
}



INT4 skygrid_copy(skygrid_t dest, const skygrid_t src)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    dest[i] = src[i];

  return i;
}



void skygrid_print(const char * comments,
                   const skygrid_t input, const char * filename)
{
  INT4 i, j;
  FILE * outfile = NULL;

  outfile = xfopen(filename, "w");

  if (comments != (char *)NULL)
    fprintf(outfile, "# %s\n", comments);

  for (i = 0; i < NUM_RA; ++i)
    {
      for (j = 0; j < NUM_DEC; ++j)
        fprintf(outfile, "% 14.8e\t", input[i*NUM_DEC + j]);
      fprintf(outfile, "\n");
    }

  xfclose(outfile);
}



void skygrid_fabs(skygrid_t absgrid, const skygrid_t input)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    absgrid[i] = fabs(input[i]);
}



void skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    sum[i] = a[i] + b[i];
}

void skygrid_subtract(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    sum[i] = a[i] - b[i];
}

void skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    result[i] = b * a[i];
}

#if 0 /* NOT USED */
/*
 * only handle 2-dimensional arrays
 */
static void make_me_an_Sarray_sequence(LALStatus *status,
                                       REAL4ArraySequence **sequence,
                                       UINT4 rows, UINT4 cols, UINT4 length)
{
  CreateArraySequenceIn params;
  UINT4Vector           dimLength;
  UINT4                 data[3];
  data[0] = 2;
  data[1] = rows;
  data[2] = cols;
  dimLength.length = 3;
  dimLength.data   = data;
  params.length    = length;
  params.dimLength = &dimLength;

  LALSCreateArraySequence(status, sequence, &params);
}
#endif



#if 0 /* NOT USED */
static void print_diagnostics(LALStatus * const status,
                              REAL4ArraySequence *sequence)
{
  UINT4 i;

  printf("statusCode = %d\n", status->statusCode);
  printf("sequence->length   = %u\n", sequence->length);
  printf("sequence->dimLength->length = %u\n", sequence->dimLength->length);
  for (i = 0; i < sequence->dimLength->length; ++i)
    {
      printf("sequence->dimLength->data[%d] = % 5d\n", i,
             sequence->dimLength->data[i]);
    }
}
#endif


FILE *xfopen(const char *path, const char *mode)
{
  FILE *f = NULL;

  f = fopen(path, mode);

  if (f == NULL)
    {
      fprintf(stderr, "%s: ", path);
      perror("could not open file");
      exit(errno);
    }
  return f;
}

int xfclose(FILE *stream)
{
  if (stream != (FILE *)NULL)
    return fclose(stream);
  else
    return 0;
}

REAL4 resp_local(REAL8 psi, REAL8 theta, REAL8 phi, GWPolarization pol)
{
  REAL8 cos_theta, cos_sq_theta, cos_2_psi, sin_2_psi;
  REAL8 half_cos_sq_theta_p1_cos_2_phi, cos_theta_sin_2_phi;
  REAL8 retval;

  cos_theta = cos(theta);
  cos_sq_theta = cos_theta * cos_theta;
  cos_2_psi = cos(2. * psi);
  sin_2_psi = sin(2. * psi);
  cos_theta_sin_2_phi = cos_theta * sin(2. * phi);
  half_cos_sq_theta_p1_cos_2_phi = (cos_sq_theta + 1.) * cos(2. * phi) / 2.;


  if (pol == gwpol_plus)
    retval = half_cos_sq_theta_p1_cos_2_phi * cos_2_psi
      - cos_theta_sin_2_phi * sin_2_psi;
  else if (pol == gwpol_cross)
    retval = half_cos_sq_theta_p1_cos_2_phi * sin_2_psi
      + cos_theta_sin_2_phi * cos_2_psi;
  else
    retval = 0.;

  return (REAL4)retval;
}



void setup_global_detectors(LALStatus *status)
{
  /*
    LALDetector det_north_pole;
    LALDetector det_south_pole;
    LALDetector det_green_equator;
    LALDetector det_green_tropic_of_cancer;
    LALDetector det_foo_tropic_of_cancer;
  */
  static BOOLEAN global_detectors_set_p = FALSE;
  LALFrDetector frdet;

  if (verbose_level & 4)
    printf("global_detectors_set_p = %d\n", global_detectors_set_p);


  if (global_detectors_set_p == FALSE)
    {
      /* Det. @ North Pole */
      (void)laldr_strlcpy(frdet.name, "North Pole", LALNameLength);
      frdet.vertexLongitudeRadians = deg_to_rad(0.);
      frdet.vertexLatitudeRadians  = deg_to_rad(90.);
      frdet.vertexElevation        = 0.;
      frdet.xArmAltitudeRadians    = 0.;
      frdet.yArmAltitudeRadians    = 0.;
      frdet.xArmAzimuthRadians     = deg_to_rad(180.);
      frdet.yArmAzimuthRadians     = deg_to_rad(90.);

      LALCreateDetector(status, &det_north_pole, &frdet,
                        LALDETECTORTYPE_IFODIFF);

      if (verbose_level & 4)
        {
          PrintLALDetector(&det_north_pole);
          (void)print_small_separator_maybe();
        }

      /* Det. @ South Pole */
      (void)laldr_strlcpy(frdet.name, "South Pole", LALNameLength);
      frdet.vertexLongitudeRadians = deg_to_rad(0.);
      frdet.vertexLatitudeRadians  = deg_to_rad(-90.);
      frdet.vertexElevation        = 0.;
      frdet.xArmAltitudeRadians    = 0.;
      frdet.yArmAltitudeRadians    = 0.;
      frdet.xArmAzimuthRadians     = deg_to_rad(180.);
      frdet.yArmAzimuthRadians     = deg_to_rad(90.);

      LALCreateDetector(status, &det_south_pole, &frdet,
                        LALDETECTORTYPE_IFODIFF);

      if (verbose_level & 4)
        {
          PrintLALDetector(&det_south_pole);
          (void)print_small_separator_maybe();
        }

      /* Det. @ (0, 0) on Earth */
      (void)laldr_strlcpy(frdet.name, "Greenwich-Equator", LALNameLength);
      frdet.vertexLongitudeRadians = deg_to_rad(0.);
      frdet.vertexLatitudeRadians  = deg_to_rad(0.);
      frdet.vertexElevation        = 0.;
      frdet.xArmAltitudeRadians    = 0.;
      frdet.yArmAltitudeRadians    = 0.;
      frdet.xArmAzimuthRadians     = deg_to_rad(180.);
      frdet.yArmAzimuthRadians     = deg_to_rad(90.);

      LALCreateDetector(status, &det_green_equator, &frdet,
                        LALDETECTORTYPE_IFODIFF);

      if (verbose_level & 4)
        {
          PrintLALDetector(&det_green_equator);
          (void)print_small_separator_maybe();
        }

      /* Det. @ (0, 23.5) on Earth */
      (void)laldr_strlcpy(frdet.name, "Greenwich-Tropic of Cancer", LALNameLength);
      frdet.vertexLongitudeRadians = deg_to_rad(0.);
      frdet.vertexLatitudeRadians  = deg_to_rad(23.5);
      frdet.vertexElevation        = 0.;
      frdet.xArmAltitudeRadians    = 0.;
      frdet.yArmAltitudeRadians    = 0.;
      frdet.xArmAzimuthRadians     = deg_to_rad(180.);
      frdet.yArmAzimuthRadians     = deg_to_rad(90.);

      LALCreateDetector(status, &det_green_tropic_of_cancer, &frdet,
                        LALDETECTORTYPE_IFODIFF);

      if (verbose_level & 4)
        {
          PrintLALDetector(&det_green_tropic_of_cancer);
          (void)print_small_separator_maybe();
        }

      /* Det. @ (38.4, 23.5) on Earth */
      (void)laldr_strlcpy(frdet.name, "Foo-Tropic of Cancer", LALNameLength);
      frdet.vertexLongitudeRadians = deg_to_rad(38.4);
      frdet.vertexLatitudeRadians  = deg_to_rad(23.5);
      frdet.vertexElevation        = 0.;
      frdet.xArmAltitudeRadians    = 0.;
      frdet.yArmAltitudeRadians    = 0.;
      frdet.xArmAzimuthRadians     = deg_to_rad(180.);
      frdet.yArmAzimuthRadians     = deg_to_rad(90.);

      LALCreateDetector(status, &det_foo_tropic_of_cancer, &frdet,
                        LALDETECTORTYPE_IFODIFF);

      if (verbose_level & 4)
        {
          PrintLALDetector(&det_foo_tropic_of_cancer);
          print_small_separator_maybe();
        }

      global_detectors_set_p = TRUE;
    }

  return;
} /* END: setup_global_detectors() */



BOOLEAN passed_special_locations_tests_p(LALStatus *status)
{
  status = NULL;
  return TRUE;
}



char *laldr_strlcpy(char *dst, const char *src, size_t len)
{
  char *retval = strncpy(dst, src, len);
  if (verbose_level & 8)
    {
      printf("sizeof(dst) = %lu\n", sizeof(dst));
      printf("strlen(src) = %lu\n", strlen(src));
      printf("TESTDR_MIN(len, strlen(src)+1) = %lu\n",
             TESTDR_MIN(len, strlen(src)+1));
    }
  dst[TESTDR_MIN(len, strlen(src) + 1) - 1] = '\0';
  return retval;
}




void fudge_factor_test(LALStatus *status)
{
  /* compute the response using a local horizon coordinate system */
  LALGPSandAcc      gps_and_acc;
  LALSource         pulsar;
  LALDetAndSource   det_and_pulsar = { (LALDetector *)NULL,
                                       (LALSource *)NULL} ;
  LALMSTUnitsAndAcc uandacc;
  LALDetAMResponse  am_response;

  REAL8 gmst1 = 0.;
  REAL8 altitude = 0.;
  REAL8 azimuth = 0.;
  INT4  cnt   = 0;
  INT4  i, j, k;
  skygrid_t plus;
  skygrid_t cross;
  skygrid_t resp_plus;
  skygrid_t resp_cros;
  skygrid_t tmp1, tmp2, tmp3;
  REAL8     fudge_factor = 0.;
  FILE *rms_diff_plus_file = (FILE *)NULL;
  FILE *rms_diff_cros_file = (FILE *)NULL;

  if (verbose_level & 4)
    printf("\n\nSTARTING FUDGE FACTOR TEST NOW...\n");

  rms_diff_plus_file = xfopen("ff_rms_diff_plus_vs_fudge.txt", "w");
  rms_diff_cros_file = xfopen("ff_rms_diff_cros_vs_fudge.txt", "w");

  /* compute LAL response at pole... */
  gps_and_acc.gps.gpsSeconds     =  13675020;
  gps_and_acc.gps.gpsNanoSeconds = 943728537;
  gps_and_acc.accuracy           = LALLEAPSEC_STRICT;
  uandacc.units    = MST_RAD;
  uandacc.accuracy = gps_and_acc.accuracy;
  LALGPStoGMST1(status, &gmst1, &(gps_and_acc.gps), &uandacc);

  if (verbose_level & 4)
    printf("gmst1 = % 20.14e\n", gmst1);

  set_source_params(&pulsar, "FOOBAR", 0., -LAL_PI_2, (REAL8)LAL_PI_2);

  det_and_pulsar.pDetector = &det_north_pole;
  det_and_pulsar.pSource   = &pulsar;

  if (verbose_level & 4)
    PrintLALDetector(det_and_pulsar.pDetector);
  print_source_maybe(det_and_pulsar.pSource);

  for (j = 0; j < NUM_RA; ++j)
    {
      det_and_pulsar.pSource->equatorialCoords.longitude =
        (REAL8)j/(REAL8)NUM_RA * LAL_TWOPI; /* RA */

      for (i = -declim; i <= declim; ++i)
        {
          cnt = j*NUM_DEC + i + declim;

          det_and_pulsar.pSource->equatorialCoords.latitude =
            asin((REAL8)i/(REAL8)declim);

          LALComputeDetAMResponse(status, &am_response,
                                  &det_and_pulsar, &gps_and_acc);

          plus[cnt]  = am_response.plus;
          cross[cnt] = am_response.cross;
        }
    }

  /* loop over fudge factors */
  if (verbose_p)
    printf("   Starting to loop over fudge_factor...\n");

  for (k = -128; k < 129; ++k)
    {
      fudge_factor = ((double)k) / 256. * (double)(LAL_PI_2)/16.;

      if (verbose_level & 8)
        {
          printf("k = %d\n", k);
          printf("fudge_factor = % 20.14e\n", fudge_factor);
        }

      for (j = 0; j < NUM_RA; ++j)
        {
          azimuth = (REAL8)j/(REAL8)NUM_RA * (REAL8)LAL_TWOPI + fudge_factor;

          for (i = -declim; i <= declim; ++i)
            {
              cnt = j*NUM_DEC + i + declim;
              altitude = asin(((REAL8)i)/((REAL8)declim));

              if (verbose_level & 8)
                {
                  printf("azimuth  = % 14.9e deg\n", rad_to_deg(azimuth));
                  printf("altitude = % 14.9e deg\n", rad_to_deg(altitude));
                }

              resp_plus[cnt] = resp_local(0., LAL_PI_2 - altitude, azimuth,
                                          gwpol_plus);
              resp_cros[cnt] = resp_local(0., LAL_PI_2 - altitude, azimuth,
                                          gwpol_cross);
            }
        }


      if (k == 0)
        {
          skygrid_print("GRASP-computed response to plus",
                        resp_plus, "ff_local_plus.txt");
          skygrid_print("GRASP-computed response to cross",
                        resp_cros, "ff_local_cros.txt");

          skygrid_print("LAL-computed response to plus",
                        plus, "ff_lal_plus.txt");
          skygrid_print("LAL-computed response to cross",
                        cross, "ff_lal_cros.txt");

          skygrid_subtract(tmp3, resp_plus, plus);
          skygrid_fabs(tmp3, tmp3);
          skygrid_print("Abs. difference between GRASP and LAL for plus",
                        tmp3, "ff_diff_plus.txt");

          skygrid_subtract(tmp3, resp_cros, cross);
          skygrid_fabs(tmp3, tmp3);
          skygrid_print("Abs. difference between GRASP and LAL for cross",
                        tmp3, "ff_diff_cros.txt");
        }

      skygrid_square(tmp1, resp_plus);
      skygrid_square(tmp2, resp_cros);
      skygrid_add(tmp1, tmp1, tmp2);

      if (verbose_p >= 4)
        printf("avg(F+^2 + Fx^2) using resp_local() = % 14.9e\n",
               skygrid_avg(tmp1));

      skygrid_subtract(tmp3, resp_plus, plus);
      if (verbose_level & 8)
        printf("RMS difference for plus  = % 20.14e\n", skygrid_rms(tmp3));
      fprintf(rms_diff_plus_file, "% 20.14e    % 20.14e\n",
              fudge_factor, skygrid_rms(tmp3));
      if (verbose_level & 8)
        printf("RMS difference for cross = % 20.14e\n", skygrid_rms(tmp3));
      fprintf(rms_diff_cros_file, "% 20.14e    % 20.14e\n",
              fudge_factor, skygrid_rms(tmp3));
    } /* for (k = -128; ..) */
  if (verbose_level & 4)
    printf("... Done with looping over fudge_factor.\n");

  xfclose(rms_diff_plus_file);
  xfclose(rms_diff_cros_file);

  return;
}



void set_source_params(LALSource * source, const char *name, REAL8 ra_rad,
                       REAL8 dec_rad, REAL8 orien_rad)
{
  (void)laldr_strlcpy(source->name, name, LALNameLength);
  source->equatorialCoords.longitude = ra_rad;
  source->equatorialCoords.latitude  = dec_rad;
  source->equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  source->orientation                = orien_rad;
} /* END: set_source_params() */



void setup_global_sources(void)
{
  /*
   * LALSource   src_0_0_p;
   * LALSource   src_0_0_c;
   * LALSource   src_0_90_p;
   * LALSource   src_0_90_c;
   * LALSource   src_0_45_p;
   * LALSource   src_0_45_c;
   */
  static BOOLEAN global_sources_set_p = FALSE;

  if (global_sources_set_p == FALSE)
    {
      set_source_params(&src_0_0_p, "RA=0deg, Dec=0deg, plus",
                        0., 0., deg_to_rad(-90.));
      set_source_params(&src_0_0_c, "RA=0deg, Dec=0deg, cros",
                        0., 0., deg_to_rad(-45.));

      print_source_maybe(&src_0_0_p);
      print_source_maybe(&src_0_0_c);

      set_source_params(&src_0_90_p, "RA=0deg, Dec=90deg, plus",
                        0., deg_to_rad(90.), deg_to_rad(-90.));
      set_source_params(&src_0_90_c, "RA=0deg, Dec=90deg, cros",
                        0., deg_to_rad(90.), deg_to_rad(-45.));

      print_source_maybe(&src_0_90_p);
      print_source_maybe(&src_0_90_c);

      set_source_params(&src_0_45_p, "RA=0deg, Dec=45deg, plus",
                        0., deg_to_rad(45.), deg_to_rad(-90.));
      set_source_params(&src_0_45_c, "RA=0deg, Dec=45deg, cros",
                        0., deg_to_rad(45.), deg_to_rad(-45.));

      print_source_maybe(&src_0_45_p);
      print_source_maybe(&src_0_45_c);

      global_sources_set_p = TRUE;
    }

  return;
} /* END: setup_global_sources */


void find_zero_gmst(LALStatus * status)
{
  REAL8             gmst1;
  LIGOTimeGPS       gps;
  LALGPSandAcc      gps_and_acc;
  LALMSTUnitsAndAcc tmp_uandacc;
  LALTimeInterval   interval;
  INT4              k;

  gmst1 = 0.;
  gps.gpsSeconds     =  13675020;
  gps.gpsNanoSeconds = 943728500;
  gps_and_acc.gps = gps;
  gps_and_acc.accuracy = LALLEAPSEC_STRICT;
  tmp_uandacc.units = MST_RAD;
  tmp_uandacc.accuracy = gps_and_acc.accuracy;
  interval.seconds = 0;
  interval.nanoSeconds =   1;

  printf("2*Pi = % 22.14e\n", 2. * LAL_PI);

  for (k = 0; k < 4096; ++k)
    {
      /*  to avoid printing out all the LAL INFO messages */
      lalDebugLevel = 0;
      LALGPStoGMST1(status, &gmst1, &(gps_and_acc.gps),
                    &tmp_uandacc);

      printf("k = %9d; GPS = %d:%d;\t\tgmst1 = % 22.14e; gmst1-2*Pi = % 20.14e\n",
             k, gps_and_acc.gps.gpsSeconds, gps_and_acc.gps.gpsNanoSeconds,
             gmst1, (gmst1-2.*(double)LAL_PI));

      /* increment observation time */
      LALIncrementGPS(status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                      &interval);
    }
  printf("AUF WIEDERSEHEN!\n");
}




void print_source_maybe(const LALSource * source)
{
  if (verbose_level & 4)
    {
      printf("Source name: %s\n", source->name);
      printf("RA:          % 14.20e rad\n",
             source->equatorialCoords.longitude);
      printf("Dec:         % 14.20e rad\n",
             source->equatorialCoords.latitude);
      printf("Orientation: % 14.20e rad\n", source->orientation);
    }
}



void crab_pulsar_test(LALStatus * status)
{
  LALSource        crab_pulsar;
  LALDetector      detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetAndSource  det_and_pulsar;
  LALDetAMResponseSeries am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries  plus_series, cross_series, scalar_series;
  LALTimeIntervalAndNSample  time_info;

  printf("MAHALO!\n");
  fflush(stdout);

  if (verbose_level & 4)
    {
      PrintLALDetector(&detector);
    }

  set_source_params(&crab_pulsar, "Crab Pulsar",
                    deg_to_rad((5. + 34./60. + 32./3600.) * 15.),
                    deg_to_rad(22. + 0. + 52./3600.), 0.);
  print_source_maybe(&crab_pulsar);

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &crab_pulsar;

  plus_series.data   = NULL;
  cross_series.data  = NULL;
  scalar_series.data = NULL;

  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(status, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(status, &(am_response_series.pCross->data), 1);
  LALSCreateVector(status, &(am_response_series.pScalar->data), 1);

  time_info.epoch.gpsSeconds     = 709398013;
  time_info.epoch.gpsNanoSeconds =         0;
  time_info.deltaT               =       864;
  time_info.nSample              =       100;
  time_info.accuracy             = LALLEAPSEC_STRICT;

  LALComputeDetAMResponseSeries(status, &am_response_series,
                                &det_and_pulsar,
                                &time_info);

  LALSPrintTimeSeries(am_response_series.pPlus, "crab_plus.txt");
  LALSPrintTimeSeries(am_response_series.pCross, "crab_cross.txt");

  LALSDestroyVector(status, &(am_response_series.pPlus->data));
  LALSDestroyVector(status, &(am_response_series.pCross->data));
  LALSDestroyVector(status, &(am_response_series.pScalar->data));
}


BOOLEAN passed_matrix_test_p(void)
{
  BOOLEAN retval = TRUE;
  LALDR_33Matrix    A;    /* test matrices */
  LALDR_33Matrix    B;
  LALDR_33Matrix    C;
  LALDR_33Matrix    D;
  LALDR_3Vector     u;    /* test vector */
  LALDR_3Vector     v;
  LALDR_3Vector     w;
  LALDR_3Vector     x;
  REAL8             c;    /* test scalar */
  REAL8             d;
  unsigned int      iter, jter;


  if (verbose_p)
  {
    printf("TEST OF MATRIX AND VECTOR FUNCTIONS\n");
    printf("-----------------------------------\n");
  }

  /* Print33Matrix */
  A[0][0] = 0.;    A[0][1] = 0.;   A[0][2] = 0.;
  A[1][0] = 0.;    A[1][1] = 0.;   A[1][2] = 0.;
  A[2][0] = 0.;    A[2][1] = 0.;   A[2][2] = 0.;

  if (verbose_p)
  {
    printf("Print33Matrix output:\n");
    LALDR_Print33Matrix(&A, "A", 0, stdout, "");
    printf("\n");
  }

  if (verbose_p)
  {
    printf("Expected output:\n");
    printf("A:\n");

    for (iter = 0; iter < 3; ++iter)
    {
      for (jter = 0; jter < 3; ++jter)
      {
          printf("% 20.14e\t", A[iter][jter]);
      }
      printf("\n");
    }
  }

  if (verbose_p)
    printf("\n- - - - -\n\n");

  A[0][0] = -4.;    A[0][1] = -3.;   A[0][2] = -2.;
  A[1][0] = -1.;    A[1][1] =  0.;   A[1][2] =  1.;
  A[2][0] =  2.;    A[2][1] =  3.;   A[2][2] =  4.;

  if (verbose_p)
  {
    printf("Print33Matrix output:\n");
    LALDR_Print33Matrix(&A, "A", 0, stdout, "");
    printf("\n");
  }

  if (verbose_p)
  {
    printf("Expected output:\n");
    printf("A:\n");

    for (iter = 0; iter < 3; ++iter)
    {
      for (jter = 0; jter < 3; ++jter)
      {
          printf("% 20.14e\t", A[iter][jter]);
      }
      printf("\n");
    }
  }

  print_separator_maybe();


  /* Zero33Matrix */
  LALDR_Zero33Matrix(&A);
  B[0][0] = 0.;    B[0][1] = 0.;   B[0][2] = 0.;
  B[1][0] = 0.;    B[1][1] = 0.;   B[1][2] = 0.;
  B[2][0] = 0.;    B[2][1] = 0.;   B[2][2] = 0.;

  print_m_results_maybe("Zero33Matrix test", &A, &B);

  if (!matrix_ok_p(&A, &B, zero_tolerance))
    {
      fprintf(stderr, "ERROR: A != B\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  /* Set33Matrix */
  LALDR_Set33Matrix(&A,
                    4., 3., 2.,
                    1., 0., -1.,
                    -2., -3., -4.);
  B[0][0] = 4.;     B[0][1] = 3.;    B[0][2] = 2.;
  B[1][0] = 1.;     B[1][1] = 0.;    B[1][2] = -1.;
  B[2][0] = -2.;    B[2][1] = -3.;   B[2][2] = -4.;

  print_m_results_maybe("Set33Matrix test", &A, &B);

  if (!matrix_ok_p(&A, &B, zero_tolerance))
    {
      fprintf(stderr, "ERROR: A != B\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* Matrix addition */
  LALDR_Add33Matrix(&C, &A, &B);
  LALDR_Set33Matrix(&D,
                    8., 6., 4.,
                    2., 0., -2.,
                    -4., -6., -8.);

  print_m_results_maybe("Add33Matrix test", &C, &D);

  if (!matrix_ok_p(&C, &D, zero_tolerance))
    {
      fprintf(stderr, "ERROR: C != D\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* Matrix multiplication */
  LALDR_Set33Matrix(&A,
                    4., 3., 2.,
                    1., 0., -4.,
                    -3., -2., -1.);

  LALDR_Set33Matrix(&B,
                    8., 17., 5.,
                    7., -23., -9.,
                    -11., 13., 31.);

  LALDR_Set33Matrix(&D,
                    31., 25., 55.,
                    52., -35., -119.,
                    -27., -18., -28.);

  LALDR_Multiply33Matrix(&C, &A, &B);

  print_m_results_maybe("Multiply33Matrix test", &C, &D);

  if (!matrix_ok_p(&C, &D, zero_tolerance))
    {
      fprintf(stderr, "ERROR: C != D\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* Scalar * matrix */
  c = 1.69e-3;
  LALDR_Set33Matrix(&D,
                    1.352e-2, 2.873e-2, 8.45e-3,
                    1.183e-2, -3.887e-2, -1.521e-2,
                    -1.859e-2, 2.197e-2, 5.239e-2);

  LALDR_ScalarMult33Matrix(&C, c, &B);

  print_m_results_maybe("ScalarMult33Matrix test", &C, &D);

  if (!matrix_ok_p(&C, &D, zero_tolerance))
    {
      fprintf(stderr, "ERROR: C != D\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* scalar product of matrices */
  d = 112.;

  c = LALDR_DotProd33Matrix(&A, &B);

  print_s_results_maybe("DotProd33Matrix test", c, d);

  if (!almost_equal_real8_p(c, d, zero_tolerance))
    {
      fprintf(stderr, "ERROR: c != d\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  /* Transpose */
  LALDR_Set33Matrix(&A,
                    4., 3., 2.,
                    1., 0., -1.,
                    -2., -3., -4.);
  LALDR_Transpose33Matrix(&B, &A);

  LALDR_Set33Matrix(&C,
                    4., 1., -2.,
                    3., 0., -3.,
                    2., -1., -4.);

  print_m_results_maybe("Transpose33Matrix test", &B, &C);

  if (!matrix_ok_p(&B, &C, zero_tolerance))
    {
      fprintf(stderr, "ERROR: B != C\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  /* Set3Vector test */
  v[0] = 1.;  v[1] = -2.;  v[2] = 3.;
  LALDR_Set3Vector(&u, 1., -2., 3.);

  print_v_results_maybe("Set3Vector test", &u, &v);

  if (!vector_ok_p(&u, &v, zero_tolerance))
    {
      fprintf(stderr, "ERROR: u != v\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* DotProd3Vector test */
  LALDR_Set3Vector(&v, 13., -11., 51.);
  c = LALDR_DotProd3Vector(&u, &v);

  d = 188.;

  print_s_results_maybe("DotProd3Vector test", c, d);

  if (!almost_equal_real8_p(c, d, zero_tolerance))
    {
      fprintf(stderr, "ERROR: c != d\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  /* CrossProd3Vector test */
  LALDR_Set3Vector(&u, 4., -11., 2.);
  LALDR_Set3Vector(&v, -3., 17., 5.);
  LALDR_CrossProd3Vector(&w, &u, &v);

  LALDR_Set3Vector(&x, -89., -26., 35.);

  print_v_results_maybe("CrossProd3Vector test", &w, &x);

  if (!vector_ok_p(&w, &x, zero_tolerance))
    {
      fprintf(stderr, "ERROR: w != x\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();


  /* OuterProd3Vector test */
  LALDR_OuterProd3Vector(&A, &u, &v);

  LALDR_Set33Matrix(&B,
                    -12., 68., 20.,
                    33., -187., -55.,
                    -6., 34., 10.);

  print_m_results_maybe("OuterProd3Vector test", &A, &B);

  if (!matrix_ok_p(&A, &B, zero_tolerance))
    {
      fprintf(stderr, "ERROR: A != B\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  return retval;
} /* END: pass_matrix_test_p() */



BOOLEAN passed_almost_equal_tests_p(void)
{
  /*
   * TEST -1: Test of almost_equal_real[48]_p() functions
   */

  REAL4 foo4, bar4;
  REAL8 foo8, bar8;

  if (verbose_p)
    {
      printf("TEST OF almost_equal_real[48]_p() functions\n");
      printf("-------------------------------------------\n");
    }

  foo4 = 0.;
  bar4 = real4_tolerance/2.;

  if (!almost_equal_real4_p(foo4, bar4, real4_tolerance))
    {
      fprintf(stderr, "ERROR: almost_equal_real4_p() failed test 1\n");

      return FALSE;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real4_p() test 1\n");
      printf("\n- - - - -\n");
    }


  bar4 = 0.;

  if (!almost_equal_real4_p(foo4, bar4, 0.))
    {
      fprintf(stderr, "ERROR: almost_equal_real4_p() failed test 2\n");

      return FALSE;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real4_p() test 2\n");
      printf("\n- - - - -\n");
    }

  bar4 = 2.*real4_tolerance;

  if (almost_equal_real4_p(foo4, bar4, real4_tolerance))
    {
      fprintf(stderr, "ERROR: almost_equal_real4_p() failed test 3\n");

      return FALSE;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real4_p() test 3\n");
      printf("\n- - - - -\n");
    }

  foo8 = 0.;
  bar8 = real8_tolerance/2.;

  if (!almost_equal_real8_p(foo8, bar8, real8_tolerance))
    {
      fprintf(stderr, "ERROR: almost_equal_real8_p() failed test 1\n");

      return 1;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real8_p() test 1\n");
      printf("\n- - - - -\n");
    }

  bar8 = 0.;

  if (!almost_equal_real8_p(foo8, bar8, 0.))
    {
      fprintf(stderr, "ERROR: almost_equal_real8_p() failed test 2\n");

      return FALSE;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real8_p() test 2\n");
      printf("\n- - - - -\n");
    }

  bar8 = 2.*real8_tolerance;

  if (almost_equal_real8_p(foo8, bar8, real8_tolerance))
    {
      fprintf(stderr, "ERROR: almost_equal_real8_p() failed test 3\n");

      return FALSE;
    }

  if (verbose_p)
    {
      printf("PASS: almost_equal_real8_p() test 3\n");
    }

  print_separator_maybe();

  return TRUE;
}  /* END: passed_almost_equal_tests_p() */



#if 0
static int local_strncasecmp(const char * a, const char * b, size_t maxlen)
{
  char tmp_a[LALNameLength];
  char tmp_b[LALNameLength];
  size_t i;

  strncpy(tmp_a, a, maxlen);
  strncpy(tmp_b, b, maxlen);

  for (i = 0; i < maxlen; ++i)
    {
      tmp_a[i] = tolower(tmp_a[i]);
      tmp_b[i] = tolower(tmp_b[i]);
    }

  return strncmp(tmp_a, tmp_b, maxlen);
}
#endif
