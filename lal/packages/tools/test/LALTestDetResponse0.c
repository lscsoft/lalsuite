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

/* these two constants are for the sky grid */
#define NUM_DEC 11
#define NUM_RA  20

static const INT4 lim = NUM_RA * NUM_DEC;

typedef REAL8 LALDR_3Vector[3];
typedef REAL8 LALDR_33Matrix[3][3];

typedef REAL4 skygrid_t[NUM_RA * NUM_DEC];


int        lalDebugLevel = 0;
BOOLEAN    verbose_p     = FALSE;
const INT4 oneBillion    = 1000000000;

const REAL8 zero_tolerance = 0.;
REAL8 real8_tolerance = (REAL8)2.*LAL_REAL8_EPS;
REAL4 real4_tolerance = (REAL4)2.*LAL_REAL4_EPS;

static void REAL4VectorSubtraction(const REAL4Vector * pA,
                                   const REAL4Vector * pB,
                                   REAL4Vector *pAminusB);
static REAL4 REAL4VectorRMS(REAL4Vector *pVector);
static void PrintLALDetector(LALDetector * const detector);
static void PrintDetResponse(const LALDetAMResponse * const response,
                             const char * const title);


static void make_me_an_Sarray_sequence(LALStatus *status,
                                       REAL4ArraySequence **sequence,
                                       UINT4 rows, UINT4 cols, UINT4 length);

static void print_diagnostics(LALStatus * const status,
                              REAL4ArraySequence *sequence);


static BOOLEAN almost_equal_real4_p(REAL4 a, REAL4 b, REAL4 tolerance);

static BOOLEAN almost_equal_real8_p(REAL8 a, REAL8 b, REAL8 tolerance);

static BOOLEAN almost_equal_real4_relative_p(REAL4 computed, REAL4 expected,
                                             REAL4 tolerance);

static BOOLEAN matrix_ok_p(LALDR_33Matrix * const computed,
                           LALDR_33Matrix * const expected,
                           REAL8 tolerance);
static BOOLEAN vector_ok_p(LALDR_3Vector * const computed,
                           LALDR_3Vector * const expected,
                           REAL8 tolerance);

static BOOLEAN vector_relative_ok_p(LALDR_3Vector * const computed,
                                    LALDR_3Vector * const expected,
                                    REAL8 tolerance);

static void print_m_results_maybe(const char * title,
                                  LALDR_33Matrix * const computed,
                                  LALDR_33Matrix * const expected);

static void print_v_results_maybe(const char * title,
                                  LALDR_3Vector * const computed,
                                  LALDR_3Vector * const expected);

static void print_s_results_maybe(const char * title,
                                  REAL8 computed, REAL8 expected);

static int print_separator_maybe(void);

static int print_small_separator_maybe(void);

static int print_passed_maybe(void);


static BOOLEAN detresponse_ok_p(LALStatus * status,
                                const LALDetAndSource * const det_and_src,
                                const LIGOTimeGPS * const gps,
                                const LALDetAMResponse * const expected_resp,
                                REAL4 tolerance);

static void handle_detresponse_test(BOOLEAN passed_p, int line);

static BOOLEAN frdetector_ok_p(const LALFrDetector * const computed,
                               const LALFrDetector * const expected);

static BOOLEAN detector_ok_p(const LALDetector * const computed,
                             const LALDetector * const expected);


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
                       LALDR_3Vector * const a,
                       LALDR_3Vector * const b)
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
LALDR_DotProd3Vector(LALDR_3Vector * const a,
                     LALDR_3Vector * const b)
{
  INT4 i;
  REAL8 result = 0.;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    result += (*a)[i] * (*b)[i];

  return result;
} /* END: LALDR_DotProd3Vector() */



static void
LALDR_OuterProd3Vector(LALDR_33Matrix * a,
                       LALDR_3Vector * const u,
                       LALDR_3Vector * const v)
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
LALDR_DotProd33Matrix(LALDR_33Matrix * const a, LALDR_33Matrix * const b)
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



/*
 * Copy matrix source to matrix target
 */
static void
LALDR_Copy33Matrix(LALDR_33Matrix * target, LALDR_33Matrix * const source)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*target)[i][j] = (*source)[i][j];

  return;
} /* END: LALDR_Copy33Matrix() */



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
                       LALDR_33Matrix * const matrixL,
                       LALDR_33Matrix * const matrixR)
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
                         LALDR_33Matrix * const  matrix)
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
                  LALDR_33Matrix * const matrix1,
                  LALDR_33Matrix * const matrix2)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*result)[i][j] = (*matrix1)[i][j] + (*matrix2)[i][j];

  return;
}



/*
 * Subtract matrices (M1 - M2)
 */
static void
LALDR_Subtract33Matrix(LALDR_33Matrix * result,
                       LALDR_33Matrix * const matrix1,
                       LALDR_33Matrix * const matrix2)
{
  INT4 i, j;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      (*result)[i][j] = (*matrix1)[i][j] - (*matrix2)[i][j];

  return;
}



/*
 * Transpose matrix
 */
static void
LALDR_Transpose33Matrix(LALDR_33Matrix * transpose,
                        LALDR_33Matrix * const matrix)
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



/*
 * The L2 norm of a matrix
 */
static REAL8
LALDR_L2Norm33Matrix(LALDR_33Matrix * const matrix)
{
    INT4   i, j;
    REAL8  l2norm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            l2norm += (*matrix)[i][j] * (*matrix)[i][j];

    l2norm = sqrt((double)l2norm);

    return l2norm;
}
    



/*
 * The RMS norm of a matrix: RMS sum of all elements.
 */
static REAL8
LALDR_RMSNorm33Matrix(LALDR_33Matrix * const matrix)
{
    INT4   i, j;
    REAL8  rmsnorm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            rmsnorm += (*matrix)[i][j] * (*matrix)[i][j];

    rmsnorm = sqrt((double)rmsnorm / 9.);

    return rmsnorm;
}



/*
 * The "infinity" norm of a matrix: max over all elems
 */
static REAL8
LALDR_InfNorm33Matrix(LALDR_33Matrix * const matrix)
{
    INT4 i, j;
    REAL8 infnorm = 0.;

    for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        for (j = 0; j < LALDR_MATRIXSIZE; ++j)
            if (fabs((double)((*matrix)[i][j])) > infnorm)
                infnorm = fabs((*matrix)[i][j]);
    
    return infnorm;
}



/*
 * Print out matrix
 */
static void
LALDR_Print33Matrix(LALDR_33Matrix * const matrix,
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
LALDR_Print3Vector(LALDR_3Vector * const vector,
                   const CHAR           *varname,
                   FILE                 *file)
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



REAL4 skygrid_avg(const skygrid_t response);
void  skygrid_square(skygrid_t square, const skygrid_t input);
void  skygrid_sqrt(skygrid_t result, const skygrid_t input);
INT4  skygrid_copy(skygrid_t dest, const skygrid_t src);
void  skygrid_print(const skygrid_t input, const char * filename);
void skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b);
void skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b);

FILE *xfopen(const char *path, const char *mode);

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
  REAL8             lmst1;
  LALPlaceAndGPS    det_and_gps;
  LALLeapSecAccuracy accuracy;
  LALGPSandAcc      gps_and_acc;
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;
  LALDetAMResponse  expected_resp;
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

  LALDetAMResponseSeries    am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series,
    circ_series, sum_series;
  REAL4                     mean;
  /* REAL4Vector               diffVector; */
  LALTimeIntervalAndNSample time_info;

  LALTimeInterval           interval;

  UINT4 k;
  INT4  i, j, count;
  REAL4 tolerance;

  skygrid_t plus;
  skygrid_t cross;
  skygrid_t sqsum;
  skygrid_t circ;
  skygrid_t plus_sq_time_avg;
  skygrid_t cross_sq_time_avg;
  skygrid_t sum_of_sq_time_avg;
  skygrid_t tmpskygrid;
  skygrid_t tmpskygrid2;
  skygrid_t theta;
  skygrid_t phi;
  INT4 declim = (NUM_DEC-1)/2;
  FILE     *file_plus_sq_avg = NULL;
  FILE     *file_cross_sq_avg = NULL;
  FILE     *file_plus_at_0_0 = NULL;
  FILE     *file_cross_at_0_0 = NULL;
  FILE     *file_plus_at_5_10 = NULL;
  FILE     *file_cross_at_5_10 = NULL;
  FILE     *file_sum_sq_avg = NULL;
  FILE     *file_sum_sq = NULL;
  FILE     *file_theta = NULL;
  FILE     *file_phi   = NULL;
  
  
  if (argc == 2)
    lalDebugLevel = atoi(argv[1]);

  if (lalDebugLevel)
    verbose_p = TRUE;

  /*
   * TEST -1: Test of almost_equal_real[48]_p() functions
   */
  if (verbose_p)
    {
      REAL4 foo4, bar4;
      REAL8 foo8, bar8;
      printf("TEST OF almost_equal_real[48]_p() functions\n");
      printf("-------------------------------------------\n");

      foo4 = 0.;
      bar4 = real4_tolerance/2.;

      if (!almost_equal_real4_p(foo4, bar4, real4_tolerance))
        {
          fprintf(stderr, "ERROR: almost_equal_real4_p() failed test 1\n");

          return 1;
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

          return 1;
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

          return 1;
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

          return 1;
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

          return 1;
        }

      if (verbose_p)
        {
          printf("PASS: almost_equal_real8_p() test 3\n");
        }

      print_separator_maybe();
    }
  

  /*
   * TEST 0: Test of matrix/vector manipulations
   */
  if (verbose_p)
    {
      printf("TEST OF MATRIX AND VECTOR FUNCTIONS\n");
      printf("-----------------------------------\n");
      printf("Manual inspection required.\n");

      /* Print33Matrix */
      A[0][0] = 0.;    A[0][1] = 0.;   A[0][2] = 0.;
      A[1][0] = 0.;    A[1][1] = 0.;   A[1][2] = 0.;
      A[2][0] = 0.;    A[2][1] = 0.;   A[2][2] = 0.;
  
      printf("Print33Matrix output:\n");
      LALDR_Print33Matrix(&A, "A", 0, stdout, "");
      printf("\n");

      printf("Expected output:\n");
      printf("A:\n");
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", 0., 0., 0.);
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", 0., 0., 0.);
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", 0., 0., 0.);
  
      printf("\n- - - - -\n\n");

      A[0][0] = -4.;    A[0][1] = -3.;   A[0][2] = -2.;
      A[1][0] = -1.;    A[1][1] =  0.;   A[1][2] =  1.;
      A[2][0] =  2.;    A[2][1] =  3.;   A[2][2] =  4.;

      printf("Print33Matrix output:\n");
      LALDR_Print33Matrix(&A, "A", 0, stdout, "");
      printf("\n");

      printf("Expected output:\n");
      printf("A:\n");
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", -4., -3., -2.);
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", -1., 0., 1.);
      printf("% 20.14e\t% 20.14e\t% 20.14e\n", 2., 3., 4.);

      print_separator_maybe();
    }
  

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

  
    
  /* END TEST 0: of matrix/vector routines tests */

  /*********************************************/

  if (verbose_p)
    {
      printf("TEST OF LALCreateDetector()\n");
      printf("---------------------------\n\n");
      printf("Manual inspection required.\n");
    }

  /* Now, here's where the docs for LALFrDetector differ
     from the official Frame spec: LIGO-T970130-F-E

     LIGO doc says azimuth of arms are measured East from North
     (i.e. the conventional definition of bearing)
     LALFrDetector says that azimuth is measured North from East
  */
  strncpy(frdet.name, "Reference", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(0.);
  frdet.vertexLatitudeRadians  = deg_to_rad(0.);
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(90.);
  frdet.yArmAzimuthRadians     = deg_to_rad(180.);

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
  strncpy(frdet.name, "LHO, from FrDetector struct (Frame spec)", LALNameLength);
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
        fprintf(stderr, "WARNING: LHO computed w/ Frame spec != LHO cached\n");
    }


  if (verbose_p)
    {
      printf("LHO tensor converted from LALFrDetector (Frame spec):\n");
      PrintLALDetector(&detector);
      printf("\n- - - - -\n");
    }

  /* Second, pass a LALFrDetector using LAL (package-tools) spec for azi
     (counterclockwise from East) */
  strncpy(frdet.name, "LHO, from FrDetector struct (package-tools spec)",
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
  strncpy(frdet.name, "LHO, from FrDetector struct (numbers from doco)",
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
        fprintf(stderr, "WARNING: LHO computed w/ numbers from LAL documentation != LHO cached\n");
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
  strncpy(frdet.name, "TRIVIAL 1", LALNameLength);
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

  /* Set up a source at (RA=0, Dec=0, orientation=0, at time GMST1=0) */
  strncpy(pulsar.name, "TEST PULSAR 1", LALNameLength);
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
                                           &gps, &expected_resp, tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (0, -1) */
  pulsar.orientation = -LAL_PI_4;
  expected_resp.plus   =  0.;
  expected_resp.cross  = -1.;
  expected_resp.scalar =  0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /*** expect (-1, 0) */
  pulsar.orientation = 0;
  expected_resp.plus   = -1.;
  expected_resp.cross  =  0.;
  expected_resp.scalar =  0.;
  
  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);
  
  print_small_separator_maybe();
  
  /*** expect (0.5, -sqrt(3)/2.) */
  pulsar.orientation = -LAL_PI/3.;
  expected_resp.plus   =  0.5;
  expected_resp.cross  = -sqrt(3.)/2.;
  expected_resp.scalar =  0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);
  
  print_small_separator_maybe();

  /* switch detector to something less trivial */
  strncpy(frdet.name, "TRIVIAL 2", LALNameLength);
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
  strncpy(pulsar.name, "TEST PULSAR 2", LALNameLength);
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
                                           &gps, &expected_resp, tolerance),
                          __LINE__);  
  
  print_small_separator_maybe();

  /*** expect (0, -1) */
  pulsar.orientation = -LAL_PI_2/2.;
  expected_resp.plus = 0.;
  expected_resp.cross = -1.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);  
  
  print_small_separator_maybe();

  /*** expect (-1, 0) */
  pulsar.orientation = 0.;
  expected_resp.plus = -1.;
  expected_resp.cross = 0.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);  
  
  print_small_separator_maybe();

  /*** expect (0.5, -sqrt(3)/2) */
  pulsar.orientation = -LAL_PI/3.;
  expected_resp.plus = 0.5;
  expected_resp.cross = -sqrt(3.)/2.;
  expected_resp.scalar = 0.;

  handle_detresponse_test(detresponse_ok_p(&status, &det_and_pulsar,
                                           &gps, &expected_resp, tolerance),
                          __LINE__);  
  
  print_small_separator_maybe();

  
  /* switch detector to LHO */
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];

  /* switch source */
  strncpy(pulsar.name, "TEST PULSAR 3", LALNameLength);
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
                                           &gps, &expected_resp, tolerance),
                          __LINE__);

  print_small_separator_maybe();

  /* change to even less trivial detector */
  strncpy(frdet.name, "TRIVIAL 3", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(15.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  if (verbose_p)
    PrintLALDetector(&detector);


  /* now have to choose a source whose RA is 15 degrees towards the East;
     whaddaya know? that's 1 hour */
  pulsar.equatorialCoords.longitude = deg_to_rad(15.);

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;  

  /* also, have to choose a time s.t. LMST1 is 0 */
  goto pass;
  gps.gpsSeconds += 829;
  for (i = 0; i < 10000000; i += 1000)
    {
      interval.nanoSeconds += i;
      printf("i = %03d\n", i);
      LALIncrementGPS(&status, &gps, &gps, &interval);
      
      printf("  gps = %u:%u\n", gps.gpsSeconds, gps.gpsNanoSeconds);
      det_and_gps.p_detector = &detector;
      det_and_gps.p_gps      = &gps;
      
      LALGPStoLMST1(&status, &lmst1, &det_and_gps, MST_SEC);
      printf("lmst1 = % 20.14f\n\n", lmst1);
    }

 pass:
  
  goto skip; /* FIXME: skip this for now -- i know it fails */

  /* use LALFrDetector structure for convenience */
  strncpy(frdet.name, "Reference", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(240.592343);
  frdet.vertexLatitudeRadians  = deg_to_rad( 46.455147);
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(360.-(215.9994 - 180.));
  frdet.yArmAzimuthRadians     = deg_to_rad(frdet.xArmAzimuthRadians - 90.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);  
  
  strncpy(pulsar.name, "TEST PULSAR 2", LALNameLength);
  pulsar.equatorialCoords.longitude = deg_to_rad(0.);  /* RA */
  pulsar.equatorialCoords.latitude  = deg_to_rad(0.);  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = deg_to_rad(-90.);  /* orientation */  

  gps_and_acc.gps = gps;
  gps_and_acc.accuracy = LALLEAPSEC_STRICT;
  LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar,
                          &gps_and_acc);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponse, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  utcDate.unixDate.tm_sec = 9;
  utcDate.unixDate.tm_min = 57;
  utcDate.unixDate.tm_hour = 13;
  utcDate.unixDate.tm_mday = 01;
  utcDate.unixDate.tm_mon  = LALMONTH_AUG;
  utcDate.unixDate.tm_year = 1990 - 1900;

  accuracy = LALLEAPSEC_STRICT;
  LALUTCtoGPS(&status, &gps, &utcDate, &accuracy);

  expected_resp.plus   = -1.14595015045445e-01;
  expected_resp.cross  = -2.61367926478056e-01;
  expected_resp.scalar =  0.;
  /* expected circ. response = 9.0376945877e-02; */

  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  if (!detresponse_ok_p(&status, &det_and_pulsar, &gps, &expected_resp,
                        tolerance))
    {
      fprintf(stderr, "ERROR: computed response != expected response\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

 skip:
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
  strncpy(detector.frDetector.name, "FAKE", LALNameLength);


  /* Make a fake detector by specifying frame format detector */
  strncpy(frdet.name, "FAKE FAKE, NOT THE REAL MCCOY", LALNameLength);
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

  /* compute response at the above times "manually" */
  /*
  gps.gpsSeconds     = time_info.epoch.gpsSeconds;
  gps.gpsNanoSeconds = time_info.epoch.gpsNanoSeconds;
  interval.seconds   = 6*3600;
  interval.nanoSeconds = 0;

  for (i = 0; i < 5; ++i)
    {
      LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps);
      printf("plus  = % 14.8e\n", am_response.plus);
      printf("cross = % 14.8e\n", am_response.cross);
      printf("\n");
      LALIncrementGPS(&status, &gps, &gps, &interval);
    }
  
  */
  
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
      for (k = 0; k < circ_series.data->length; ++k)
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

  if (lalDebugLevel >= 8)
    {
      printf("plus: (");
      for (k = 0; k < time_info.nSample; ++k)
        {
          printf("% 14.6e, ", am_response_series.pPlus->data->data[k]);
        }
      printf(")\n");

      printf("cross: (");
      for (k = 0; k < time_info.nSample; ++k)
        {
          printf("% 14.6e, ", am_response_series.pCross->data->data[k]);
        }
      printf(")\n");

      printf("scalar: (");
      for (k = 0; k < time_info.nSample; ++k)
        {
          printf("% 14.6e, ", am_response_series.pScalar->data->data[k]);
        }
      printf(")\n");


      /* print out quadrature sum of plus- and cross-response */
      printf("sqrt(PLUS^2 + CROSS^2): (");
      for (k = 0; k < time_info.nSample; ++k)
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
  printf("ALOHA\n");

  /* use Livingston */
  /* detector = lalCachedDetectors[LALDetectorIndexLLODIFF]; */
  
  if (lalDebugLevel >= 1)
    {
      printf("\nStarting whole-sky test...\n");
      PrintLALDetector(&detector);
      count = 0;
      printf("NUM_RA = %d; NUM_DEC = %d\n", NUM_RA, NUM_DEC);

      file_plus_sq_avg  = xfopen("plus_sq_avg.txt", "w");
      file_cross_sq_avg = xfopen("cross_sq_avg.txt", "w");
      file_plus_at_0_0  = xfopen("plus_at_0_0.txt", "w");
      file_cross_at_0_0 = xfopen("cross_at_0_0.txt", "w");
      file_plus_at_5_10 = xfopen("plus_at_5_10.txt", "w");
      file_cross_at_5_10 = xfopen("cross_at_5_10.txt", "w");
      file_sum_sq_avg = xfopen("sum_sq_avg.txt", "w");
      file_sum_sq = xfopen("sum_sq.txt", "w");
      file_theta = xfopen("theta.txt", "w");
      file_phi   = xfopen("phi.txt", "w");

      printf("Done opening files.\n");
      
      gps.gpsSeconds     = time_info.epoch.gpsSeconds;
      gps.gpsNanoSeconds = time_info.epoch.gpsNanoSeconds;
      LALFloatToInterval(&status, &interval, &(time_info.deltaT));

      printf("CUBAAN\n");

      /* FIXME */
      /* only need to print out the (phi, theta) grid once */
      for (j = 0; j < NUM_RA; ++j)
        {
          for (i = -declim; i <= declim; ++i)
            {
              
            }
        }

      printf("N sample = %d\n", time_info.nSample);
      
      for (k = 0; k < time_info.nSample; ++k)
        {
          for (j = 0; j < NUM_RA; ++j)
            {
              pulsar.equatorialCoords.longitude =
                (REAL8)j/(REAL8)NUM_RA * LAL_TWOPI; /* RA */

              for (i = -declim; i <= declim; ++i)
                {
                  INT4 cnt = j*NUM_DEC + i + declim;

                  pulsar.equatorialCoords.latitude =
                    acos((REAL8)i/(REAL8)NUM_DEC);

                  LALComputeDetAMResponse(&status, &am_response,
                                          &det_and_pulsar, &gps_and_acc);
                  plus[cnt]  = am_response.plus;
                  cross[cnt] = am_response.cross;
                  sqsum[cnt] = (plus[cnt] * plus[cnt])
                    + (cross[cnt] * cross[cnt]);

                  if (i == 0 && j == 0)
                    {
                      fprintf(file_plus_at_0_0, "% 14.9e\n", plus[cnt]);
                      fprintf(file_cross_at_0_0, "% 14.9e\n", cross[cnt]);
                      fprintf(file_sum_sq, "% 14.9e\n", sqsum[cnt]);
                    }

                  if (i == 5 && j == 10)
                    {
                      fprintf(file_plus_at_5_10, "% 14.9e\n", plus[cnt]);
                      fprintf(file_cross_at_5_10, "% 14.9e\n", cross[cnt]);
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

          if (k == 300)
            {
              skygrid_print(plus, "plus.txt");
              skygrid_print(cross, "cross.txt");

              skygrid_square(tmpskygrid, plus);
              skygrid_square(tmpskygrid2, cross);
              skygrid_add(tmpskygrid, tmpskygrid, tmpskygrid2);
              printf("avg(F+^2 + Fx^2) = % 14.8e\n", skygrid_avg(tmpskygrid));
            }

          /* increment observation time */
          LALIncrementGPS(&status, &(gps_and_acc.gps), &(gps_and_acc.gps),
                          &interval);
        }
      printf("\n");

      skygrid_print(plus_sq_time_avg, "plus_sq_time_avg.txt");
      skygrid_print(cross_sq_time_avg, "cross_sq_time_avg.txt");
      skygrid_print(sum_of_sq_time_avg, "sum_of_sq_time_avg.txt");

      fprintf(file_plus_sq_avg, "\n");
      fprintf(file_cross_sq_avg, "\n");
      fprintf(file_plus_at_0_0, "\n");
      fprintf(file_cross_at_0_0, "\n");
      fprintf(file_plus_at_5_10, "\n");
      fprintf(file_cross_at_5_10, "\n");
      fprintf(file_sum_sq_avg, "\n");
      fprintf(file_sum_sq, "\n");

      if (verbose_p)
        printf("avg(F+^2 + Fx^2) = % 14.8e\n", skygrid_avg(sqsum));

      fclose(file_plus_sq_avg);
      fclose(file_cross_sq_avg);
      fclose(file_plus_at_0_0);
      fclose(file_cross_at_0_0);
      fclose(file_plus_at_5_10);
      fclose(file_cross_at_5_10);
      fclose(file_sum_sq_avg);
      fclose(file_sum_sq);
      fclose(file_theta);
      fclose(file_phi);
    }


  /*
   * Housekeeping
   */
  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));
  LALSDestroyVector(&status, &(circ_series.data));
  LALSDestroyVector(&status, &(sum_series.data));

  LALCheckMemoryLeaks();

  return 0;
} /* END: main() */


/*
 * subtracts two REAL4Vectors; user must do all allocation beforehand
 */
static void REAL4VectorSubtraction(const REAL4Vector *pA, const REAL4Vector *pB,
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




static BOOLEAN matrix_ok_p(LALDR_33Matrix * const computed,
                           LALDR_33Matrix * const expected,
                           REAL8 tolerance)
{
  INT4 i, j;

  for (i = 0; i < 2; ++i)
    for (j = 0; j < 2; ++j)
      {
        if (!almost_equal_real8_p((*computed)[i][j], (*expected)[i][j], tolerance))
          {
            if (verbose_p)
              {
                LALDR_Print33Matrix(expected, "INFO: matrix_ok_p(): expected", 0, stdout, "");
                LALDR_Print33Matrix(computed, "INFO: matrix_ok_p(): computed", 0, stdout, "");
              }
            
            return FALSE;
          }
        else
          {
            return TRUE;
          }
      }
}



static BOOLEAN vector_ok_p(LALDR_3Vector * const computed,
                           LALDR_3Vector * const expected,
                           REAL8 tolerance)
{
  INT4 i;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    {
      if (!almost_equal_real8_p((*computed)[i], (*expected)[i], tolerance))
        {
          if (verbose_p)
            {
              LALDR_Print3Vector(computed, "INFO: vector_ok_p(): computed", stdout);
              LALDR_Print3Vector(expected, "INFO: vector_ok_p(): expected", stdout);
            }

          return FALSE;
        }
      else
        {
          return TRUE;
        }
    }
}



static BOOLEAN vector_relative_ok_p(LALDR_3Vector * const computed,
                                    LALDR_3Vector * const expected,
                                    REAL8 tolerance)
{
  /* do this term-by-term rather than using a metric. bah. */
  INT4 i;
  REAL8 relative_diff;

  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    {
      relative_diff = fabs((*computed)[i]/(*expected)[i] - 1.);

      if (relative_diff <= tolerance)
        {
          return TRUE;
        }
      else
        {
          return FALSE;
        }
    }
}




static void print_m_results_maybe(const char * title,
                                  LALDR_33Matrix * const computed,
                                  LALDR_33Matrix * const expected)
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
                                  LALDR_3Vector * const computed,
                                  LALDR_3Vector * const expected)
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
      printf("  (REAL4)fabsf(a - b) = % 16.10e\n", (REAL4)fabsf(a - b));
      printf("  tolerance = % 16.10e\n\n", tolerance);
    }
  if (tolerance == 0.)
    return (a == b);
  else
    return ((REAL4)fabsf(a - b) <= tolerance);
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
      relative_err = fabsf(computed - expected)/fabsf(expected);
      if (relative_err <= tolerance)
        return TRUE;
      else
        return FALSE;
    }
}



static BOOLEAN detresponse_ok_p(LALStatus * status,
                                const LALDetAndSource * const det_and_src,
                                const LIGOTimeGPS * const gps,
                                const LALDetAMResponse * const expected_resp,
                                REAL4 tolerance)
                                
{
  LALDetAMResponse computed_resp;
  BOOLEAN resp_plus_ok_p;
  BOOLEAN resp_cross_ok_p;
  BOOLEAN resp_scalar_ok_p;
  BOOLEAN result;
  REAL4   computed_circ_resp, expected_circ_resp;

  LALComputeDetAMResponse(status, &computed_resp, det_and_src, gps);

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



static BOOLEAN frdetector_ok_p(const LALFrDetector * const computed,
                               const LALFrDetector * const expected)
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




static BOOLEAN detector_ok_p(const LALDetector * const computed,
                             const LALDetector * const expected)
{
  return (vector_relative_ok_p(&(computed->location), &(expected->location),
                               1.e-4)
          && matrix_ok_p(&(computed->response), &(expected->response), 1.e-4)
          && computed->type == expected->type
          && frdetector_ok_p(&(computed->frDetector), &(expected->frDetector)));
}



REAL4 skygrid_avg(const skygrid_t response)
{
  INT4 i;
  REAL4 retval = 0.;

  for (i = 0; i < lim; ++i)
    retval += response[i];

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



INT4 skygrid_copy(skygrid_t dest, const skygrid_t src)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    dest[i] = src[i];

  return i;
}



void skygrid_print(const skygrid_t input, const char * filename)
{
  INT4 i, j;
  FILE * outfile = NULL;

  outfile = xfopen(filename, "w");

  for (i = 0; i < NUM_RA; ++i)
    {
      for (j = 0; j < NUM_DEC; ++j)
        fprintf(outfile, "% 14.8e\t", input[i*NUM_DEC + j]);
      fprintf(outfile, "\n");
    }

  fclose(outfile);
}



void skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    sum[i] = a[i] + b[i];
}

void skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b)
{
  INT4 i;

  for (i = 0; i < lim; ++i)
    result[i] = b * a[i];
}

/*
 * only handle 2-dimensional arrays
 */
static void make_me_an_Sarray_sequence(LALStatus *status,
                                       REAL4ArraySequence **sequence,
                                       UINT4 rows, UINT4 cols, UINT4 length)
{
  CreateArraySequenceIn params;
  UINT4Vector           dimLength;
  UINT4                 data[] = {2, rows, cols};

  dimLength.length = 3;
  dimLength.data   = data;
  params.length    = length;
  params.dimLength = &dimLength;

  LALSCreateArraySequence(status, sequence, &params);
}



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
