/* <lalVerbatim file="LALTestDetResponse0CV">
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
#include <lal/LALConfig.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>

#include <lal/PrintFTSeries.h>

NRCSID( LALTESTDETRESPONSE0C, "$Id$" );

#define FALSE 0
#define TRUE  1
#define LALDR_MATRIXSIZE 3

typedef REAL8 LALDR_3Vector[3];
typedef REAL8 LALDR_33Matrix[3][3];

int        lalDebugLevel = 0;
BOOLEAN    verbose_p     = FALSE;
const INT4 oneBillion    = 1000000000;

REAL8 default_tolerance = (REAL8)0.;

static void REAL4VectorSubtraction(const REAL4Vector * pA,
                                   const REAL4Vector * pB,
                                   REAL4Vector *pAminusB);
static REAL4 REAL4VectorRMS(REAL4Vector *pVector);
static void PrintLALDetector(LALDetector * const detector);
static void PrintDetResponse(const LALDetAMResponse * const response,
                             const char * const title);

static BOOLEAN matrix_ok_p(LALDR_33Matrix * const computed,
                           LALDR_33Matrix * const expected,
                           REAL8 tolerance);
static BOOLEAN vector_ok_p(LALDR_3Vector * const computed,
                           LALDR_3Vector * const expected,
                           REAL8 tolerance);
static BOOLEAN scalar_ok_p(REAL8 computed, REAL8 expected, REAL8 tolerance);

static void print_m_results_maybe(const char * title,
                                  LALDR_33Matrix * const computed,
                                  LALDR_33Matrix * const expected);

static void print_v_results_maybe(const char * title,
                                  LALDR_3Vector * const computed,
                                  LALDR_3Vector * const expected);

static void print_s_results_maybe(const char * title,
                                  REAL8 computed, REAL8 expected);

static int print_separator_maybe(void);

static int print_passed_maybe(void);


static BOOLEAN detresponse_ok_p(LALStatus * status,
                                const LALDetAndSource * const det_and_src,
                                const LIGOTimeGPS * const gps,
                                const LALDetAMResponse * const expected_resp);

static BOOLEAN frdetector_ok_p(const LALFrDetector * const computed,
                               const LALFrDetector * const expected);

static BOOLEAN detector_ok_p(const LALDetector * const computed,
                             const LALDetector * const expected);

/*
 * Private functions
 */
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
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  /* REAL4Vector               diffVector; */
  LALTimeIntervalAndNSample time_info;

  UINT4 i;
  
  if (argc == 2)
    lalDebugLevel = atoi(argv[1]);

  if (lalDebugLevel)
    verbose_p = TRUE;

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

  if (!matrix_ok_p(&A, &B, default_tolerance))
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
  
  if (!matrix_ok_p(&A, &B, default_tolerance))
    {
      fprintf(stderr, "ERROR: A != B\n");
      
      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  if (verbose_p)
    printf("\n* * * * * * * * * *\n\n");

  
  /* Matrix addition */
  LALDR_Add33Matrix(&C, &A, &B);
  LALDR_Set33Matrix(&D,
                    8., 6., 4.,
                    2., 0., -2.,
                    -4., -6., -8.);

  print_m_results_maybe("Add33Matrix test", &C, &D);
  
  if (!matrix_ok_p(&C, &D, default_tolerance))
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

  if (!matrix_ok_p(&C, &D, default_tolerance))
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

  if (!matrix_ok_p(&C, &D, default_tolerance))
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

  if (!scalar_ok_p(c, d, default_tolerance))
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

  if (!matrix_ok_p(&B, &C, default_tolerance))
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

  if (!vector_ok_p(&u, &v, default_tolerance))
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

  if (!scalar_ok_p(c, d, default_tolerance))
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

  if (!vector_ok_p(&w, &x, default_tolerance))
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

  if (!matrix_ok_p(&A, &B, default_tolerance))
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
     LALFrDetector says that azimuth is measure North from East
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
   * Once again, I wish people would stop using jargon.  This document
   * claims, in a footnote, that the arm azimuths (is this the correct
   * plural?) are in Northing and Easting, BUT Table 4 for LHO and Table
   * 18 for LLO specify azimuths in, variously, Northing and Southing
   * (i.e. deviations towards North and South respectively)
   *
   * Regardless, convert to plain-old bearings, which are measured
   * Clockwise (East) from North, in radians, here are the azimuths for:
   * 
   * 1. LHO -- x-arm azimuth = 305.9994deg = 5.340697rad
   *           y-arm azimuth = 215.9994deg = 3.769901rad
   *    (check: x-arm azi - y-arm azi = 90.0000deg)
   *
   * 2. LLO -- x-arm azimuth = 197.7165deg
   *           y-arm azimuth = 107.7165deg
   *    (check: x-arm azi - y-arm azi = 90.0000deg)
   *
   */

  /* First, pass a LALFrDetector using Frame spec for azimuths */
  strncpy(frdet.name, "LHO, from FrDetector struct (Frame spec)", LALNameLength);
  frdet.vertexLongitudeRadians = (REAL8)deg_to_rad(-119. - 25./60. - 27.5657/3600.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(46. + 27./60. + 18.528/3600.);
  frdet.vertexElevation        = 142.554;
  frdet.xArmAltitudeRadians    = -6.195e-04;
  frdet.yArmAltitudeRadians    =  1.250e-05;
  frdet.xArmAzimuthRadians     = deg_to_rad(305.9994);
  frdet.yArmAzimuthRadians     = deg_to_rad(215.9994);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

  /*
  if (!detector_ok_p(&detector, &(lalCachedDetectors[LALDetectorIndexLHODIFF])))
    {
      if (verbose_p)
        fprintf(stderr, "LHO computed w/ Frame spec != LHO cached\n");
      
      return 1;
    }
  */

  if (verbose_p)
    {
      printf("LHO tensor converted from LALFrDetector (Frame spec):\n");
      PrintLALDetector(&detector);
      printf("\n- - -\n");
    }

  /* Second, pass a LALFrDetector using LAL (package-tools) spec for azi */
  strncpy(frdet.name, "LHO, from FrDetector struct (package-tools spec)",
          LALNameLength);
  frdet.vertexLongitudeRadians  = (REAL8)deg_to_rad(-119. - 25./60. - 27.5657/3600.);
  frdet.vertexLatitudeRadians  = (REAL8)deg_to_rad(46. + 27./60. + 18.528/3600.);
  frdet.vertexElevation        = 142.554;
  frdet.xArmAltitudeRadians    = -6.195e-04;
  frdet.yArmAltitudeRadians    =  1.250e-05;
  frdet.xArmAzimuthRadians     = deg_to_rad( 54.0006);
  frdet.yArmAzimuthRadians     = deg_to_rad(144.0006);
  /* check: y-arm azi - x-arm azi = 90deg */

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);

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

  if (verbose_p)
    {
      printf("TEST OF LALDetAMResponse()\n");
      printf("--------------------------\n\n");
      printf("Tolerance: LAL_REAL4_EPS = % 15.8e\n\n", LAL_REAL4_EPS);
    }

  /* TEST 1: test LALDetAMResponse() at particular points */
  
  /*
   * As per John Whelan's suggestion, directly create a LALDetector 
   * structure rather than starting with a LALFrDetector structure.
   */
  detector.location[0] = LAL_AWGS84_SI;
  detector.location[1] = 0.;
  detector.location[2] = 0.;
  detector.response[0][0] = 0.;
  detector.response[1][1] = 0.5;
  detector.response[2][2] = -0.5;
  detector.response[0][1] = detector.response[1][0] = 0.;
  detector.response[0][2] = detector.response[2][0] = 0.;
  detector.response[1][2] = detector.response[2][1] = 0.;
  detector.type = LALDETECTORTYPE_ABSENT;

  if (status.statusCode && lalDebugLevel)
    {
      fprintf(stderr,
              "LALTestDetResponse0: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /*
   * Set a GPS time that's close to 0h GMST1. (Found this by trial and
   * error.) 
   */
  gps.gpsSeconds     =     61094;
  gps.gpsNanoSeconds = 640000000;

  /*
   * Set up a source at (RA=0, Dec=0, orientation=0, at time GMST1=0)
   */
  strncpy(pulsar.name, "TEST PULSAR 1", LALNameLength);
  pulsar.equatorialCoords.longitude = 0.;  /* RA */
  pulsar.equatorialCoords.latitude  = 0.;  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = 0.;  /* orientation */


  /*
   * Stuff Detector and Source into structure
   */
  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*
   * Compute instantaneous AM response
   */
  LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponse, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }
  
  expected_resp.plus = 1.;
  expected_resp.cross = 0.;
  expected_resp.scalar = 0.;

  if (!detresponse_ok_p(&status, &det_and_pulsar, &gps, &expected_resp))
    {
      fprintf(stderr, "ERROR: computed response != expected response\n");

      return 1;
    }
  else
    {
      print_passed_maybe();
    }

  print_separator_maybe();

  /* another point */

  goto skip; /* skip this for now -- i know it fails */

  /* use LALFrDetector structure for convenience */
  strncpy(frdet.name, "Reference", LALNameLength);
  frdet.vertexLongitudeRadians = deg_to_rad(-71.);
  frdet.vertexLatitudeRadians  = deg_to_rad(85.);
  frdet.xArmAltitudeRadians    = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(45.);
  frdet.yArmAzimuthRadians     = deg_to_rad(135.);

  LALCreateDetector(&status, &detector, &frdet, LALDETECTORTYPE_IFODIFF);  
  
  strncpy(pulsar.name, "TEST PULSAR 2", LALNameLength);
  pulsar.equatorialCoords.longitude = deg_to_rad(2.64721);  /* RA */
  pulsar.equatorialCoords.latitude  = deg_to_rad(27.3022);  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = 0.;  /* orientation */  

  LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponse, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  expected_resp.plus   =  8.7130418971e-02;
  expected_resp.cross  = -2.4005883737e-02;
  expected_resp.scalar =  9.0376945877e-02;

  if (!detresponse_ok_p(&status, &det_and_pulsar, &gps, &expected_resp))
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
  detector.location[0] = LAL_AWGS84_SI;
  detector.location[1] = 0.;
  detector.location[2] = 0.;
  detector.response[0][0] = 0.;
  detector.response[1][1] = 0.5;
  detector.response[2][2] = -0.5;
  detector.response[0][1] = detector.response[1][0] = 0.;
  detector.response[0][2] = detector.response[2][0] = 0.;
  detector.response[1][2] = detector.response[2][1] = 0.;
  detector.type = LALDETECTORTYPE_ABSENT;

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;
  
  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(&status, &(am_response_series.pCross->data), 1);
  LALSCreateVector(&status, &(am_response_series.pScalar->data), 1);

  if (lalDebugLevel > 0)
    {
      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCros->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
    }
    
  time_info.epoch.gpsSeconds     = 61094;
  time_info.epoch.gpsNanoSeconds = 640000000;
  time_info.deltaT               = 1800;
  time_info.nSample              = 25;

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
    }

  if (lalDebugLevel >= 8)
    {
      printf("plus: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%14.6e, ", am_response_series.pPlus->data->data[i]);
        }
      printf(")\n");

      printf("cross: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%14.6e, ", am_response_series.pCross->data->data[i]);
        }
      printf(")\n");

      printf("scalar: (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%14.6e, ", am_response_series.pScalar->data->data[i]);
        }
      printf(")\n");


      /* print out quadrature sum of plus- and cross-response */
      printf("sqrt(PLUS^2 + CROSS^2): (");
      for (i = 0; i < time_info.nSample; ++i)
        {
          printf("%1.6e, ", sqrt(am_response_series.pPlus->data->data[i] *
                                 am_response_series.pPlus->data->data[i] +
                                 am_response_series.pCross->data->data[i] *
                                 am_response_series.pCross->data->data[i]));
        }
      printf(")\n");
    }


  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));

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

  if (tolerance == 0.)
    {
      for (i = 0; i < 2; ++i)
        for (j = 0; j < 2; ++j)
          {
            if ((*computed)[i][j] != (*expected)[i][j])
              {
                if (verbose_p)
                  {
                    LALDR_Print33Matrix(expected, "expected", 0, stdout, "");
                    LALDR_Print33Matrix(computed, "computed", 0, stdout, "");
                  }

                return FALSE;
              }
            else
              {
                return TRUE;
              }
          }
    }
  else
    {
      for (i = 0; i < 2; ++i)
        for (j = 0; j < 2; ++j)
          {
            if (fabs((*computed)[i][j] - (*expected)[i][j]) > tolerance)
              {
                if (verbose_p)
                  {
                    LALDR_Print33Matrix(expected, "expected", 0, stdout, "");
                    LALDR_Print33Matrix(computed, "computed", 0, stdout, "");
                  }

                return FALSE;
              }
            else
              {
                return TRUE;
              }
          }
    }
}



static BOOLEAN vector_ok_p(LALDR_3Vector * const computed,
                           LALDR_3Vector * const expected,
                           REAL8 tolerance)
{
  INT4 i;

  if (tolerance == 0.)
    {
      for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        {
          if ((*computed)[i] != (*expected)[i])
            {
              if (verbose_p)
                {
                  LALDR_Print3Vector(computed, "computed", stdout);
                  LALDR_Print3Vector(expected, "expected", stdout);
                }

              return FALSE;
            }
          else
            {
              return TRUE;
            }
        }
    }
  else
    {
      for (i = 0; i < LALDR_MATRIXSIZE; ++i)
        {
          if (fabs((*computed)[i] - (*expected)[i]) > tolerance)
            {
              if (verbose_p)
                {
                  LALDR_Print3Vector(computed, "computed", stdout);
                  LALDR_Print3Vector(expected, "expected", stdout);
                }

              return FALSE;
            }
          else
            {
              return TRUE;
            }
        }
    }
}



static BOOLEAN scalar_ok_p(REAL8 computed, REAL8 expected, REAL8 tolerance)
{
  if (tolerance == 0.)
    {
      if (computed != expected)
        {
          if (verbose_p)
            {
              printf("computed = %22.14e\n", computed);
              printf("expected = %22.14e\n", expected);
            }

          return FALSE;
        }
      else
        {
          return TRUE;
        }
    }
  else
    {
      if (fabs(computed - expected) > tolerance)
        {
          if (verbose_p)
            {
              printf("computed = %22.14e\n", computed);
              printf("expected = %22.14e\n", expected);
            }

          return FALSE;
        }
      else
        {
          return TRUE;
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



static int print_passed_maybe(void)
{
  if (verbose_p)
    return printf("\nOK\n");
  else
    return 0;
}



static BOOLEAN detresponse_ok_p(LALStatus * status,
                                const LALDetAndSource * const det_and_src,
                                const LIGOTimeGPS * const gps,
                                const LALDetAMResponse * const expected_resp)
                                
{
  LALDetAMResponse computed_resp;
  BOOLEAN resp_plus_ok_p;
  BOOLEAN resp_cross_ok_p;
  BOOLEAN resp_scalar_ok_p;

  LALComputeDetAMResponse(status, &computed_resp, det_and_src, gps);

  if (verbose_p)
    {
      PrintDetResponse(&computed_resp, "Computed");
      PrintDetResponse(expected_resp, "Expected");
    }

  resp_plus_ok_p = (REAL4)fabs((double)(computed_resp.plus - expected_resp->plus)) < LAL_REAL4_EPS;

  resp_cross_ok_p = (REAL4)fabs((double)(computed_resp.cross - expected_resp->cross)) < LAL_REAL4_EPS;

  resp_scalar_ok_p = (REAL4)fabs((double)(computed_resp.scalar - expected_resp->scalar)) < LAL_REAL4_EPS;

  return (resp_plus_ok_p && resp_cross_ok_p && resp_scalar_ok_p);
}




static void PrintDetResponse(const LALDetAMResponse * const response,
                             const char * const title)
{
  printf("%s:\n", title);
  printf("      plus = % 15.8e\n", response->plus);
  printf("     cross = % 15.8e\n", response->cross);
  printf("    scalar = % 15.8e\n", response->scalar);
}



static BOOLEAN frdetector_ok_p(const LALFrDetector * const computed,
                               const LALFrDetector * const expected)
{
  /* let's not bother with comparing the name */

  return (computed->vertexLongitudeRadians == expected->vertexLongitudeRadians
          && computed->vertexLatitudeRadians == expected->vertexLatitudeRadians
          && computed->vertexElevation == expected->vertexElevation
          && computed->xArmAltitudeRadians == expected->xArmAltitudeRadians
          && computed->xArmAzimuthRadians == expected->xArmAzimuthRadians
          && computed->yArmAltitudeRadians == expected->yArmAltitudeRadians
          && computed->yArmAzimuthRadians == expected->yArmAzimuthRadians);
}




static BOOLEAN detector_ok_p(const LALDetector * const computed,
                             const LALDetector * const expected)
{
  return (vector_ok_p(&(computed->location), &(expected->location), default_tolerance) &&
          matrix_ok_p(&(computed->response), &(expected->response), default_tolerance) &&
          computed->type == expected->type &&
          frdetector_ok_p(&(computed->frDetector), &(expected->frDetector)));
}
