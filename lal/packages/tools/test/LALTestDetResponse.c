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

/* <lalVerbatim file="LALTestDetResponseCV">
   $Id$
   </lalVerbatim
*/

/*
<lalLaTeX>

\subsection{Program {\texttt{LALTestDetResponse.c}}}
\label{ss:LALTestDetResponse.c}

\subsubsection*{Usage}

\subsubsection*{Description}

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

NRCSID( LALTESTDETRESPONSEC, "$Id$" );


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

/*
 * Cross product of two 3-vectors:
 *    result = a x b
 */
static void
LALDR_CrossProd3Vector(REAL8 result[3],
                       const REAL8 a[3], const REAL8 b[3])
{
  result[0] =  a[1]*b[2] - a[2]*b[1];
  result[1] = -a[0]*b[2] + a[2]*b[0];
  result[2] =  a[0]*b[1] - a[1]*b[0];

  return;
} /* END: LALDR_CrossProd3Vector() */



/*
 * Dot product of two 3-vectors
 */
static REAL8
LALDR_DotProd3Vector(REAL8 a[3], REAL8 b[3])
{
  INT4 i;
  REAL8 result = 0.;

  for (i = 0; i < 3; ++i)
    result += a[i] * b[i];

  return result;
} /* END: LALDR_DotProd3Vector() */



/*
 * Scalar product of two 3x3 matrices
 */
static REAL8
LALDR_DotProd33Matrix(REAL8 a[3][3], REAL8 b[3][3])
{
  INT4 i, j;
  REAL8 result = 0.;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      result += a[i][j] * b[i][j];

  return result;
} /* END: LALDR_DotProd33Matrix() */




/*
 * Sets all elements of a 3x3 matrix
 */
static void
LALDR_Set33Matrix(REAL8 matrix[3][3],
                  REAL8 a11, REAL8 a12, REAL8 a13,
                  REAL8 a21, REAL8 a22, REAL8 a23,
                  REAL8 a31, REAL8 a32, REAL8 a33)
{
  matrix[0][0] = a11;  matrix[0][1] = a12;  matrix[0][2] = a13;

  matrix[1][0] = a21;  matrix[1][1] = a22;  matrix[1][2] = a23;

  matrix[2][0] = a31;  matrix[2][1] = a32;  matrix[2][2] = a33;

  return;
} /* END: LALDR_Set33Matrix() */



/*
 * Copy matrix source to matrix target
 */
static void
LALDR_Copy33Matrix(REAL8 target[3][3], const REAL8 source[3][3])
{
  INT4 i, j;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      target[i][j] = source[i][j];

  return;
} /* END: LALDR_Copy33Matrix() */



/*
 * Zero matrix
 */
static void
LALDR_Zero33Matrix(REAL8 matrix[3][3])
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
LALDR_Multiply33Matrix(REAL8 product[3][3],
                       REAL8 matrixL[3][3],
                       REAL8 matrixR[3][3])
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
  for (i = 0; i < 3; ++i)
    for (k = 0; k < 3; ++k)
      for (j = 0; j < 3; ++j)
        product[i][k] += matrixL[i][j] * matrixR[j][k];

  return;
}



/*
 * Scalar multiply
 */
static void
LALDR_ScalarMult33Matrix(REAL8 result[3][3],
                         REAL8 coefficient,
                         const REAL8 matrix[3][3])
{
  INT4 i, j;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      result[i][j] = coefficient * matrix[i][j];

  return;
}



/*
 * Add matrix
 */
static void
LALDR_Add33Matrix(REAL8 result[3][3],
                  const REAL8 matrix1[3][3],
                  const REAL8 matrix2[3][3])
{
  INT4 i, j;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      result[i][j] = matrix1[i][j] + matrix2[i][j];

  return;
}



/*
 * Subtract matrices (M1 - M2)
 */
static void
LALDR_Subtract33Matrix(REAL8 result[3][3],
                  const REAL8 matrix1[3][3],
                  const REAL8 matrix2[3][3])
{
  INT4 i, j;

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      result[i][j] = matrix1[i][j] - matrix2[i][j];

  return;
}



/*
 * Transpose matrix
 */
static void
LALDR_Transpose33Matrix(REAL8 transpose[3][3],
                        REAL8 matrix[3][3])
{
  INT4 i, j;

  /*
   * Zero out output matrix
   */
  LALDR_Set33Matrix(transpose,
                    0., 0., 0.,
                    0., 0., 0.,
                    0., 0., 0.);

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      transpose[i][j] = matrix[j][i];

  return;
}



/*
 * The L2 norm of a matrix
 */
static REAL8
LALDR_L2Norm33Matrix(const REAL8 matrix[3][3])
{
    INT4   i, j;
    REAL8  l2norm = 0.;

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            l2norm += matrix[i][j] * matrix[i][j];

    l2norm = sqrt((double)l2norm);

    return l2norm;
}




/*
 * The RMS norm of a matrix: RMS sum of all elements.
 */
static REAL8
LALDR_RMSNorm33Matrix(const REAL8 matrix[3][3])
{
    INT4   i, j;
    REAL8  rmsnorm = 0.;

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            rmsnorm += matrix[i][j] * matrix[i][j];

    rmsnorm = sqrt((double)rmsnorm / 9.);

    return rmsnorm;
}



/*
 * The "infinity" norm of a matrix: max over all elems
 */
static REAL8
LALDR_InfNorm33Matrix(const REAL8 matrix[3][3])
{
    INT4 i, j;
    REAL8 infnorm = 0.;

    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            if (fabs((double)(matrix[i][j])) > infnorm)
                infnorm = fabs(matrix[i][j]);

    return infnorm;
}



/*
 * Print out matrix
 */
static void
LALDR_Print33Matrix(const REAL8 matrix[3][3],
                    const CHAR *varname,
                    UINT4       format,
                    FILE       *file,
                    const CHAR *graph_title)
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
      for (i = 0; i < 3; ++i)
        {
          for (j = 0; j < 3; ++j)
            fprintf(file, "%20.14e\t", matrix[i][j]);

          fprintf(file, "\n");
        }
    }


  /*
     * Maple V.3 format: makes a patchcontour plot
     */
  if (format == 1)
    {
      max = 0.;
      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          if (fabs(matrix[i][j]) > max)
            max = fabs(matrix[i][j]);


      /* load linalg and plots packages for Maple; plot options */
      fprintf(file, "with(linalg): with(plots):\n");
      fprintf(file, "setoptions3d(shading=ZHUE,style=PATCHCONTOUR,");
      fprintf(file, "projection=0.5,axes=BOXED,");
      fprintf(file, "title=`%s (abs. max = %10.5e)`):\n",
              graph_title, fabs(max));

      /* input "varname" is name of variable */
      fprintf(file, "%s := array(1 .. 3, 1 .. 3,[",
              varname);

      for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
          {
            /* last entry doesn't have trailing comma */
            if (i == 2 && j == 2)
              fprintf(file, "(3, 3)=%20.14e)",
                      fabs(matrix[i][j]));
            else
              fprintf(file, "(%d, %d)=%20.14e,",
                      j, i, fabs(matrix[i][j]));
          }

      fprintf(file, "]):\n");
      fprintf(file, "matrixplot(%s);\n", varname);
    }

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
LALDR_EulerRotation(REAL8        rotationMatrix[3][3],
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

    REAL8 e_xarm[3];
    REAL8 e_yarm[3];

    REAL8 normal[3];

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



/*
 * Make epsilon larger: GMST1 routine is only good to about 1.2e-5
 */
#define LALDR_REAL4_EPS 1.5e-05

void TestOnePoint(LALStatus              *status,
                  INT4                   *retval,
                  const LALFrDetector    *p_frDetector,
                  const LALSource        *p_source,
                  const LIGOTimeGPS      *p_gps,
                  const LALDetAMResponse *p_expectedResponse,
                  INT4                    use_LHO_p);

static REAL8 deg_to_rad(REAL8 deg);

int lalDebugLevel = 0;
int verbose_p       = 1;


int main(int argc, char *argv[])
{
  static LALStatus status;
  LALSource        pulsar;
  LALFrDetector    frdet;
  LIGOTimeGPS      gps;
  LALDate          date;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
  LALDetAMResponse expectedResponse;
  INT4             passed;
  INT4             testNumber;

  if (argc == 2)
    {
      lalDebugLevel = atoi(argv[1]);
      if (lalDebugLevel >= 8)
        verbose_p = 1;
      else
        verbose_p = 0;
    }

  if (verbose_p > 0)
    printf("LALDR_REAL4_EPS = %2.9e\n", LALDR_REAL4_EPS);

  /* Set source coord system to equatorial */
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /*****************
   * Test number 1
   *****************/
  testNumber = 1;

  /* source: RA=0, Dec=0, psi=0 */
  strcpy(pulsar.name, "TEST1: 0,0,0");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 0.;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 2
   ****************/
  testNumber = 2;

  /* source: RA=0, Dec=0, psi=pi/2 */
  strcpy(pulsar.name, "TEST2: 0,0,pi/2");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI_2;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)-1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 3
   ****************/
  testNumber = 3;

  /* source: RA=0, Dec=0, psi=pi */
  strcpy(pulsar.name, "TEST3: 0,0,pi");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 4
   ****************/
  testNumber = 4;

  /* source: RA=0, Dec=0, psi=3*pi/2 */
  strcpy(pulsar.name, "TEST4: 0,0,3*pi/2");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI + LAL_PI_2;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)-1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 5
   ****************/
  testNumber = 5;

  /* source: RA=0, Dec=0, psi=pi/4 */
  strcpy(pulsar.name, "TEST5: 0,0,pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)-1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 6
   ****************/
  testNumber = 6;

  /* source: RA=0, Dec=0, psi=pi/4 */
  strcpy(pulsar.name, "TEST6: 0,0,3*pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 3.*LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 7
   ****************/
  testNumber = 7;

  /* source: RA=0, Dec=0, psi=5*pi/4 */
  strcpy(pulsar.name, "TEST7: 0,0,5*pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 5.*LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)-1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 8
   ****************/
  testNumber = 8;

  /* source: RA=36deg, dec=0, orientation=0 */
  strcpy(pulsar.name, "TEST8: 36deg, 0deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(36.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(0.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST3: (36, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 36.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 9
   ****************/
  testNumber = 9;

  /* source: RA=0deg, dec=60deg, orientation=0 */
  strcpy(pulsar.name, "TEST9: 0deg, 60deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(0.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (0E, 60N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST3: (0, 60), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 10
   ****************/
  testNumber = 10;

  /* source: RA=30deg, dec=60, orientation=0 */
  strcpy(pulsar.name, "TEST10: 30deg, 60deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST4: (30, 60), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 11
   ****************/
  testNumber = 11;

  /* source: RA=30deg, dec=60deg, orientation=-30deg */
  strcpy(pulsar.name, "TEST11: 30deg, 60deg, -30deg");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(-30.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST5: (30, 60), (30, 120)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = LAL_PI / 6.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + LAL_PI / 6.;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 12
   ****************/
  testNumber = 12;

  /* source: RA=30deg, dec=60deg, orientation=-75deg */
  strcpy(pulsar.name, "TEST12: 30deg, 60deg, -75deg");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(-75.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST5: (30, 60), (30, 120)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = LAL_PI / 6.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + LAL_PI / 6.;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 13
   ****************/
  testNumber = 13;

  /* source: RA=27deg, dec=53.8deg, orientation=0deg */
  strcpy(pulsar.name, "TEST13: 27deg, -53.8deg, 0.");
  pulsar.equatorialCoords.longitude = deg_to_rad(27.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(-53.8);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (27E, 53.8N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST6: (27, 53.8), (0,90)");
  frdet.vertexLongitudeDegrees = 27.;
  frdet.vertexLatitudeDegrees  = 53.8;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 14
   ****************/
  testNumber = 14;

  /* source: RA=27deg, dec=53.8deg, orientation=73.19deg */
  strcpy(pulsar.name, "TEST14: 27deg, -53.8deg, 73.19");
  pulsar.equatorialCoords.longitude = deg_to_rad(27.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(-53.8);
  pulsar.orientation                = deg_to_rad(73.19);

  /* detector: (27E, 53.8N) x-arm: bearing 16.81; y-arm: bearing 286.81 */
  strcpy(frdet.name, "TEST7: (27, 53.8), (16.81,163.19)");
  frdet.vertexLongitudeDegrees = 27.;
  frdet.vertexLatitudeDegrees  = 53.8;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(73.19);  /* measured from E */
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + deg_to_rad(73.19);

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 15
   ****************/
  testNumber = 15;

  /* source */
  strcpy(pulsar.name, "TEST15: 90deg, 0deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(90.);
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 0.;

  /* detector */
  strcpy(frdet.name, "TEST7");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4).5;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 16
   ****************/
  testNumber = 16;

  /* source */
  strcpy(pulsar.name, "TEST16: 0deg, 0deg, -30deg");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = deg_to_rad(-30.);

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST8");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)0.5;
  expectedResponse.cross  = (REAL4)sqrt(3.)/2.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }




  return 0;
}


/* returns 1 if computed values are as expected, 0 otherwise */
/* only tests differential mode IFOs, and only plus and cross components */
void TestOnePoint(LALStatus              *status,
                  INT4                   *retval,
                  const LALFrDetector    *p_frDetector,
                  const LALSource        *p_source,
                  const LIGOTimeGPS      *p_gps,
                  const LALDetAMResponse *p_expectedResponse,
                  INT4                    use_LHO_p)
{
  LALDetector       detector;
  LALDetAndSource   detAndSource;
  LALDetAMResponse  computedResponse;
  char              infostr[1024];

  INITSTATUS(status, "TestOnePoint", LALTESTDETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  if (use_LHO_p == 0)
    {
      TRY( LALCreateDetector(status->statusPtr, &detector, p_frDetector,
                             LALDETECTORTYPE_IFODIFF), status );

      detAndSource.pDetector = &detector;
    }
  else
    {
      detAndSource.pDetector = &(lalCachedDetectors[LALDetectorIndexLHODIFF]);
    }

  detAndSource.pSource   = p_source;

  if (lalDebugLevel > 0 || verbose_p > 0)
    {
      sprintf(infostr,
              "Detector: (%2.9e deg. E, %2.9e deg. N)\n          (x-azi=%2.9e rad., y-azi=%2.9e rad.)",
              p_frDetector->vertexLongitudeDegrees,
              p_frDetector->vertexLatitudeDegrees,
              p_frDetector->xArmAzimuthRadians,
              p_frDetector->yArmAzimuthRadians);
      LALInfo(status, infostr);

      sprintf(infostr,
              "Source: (%2.9e RA, %2.9e Dec, %2.9e orientation) radians",
              p_source->equatorialCoords.longitude,
              p_source->equatorialCoords.latitude,
              p_source->orientation);
      LALInfo(status, infostr);
    }

  TRY( LALComputeDetAMResponse(status->statusPtr, &computedResponse,
                               &detAndSource, p_gps), status);

  if (lalDebugLevel > 0 || verbose_p > 0)
    {
      printf("expected: plus=%2.9e, cross=%2.9e  +/- %2.9e\n",
             p_expectedResponse->plus, p_expectedResponse->cross,
             LALDR_REAL4_EPS);
      printf("computed: plus=%2.9e, cross=%2.9e\n",
             computedResponse.plus, computedResponse.cross);
      printf("diff.:    plus=%2.9e, cross=%2.9e\n",
             (computedResponse.plus - p_expectedResponse->plus),
             (computedResponse.cross - p_expectedResponse->cross));
    }

  if ((REAL4)fabs((double)(computedResponse.plus - p_expectedResponse->plus))
      > LALDR_REAL4_EPS ||
      (REAL4)fabs((double)(computedResponse.cross - p_expectedResponse->cross))
      > LALDR_REAL4_EPS)
    {
      *retval = 0;
    }
  else
    {
      *retval = 1;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* Computes azimuth given dec, lat, alt (all in radians)
 *   range [0, Pi] */
REAL8 AM_hor2azimuth(REAL8 dec,REAL8 lat, REAL8 alt)
{
    REAL8 tmp, tmp2;
    REAL8 retval;

    tmp = 0.;
    tmp2 = 0.;
    retval = 0.;


    /*
     * Floating point checks for equality are kosher using IEEE arithmetic
     */
    if (lat == (REAL8)LAL_PI_2)
        return (REAL8)LAL_PI;  /* North Pole */
    else if (lat == -(REAL8)LAL_PI_2)
        return (REAL8)0.;    /* South Pole */
    else
    {
        if (alt == (REAL8)LAL_PI_2)
            return (REAL8)0.; /* Zenith or Nadir */
        else
        {
            tmp = sin(dec) / (cos(alt) * cos(lat))
                - tan(alt) * tan(lat);

            if (fabs(tmp) > 1.)
            {
                if (tmp < 0.)
                    tmp = -1.;

                if (tmp > 0.)
                    tmp = 1.;
            }

            if (lat == 0.)
            {
                tmp2 = alt + dec;

                if (dec < 0. && tmp2 == -(REAL8)LAL_PI_2)
                    tmp = (double)-1.;

                if (dec > 0. && tmp2 == (REAL8)LAL_PI_2)
                    tmp = (double)1.;
            }

            /* printf("Hor2Azi: azi = %19.15e\n", rad_to_deg(acos(tmp))); */

            return (REAL8)acos(tmp);
        }
    }
}

/* Computes altitude given dec (rad), lat (rad), ha (sid. sec) */
REAL8 AM_hor2altitude(REAL8 dec, REAL8 lat, REAL8 ha)
{
    REAL8 ha_r;

    ha_r = sidsec_to_rad(ha);

    if (lat == 0.)
    {
        if (ha == 21600. || ha == 64800.)
            return (REAL8)0.;
        else if (ha == 0.)
            return (REAL8)(LAL_PI_2 - fabs(dec));
        else if (ha == 43200.)
            return (REAL8)(fabs(dec) - LAL_PI_2);
        else
            return rad_to_deg(asin(cos(dec) * cos(ha)));
    }
    else if (dec == 0.)
    {
        if (ha == 21600. || ha == 64800.)
            return (REAL8)0.;
        else if (ha == 0.)
            return (REAL8)((REAL8)LAL_PI_2 - fabs(lat));
        else if (ha == 43200.)
            return (REAL8)(fabs(lat) - (REAL8)LAL_PI_2);
        else
            return rad_to_deg(asin(cos(lat) * cos(ha)));
    }
    else
        return rad_to_deg(asin(sin(dec) * sin(lat)
                               + cos(dec) * cos(lat) * cos(ha)));
}

/*
 * Converts src coords in Equatorial to Horizon coords.
 * lst  -  Local Sidereal Time in seconds
 */
void AM_Eq2Hor(Status                 *status,
               SkyPosition            *srcHorizon,
               const SkyPosition      *srcEquatorial,
               const EarthPosition    *observer,
               REAL8                   lst)
{
    REAL8 ha; /* in seconds */

    INITSTATUS (status, "Eq2Hor", LALTESTDETRESPONSEC);

    /* compute Hour Angle of source in seconds */
    ha = lst - rad_to_deg(srcEquatorial->longitude) * (REAL8)3600.;

    if (ha > 86400.)
        ha -= 86400.;

    srcHorizon->latitude = AM_hor2altitude(srcEquatorial->latitude,
                                           observer->geodetic.latitude,
                                           ha);

    srcHorizon->longitude = AM_hor2azimuth(srcEquatorial->latitude,
                                           observer->geodetic.latitude,
                                           srcHorizon->latitude);

    /* take care of azi of source based on HA */
    if (ha < 0. || ha > 43200.)
        srcHorizon->longitude = (REAL8)LAL_2_PI - srcHorizon->longitude;

    if (srcHorizon->longitude == (REAL8)LAL_2_PI)
        srcHorizon->longitude = 0.;

    RETURN (status);
}


/*
 * st_response() -- response according to Schutz and Tinto
 */
void
st_response(LALStatus             *status,
            LALDetAMResponse      *pResponse,
            const LALDetAndSource *pDetAndSrc,
            const LIGOTimeGPS     *pGPS)
{
  char infostr[1024];

  /* Euler rotations to transfrom from src to Earth-fixed */
  REAL8 r1[3][3];
  REAL8 r2[3][3];
  REAL8 r3[3][3];
  REAL8 tmpmat[3][3];
  REAL8 rtot[3][3];
  REAL8 rtot_T[3][3];

  REAL8 xi, theta, phi, psi, open;
  REAL8 sin2open;



  INITSTATUS( status, "st_response", LALTESTDETRESPONSEC );

  /*
   * Convert source coords to alt-azi
   */


  xi = deg_to_rad(pDetAndSrc->pSource->orientation) + LAL_PI_2;





