/*<lalVerbatim file="DetResponseCV">

Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim> */

/*
<lalLaTeX>

\subsection{Module \texttt{DetResponse.c}}
\label{ss:DetResponse.c}

Computes the response of a detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DetResponseCP}
\idx{LALComputeDetAMResponse()}
\idx{LALComputeDetAMResponseSeries()}

\subsubsection*{Description}

These routines compute the antenna beam pattern for all supported detector
types.  \texttt{LALComputeDetAMResponse()} computes the response at one
instance in time, and \texttt{LALComputeDetAMResponseSeries()} computes a
vector of response for some length of time. 

\subsubsection*{Algorithm}

This code is a translation of the algorithm in the Maple worksheet by
Anderson, \textit{et al.}~\cite{ABCF:2000}.  We compute the $h$-tensors for
$+$- and $\times$-polarized in the Earth-fixed frame, and then contract
them (take the scalar product) with the detector response tensors as
described in the \texttt{DetectorSite.h} section of the \texttt{tools}
package. This is a translation of the Maple worksheet by
Anderson,~\textit{et al.}~\cite{ABCF:2000}.

\texttt{DetectorSite.h} in the \texttt{tools} package  provides predefined
\texttt{LALDetector} structures representing most current detectors,
including LIGO (Hanford and Livingston), and GEO.

\subsubsection*{Uses}
\texttt{LALGPStoGMST1()}

\subsubsection*{Notes}

For examples of usage, please see the test programs in the \texttt{test}
directory.

\vfill{\footnotesize\input{DetResponseCV}}

</lalLaTeX>
*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALError.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetResponse.h>

NRCSID( DETRESPONSEC, "$Id$" );


/* Matrix manipulation routines, and other stuff */
/* #include "DR_Util.h" */



/*
 * These are private utility functions for the DetResponse package: vector
 * and array manipulations, etc.
 * Author: David Chin <dwchin@umich.edu> +1-734-730-1274
 * $Id$
 */

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
#if 0
static void
LALDR_Print33Matrix(const REAL8 matrix[3][3],
                    const CHAR *varname,
                    UINT4       format,
                    FILE *file,
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
#endif





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
#if 0
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
#endif



/*
 * Another way of computing detector response, making use
 * of the detector response tensor stored in LALDetector
 */

/* <lalVerbatim file="DetResponseCP"> */
void
LALComputeDetAMResponse( LALStatus             *status,
                         LALDetAMResponse      *pResponse,
                         const LALDetAndSource *pDetAndSrc,
                         const LIGOTimeGPS     *pGPS)
{ /* </lalVerbatim> */

  char infostr[1024];
  
  /* We need to express the metric perturbation tensor in the Earth-fixed 
   * frame, and then take its dot product with the detector response tensor
   * as stored in the LALDetector structure. */

  /* "unit" plus perturbation */
  REAL8 e_plus[3][3];
  REAL8 e_plus_E[3][3];

  /* "unit" cross perturbation */
  REAL8 e_cros[3][3]; 
  REAL8 e_cros_E[3][3];

  /* Euler angles to transform from src to Earth-fixed */
  REAL8 psi;      /* orientation of src: ccw "N of W" */
  REAL8 theta;    /* Pi/2 + Declination */
  REAL8 phi;      /* RA - GMST1 */
  REAL8 gmst;     /* GMST1, in RADIANS */

  /* Euler rotations to transfrom from src to Earth-fixed */
  REAL8 r1[3][3];
  REAL8 r2[3][3];
  REAL8 r3[3][3];
  REAL8 tmpmat[3][3];
  REAL8 rtot[3][3];
  REAL8 rtot_T[3][3];

  /* Since LALDetector stores response in REAL4, and we do everything
   * here in REAL8, I need a temporary one to store response */
  REAL8 response[3][3];
  INT4 i, j;

   
  INITSTATUS( status, "LALComputeDetAMResponse", DETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pResponse != (LALDetAMResponse *)NULL, status,
         DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

  ASSERT(pDetAndSrc != NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  ASSERT(pGPS != (LIGOTimeGPS *)NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  /* source coordinates must be in equatorial system */
  ASSERT(pDetAndSrc->pSource->equatorialCoords.system == COORDINATESYSTEM_EQUATORIAL, status,
         DETRESPONSEH_ESRCNOTEQUATORIAL, DETRESPONSEH_MSGESRCNOTEQUATORIAL);


  /*
   * Initialize e_plus and e_cros (in the obvious frame)
   */
  LALDR_Set33Matrix(e_plus, 1.,  0., 0.,
                    0., -1., 0.,
                    0.,  0., 0.);

  LALDR_Set33Matrix(e_cros, 0., 1., 0.,
                    1., 0., 0.,
                    0., 0., 0.);

  /*
   * Here are the rotations required to get the perturbation tensors
   * into Earth-fixed basis
   * Extract the angles (for convenience):
   */
  /* compute phi */
  TRY( LALGPStoGMST1(status->statusPtr, &gmst, pGPS, MST_RAD), status );

  if (lalDebugLevel > 0)
    {
      sprintf(infostr, "gmst = %2.9e radians", gmst);
      LALInfo(status, infostr);
    }
   
  phi = pDetAndSrc->pSource->equatorialCoords.longitude - gmst
    - (REAL8)LAL_PI_2;

  /* compute theta */
  theta = pDetAndSrc->pSource->equatorialCoords.latitude + (REAL8)LAL_PI_2;

  /* psi */
  psi = pDetAndSrc->pSource->orientation;

  /* Now, form the Euler rotations, and compute total rotation */
  LALDR_EulerRotation(r1, -psi,   zAxis);
  LALDR_EulerRotation(r2, -theta, xAxis);
  LALDR_EulerRotation(r3, -phi,   zAxis);

  if (lalDebugLevel >= 8)
    {
      sprintf(infostr,
              "Rz(-psi) = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              r1[0][0], r1[0][1], r1[0][2],
              r1[1][0], r1[1][1], r1[1][2],
              r1[2][0], r1[2][1], r1[2][2]); 

      LALInfo(status, infostr);

      sprintf(infostr,
              "Rx(-theta) = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              r2[0][0], r2[0][1], r2[0][2],
              r2[1][0], r2[1][1], r2[1][2],
              r2[2][0], r2[2][1], r2[2][2]);

      LALInfo(status, infostr);

      sprintf(infostr,
              "Rz(-phi) = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              r3[0][0], r3[0][1], r3[0][2],
              r3[1][0], r3[1][1], r3[1][2],
              r3[2][0], r3[2][1], r3[2][2]);

      LALInfo(status, infostr);
    }

  LALDR_Multiply33Matrix(tmpmat, r2, r1);
  LALDR_Multiply33Matrix(rtot, r3, tmpmat);

  if (lalDebugLevel >= 8)
    {
      sprintf(infostr,
              "Rtot = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              rtot[0][0], rtot[0][1], rtot[0][2],
              rtot[1][0], rtot[1][1], rtot[1][2],
              rtot[2][0], rtot[2][1], rtot[2][2]);

      LALInfo(status, infostr);      
    }
  

  /*
   * Now, get the "unit" perturbation tensors in the Earth fixed frame.
   *
   *    e_plus_Earth = rtot &* e_plus &* transpose(rtot)
   *    e_cros_Earth = rtot &* e_cros &* transpose(rtot)
   */
  LALDR_Transpose33Matrix(rtot_T, rtot);
  LALDR_Multiply33Matrix(tmpmat, e_plus, rtot_T);
  LALDR_Multiply33Matrix(e_plus_E, rtot, tmpmat);

  LALDR_Transpose33Matrix(rtot_T, rtot);
  LALDR_Multiply33Matrix(tmpmat, e_cros, rtot_T);
  LALDR_Multiply33Matrix(e_cros_E, rtot, tmpmat);

  /*
   * Then, F_plus = det_response &. e_plus_E,
   *       F_cros = det_response &. e_cros_E
   */
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      response[i][j] = (REAL8)(pDetAndSrc->pDetector->response[i][j]);

  if (lalDebugLevel >= 8)
    {

      sprintf(infostr,
              "LAL detector response tensor = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              pDetAndSrc->pDetector->response[0][0],
              pDetAndSrc->pDetector->response[0][1],
              pDetAndSrc->pDetector->response[0][2],
              pDetAndSrc->pDetector->response[1][0],
              pDetAndSrc->pDetector->response[1][1],
              pDetAndSrc->pDetector->response[1][2],
              pDetAndSrc->pDetector->response[2][0],
              pDetAndSrc->pDetector->response[2][1],
              pDetAndSrc->pDetector->response[2][2]);

      LALInfo(status, infostr);
       
      sprintf(infostr,
              "e_Plus_E = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              e_plus_E[0][0], e_plus_E[0][1], e_plus_E[0][2],
              e_plus_E[1][0], e_plus_E[1][1], e_plus_E[1][2],
              e_plus_E[2][0], e_plus_E[2][1], e_plus_E[2][2]);

      LALInfo(status, infostr);

      sprintf(infostr,
              "e_cros_E = \n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n\t%2.9e   %2.9e   %2.9e\n",
              e_cros_E[0][0], e_cros_E[0][1], e_cros_E[0][2],
              e_cros_E[1][0], e_cros_E[1][1], e_cros_E[1][2],
              e_cros_E[2][0], e_cros_E[2][1], e_cros_E[2][2]);

      LALInfo(status, infostr);      
    }

  pResponse->plus  = (REAL4)LALDR_DotProd33Matrix(response, e_plus_E);
  pResponse->cross = (REAL4)LALDR_DotProd33Matrix(response, e_cros_E);

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALDetResponse() */



/*
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 */
/* <lalVerbatim file="DetResponseCP"> */
void
LALComputeDetAMResponseSeries( LALStatus              *status,
                               LALDetAMResponseSeries *pResponseSeries,
                               const LALDetAndSource  *pDetAndSource,
                               const LALTimeIntervalAndNSample *pTimeInfo)
{ /* </lalVerbatim> */
  /* Want to loop over the time and call LALComputeDetAMResponse() */
  LALDetAMResponse  instResponse;
  LIGOTimeGPS       gps;
  LALTimeInterval   dt;
  UINT4             i;
  char              infostr[128];

  INITSTATUS(status, "LALComputeDetAMResponseSeries", DETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  if (lalDebugLevel >= 8)
    {
      sprintf(infostr,
              "pResponseSeries->pPlus->data->length = %d\npTimeInfo->nSample = %d\n",
              pResponseSeries->pPlus->data->length, pTimeInfo->nSample);
      LALInfo(status, infostr);
    }

  /*
   * Error-checking assertions
   */
  ASSERT(pResponseSeries != (LALDetAMResponseSeries *)NULL, status,
         DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

  ASSERT(pDetAndSource != (LALDetAndSource *)NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  ASSERT(pTimeInfo != (LALTimeIntervalAndNSample *)NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  /*
   * Set sampling parameters
   */
  pResponseSeries->pPlus->epoch       = pTimeInfo->epoch;
  pResponseSeries->pPlus->deltaT      = pTimeInfo->deltaT;
  pResponseSeries->pPlus->f0          = 0.;
  pResponseSeries->pPlus->sampleUnits = lalDimensionlessUnit;

  pResponseSeries->pCross->epoch       = pTimeInfo->epoch;
  pResponseSeries->pCross->deltaT      = pTimeInfo->deltaT;
  pResponseSeries->pCross->f0          = 0.;
  pResponseSeries->pCross->sampleUnits = lalDimensionlessUnit;

  pResponseSeries->pScalar->epoch       = pTimeInfo->epoch;
  pResponseSeries->pScalar->deltaT      = pTimeInfo->deltaT;
  pResponseSeries->pScalar->f0          = 0.;
  pResponseSeries->pScalar->sampleUnits = lalDimensionlessUnit;
  
  /*
   * Ensure enough memory for requested vectors
   */
  if (pResponseSeries->pPlus->data->length < pTimeInfo->nSample)
    {
      if (lalDebugLevel >= 8)
        LALInfo(status, "plus sequence too short -- reallocating");

      TRY( LALSDestroyVector(status->statusPtr,
                             &(pResponseSeries->pPlus->data)), status );

      TRY( LALSCreateVector(status->statusPtr,
                            &(pResponseSeries->pPlus->data),
                            pTimeInfo->nSample), status );

      if (lalDebugLevel > 0)
        printf("pResponseSeries->pPlus->data->length = %d\n",
               pResponseSeries->pPlus->data->length);
      
    }

  if (pResponseSeries->pCross->data->length < pTimeInfo->nSample)
    {
      if (lalDebugLevel >= 8)      
        LALInfo(status, "cross sequence too short -- reallocating");

      TRY( LALSDestroyVector(status->statusPtr,
                             &(pResponseSeries->pCross->data)), status );

      TRY( LALSCreateVector(status->statusPtr,
                            &(pResponseSeries->pCross->data),
                            pTimeInfo->nSample), status );

    }

  if (pResponseSeries->pScalar->data->length < pTimeInfo->nSample)
    {
      if (lalDebugLevel >= 8)      
        LALInfo(status, "scalar sequence too short -- reallocating");

      TRY( LALSDestroyVector(status->statusPtr,
                             &(pResponseSeries->pScalar->data)), status );

      TRY( LALSCreateVector(status->statusPtr,
                            &(pResponseSeries->pScalar->data),
                            pTimeInfo->nSample), status );
    }
  

  /*
   * Loop to compute each element in time series
   */
  gps = pTimeInfo->epoch;
  dt.seconds     = (INT4)(pTimeInfo->deltaT);
  dt.nanoSeconds = (INT4)(pTimeInfo->deltaT - (REAL8)(dt.seconds)) *
    1000000000;

  for (i = 0; i < pTimeInfo->nSample; ++i)
    {
      if (lalDebugLevel >= 8)
        {
          sprintf(infostr, "LALComputeDetAMResponseSeries: i = %d\n", i);
          LALInfo(status, infostr);
        }
      
      TRY( LALComputeDetAMResponse(status->statusPtr, &instResponse,
                                   pDetAndSource, &gps), status );

      pResponseSeries->pPlus->data->data[i]   = instResponse.plus;
      pResponseSeries->pCross->data->data[i]  = instResponse.cross;
      pResponseSeries->pScalar->data->data[i] = instResponse.scalar;
      
      gps.gpsSeconds     += dt.seconds;
      gps.gpsNanoSeconds += dt.nanoSeconds;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALComputeDetAMResponseSeries() */
