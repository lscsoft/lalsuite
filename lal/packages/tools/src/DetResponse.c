/*<lalVerbatim file="DetResponseCV">

Author: David Chin <dwchin@umich.edu> +1-734-709-9119
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
Anderson, \textit{et al.}~\cite{tools:Anderson:2000}.  We compute the $h$-tensors for
$+$- and $\times$-polarized in the Earth-fixed frame, and then contract
them (take the scalar product) with the detector response tensors as
described in the \texttt{DetectorSite.h} section of the \texttt{tools}
package. 

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

/* static const INT4 oneBillion = 1000000000; */

#define LALDR_MATRIXSIZE 3

/* I love C, really I do.
 * These typedefs are to get rid of type mismatch warnings. (Lint is our
 * friend.)  */
typedef REAL8 LALDR_3Vector[LALDR_MATRIXSIZE];
typedef REAL8 LALDR_33Matrix[LALDR_MATRIXSIZE][LALDR_MATRIXSIZE];

/*
 * Private functions
 */
#if 0
static REAL8 deg_to_rad(REAL8 degrees)
{
  return degrees * (REAL8)LAL_PI / (REAL8)180.;
}

static REAL8 rad_to_deg(REAL8 radians)
{
  return radians * (REAL8)180. / (REAL8)LAL_PI;
}
#endif


/* axis for LALDR_EulerRotation() */
typedef enum { xAxis = 1, yAxis = 2, zAxis = 3 } LALDR_Axis_t;

/* This #if is so I can hide this block in Emacs's hide-ifdef-mode */
#if 1
#if 0
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
#endif



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



/*
 * Copy matrix source to matrix target
 */
#if 0
static void
LALDR_Copy33Matrix(LALDR_33Matrix * target, LALDR_33Matrix * source)
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
#endif



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
#if 0
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



/*
 * The L2 norm of a matrix
 */
#if 0
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
#endif
     





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





/* YES, THIS IS MEANT TO BE #if'ed out */
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


/*
 * Another way of computing detector response, making use
 * of the detector response tensor stored in LALDetector
 */

/* <lalVerbatim file="DetResponseCP"> */
void
LALComputeDetAMResponse( LALStatus             *status,
                         LALDetAMResponse      *pResponse,
                         const LALDetAndSource *pDetAndSrc,
                         const LALGPSandAcc    *pGPSandAcc)
{ /* </lalVerbatim> */

  char infostr[1024];
  LALMSTUnitsAndAcc units_and_acc;
  
  /* We need to express the metric perturbation tensor in the Earth-fixed 
   * frame, and then take its dot product with the detector response tensor
   * as stored in the LALDetector structure. */

  /* "unit" plus perturbation */
  LALDR_33Matrix e_plus;
  LALDR_33Matrix e_plus_E;

  /* "unit" cross perturbation */
  LALDR_33Matrix e_cros; 
  LALDR_33Matrix e_cros_E;

  /* Euler angles to transform from src to Earth-fixed */
  REAL8 psi;      /* orientation of src: ccw "N of W" */
  REAL8 theta;    /* Pi/2 + Declination */
  REAL8 phi;      /* RA - GMST1 */
  REAL8 gmst;     /* GMST1, in RADIANS */

  /* Euler rotations to transfrom from src to Earth-fixed */
  LALDR_33Matrix r1;
  LALDR_33Matrix r2;
  LALDR_33Matrix r3;
  LALDR_33Matrix tmpmat;
  LALDR_33Matrix rtot;
  LALDR_33Matrix rtot_T;

  /* Since LALDetector stores response in REAL4, and we do everything
   * here in REAL8, I need a temporary one to store response */
  LALDR_33Matrix response;
  INT4 i, j;

   
  INITSTATUS( status, "LALComputeDetAMResponse", DETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pResponse != (LALDetAMResponse *)NULL, status,
         DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

  ASSERT(pDetAndSrc != NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  ASSERT(pGPSandAcc != (LALGPSandAcc *)NULL, status,
         DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

  /* source coordinates must be in equatorial system */
  ASSERT(pDetAndSrc->pSource->equatorialCoords.system ==
         COORDINATESYSTEM_EQUATORIAL, status,
         DETRESPONSEH_ESRCNOTEQUATORIAL, DETRESPONSEH_MSGESRCNOTEQUATORIAL);


  /*
   * Initialize e_plus and e_cros (in the obvious frame)
   */
  LALDR_Set33Matrix(&e_plus, 1.,  0., 0.,
                    0., -1., 0.,
                    0.,  0., 0.);

  LALDR_Set33Matrix(&e_cros, 0., 1., 0.,
                    1., 0., 0.,
                    0., 0., 0.);

  /*
   * Here are the rotations required to get the perturbation tensors
   * into Earth-fixed basis
   * Extract the angles (for convenience):
   */
  units_and_acc.units = MST_RAD;
  units_and_acc.accuracy = pGPSandAcc->accuracy;
  /* compute phi */
  TRY( LALGPStoGMST1(status->statusPtr, &gmst, &(pGPSandAcc->gps),
                     &units_and_acc), status );

  if (lalDebugLevel > 4)
    {
      sprintf(infostr, "gps  = %9d:%9d\n", pGPSandAcc->gps.gpsSeconds,
              pGPSandAcc->gps.gpsNanoSeconds);
      sprintf(infostr, "gmst = %g radians", gmst);
      LALInfo(status, infostr);
    }
   
  phi = pDetAndSrc->pSource->equatorialCoords.longitude - gmst
    - (REAL8)LAL_PI_2;

  /* compute theta */
  theta = pDetAndSrc->pSource->equatorialCoords.latitude + (REAL8)LAL_PI_2;

  /* psi */
  psi = pDetAndSrc->pSource->orientation;

  /* Now, form the Euler rotations, and compute total rotation */
  LALDR_EulerRotation(&r1, -psi,   zAxis);
  LALDR_EulerRotation(&r2, -theta, xAxis);
  LALDR_EulerRotation(&r3, -phi,   zAxis);

#ifdef DEBUG  
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
#endif  

  LALDR_Multiply33Matrix(&tmpmat, &r2, &r1);
  LALDR_Multiply33Matrix(&rtot, &r3, &tmpmat);

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
  LALDR_Transpose33Matrix(&rtot_T, &rtot);
  LALDR_Multiply33Matrix(&tmpmat, &e_plus, &rtot_T);
  LALDR_Multiply33Matrix(&e_plus_E, &rtot, &tmpmat);

  LALDR_Transpose33Matrix(&rtot_T, &rtot);
  LALDR_Multiply33Matrix(&tmpmat, &e_cros, &rtot_T);
  LALDR_Multiply33Matrix(&e_cros_E, &rtot, &tmpmat);

  /*
   * Then, F_plus = det_response &. e_plus_E,
   *       F_cros = det_response &. e_cros_E
   */
  for (i = 0; i < LALDR_MATRIXSIZE; ++i)
    for (j = 0; j < LALDR_MATRIXSIZE; ++j)
      response[i][j] = (REAL8)(pDetAndSrc->pDetector->response[i][j]);

#ifdef DEBUG  
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
#endif  

  pResponse->plus  = (REAL4)LALDR_DotProd33Matrix(&response, &e_plus_E);
  pResponse->cross = (REAL4)LALDR_DotProd33Matrix(&response, &e_cros_E);
  /* FIXME: scalar response not implemented, yet.  Will have to
     read Waggoner's paper to do this. */
  pResponse->scalar = 0.;

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
  LIGOTimeGPS       tmpgps;
  LALGPSandAcc      gps_and_acc;
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
   * Set names
   */
  pResponseSeries->pPlus->name[0]   = '\0';
  pResponseSeries->pCross->name[0]  = '\0';
  pResponseSeries->pScalar->name[0] = '\0';

  strncpy(pResponseSeries->pPlus->name, "plus", LALNameLength);
  strncpy(pResponseSeries->pCross->name, "cross", LALNameLength);
  strncpy(pResponseSeries->pScalar->name, "scalar", LALNameLength);

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
      if (lalDebugLevel & 0x08)
        LALInfo(status, "scalar sequence too short -- reallocating");

      TRY( LALSDestroyVector(status->statusPtr,
                             &(pResponseSeries->pScalar->data)), status );

      TRY( LALSCreateVector(status->statusPtr,
                            &(pResponseSeries->pScalar->data),
                            pTimeInfo->nSample), status );
    }
  

  /*
   * Loop to compute each element in time series; rint(3) is a std C
   * function that rounds floating point numbers properly.
   */
  gps = pTimeInfo->epoch;

  TRY(LALFloatToInterval(status->statusPtr, &dt, &(pTimeInfo->deltaT)),
      status);

  for (i = 0; i < pTimeInfo->nSample; ++i)
    {
      if (lalDebugLevel >= 8)
        {
          sprintf(infostr, "LALComputeDetAMResponseSeries: i = %d\n", i);
          LALInfo(status, infostr);
        }

      gps_and_acc.gps = gps;
      gps_and_acc.accuracy = pTimeInfo->accuracy;
      TRY( LALComputeDetAMResponse(status->statusPtr, &instResponse,
                                   pDetAndSource, &gps_and_acc), status );

      pResponseSeries->pPlus->data->data[i]   = instResponse.plus;
      pResponseSeries->pCross->data->data[i]  = instResponse.cross;
      pResponseSeries->pScalar->data->data[i] = instResponse.scalar;

      tmpgps = gps;
      
      TRY(LALIncrementGPS(status->statusPtr, &gps, &tmpgps, &dt), status);

    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALComputeDetAMResponseSeries() */
