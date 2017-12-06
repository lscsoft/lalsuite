#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_deriv.h>

#include "LALSimIMRSpinEOB.h"
#include "LALSimIMRSpinEOBHcapNumericalDerivative.c"
#include "LALSimIMRSpinEOBHcapNumericalDerivativePrec.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMREOBFactorizedWaveform.c"

#ifndef _LALSIMIMRSPINPRECEOBINITIALCONDITIONS_C
#define _LALSIMIMRSPINPRECEOBINITIALCONDITIONS_C

/* The following values are used to scale the inputs to the equations being
 * solved to get the initial spherical orbit. By scaling we bring different
 * inputs to the same scale, and drastically increase the rate of convergence.
 * */
REAL8 scale1 = 1, scale2 = 2, scale3 = 200;

/**
 *  Calculate the dot product of two vectors
 */
static		inline
		REAL8
CalculateDotProductPrec(
                        const REAL8 a[], /**<< first vector */
                        const REAL8 b[] /**<< second vector */)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


/**
 * Calculates the ith component of the cross product of two vectors,
 * where i is in the range 0-2.
 */
static		inline
		REAL8
CalculateCrossProductPrec(
                          const INT4 i, /**<< i-th component of cross product */
                          const REAL8 a[],  /**<< first vector */
                          const REAL8 b[] /**<< second vector */)
{
	return a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3];
}


/**
 * Normalizes the given vector
 */
static		inline
int
NormalizeVectorPrec(REAL8 a[])
{
	REAL8		norm = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	a[0] /= norm;
	a[1] /= norm;
	a[2] /= norm;

	return XLAL_SUCCESS;
}

/**
 * Calculate the rotation matrix and its inverse.
 * Rotate the ex-ey-ez frame to the r-v-L frame.
 * This static function is called only by XLALSimIMRSpinEOBInitialConditionsPrec
 */
static int
CalculateRotationMatrixPrec(
			gsl_matrix * rotMatrix,	/**< OUTPUT, rotation matrix */
			gsl_matrix * rotInverse,	/**< OUTPUT, rotation matrix inversed */
			REAL8 r[],	/**< position vector */
			REAL8 v[],	/**< velocity vector */
			REAL8 L[]	/**< orbital angular momentum */
)
{

	/** a, b, g are the angles alpha, beta and gamma */
	/* Use a, b and g to avoid shadowing gamma and getting a warning */
	REAL8		a       , b, g;
	REAL8		cosa    , sina, cosb, sinb, cosg, sing;

	/* Calculate the Euclidean Euler angles */

	/* We need to avoid a singularity if L is along z axis */
	if (L[2] > 0.9999) {
		a = b = g = 0.0;
	} else {
		a = atan2(L[0], -L[1]);
		b = atan2(sqrt(L[0] * L[0] + L[1] * L[1]), L[2]);
		g = atan2(r[2], v[2]);
	}

	if ((cosa = cos(a)) < 1.0e-16)
		cosa = 0.0;
	if ((sina = sin(a)) < 1.0e-16)
		sina = 0.0;
	if ((cosb = cos(b)) < 1.0e-16)
		cosb = 0.0;
	if ((sinb = sin(b)) < 1.0e-16)
		sinb = 0.0;
	if ((cosg = cos(g)) < 1.0e-16)
		cosg = 0.0;
	if ((sing = sin(g)) < 1.0e-16)
		sing = 0.0;

	/**
	 * Implement the Rotation Matrix following the "x-convention"
	 * 1. rotate about the z-axis by an angle a, rotation matrix Rz(a);
	 * 2. rotate about the former x-axis (now x') by an angle b, rotation matrix Rx(b);
	 * 3. rotate about the former z-axis (now z') by an algle g, rotation matrix Rz(g);
	 */
	/* populate the matrix */
	/*
	 * gsl_matrix_set( rotMatrix, 0, 0, cosg*cosa - cosb*sina*sing );
	 * gsl_matrix_set( rotMatrix, 0, 1, cosg*sina + cosb*cosa*sing );
	 * gsl_matrix_set( rotMatrix, 0, 2, sing*sinb ); gsl_matrix_set(
	 * rotMatrix, 1, 0, -sing*cosa - cosb*sina*cosg ); gsl_matrix_set(
	 * rotMatrix, 1, 1, -sing*sina + cosb*cosa*cosg ); gsl_matrix_set(
	 * rotMatrix, 1, 2, cosg*sinb ); gsl_matrix_set( rotMatrix, 2, 0,
	 * sinb*sina ); gsl_matrix_set( rotMatrix, 2, 1, -sinb*cosa );
	 * gsl_matrix_set( rotMatrix, 2, 2, cosb );
	 */
	//Better to avoid angle singularity
		gsl_matrix_set(rotMatrix, 0, 0, r[0]);
	gsl_matrix_set(rotMatrix, 0, 1, r[1]);
	gsl_matrix_set(rotMatrix, 0, 2, r[2]);
	gsl_matrix_set(rotMatrix, 1, 0, v[0]);
	gsl_matrix_set(rotMatrix, 1, 1, v[1]);
	gsl_matrix_set(rotMatrix, 1, 2, v[2]);
	gsl_matrix_set(rotMatrix, 2, 0, L[0]);
	gsl_matrix_set(rotMatrix, 2, 1, L[1]);
	gsl_matrix_set(rotMatrix, 2, 2, L[2]);


	/* Now populate the transpose (which should be the inverse) */
	gsl_matrix_transpose_memcpy(rotInverse, rotMatrix);

	/* Test that the code does what it should do */

	/*
	 * gsl_matrix *ab = gsl_matrix_alloc( 3, 3 ); gsl_blas_dgemm(
	 * CblasNoTrans, CblasNoTrans, 1., rotMatrix, rotInverse, 0., ab );
	 */

	/*
	 * XLAL_PRINT_INFO( "The generated rotation matrix: \n" ); for ( int i = 0; i
	 * < 3; i++ ) { for ( int j = 0; j < 3; j++ ) { XLAL_PRINT_INFO( "%.16e\t",
	 * gsl_matrix_get( rotMatrix, i, j ) ); } XLAL_PRINT_INFO( "\n" ); }
	 */

	/*
	 * XLAL_PRINT_INFO( "Is the transpose of the rotation matrix the inverse?\n"
	 * ); for ( int i = 0; i < 3; i++ ) { for ( int j = 0; j < 3; j++ ) {
	 * XLAL_PRINT_INFO( "%.16e\t", gsl_matrix_get( ab, i, j ) ); } XLAL_PRINT_INFO( "\n" );
	 * }
	 */
	/* gsl_matrix_free( ab ); */

	return XLAL_SUCCESS;
}


/**
 * Applies a rotation matrix to a given vector
 */
static int
ApplyRotationMatrixPrec(
		    gsl_matrix * rotMatrix,	/**< rotation matrix */
		    REAL8 a[]	/**< OUTPUT, vector rotated */
)
{

	double b[3];
	gsl_vector_view tmpView1 = gsl_vector_view_array(a, 3);
	gsl_vector_view tmpView2 = gsl_vector_view_array(b, 3);
	gsl_vector     *tmpVec1 = &tmpView1.vector;
	gsl_vector     *tmpVec2 = &tmpView2.vector;

	gsl_blas_dgemv(CblasNoTrans, 1.0, rotMatrix, tmpVec1, 0.0, tmpVec2);

	gsl_vector_memcpy(tmpVec1, tmpVec2);

	return XLAL_SUCCESS;
}

/**
 * Performs a co-ordinate transformation from spherical to Cartesian co-ordinates.
 * In the code from Tyson Littenberg, this was only applicable to the special theta=pi/2, phi=0 case.
 * This special transformation is a static function called only by
 * GSLSpinHamiltonianDerivWrapperPrec and XLALSimIMRSpinEOBInitialConditionsPrec
 */
static int
SphericalToCartesianPrec(
		     REAL8 qCart[],	/**<< OUTPUT, position vector in Cartesean coordinates */
		     REAL8 pCart[],	/**<< OUTPUT, momentum vector in Cartesean coordinates */
		     const REAL8 qSph[],	/**<< position vector in spherical coordinates */
		     const REAL8 pSph[]	/**<< momentum vector in spherical coordinates */
)
{

	REAL8		r;
	REAL8		pr      , pTheta, pPhi;

	REAL8		x       , y, z;
	REAL8		pX      , pY, pZ;

	r = qSph[0];
	pr = pSph[0];
	pTheta = pSph[1];
	pPhi = pSph[2];

	x = r;
	y = 0.0;
	z = 0.0;

	pX = pr;
	pY = pPhi / r;
	pZ = -pTheta / r;

	/* Copy these values to the output vectors */
	qCart[0] = x;
	qCart[1] = y;
	qCart[2] = z;
	pCart[0] = pX;
	pCart[1] = pY;
	pCart[2] = pZ;

	return XLAL_SUCCESS;
}


/**
 * Perform a co-ordinate transformation from Cartesian to spherical co-ordinates.
 * In the code from Tyson Littenberg, this was only applicable to the special theta=pi/2, phi=0 case.
 * This special transformation is a static function called only by XLALSimIMRSpinEOBInitialConditionsPrec
 */
static int
CartesianToSphericalPrec(
		     REAL8 qSph[],	/**<< OUTPUT, position vector in spherical coordinates */
		     REAL8 pSph[],	/**<< OUTPUT, momentum vector in Cartesean coordinates */
		     const REAL8 qCart[],	/**<< position vector in spherical coordinates */
		     const REAL8 pCart[]	/**<< momentum vector in Cartesean coordinates */
)
{

	REAL8		r;
	REAL8		pr      , pTheta, pPhi;

	REAL8		x;
	//, y, z;
	REAL8		pX      , pY, pZ;

	x = qCart[0];
	//y = qCart[1];
	//z = qCart[2];
	pX = pCart[0];
	pY = pCart[1];
	pZ = pCart[2];

	r = x;
	pr = pX;
	pTheta = -r * pZ;
	pPhi = r * pY;

	/* Fill the output vectors */
	qSph[0] = r;
	qSph[1] = LAL_PI_2;
	qSph[2] = 0.0;
	pSph[0] = pr;
	pSph[1] = pTheta;
	pSph[2] = pPhi;

	return XLAL_SUCCESS;
}


/**
 * Root function for gsl root finders.
 * The gsl root finder will look for gsl_vector *x position in parameter space
 * where the returned values in gsl_vector *f are zero.
 * The functions on which we search for roots are:
 * dH/dr, dH/dptheta and dH/dpphi-omega,
 * namely, Eqs. 4.8 and 4.9 of Buonanno et al. PRD 74, 104005 (2006).
 */
static int
XLALFindSphericalOrbitPrec(
          const gsl_vector * x,	/**<< Parameters requested by gsl root finder */
		      void *params,	        /**<< Spin EOB parameters */
		      gsl_vector * f	      /**<< Function values for the given parameters*/
)
{
	int		debugPK = 0;
	SEOBRootParams *rootParams = (SEOBRootParams *) params;
	REAL8		mTotal = rootParams->params->eobParams->m1 + rootParams->params->eobParams->m2;
	if (debugPK)
		XLAL_PRINT_INFO("\n\nmtotal in XLALFindSphericalOrbitPrec = %f\n", (float)mTotal);

	REAL8		py      , pz, r, ptheta, pphi;

	/* Numerical derivative of Hamiltonian wrt given value */
	REAL8		dHdx    , dHdpy, dHdpz;
	REAL8		dHdr    , dHdptheta, dHdpphi;

	/* Populate the appropriate values */
	/* In the special theta=pi/2 phi=0 case, r is x */
	rootParams->values[0] = r =  gsl_vector_get(x, 0)/scale1;
	rootParams->values[4] = py = gsl_vector_get(x, 1)/scale2;
	rootParams->values[5] = pz = gsl_vector_get(x, 2)/scale3;

    if(isnan(rootParams->values[0])) {
        rootParams->values[0] = 100.;
    }
    if(isnan(rootParams->values[4])) {
        rootParams->values[4] = 0.1;
    }
    if(isnan(rootParams->values[5])) {
        rootParams->values[5] = 0.01;
    }
	if (debugPK) {
    XLAL_PRINT_INFO("%3.10f %3.10f %3.10f %3.10f %3.10f %3.10f\n",
      rootParams->values[0], rootParams->values[1], rootParams->values[2],
      rootParams->values[3], rootParams->values[4], rootParams->values[5]);
    XLAL_PRINT_INFO("Values r = %.16e, py = %.16e, pz = %.16e\n", r, py, pz);
    fflush(NULL);
  }

	ptheta = -r * pz;
	pphi = r * py;

	if (debugPK)
		XLAL_PRINT_INFO("Input Values r = %.16e, py = %.16e, pz = %.16e\n pthetha = %.16e pphi = %.16e\n", r, py, pz, ptheta, pphi);

	/* dH by dR and dP */
	REAL8		tmpDValues[14];
	int UNUSED	status;
	for (int i = 0; i < 3; i++) {
		rootParams->values[i + 6] /= mTotal * mTotal;
		rootParams->values[i + 9] /= mTotal * mTotal;
	}
    UINT4 oldignoreflux = rootParams->params->ignoreflux;
    rootParams->params->ignoreflux = 1;
	status = XLALSpinPrecHcapNumericalDerivative(0, rootParams->values, tmpDValues, rootParams->params);
    rootParams->params->ignoreflux = oldignoreflux;
	for (int i = 0; i < 3; i++) {
		rootParams->values[i + 6] *= mTotal * mTotal;
		rootParams->values[i + 9] *= mTotal * mTotal;
	}
	REAL8	rvec[3] =
        {rootParams->values[0], rootParams->values[1], rootParams->values[2]};
	REAL8	pvec[3] =
        {rootParams->values[3], rootParams->values[4], rootParams->values[5]};
	REAL8	chi1vec[3] =
        {rootParams->values[6], rootParams->values[7], rootParams->values[8]};
	REAL8	chi2vec[3] =
        {rootParams->values[9], rootParams->values[10], rootParams->values[11]};
	REAL8	Lvec[3] = {CalculateCrossProductPrec(0, rvec, pvec),
    CalculateCrossProductPrec(1, rvec, pvec), CalculateCrossProductPrec(2, rvec, pvec)};
	REAL8	theta1 = acos(CalculateDotProductPrec(chi1vec, Lvec)
                / sqrt(CalculateDotProductPrec(chi1vec, chi1vec))
                / sqrt(CalculateDotProductPrec(Lvec, Lvec)));
	REAL8	theta2 = acos(CalculateDotProductPrec(chi2vec, Lvec)
                / sqrt(CalculateDotProductPrec(chi2vec, chi2vec))
                / sqrt(CalculateDotProductPrec(Lvec, Lvec)));

  if (debugPK) {
		XLAL_PRINT_INFO("rvec = %.16e %.16e %.16e\n", rvec[0], rvec[1], rvec[2]);
		XLAL_PRINT_INFO("pvec = %.16e %.16e %.16e\n", pvec[0], pvec[1], pvec[2]);
		XLAL_PRINT_INFO("theta1 = %.16e\n", theta1);
		XLAL_PRINT_INFO("theta2 = %.16e\n", theta2);
	}

  if (theta1 > 1.0e-6 && theta2 >= 1.0e-6) {
		dHdx = -tmpDValues[3];
		dHdpy = tmpDValues[1];
		dHdpz = tmpDValues[2];
	} else {
		//rootParams->values[5] = 0.;
		//rootParams->values[6] = 0.;
		//rootParams->values[7] = 0.;
		//rootParams->values[8] = sqrt(CalculateDotProductPrec(chi1vec, chi1vec));
    //rootParams->values[9] = 0.;
		//rootParams->values[10]= 0.;
		//rootParams->values[11]= sqrt(CalculateDotProductPrec(chi2vec, chi2vec));

		dHdx = XLALSpinPrecHcapNumDerivWRTParam(0,
                  rootParams->values, rootParams->params);
		dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4,
                  rootParams->values, rootParams->params);
		dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5,
                  rootParams->values, rootParams->params);
	}

	if (XLAL_IS_REAL8_FAIL_NAN(dHdx)) { XLAL_ERROR(XLAL_EDOM); }
	if (XLAL_IS_REAL8_FAIL_NAN(dHdpy)) { XLAL_ERROR(XLAL_EDOM); }
	if (XLAL_IS_REAL8_FAIL_NAN(dHdpz)) { XLAL_ERROR(XLAL_EDOM); }
	if (debugPK)
		XLAL_PRINT_INFO("dHdx = %.16e, dHdpy = %.16e, dHdpz = %.16e\n", dHdx, dHdpy, dHdpz);

	/* Now convert to spherical polars */
	dHdr      = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
	dHdptheta = -dHdpz / r;
	dHdpphi   = dHdpy / r;

	if (debugPK)
		XLAL_PRINT_INFO("dHdr = %.16e dHdptheta = %.16e dHdpphi = %.16e\n",
            dHdr, dHdptheta, dHdpphi);
	/* populate the function vector */
	gsl_vector_set(f, 0, dHdr);
	gsl_vector_set(f, 1, dHdptheta);
	gsl_vector_set(f, 2, dHdpphi - rootParams->omega);

	if (debugPK)
		XLAL_PRINT_INFO("Current funcvals = %.16e %.16e %.16e\n",
        gsl_vector_get(f, 0), gsl_vector_get(f, 1), gsl_vector_get(f, 2) );

  /* Rescale back */
  rootParams->values[0] *= scale1;
  rootParams->values[4] *= scale2;
  rootParams->values[5] *= scale3;

	return XLAL_SUCCESS;
}


/**
 * Wrapper for calculating specific derivative of the Hamiltonian in spherical co-ordinates,
 * dH/dr, dH/dptheta and dH/dpphi.
 * It only works for the specific co-ord system we use here
 */
static double
GSLSpinHamiltonianDerivWrapperPrec(double x,	/**<< Derivative at x */
			       void *params /**<< Function parameters */ )
{
	int		debugPK = 0;
	HcapSphDeriv2Params *dParams = (HcapSphDeriv2Params *) params;
	REAL8		mTotal = dParams->params->eobParams->m1 + dParams->params->eobParams->m2;

	REAL8		sphValues[12];
	REAL8		cartValues[12];

	REAL8		dHdr    , dHdx, dHdpy, dHdpz;
	REAL8		r       , ptheta, pphi;

	memcpy(sphValues, dParams->sphValues, sizeof(sphValues));
	sphValues[dParams->varyParam1] = x;

	SphericalToCartesianPrec(cartValues, cartValues + 3, sphValues, sphValues + 3);
	memcpy(cartValues + 6, sphValues + 6, 6 * sizeof(REAL8));

	r = sphValues[0];
	ptheta = sphValues[4];
	pphi = sphValues[5];

	/* New code to compute derivatives w.r.t. cartesian variables */
	REAL8		tmpDValues[14];
	int UNUSED	status;
	for (int i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT4 oldignoreflux = dParams->params->ignoreflux;
    dParams->params->ignoreflux = 1;
	status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, dParams->params);
    dParams->params->ignoreflux = oldignoreflux;
	for (int i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	double		rvec    [3] = {cartValues[0], cartValues[1], cartValues[2]};
	double		pvec    [3] = {cartValues[3], cartValues[4], cartValues[5]};
	double		chi1vec [3] = {cartValues[6], cartValues[7], cartValues[8]};
	double		chi2vec [3] = {cartValues[9], cartValues[10], cartValues[11]};
	double		Lvec    [3] = {CalculateCrossProductPrec(0, rvec, pvec), CalculateCrossProductPrec(1, rvec, pvec), CalculateCrossProductPrec(2, rvec, pvec)};
	double		theta1 = acos(CalculateDotProductPrec(chi1vec, Lvec) / sqrt(CalculateDotProductPrec(chi1vec, chi1vec)) / sqrt(CalculateDotProductPrec(Lvec, Lvec)));
	double		theta2 = acos(CalculateDotProductPrec(chi2vec, Lvec) / sqrt(CalculateDotProductPrec(chi2vec, chi2vec)) / sqrt(CalculateDotProductPrec(Lvec, Lvec)));
	if (debugPK) {
		XLAL_PRINT_INFO("In GSLSpinHamiltonianDerivWrapperPrec");
		XLAL_PRINT_INFO("rvec = %.16e %.16e %.16e\n", rvec[0], rvec[1], rvec[2]);
		XLAL_PRINT_INFO("pvec = %.16e %.16e %.16e\n", pvec[0], pvec[1], pvec[2]);
		XLAL_PRINT_INFO("theta1 = %.16e\n", theta1);
		XLAL_PRINT_INFO("theta2 = %.16e\n", theta2);
	}
	if (theta1 > 1.0e-6 && theta2 >= 1.0e-6) {
		switch (dParams->varyParam2) {
		case 0:
			/* dHdr */
			dHdx = -tmpDValues[3];
			//XLALSpinPrecHcapNumDerivWRTParam(0, cartValues, dParams->params);
			dHdpy = tmpDValues[1];
			//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			dHdpz = tmpDValues[2];
			//XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);

			dHdr = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
			//XLAL_PRINT_INFO("dHdr = %.16e\n", dHdr);
			return dHdr;

			break;
		case 4:
			/* dHdptheta */
			dHdpz = tmpDValues[2];
			//XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);
			return -dHdpz / r;
			break;
		case 5:
			/* dHdpphi */
			dHdpy = tmpDValues[1];
			//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			return dHdpy / r;
			break;
		default:
			XLALPrintError("This option is not supported in the second derivative function!\n");
			XLAL_ERROR_REAL8(XLAL_EINVAL);
			break;
		}
	} else {
		switch (dParams->varyParam2) {
		case 0:
			/* dHdr */
			dHdx = XLALSpinPrecHcapNumDerivWRTParam(0, cartValues, dParams->params);
			dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);

			dHdr = dHdx - dHdpy * pphi / (r * r) + dHdpz * ptheta / (r * r);
			//XLAL_PRINT_INFO("dHdr = %.16e\n", dHdr);
			return dHdr;

			break;
		case 4:
			/* dHdptheta */
			dHdpz = XLALSpinPrecHcapNumDerivWRTParam(5, cartValues, dParams->params);
			return -dHdpz / r;
			break;
		case 5:
			/* dHdpphi */
			dHdpy = XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, dParams->params);
			return dHdpy / r;
			break;
		default:
			XLALPrintError("This option is not supported in the second derivative function!\n");
			XLAL_ERROR_REAL8(XLAL_EINVAL);
			break;
		}
	}
}


/**
 * Function to calculate the second derivative of the Hamiltonian.
* The derivatives are taken with respect to indices idx1, idx2    */
static REAL8
XLALCalculateSphHamiltonianDeriv2Prec(
				  const int idx1,	/**<< Derivative w.r.t. index 1 */
				  const int idx2,	/**<< Derivative w.r.t. index 2 */
				  const REAL8 values[],	/**<< Dynamical variables in spherical coordinates */
				  SpinEOBParams * params	/**<< Spin EOB Parameters */
)
{

	static const REAL8 STEP_SIZE = 1.0e-5;

	REAL8		result;
	REAL8 UNUSED	absErr;

	HcapSphDeriv2Params dParams;

	gsl_function	F;
	INT4 UNUSED	gslStatus;

	dParams.sphValues = values;
	dParams.varyParam1 = idx1;
	dParams.varyParam2 = idx2;
	dParams.params = params;

	/*
	 * XLAL_PRINT_INFO( " In second deriv function: values\n" ); for ( int i = 0;
	 * i < 12; i++ ) { XLAL_PRINT_INFO( "%.16e ", values[i] ); } XLAL_PRINT_INFO( "\n" );
	 */
	F.function = GSLSpinHamiltonianDerivWrapperPrec;
	F.params = &dParams;

	/* GSL seemed to give weird answers - try my own code */
	/*
	 * result = GSLSpinHamiltonianDerivWrapperPrec( values[idx1] + STEP_SIZE,
	 * &dParams ) - GSLSpinHamiltonianDerivWrapperPrec( values[idx1] -
	 * STEP_SIZE, &dParams ); XLAL_PRINT_INFO( "%.16e - %.16e / 2h\n",
	 * GSLSpinHamiltonianDerivWrapperPrec( values[idx1] + STEP_SIZE, &dParams
	 * ), GSLSpinHamiltonianDerivWrapperPrec( values[idx1] - STEP_SIZE,
	 * &dParams ) );
	 *
	 * result = result / ( 2.*STEP_SIZE );
	 */

	XLAL_CALLGSL(gslStatus = gsl_deriv_central(&F, values[idx1],
					      STEP_SIZE, &result, &absErr));

	if (gslStatus != GSL_SUCCESS) {
		XLALPrintError("XLAL Error %s - Failure in GSL function\n", __func__);
		XLAL_ERROR_REAL8(XLAL_EDOM);
	}
	//XLAL_PRINT_INFO("Second deriv abs err = %.16e\n", absErr);

	//XLAL_PRINT_INFO("RESULT = %.16e\n", result);
	return result;
}

/**
 * Main function for calculating the spinning EOB initial conditions, following the
 * quasi-spherical initial conditions described in Sec. IV A of
 * Buonanno, Chen & Damour PRD 74, 104005 (2006).
 * All equation numbers in the comments of this function refer to this paper.
 *
 * Inputs are binary parameters (masses, spin vectors and inclination),
 * EOB model parameters and initial frequency.
 *
 * Outputs are initial dynamical variables in a vector
 * (x, y, z, px, py, pz, S1x, S1y, S1z, S2x, S2y, S2z).
 *
 * In the paper, the initial conditions are solved for a given initial radius,
 * while in this function, they are solved for a given inital orbital frequency.
 *
 * This difference is reflected in solving Eq. (4.8).
 * The calculation is done in 5 steps:
 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame, where
 * LNhat0 and N0 are initial normal to orbital plane and initial orbital separation;
 * STEP 2) After rotation in STEP 1, in spherical coordinates,
 * phi0 and theta0 are given directly in Eq. (4.7),
 * r0, pr0, ptheta0 and pphi0 are obtained by solving Eqs. (4.8) and (4.9)
 * (using gsl_multiroot_fsolver).
 * At this step, we find initial conditions for a spherical orbit without
 * radiation reaction.
 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
 * L0 is the initial orbital angular momentum and
 * L0 is calculated using initial position and linear momentum obtained in STEP 2.
 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq. (4.14),
 * then initial dr/dt using Eq. (4.10), and finally pr0 using Eq. (4.15).
 * STEP 5) Rotate back to the original inertial frame by inverting the rotation of STEP 3 and
 * then inverting the rotation of STEP 1.
 */

static int
XLALSimIMRSpinEOBInitialConditionsPrec(
				   REAL8Vector * initConds,	/**<< OUTPUT, Initial dynamical variables */
				   const REAL8 mass1,	/**<< mass 1 */
				   const REAL8 mass2,	/**<< mass 2 */
				   const REAL8 fMin,	/**<< Initial frequency (given) */
				   const REAL8 inc,	/**<< Inclination */
				   const REAL8 spin1[],	/**<< Initial spin vector 1 */
				   const REAL8 spin2[],	/**<< Initial spin vector 2 */
				   SpinEOBParams * params	/**<< Spin EOB parameters */
)
{

#ifndef LAL_NDEBUG
	if (!initConds) {
		XLAL_ERROR(XLAL_EINVAL);
	}
#endif

	int	debugPK = 0; int printPK = 0;
  FILE* UNUSED out = NULL;

	if (printPK) {
		XLAL_PRINT_INFO("Inside the XLALSimIMRSpinEOBInitialConditionsPrec function!\n");
		XLAL_PRINT_INFO(
    "Inputs: m1 = %.16e, m2 = %.16e, fMin = %.16e, inclination = %.16e\n",
      mass1, mass2, (double)fMin, (double)inc);
		XLAL_PRINT_INFO("Inputs: mSpin1 = {%.16e, %.16e, %.16e}\n",
      spin1[0], spin1[1], spin1[2]);
		XLAL_PRINT_INFO("Inputs: mSpin2 = {%.16e, %.16e, %.16e}\n",
      spin2[0], spin2[1], spin2[2]);
		fflush(NULL);
	}
	static const int UNUSED lMax = 8;

	int		i;

	/* Variable to keep track of whether the user requested the tortoise */
	int		tmpTortoise;

	UINT4		SpinAlignedEOBversion;

	REAL8		mTotal;
	REAL8		eta;
	REAL8		omega   , v0;	/* Initial velocity and angular
					 * frequency */

	REAL8		ham;	/* Hamiltonian */

	REAL8		LnHat    [3];	/* Initial orientation of angular
					 * momentum */
	REAL8		rHat     [3];	/* Initial orientation of radial
					 * vector */
	REAL8		vHat     [3];	/* Initial orientation of velocity
					 * vector */
	REAL8		Lhat     [3];	/* Direction of relativistic ang mom */
	REAL8		qHat     [3];
	REAL8		pHat     [3];

	/* q and p vectors in Cartesian and spherical coords */
	REAL8		qCart    [3], pCart[3];
	REAL8		qSph     [3], pSph[3];

	/* We will need to manipulate the spin vectors */
	/* We will use temporary vectors to do this */
	REAL8		tmpS1    [3];
	REAL8		tmpS2    [3];
	REAL8		tmpS1Norm[3];
	REAL8		tmpS2Norm[3];

	REAL8Vector	qCartVec, pCartVec;
	REAL8Vector	s1Vec, s2Vec, s1VecNorm, s2VecNorm;
	REAL8Vector	sKerr, sStar;
	REAL8		sKerrData[3], sStarData[3];
	REAL8		a = 0.;
	//, chiS, chiA;
	//REAL8 chi1, chi2;

	/*
	 * We will need a full values vector for calculating derivs of
	 * Hamiltonian
	 */
	REAL8		sphValues[12];
	REAL8		cartValues[12];

	/* Matrices for rotating to the new basis set. */
	/* It is more convenient to calculate the ICs in a simpler basis */
	gsl_matrix     *rotMatrix = NULL;
	gsl_matrix     *invMatrix = NULL;
	gsl_matrix     *rotMatrix2 = NULL;
	gsl_matrix     *invMatrix2 = NULL;

	/* Root finding stuff for finding the spherical orbit */
	SEOBRootParams	rootParams;
	const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid;
	gsl_multiroot_fsolver *rootSolver = NULL;

	gsl_multiroot_function rootFunction;
	gsl_vector     *initValues = NULL;
	gsl_vector     *finalValues = NULL;
	INT4 gslStatus;
        INT4 cntGslNoProgress = 0, MAXcntGslNoProgress = 5;
        //INT4 cntGslNoProgress = 0, MAXcntGslNoProgress = 50;
        REAL8 multFacGslNoProgress = 3./5.;
	//const int	maxIter = 2000;
	const int	maxIter = 10000;

	memset(&rootParams, 0, sizeof(rootParams));

	mTotal = mass1 + mass2;
	eta = mass1 * mass2 / (mTotal * mTotal);
	memcpy(tmpS1, spin1, sizeof(tmpS1));
	memcpy(tmpS2, spin2, sizeof(tmpS2));
	memcpy(tmpS1Norm, spin1, sizeof(tmpS1Norm));
	memcpy(tmpS2Norm, spin2, sizeof(tmpS2Norm));
	for (i = 0; i < 3; i++) {
		tmpS1Norm[i] /= mTotal * mTotal;
		tmpS2Norm[i] /= mTotal * mTotal;
	}
	SpinAlignedEOBversion = params->seobCoeffs->SpinAlignedEOBversion;
	/* We compute the ICs for the non-tortoise p, and convert at the end */
	tmpTortoise = params->tortoise;
	params->tortoise = 0;

	EOBNonQCCoeffs *nqcCoeffs = NULL;
	nqcCoeffs = params->nqcCoeffs;

	/*
	 * STEP 1) Rotate to LNhat0 along z-axis and N0 along x-axis frame,
	 * where LNhat0 and N0 are initial normal to orbital plane and
	 * initial orbital separation;
	 */

	/* Set the initial orbital ang mom direction. Taken from STPN code */
	LnHat[0] = sin(inc);
	LnHat[1] = 0.;
	LnHat[2] = cos(inc);

	/*
	 * Set the radial direction - need to take care to avoid singularity
	 * if L is along z axis
	 */
	if (LnHat[2] > 0.9999) {
		rHat[0] = 1.;
		rHat[1] = rHat[2] = 0.;
	} else {
		REAL8		theta0 = atan(-LnHat[2] / LnHat[0]);	/* theta0 is between 0
									 * and Pi */
		rHat[0] = sin(theta0);
		rHat[1] = 0;
		rHat[2] = cos(theta0);
	}

	/* Now we can complete the triad */
	vHat[0] = CalculateCrossProductPrec(0, LnHat, rHat);
	vHat[1] = CalculateCrossProductPrec(1, LnHat, rHat);
	vHat[2] = CalculateCrossProductPrec(2, LnHat, rHat);

	NormalizeVectorPrec(vHat);

	/* Vectors BEFORE rotation */
	if (printPK) {
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
        i, LnHat[i], i, rHat[i], i, vHat[i]);

		XLAL_PRINT_INFO("\n\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);
    fflush(NULL);
	}

	/* Allocate and compute the rotation matrices */
	XLAL_CALLGSL(rotMatrix = gsl_matrix_alloc(3, 3));
	XLAL_CALLGSL(invMatrix = gsl_matrix_alloc(3, 3));
	if (!rotMatrix || !invMatrix) {
		if (rotMatrix)
			gsl_matrix_free(rotMatrix);
		if (invMatrix)
			gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	if (CalculateRotationMatrixPrec(rotMatrix, invMatrix, rHat, vHat, LnHat) == XLAL_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	/* Rotate the orbital vectors and spins */
	ApplyRotationMatrixPrec(rotMatrix, rHat);
	ApplyRotationMatrixPrec(rotMatrix, vHat);
	ApplyRotationMatrixPrec(rotMatrix, LnHat);
	ApplyRotationMatrixPrec(rotMatrix, tmpS1);
	ApplyRotationMatrixPrec(rotMatrix, tmpS2);
	ApplyRotationMatrixPrec(rotMatrix, tmpS1Norm);
	ApplyRotationMatrixPrec(rotMatrix, tmpS2Norm);

	/* See if Vectors have been rotated fine */
	if (printPK) {
		XLAL_PRINT_INFO("\nAfter applying rotation matrix:\n\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" LnHat[%d] = %.16e, rHat[%d] = %.16e, vHat[%d] = %.16e\n",
                i, LnHat[i], i, rHat[i], i, vHat[i]);

		XLAL_PRINT_INFO("\n");
		for (i = 0; i < 3; i++)
			XLAL_PRINT_INFO(" s1[%d] = %.16e, s2[%d] = %.16e\n", i, tmpS1[i], i, tmpS2[i]);

    fflush(NULL);
	}
	/*
	 * STEP 2) After rotation in STEP 1, in spherical coordinates, phi0
	 * and theta0 are given directly in Eq. (4.7), r0, pr0, ptheta0 and
	 * pphi0 are obtained by solving Eqs. (4.8) and (4.9) (using
	 * gsl_multiroot_fsolver). At this step, we find initial conditions
	 * for a spherical orbit without radiation reaction.
	 */

  /* Initialise the gsl stuff */
	XLAL_CALLGSL(rootSolver = gsl_multiroot_fsolver_alloc(T, 3));
	if (!rootSolver) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	XLAL_CALLGSL(initValues = gsl_vector_calloc(3));
	if (!initValues) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}

	rootFunction.f = XLALFindSphericalOrbitPrec;
	rootFunction.n = 3;
	rootFunction.params = &rootParams;

	/* Calculate the initial velocity from the given initial frequency */
	omega = LAL_PI * mTotal * LAL_MTSUN_SI * fMin;
	v0 = cbrt(omega);

	/* Given this, we can start to calculate the initial conditions */
	/* for spherical coords in the new basis */
	rootParams.omega = omega;
	rootParams.params = params;

	/* To start with, we will just assign Newtonian-ish ICs to the system */
	rootParams.values[0] = scale1 * 1. / (v0 * v0);	/* Initial r */
	rootParams.values[4] = scale2 * v0;	            /* Initial p */
	rootParams.values[5] = scale3 * 1e-3;
	//PK
  memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
	memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

	if (printPK) {
    XLAL_PRINT_INFO("ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
      rootParams.values[0]/scale1, rootParams.values[4]/scale2,
      rootParams.values[5]/scale3);
    fflush(NULL);
  }

	gsl_vector_set(initValues, 0, rootParams.values[0]);
	gsl_vector_set(initValues, 1, rootParams.values[4]);
  gsl_vector_set(initValues, 2, rootParams.values[5]);

	gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);

	/* We are now ready to iterate to find the solution */
	i = 0;

  if(debugPK){ out = fopen("ICIterations.dat", "w"); }
	do {
		XLAL_CALLGSL(gslStatus = gsl_multiroot_fsolver_iterate(rootSolver));
		if (debugPK) {
      fprintf( out, "%d\t", i );

      /* Write to file */
      fprintf( out, "%.16e\t%.16e\t%.16e\t",
        rootParams.values[0]/scale1, rootParams.values[4]/scale2,
        rootParams.values[5]/scale3 );

      /* Residual Function values whose roots we are trying to find */
      finalValues = gsl_multiroot_fsolver_f(rootSolver);

      /* Write to file */
      fprintf( out, "%.16e\t%.16e\t%.16e\t",
        gsl_vector_get(finalValues, 0),
        gsl_vector_get(finalValues, 1),
        gsl_vector_get(finalValues, 2) );

      /* Step sizes in each of function variables */
      finalValues = gsl_multiroot_fsolver_dx(rootSolver);

      /* Write to file */
      fprintf( out, "%.16e\t%.16e\t%.16e\t%d\n",
        gsl_vector_get(finalValues, 0)/scale1,
        gsl_vector_get(finalValues, 1)/scale2,
        gsl_vector_get(finalValues, 2)/scale3,
        gslStatus );
		}

    if (gslStatus == GSL_ENOPROG || gslStatus == GSL_ENOPROGJ) {
      XLAL_PRINT_INFO(
        "\n NO PROGRESS being made by Spherical orbit root solver\n");

      /* Print Residual Function values whose roots we are trying to find */
      finalValues = gsl_multiroot_fsolver_f(rootSolver);
      XLAL_PRINT_INFO("Function value here given by the following:\n");
      XLAL_PRINT_INFO(" F1 = %.16e, F2 = %.16e, F3 = %.16e\n",
          gsl_vector_get(finalValues, 0),
		       gsl_vector_get(finalValues, 1), gsl_vector_get(finalValues, 2));

      /* Print Step sizes in each of function variables */
      finalValues = gsl_multiroot_fsolver_dx(rootSolver);
//      XLAL_PRINT_INFO("Stepsizes in each dimension:\n");
//      XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
//          gsl_vector_get(finalValues, 0)/scale1,
//	       gsl_vector_get(finalValues, 1)/scale2,
//          gsl_vector_get(finalValues, 2)/scale3);

      /* Only allow this flag to be caught MAXcntGslNoProgress no. of times */
      cntGslNoProgress += 1;
      if (cntGslNoProgress >= MAXcntGslNoProgress) {
        cntGslNoProgress = 0;

        if(multFacGslNoProgress < 1.){ multFacGslNoProgress *= 1.02; }
        else{ multFacGslNoProgress /= 1.01; }

      } 
      /* Now that no progress is being made, we need to reset the initial guess
       * for the (r,pPhi, pTheta) and reset the integrator */
      rootParams.values[0] = scale1 * 1. / (v0 * v0);	/* Initial r */
      rootParams.values[4] = scale2 * v0;	            /* Initial p */
      if( cntGslNoProgress % 2 )
        rootParams.values[5] = scale3 * 1e-3 / multFacGslNoProgress;
      else
        rootParams.values[5] = scale3 * 1e-3 * multFacGslNoProgress;
      //PK
      memcpy(rootParams.values + 6, tmpS1, sizeof(tmpS1));
      memcpy(rootParams.values + 9, tmpS2, sizeof(tmpS2));

      if (printPK) {
        XLAL_PRINT_INFO("New ICs guess: x = %.16e, py = %.16e, pz = %.16e\n",
                rootParams.values[0]/scale1, rootParams.values[4]/scale2,
                rootParams.values[5]/scale3);
        fflush(NULL);
      }

      gsl_vector_set(initValues, 0, rootParams.values[0]);
      gsl_vector_set(initValues, 1, rootParams.values[4]);
      gsl_vector_set(initValues, 2, rootParams.values[5]);
      gsl_multiroot_fsolver_set(rootSolver, &rootFunction, initValues);
    }
    else if (gslStatus == GSL_EBADFUNC) {
      XLALPrintError(
      "Inf or Nan encountered in evaluluation of spherical orbit Eqn\n");
			gsl_multiroot_fsolver_free(rootSolver);
			gsl_vector_free(initValues);
			gsl_matrix_free(rotMatrix);
			gsl_matrix_free(invMatrix);
			XLAL_ERROR(XLAL_EDOM);
    }
		else if (gslStatus != GSL_SUCCESS) {
			XLALPrintError("Error in GSL iteration function!\n");
			gsl_multiroot_fsolver_free(rootSolver);
			gsl_vector_free(initValues);
			gsl_matrix_free(rotMatrix);
			gsl_matrix_free(invMatrix);
			XLAL_ERROR(XLAL_EDOM);
		}

    /* different ways to test convergence of the method */
		XLAL_CALLGSL(gslStatus = gsl_multiroot_test_residual(rootSolver->f, 1.0e-8));
    /*XLAL_CALLGSL(gslStatus= gsl_multiroot_test_delta(
          gsl_multiroot_fsolver_dx(rootSolver),
          gsl_multiroot_fsolver_root(rootSolver),
          1.e-8, 1.e-5));*/
		i++;
	}
	while (gslStatus == GSL_CONTINUE && i <= maxIter);

  if(debugPK) { fflush(NULL); fclose(out); }

	if (i > maxIter && gslStatus != GSL_SUCCESS) {
		gsl_multiroot_fsolver_free(rootSolver);
		gsl_vector_free(initValues);
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		//XLAL_ERROR(XLAL_EMAXITER);
		XLAL_ERROR(XLAL_EDOM);
	}
	finalValues = gsl_multiroot_fsolver_root(rootSolver);

	if (printPK) {
		XLAL_PRINT_INFO("Spherical orbit conditions here given by the following:\n");
		XLAL_PRINT_INFO(" x = %.16e, py = %.16e, pz = %.16e\n",
           gsl_vector_get(finalValues, 0)/scale1,
		       gsl_vector_get(finalValues, 1)/scale2,
           gsl_vector_get(finalValues, 2)/scale3);
	}
	memset(qCart, 0, sizeof(qCart));
	memset(pCart, 0, sizeof(pCart));

	qCart[0] = gsl_vector_get(finalValues, 0)/scale1;
	pCart[1] = gsl_vector_get(finalValues, 1)/scale2;
	pCart[2] = gsl_vector_get(finalValues, 2)/scale3;


	/* Free the GSL root finder, since we're done with it */
	gsl_multiroot_fsolver_free(rootSolver);
	gsl_vector_free(initValues);


	/*
	 * STEP 3) Rotate to L0 along z-axis and N0 along x-axis frame, where
	 * L0 is the initial orbital angular momentum and L0 is calculated
	 * using initial position and linear momentum obtained in STEP 2.
	 */

	/* Now we can calculate the relativistic L and rotate to a new basis */
	memcpy(qHat, qCart, sizeof(qCart));
	memcpy(pHat, pCart, sizeof(pCart));

	NormalizeVectorPrec(qHat);
	NormalizeVectorPrec(pHat);

	Lhat[0] = CalculateCrossProductPrec(0, qHat, pHat);
	Lhat[1] = CalculateCrossProductPrec(1, qHat, pHat);
	Lhat[2] = CalculateCrossProductPrec(2, qHat, pHat);

	NormalizeVectorPrec(Lhat);

	XLAL_CALLGSL(rotMatrix2 = gsl_matrix_alloc(3, 3));
	XLAL_CALLGSL(invMatrix2 = gsl_matrix_alloc(3, 3));

	if (CalculateRotationMatrixPrec(rotMatrix2, invMatrix2, qHat, pHat, Lhat) == XLAL_FAILURE) {
		gsl_matrix_free(rotMatrix);
		gsl_matrix_free(invMatrix);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	ApplyRotationMatrixPrec(rotMatrix2, rHat);
	ApplyRotationMatrixPrec(rotMatrix2, vHat);
	ApplyRotationMatrixPrec(rotMatrix2, LnHat);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS1);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS2);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS1Norm);
	ApplyRotationMatrixPrec(rotMatrix2, tmpS2Norm);
	ApplyRotationMatrixPrec(rotMatrix2, qCart);
	ApplyRotationMatrixPrec(rotMatrix2, pCart);

        gsl_matrix_free(rotMatrix);
        gsl_matrix_free(rotMatrix2);

        if (printPK) {
		XLAL_PRINT_INFO("qCart after rotation2 %3.10f %3.10f %3.10f\n", qCart[0], qCart[1], qCart[2]);
		XLAL_PRINT_INFO("pCart after rotation2 %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
		XLAL_PRINT_INFO("S1 after rotation2 %3.10f %3.10f %3.10f\n", tmpS1Norm[0], tmpS1Norm[1], tmpS1Norm[2]);
		XLAL_PRINT_INFO("S2 after rotation2 %3.10f %3.10f %3.10f\n", tmpS2Norm[0], tmpS2Norm[1], tmpS2Norm[2]);
	}
	/*
	 * STEP 4) In the L0-N0 frame, we calculate (dE/dr)|sph using Eq.
	 * (4.14), then initial dr/dt using Eq. (4.10), and finally pr0 using
	 * Eq. (4.15).
	 */

	/* Now we can calculate the flux. Change to spherical co-ords */
	CartesianToSphericalPrec(qSph, pSph, qCart, pCart);
	memcpy(sphValues, qSph, sizeof(qSph));
	memcpy(sphValues + 3, pSph, sizeof(pSph));
	memcpy(sphValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(sphValues + 9, tmpS2, sizeof(tmpS2));

	memcpy(cartValues, qCart, sizeof(qCart));
	memcpy(cartValues + 3, pCart, sizeof(pCart));
	memcpy(cartValues + 6, tmpS1, sizeof(tmpS1));
	memcpy(cartValues + 9, tmpS2, sizeof(tmpS2));

	REAL8		dHdpphi , d2Hdr2, d2Hdrdpphi;
	REAL8		rDot    , dHdpr, flux, dEdr;

	d2Hdr2 = XLALCalculateSphHamiltonianDeriv2Prec(0, 0, sphValues, params);
	d2Hdrdpphi = XLALCalculateSphHamiltonianDeriv2Prec(0, 5, sphValues, params);

	if (printPK)
		XLAL_PRINT_INFO("d2Hdr2 = %.16e, d2Hdrdpphi = %.16e\n", d2Hdr2, d2Hdrdpphi);

	/* New code to compute derivatives w.r.t. cartesian variables */

	REAL8		tmpDValues[14];
	int UNUSED	status;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] /= mTotal * mTotal;
		cartValues[i + 9] /= mTotal * mTotal;
	}
    UINT4 oldignoreflux = params->ignoreflux;
    params->ignoreflux = 1;
	status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
    params->ignoreflux = oldignoreflux;
	for (i = 0; i < 3; i++) {
		cartValues[i + 6] *= mTotal * mTotal;
		cartValues[i + 9] *= mTotal * mTotal;
	}

	dHdpphi = tmpDValues[1] / sqrt(cartValues[0] * cartValues[0] + cartValues[1] * cartValues[1] + cartValues[2] * cartValues[2]);
	//XLALSpinPrecHcapNumDerivWRTParam(4, cartValues, params) / sphValues[0];

	dEdr = -dHdpphi * d2Hdr2 / d2Hdrdpphi;

	if (printPK)
		XLAL_PRINT_INFO("d2Hdr2 = %.16e d2Hdrdpphi = %.16e dHdpphi = %.16e\n",
            d2Hdr2, d2Hdrdpphi, dHdpphi);

	if (d2Hdr2 != 0.0) {
		/* We will need to calculate the Hamiltonian to get the flux */
		s1Vec.length = s2Vec.length = s1VecNorm.length = s2VecNorm.length = sKerr.length = sStar.length = 3;
		s1Vec.data = tmpS1;
		s2Vec.data = tmpS2;
		s1VecNorm.data = tmpS1Norm;
		s2VecNorm.data = tmpS2Norm;
		sKerr.data = sKerrData;
		sStar.data = sStarData;

		qCartVec.length = pCartVec.length = 3;
		qCartVec.data = qCart;
		pCartVec.data = pCart;

		//chi1 = tmpS1[0] * LnHat[0] + tmpS1[1] * LnHat[1] + tmpS1[2] * LnHat[2];
		//chi2 = tmpS2[0] * LnHat[0] + tmpS2[1] * LnHat[1] + tmpS2[2] * LnHat[2];

		//if (debugPK)
			//XLAL_PRINT_INFO("magS1 = %.16e, magS2 = %.16e\n", chi1, chi2);

		//chiS = 0.5 * (chi1 / (mass1 * mass1) + chi2 / (mass2 * mass2));
		//chiA = 0.5 * (chi1 / (mass1 * mass1) - chi2 / (mass2 * mass2));

		XLALSimIMRSpinEOBCalculateSigmaKerr(&sKerr, mass1, mass2, &s1Vec, &s2Vec);
		XLALSimIMRSpinEOBCalculateSigmaStar(&sStar, mass1, mass2, &s1Vec, &s2Vec);

		/*
		 * The a in the flux has been set to zero, but not in the
		 * Hamiltonian
		 */
		a = sqrt(sKerr.data[0] * sKerr.data[0] + sKerr.data[1] * sKerr.data[1] + sKerr.data[2] * sKerr.data[2]);
		//XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(params->eobParams->hCoeffs, mass1, mass2, eta, /* a */ 0.0, chiS, chiA);
		//XLALSimIMRCalculateSpinPrecEOBHCoeffs(params->seobCoeffs, eta, a);
		ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &qCartVec, &pCartVec, &s1VecNorm, &s2VecNorm, &sKerr, &sStar, params->tortoise, params->seobCoeffs);

		if (printPK)
			XLAL_PRINT_INFO("Stas: hamiltonian in ICs at this point is %.16e\n", ham);

		/* And now, finally, the flux */
		REAL8Vector	polarDynamics, cartDynamics;
		REAL8		polarData[4], cartData[12];

		polarDynamics.length = 4;
		polarDynamics.data = polarData;

		polarData[0] = qSph[0];
		polarData[1] = 0.;
		polarData[2] = pSph[0];
		polarData[3] = pSph[2];

		cartDynamics.length = 12;
		cartDynamics.data = cartData;

		memcpy(cartData, qCart, 3 * sizeof(REAL8));
		memcpy(cartData + 3, pCart, 3 * sizeof(REAL8));
		memcpy(cartData + 6, tmpS1Norm, 3 * sizeof(REAL8));
		memcpy(cartData + 9, tmpS2Norm, 3 * sizeof(REAL8));

		//XLAL_PRINT_INFO("Stas: starting FLux calculations\n");

		flux = XLALInspiralPrecSpinFactorizedFlux(&polarDynamics, &cartDynamics, nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion);
		/*
		 * flux  = XLALInspiralSpinFactorizedFlux( &polarDynamics,
		 * nqcCoeffs, omega, params, ham, lMax, SpinAlignedEOBversion
		 * );
		 */
		//XLAL_PRINT_INFO("Stas flux = %.16e \n", flux);
		//exit(0);
		flux = flux / eta;

		rDot = -flux / dEdr;
		if (debugPK) {
			XLAL_PRINT_INFO("Stas here I am 2  \n");
		}
		/*
		 * We now need dHdpr - we take it that it is safely linear up
		 * to a pr of 1.0e-3 PK: Ideally, the pr should be of the
		 * order of other momenta, in order for its contribution to
		 * the Hamiltonian to not get buried in the numerical noise
		 * in the numerically larger momenta components
		 */
		cartValues[3] = 1.0e-3;
		for (i = 0; i < 3; i++) {
			cartValues[i + 6] /= mTotal * mTotal;
			cartValues[i + 9] /= mTotal * mTotal;
		}
        oldignoreflux = params->ignoreflux;
        params->ignoreflux = 1;
        params->seobCoeffs->updateHCoeffs = 1;
		status = XLALSpinPrecHcapNumericalDerivative(0, cartValues, tmpDValues, params);
        params->ignoreflux = oldignoreflux;
		for (i = 0; i < 3; i++) {
			cartValues[i + 6] *= mTotal * mTotal;
			cartValues[i + 9] *= mTotal * mTotal;
		}
        REAL8		csi = sqrt(XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, qSph[0], eta, a)*XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, qSph[0], eta, a)) / (qSph[0] * qSph[0] + a * a);

		dHdpr = csi*tmpDValues[0];
		//XLALSpinPrecHcapNumDerivWRTParam(3, cartValues, params);

		if (debugPK) {
			XLAL_PRINT_INFO("Ingredients going into prDot:\n");
			XLAL_PRINT_INFO("flux = %.16e, dEdr = %.16e, dHdpr = %.16e, dHdpr/pr = %.16e\n", flux, dEdr, dHdpr, dHdpr / cartValues[3]);
		}
		/*
		 * We can now calculate what pr should be taking into account
		 * the flux
		 */
		pSph[0] = rDot / (dHdpr / cartValues[3]);
	} else {
		/*
		 * Since d2Hdr2 has evaluated to zero, we cannot do the
		 * above. Just set pr to zero
		 */
		//XLAL_PRINT_INFO("d2Hdr2 is zero!\n");
		pSph[0] = 0;
	}

	/* Now we are done - convert back to cartesian coordinates ) */
	SphericalToCartesianPrec(qCart, pCart, qSph, pSph);

	/*
	 * STEP 5) Rotate back to the original inertial frame by inverting
	 * the rotation of STEP 3 and then  inverting the rotation of STEP 1.
	 */

	/* Undo rotations to get back to the original basis */
	/* Second rotation */
	ApplyRotationMatrixPrec(invMatrix2, rHat);
	ApplyRotationMatrixPrec(invMatrix2, vHat);
	ApplyRotationMatrixPrec(invMatrix2, LnHat);
	ApplyRotationMatrixPrec(invMatrix2, tmpS1);
	ApplyRotationMatrixPrec(invMatrix2, tmpS2);
	ApplyRotationMatrixPrec(invMatrix2, tmpS1Norm);
	ApplyRotationMatrixPrec(invMatrix2, tmpS2Norm);
	ApplyRotationMatrixPrec(invMatrix2, qCart);
	ApplyRotationMatrixPrec(invMatrix2, pCart);

	/* First rotation */
	ApplyRotationMatrixPrec(invMatrix, rHat);
	ApplyRotationMatrixPrec(invMatrix, vHat);
	ApplyRotationMatrixPrec(invMatrix, LnHat);
	ApplyRotationMatrixPrec(invMatrix, tmpS1);
	ApplyRotationMatrixPrec(invMatrix, tmpS2);
	ApplyRotationMatrixPrec(invMatrix, tmpS1Norm);
	ApplyRotationMatrixPrec(invMatrix, tmpS2Norm);
	ApplyRotationMatrixPrec(invMatrix, qCart);
	ApplyRotationMatrixPrec(invMatrix, pCart);

        gsl_matrix_free(invMatrix);
        gsl_matrix_free(invMatrix2);

        /* If required, apply the tortoise transform */
	if (tmpTortoise) {
		REAL8		r = sqrt(qCart[0] * qCart[0] + qCart[1] * qCart[1] + qCart[2] * qCart[2]);
		REAL8		deltaR = XLALSimIMRSpinPrecEOBHamiltonianDeltaR(params->seobCoeffs, r, eta, a);
		REAL8		deltaT = XLALSimIMRSpinPrecEOBHamiltonianDeltaT(params->seobCoeffs, r, eta, a);
		REAL8		csi = sqrt(deltaT * deltaR) / (r * r + a * a);

		REAL8		pr = (qCart[0] * pCart[0] + qCart[1] * pCart[1] + qCart[2] * pCart[2]) / r;

		params->tortoise = tmpTortoise;

		if (debugPK) {
			XLAL_PRINT_INFO("Applying the tortoise to p (csi = %.26e)\n", csi);
			XLAL_PRINT_INFO("pCart = %3.10f %3.10f %3.10f\n", pCart[0], pCart[1], pCart[2]);
		}
		for (i = 0; i < 3; i++) {
			pCart[i] = pCart[i] + qCart[i] * pr * (csi - 1.) / r;
		}
	}


    /* Now copy the initial conditions back to the return vector */
	memcpy(initConds->data, qCart, sizeof(qCart));
	memcpy(initConds->data + 3, pCart, sizeof(pCart));
	memcpy(initConds->data + 6, tmpS1Norm, sizeof(tmpS1Norm));
	memcpy(initConds->data + 9, tmpS2Norm, sizeof(tmpS2Norm));

    for (i=0; i<12; i++) {
        if (fabs(initConds->data[i]) <=1.0e-15) {
            initConds->data[i] = 0.;
        }
    }

	if (debugPK) {
		XLAL_PRINT_INFO("THE FINAL INITIAL CONDITIONS:\n");
		XLAL_PRINT_INFO(" %.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n%.16e %.16e %.16e\n", initConds->data[0], initConds->data[1], initConds->data[2],
		       initConds->data[3], initConds->data[4], initConds->data[5], initConds->data[6], initConds->data[7], initConds->data[8],
		       initConds->data[9], initConds->data[10], initConds->data[11]);
	}
	return XLAL_SUCCESS;
}

#endif				/* _LALSIMIMRSPINPRECEOBINITIALCONDITIONS_C */
