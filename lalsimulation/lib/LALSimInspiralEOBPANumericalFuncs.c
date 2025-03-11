#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimSphHarmMode.h>
#include <LALSimInspiralWaveformFlags.h>
#include <lal/LALDict.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#include "LALSimInspiralEOBPostAdiabatic.h"

/**
 * Function which calculates the 4-order first finite difference
 * derivative of a numerical function
 */
int
XLALFDDerivative1Order4(
	REAL8Vector *XVec,
    /**<< An array of X values */
	REAL8Vector *YVec,
    /**<< An array of Y values */
	REAL8Vector *derivativeVec
    /**<< OUTPUT, the derivative dY/dX */
)
{
	REAL8 fourthOrderCoeffs[5][5] = {
		{-25./12., 4., -3., 4./3., -1./4.},
		{-1./4., -5./6., 3./2., -1./2., 1./12.},
		{1./12., -2./3., 0, 2./3., -1./12.},
		{-1./12., 1./2., -3./2., 5./6., 1./4.},
		{1./4., -4./3., 3., -4., 25./12.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += (
                    fourthOrderCoeffs[0][j] * YVec->data[j]
                );
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += (
                    fourthOrderCoeffs[1][j] * YVec->data[j]
                );
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += (
                    fourthOrderCoeffs[3][j] * YVec->data[i+j-3]
                );
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += (
                    fourthOrderCoeffs[4][j] * YVec->data[i+j-4]
                );
			}
		}
		else
		{
			for (j = 0; j <= 4; j++)
			{
				derivativeVec->data[i] += (
                    fourthOrderCoeffs[2][j] * YVec->data[i+j-2]
                );
			}
		}

		derivativeVec->data[i] /= h;
	}

    return XLAL_SUCCESS;
}

/**
 * Function which calculates the 2-order first finite difference
 * derivative of a numerical function
 */
int
XLALFDDerivative1Order2(
	REAL8Vector *XVec,
    /**<< An array of X values */
	REAL8Vector *YVec,
    /**<< An array of Y values */
	REAL8Vector *derivativeVec
    /**<< OUTPUT, the derivative dY/dX */
)
{
	REAL8 secondOrderCoeffs[3][3] = {
		{-3./2., 2, -1./2.},
		{-1./2., 0, 1./2.},
		{1./2., -2., 3./2.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;


	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += (
                    secondOrderCoeffs[0][j] * YVec->data[j]
                );
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += (
                    secondOrderCoeffs[2][j] * YVec->data[i+j-2]
                );
			}
		}
		else
		{
			for (j = 0; j <= 2; j++)
			{
				derivativeVec->data[i] += (
                    secondOrderCoeffs[1][j] * YVec->data[i+j-1]
                );
			}
		}

		derivativeVec->data[i] /= h;
	}

    return XLAL_SUCCESS;
}

/**
 * Function which calculates the 6-order first finite difference
 * derivative of a numerical function
 */
int
XLALFDDerivative1Order6(
	REAL8Vector *XVec,
    /**<< An array of X values */
	REAL8Vector *YVec,
    /**<< An array of Y values */
	REAL8Vector *derivativeVec
    /**<< OUTPUT, the derivative dY/dX */
)
{
	REAL8 sixthOrderCoeffs[7][7] = {
		{-49./20., 6., -15./2., 20./3., -15./4., 6./5., -1./6.},
		{-1./6., -77./60., 5./2., -5./3., 5./6., -1./4., 1./30.},
		{1./30., -2./5., -7./12., 4./3., -1./2., 2./15., -1./60.},
		{-1./60., 3./20., -3./4., 0, 3./4., -3./20., 1./60.},
		{1./60., -2./15., 1./2., -4./3., 7./12., 2./5., -1./30.},
		{-1./30., 1./4., -5./6., 5./3., -5./2., 77./60., 1./6.},
		{1./6., -6./5., 15./4., -20./3., 15./2., -6., 49./20.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;


	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[0][j] * YVec->data[j]
                );
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[1][j] * YVec->data[j]
                );
			}
		}
		else if (i == 2)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[2][j] * YVec->data[j]
                );
			}
		}
		else if (i == vecLength-3)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[4][j] * YVec->data[i+j-4]
                );
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[5][j] * YVec->data[i+j-5]
                );
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[6][j] * YVec->data[i+j-6]
                );
			}
		}
		else
		{
			for (j = 0; j <= 6; j++)
			{
				derivativeVec->data[i] += (
                    sixthOrderCoeffs[3][j] * YVec->data[i+j-3]
                );
			}
		}

		derivativeVec->data[i] /= h;
	}

    return XLAL_SUCCESS;
}

/**
 * Function which calculates the 8-order first finite difference
 * derivative of a numerical function
 */
int
XLALFDDerivative1Order8(
	REAL8Vector *XVec,
    /**<< An array of X values */
	REAL8Vector *YVec,
    /**<< An array of Y values */
	REAL8Vector *derivativeVec
    /**<< OUTPUT, the derivative dY/dX */
)
{
	REAL8 eightOrderCoeffs[9][9] = {
		{-761./280., 8., -14., 56./3., -35./2., 56./5., -14./3., 8./7., -1./8.},
		{-1./8., -223./140., 7./2., -7./2., 35./12., -7./4., 7./10., -1./6., 1./56.},
		{1./56., -2./7., -19./20., 2., -5./4., 2./3., -1./4., 2./35., -1./168.},
		{-1./168., 1./14., -1./2., -9./20., 5./4., -1./2., 1./6., -1./28., 1./280.},
		{1./280., -4./105., 1./5., -4./5., 0, 4./5., -1./5., 4./105., -1./280.},
		{-1./280., 1./28., -1./6., 1./2., -5./4., 9./20., 1./2., -1./14., 1./168.},
		{1./168., -2./35., 1./4., -2./3., 5./4., -2., 19./20., 2./7., -1./56.},
		{-1./56., 1./6., -7./10., 7./4., -35./12., 7./2., -7./2., 223./140., 1./8.},
		{1./8., -8./7., 14./3., -56./5., 35./2., -56./3., 14., -8., 761./280.}
	};

	UINT4 vecLength;
	vecLength = XVec->length;


	REAL8 h;
	h = fabs(XVec->data[0] - XVec->data[1]);

	UINT4 i;
	UINT4 j;

	for (i = 0; i <= vecLength-1; i++)
	{
		if (i == 0)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[0][j] * YVec->data[j]
                );
			}
		}
		else if (i == 1)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[1][j] * YVec->data[j]
                );
			}
		}
		else if (i == 2)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[2][j] * YVec->data[j]
                );
			}
		}
		else if (i == 3)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[3][j] * YVec->data[j]
                );
			}
		}
		else if (i == vecLength-4)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[5][j] * YVec->data[i+j-5]
                );
			}
		}
		else if (i == vecLength-3)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[6][j] * YVec->data[i+j-6]
                );
			}
		}
		else if (i == vecLength-2)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[7][j] * YVec->data[i+j-7]
                );
			}
		}
		else if (i == vecLength-1)
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[8][j] * YVec->data[i+j-8]
                );
			}
		}
		else
		{
			for (j = 0; j <= 8; j++)
			{
				derivativeVec->data[i] += (
                    eightOrderCoeffs[4][j] * YVec->data[i+j-4]
                );
			}
		}

		derivativeVec->data[i] /= h;
	}

	return XLAL_SUCCESS;
}

/**
 * Function which calculates the 3-order cumulative derivative of a
 * numerical function
 */
int
XLALCumulativeIntegral3(
	REAL8Vector *XVec,
    /**<< An array of X values */
	REAL8Vector *YVec,
    /**<< An array of Y values */
	REAL8Vector *integralVec
    /**<< OUTPUT, the integral of Y dX */
)
{
	UINT4 vecLength;
	vecLength = XVec->length;

	REAL8Vector *XVecExt = XLALCreateREAL8Vector(vecLength+2);
	REAL8Vector *YVecExt = XLALCreateREAL8Vector(vecLength+2);

	memset(XVecExt->data, 0, XVecExt->length * sizeof(REAL8));
	memset(YVecExt->data, 0, YVecExt->length * sizeof(REAL8));

	REAL8 *X0, *X1, *X2, *X3;
	REAL8 *Y0, *Y1, *Y2, *Y3;

	REAL8 a, b, c, d, e, h, g, z;
	const REAL8 oo12 = 0.08333333333333333;

	UINT4 i;

	for (i=1; i < vecLength+1; i++)
	{
		XVecExt->data[i] = XVec->data[i-1];
		YVecExt->data[i] = YVec->data[i-1];
	}

	XVecExt->data[0] = XVec->data[3];
	XVecExt->data[vecLength+1] = XVec->data[vecLength-4];

	YVecExt->data[0] = YVec->data[3];
	YVecExt->data[vecLength+1] = YVec->data[vecLength-4];

	X0 = &XVecExt->data[0];
	X1 = &XVecExt->data[1];
	X2 = &XVecExt->data[2];
	X3 = &XVecExt->data[3];

	Y0 = &YVecExt->data[0];
	Y1 = &YVecExt->data[1];
	Y2 = &YVecExt->data[2];
	Y3 = &YVecExt->data[3];

	for (i=0; i < vecLength-1; i++)
	{
		a = X1[i] - X0[i];
		b = X2[i] - X1[i];
		c = X3[i] - X2[i];
		d = Y1[i] - Y0[i];
		e = Y2[i] - Y1[i];
		h = Y3[i] - Y2[i];
		g = 0.5 * (Y1[i]+Y2[i]);
		z = (
            b * g
            + oo12 * b * b * (
                c * b * (2 * c + b) * (c + b) * d
                - a * c * (c - a) * (2 * c + 2 * a + 3 * b) * e
                - a * b * (2 * a + b) * (a + b) * h
            ) / (a * c * (a + b) * (c + b) * (c + a + b))
        );
		integralVec->data[i+1] = integralVec->data[i] + z;
	}
	XLALDestroyREAL8Vector(XVecExt);
	XLALDestroyREAL8Vector(YVecExt);

	return XLAL_SUCCESS;
}
