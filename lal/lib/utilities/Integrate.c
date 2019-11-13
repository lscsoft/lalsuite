/*
*  Copyright (C) 2007 Jolien Creighton, Drew Keppel
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

/* ---------- see Integrate.h for doxygen documentation ---------- */
#include <math.h>
#include <string.h>
#include <lal/Integrate.h>
#include <lal/Interpolate.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALError.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static REAL8
XLALREAL8Trapezoid (
    REAL8 integral,
    REAL8 (*f)(REAL8 x, void *params),
    void *params,
    REAL8 xmin,
    REAL8 xmax,
    IntegralType type,
    int refinement
    )
{
  switch (type)
  {
    default:
      break;
  }

  if (refinement)
  {
    INT4  i;
    INT4  n   = 1 << (refinement - 1);
    REAL8 dx  = (xmax - xmin)/n;
    REAL8 x   = xmin + dx/2;
    REAL8 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL8 y;
      y = f (x, params);
      if (xlalErrno)
        XLAL_ERROR_REAL8(XLAL_EFUNC);
      sum += y;
      x   += dx;
    }
    integral += (xmax - xmin)*sum/n;
    integral /= 2;
  }
  else
  {
    REAL8 y_0;
    REAL8 y_1;
    y_0 = f (xmin, params);
    if (xlalErrno)
      XLAL_ERROR_REAL8(XLAL_EFUNC);
    y_1 = f (xmax, params);
    if (xlalErrno)
      XLAL_ERROR_REAL8(XLAL_EFUNC);
    integral = (xmax - xmin)*(y_1 + y_0)/2;
  }

  return integral;
}


static INT4
ThreePow (INT4 n)
{
  INT4 p = 1;
  while (n-- > 0)
  {
    p *= 3;
  }
  return p;
}


static REAL8
DEqualsX (REAL8 x, REAL8 UNUSED a, REAL8 UNUSED b, REAL8 *jac)
{
  *jac = 1;
  return x;
}

static REAL8
DEqualsInvX (REAL8 x, REAL8 UNUSED a, REAL8 UNUSED b, REAL8 *jac)
{
  REAL8 invx = 1/x;
  *jac = invx*invx;
  return invx;
}

static REAL8
DEqualsAPlusXSq (REAL8 x, REAL8 a, REAL8 UNUSED b, REAL8 *jac)
{
  *jac = 2*x;
  return a + x*x;
}

static REAL8
DEqualsBMinusXSq (REAL8 x, REAL8 UNUSED a, REAL8 b, REAL8 *jac)
{
  *jac = 2*x;
  return b - x*x;
}

static REAL8
DEqualsMinusLogX (REAL8 x, REAL8 UNUSED a, REAL8 UNUSED b, REAL8 *jac)
{
  *jac = 1/x;
  return -log(x);
}

static REAL8
XLALREAL8Midpoint (
    REAL8 integral,
    REAL8 (*f)(REAL8 x, void *params),
    void *params,
    REAL8 a,
    REAL8 b,
    IntegralType type,
    int refinement
    )
{
  REAL8 (*ChangeOfVariables)(REAL8, REAL8, REAL8, REAL8 *);
  REAL8 xmax;
  REAL8 xmin;

  switch (type)
  {
    case OpenInterval:
      ChangeOfVariables = DEqualsX;
      xmax = b;
      xmin = a;
      break;

    case SingularLowerLimit:
      ChangeOfVariables = DEqualsAPlusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case SingularUpperLimit:
      ChangeOfVariables = DEqualsBMinusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case InfiniteDomainPow:
      if (!((b > 0 && a > 0) || (b < 0 && a < 0)))
        XLAL_ERROR_REAL8(XLAL_EDOM);
      ChangeOfVariables = DEqualsInvX;
      xmax = 1/a;
      xmin = 1/b;
      break;

    case InfiniteDomainExp:
      ChangeOfVariables = DEqualsMinusLogX;
      xmax = exp(-a);
      xmin = 0;
      break;

    default: /* unrecognized type */
      XLAL_ERROR_REAL8(XLAL_EINVAL);
  }

  if (refinement)
  {
    INT4  i;
    INT4  n   = ThreePow (refinement - 1);
    REAL8 dx1 = (xmax - xmin)/(3*n);
    REAL8 dx2 = 2*dx1;
    REAL8 x   = xmin + dx1/2;
    REAL8 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL8 y;
      REAL8 jac;

      y = f (ChangeOfVariables (x, a, b, &jac), params);
      if (xlalErrno)
        XLAL_ERROR_REAL8(XLAL_EFUNC);
      sum += y*jac;
      x   += dx2;

      y = f (ChangeOfVariables (x, a, b, &jac), params);
      if (xlalErrno)
        XLAL_ERROR_REAL8(XLAL_EFUNC);
      sum += y*jac;
      x   += dx1;
    }
    integral += (xmax - xmin)*sum/n;
    integral /= 3;
  }
  else
  {
    REAL8 x = (xmax + xmin)/2;
    REAL8 y;
    REAL8 jac;
    y = f (ChangeOfVariables (x, a, b, &jac), params);
    if (xlalErrno)
      XLAL_ERROR_REAL8(XLAL_EFUNC);
    integral = (xmax - xmin)*y*jac;
  }

  return integral;
}



REAL8
XLALREAL8RombergIntegrate (
    REAL8 (*f)(REAL8 x, void *params),
    void *params,
    REAL8 xmin,
    REAL8 xmax,
    IntegralType type
    )
{
  const REAL8 epsilon = 1e-15;
  enum { MaxSteps     = 20 };
  enum { Order        = 4  };

  REAL8 (*Algorithm)(REAL8, REAL8 (*)(REAL8, void *), void *, REAL8, REAL8, IntegralType, int);

  REAL8          integral[MaxSteps + 1];
  REAL8          stepSize[MaxSteps + 1];
  REAL8          refineFactor;
  REAL8          temp = 0;
  int            refinement;

  if (f == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (xmax <= xmin)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  switch (type)
  {
    case ClosedInterval:
      Algorithm    = XLALREAL8Trapezoid;
      refineFactor = 0.25;
      break;

    case OpenInterval:
    case SingularLowerLimit:
    case SingularUpperLimit:
    case InfiniteDomainPow:
    case InfiniteDomainExp:
      Algorithm    = XLALREAL8Midpoint;
      refineFactor = 1.0/9.0;
      break;

    default: /* unrecognized type */
      XLAL_ERROR_REAL8(XLAL_EINVAL);
  }

  stepSize[0] = 1;
  for (refinement = 0; refinement < MaxSteps; ++refinement)
  {
    temp = Algorithm (temp, f, params, xmin, xmax, type, refinement);
    if (xlalErrno)
      XLAL_ERROR_REAL8(XLAL_EFUNC);

    integral[refinement] = temp;

    if (refinement >= Order)
    {
      REAL8 dy;
      REAL8 y;

      /* extrapolate to continuum limit (stepSize = 0) */
      dy = XLALREAL8PolynomialInterpolation (&y, 0, integral + refinement - Order, stepSize + refinement - Order, Order + 1);
      if (xlalErrno)
        XLAL_ERROR_REAL8(XLAL_EFUNC);

      if (fabs(dy) < epsilon*fabs(y))
      {
        return y;
      }
    }

    integral[refinement + 1] = temp;
    stepSize[refinement + 1] = refineFactor*stepSize[refinement];
  }

  XLAL_ERROR_REAL8(XLAL_EMAXITER);
}
