/*----------------------------------------------------------------------- 
 * 
 * File Name: Integrate.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _INTEGRATE_H
#include "Integrate.h"
#ifndef _INTEGRATE_H
#define _INTEGRATE_H
#endif
#endif

#ifndef _INTERPOLATE_H
#include "Interpolate.h"
#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H
#endif
#endif


NRCSID (INTEGRATEC, "$Id$");

typedef struct
tagIntegralState
{
  REAL4 integral;    /* value of integral at current refinement level */
  INT4  refinement;  /* levels of refinement                          */
}
IntegralState;

static void
Trapezoid (
    Status        *status,
    IntegralState *output,
    IntegrateIn   *input,
    void          *params
    )
{
  INITSTATUS (status, INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  if (output->refinement)
  {
    INT4  i;
    INT4  n   = 1 << (output->refinement - 1);
    REAL4 dx  = (input->xmax - input->xmin)/n;
    REAL4 x   = input->xmin + dx/2;
    REAL4 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL4 y;
      input->function (status->statusPtr, &y, x, params);
      CHECKSTATUSPTR (status);
      sum += y;
      x   += dx;
    }
    output->integral += (input->xmax - input->xmin)*sum/n;
    output->integral /= 2;
  }
  else
  {
    REAL4 y0;
    REAL4 y1;
    input->function (status->statusPtr, &y0, input->xmin, params);
    CHECKSTATUSPTR (status);
    input->function (status->statusPtr, &y1, input->xmax, params);
    CHECKSTATUSPTR (status);
    output->integral = (input->xmax - input->xmin)*(y1 - y0)/2;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
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


static REAL4
EqualsX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 1;
  return x;
}

static REAL4
EqualsInvX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  REAL4 invx = 1/x;
  *jac = invx*invx;
  return invx;
}

static REAL4
EqualsAPlusXSq (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 2*x;
  return a + x*x;
}

static REAL4
EqualsBMinusXSq (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 2*x;
  return b - x*x;
}

static REAL4
EqualsMinusLogX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 1/x;
  return -log(x);
}

static void
Midpoint (
    Status        *status,
    IntegralState *output,
    IntegrateIn   *input,
    void          *params
    )
{
  REAL4 (*ChangeOfVariables)(REAL4, REAL4, REAL4, REAL4 *);
  REAL4 a = input->xmin;
  REAL4 b = input->xmax;
  REAL4 xmax;
  REAL4 xmin;

  INITSTATUS (status, INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  switch (input->type)
  {
    case OpenInterval:
      ChangeOfVariables = EqualsX;
      xmax = b;
      xmin = a;
      break;

    case SingularLowerLimit:
      ChangeOfVariables = EqualsAPlusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case SingularUpperLimit:
      ChangeOfVariables = EqualsBMinusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case InfiniteDomainPow:
      ASSERT ((b > 0 && a > 0) || (b < 0 && a < 0), status,
              INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);
      ChangeOfVariables = EqualsInvX;
      xmax = 1/a;
      xmin = 1/b;
      break;

    case InfiniteDomainExp:
      ChangeOfVariables = EqualsMinusLogX;
      xmax = exp(-a);
      xmin = 0;
      break;

    default: /* unrecognized type */
      ABORT (status, INTEGRATE_ETYPE, INTEGRATE_MSGETYPE);
  }

  if (output->refinement)
  {
    INT4  i;
    INT4  n   = ThreePow (output->refinement - 1);
    REAL4 dx1 = (xmax - xmin)/(3*n);
    REAL4 dx2 = 2*dx1;
    REAL4 x   = xmin + dx1/2;
    REAL4 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL4 y;
      REAL4 jac;

      input->function (status->statusPtr, &y,
                       ChangeOfVariables (x, a, b, &jac), params);
      CHECKSTATUSPTR (status);
      sum += y*jac;
      x   += dx2;

      input->function (status->statusPtr, &y,
                       ChangeOfVariables (x, a, b, &jac), params);
      CHECKSTATUSPTR (status);
      sum += y*jac;
      x   += dx1;
    }
    output->integral += (xmax - xmin)*sum/n;
    output->integral /= 3;
  }
  else
  {
    REAL4 x = (xmax + xmin)/2;
    REAL4 y;
    REAL4 jac;
    input->function (status->statusPtr, &y,
                     ChangeOfVariables (x, a, b, &jac), params);
    CHECKSTATUSPTR (status);
    output->integral = (xmax - xmin)*y*jac;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}





void
RombergIntegrate (
    Status      *status,
    REAL4       *result,
    IntegrateIn *input,
    void        *params
    )
{
  const REAL4 epsilon = 1e-6;
  enum { MaxSteps     = 20 };
  enum { Order        = 4  };

  void (*Algorithm)(Status *, IntegralState *, IntegrateIn *, void *);

  IntegralState state;
  REAL4         integral[MaxSteps + 1];
  REAL4         stepSize[MaxSteps + 1];
  REAL4         refineFactor;

  INITSTATUS (status, INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  ASSERT (result, status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->function,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->xmax > input->xmin, status,
          INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);

  switch (input->type)
  {
    case ClosedInterval:
      Algorithm    = Trapezoid;
      refineFactor = 0.25;
      break;

    case OpenInterval:
    case SingularLowerLimit:
    case SingularUpperLimit:
    case InfiniteDomainPow:
    case InfiniteDomainExp:
      Algorithm    = Midpoint;
      refineFactor = 1.0/9.0;
      break;

    default: /* unrecognized type */
      ABORT (status, INTEGRATE_ETYPE, INTEGRATE_MSGETYPE);
  }

  stepSize[0] = 1;
  for (state.refinement = 0; state.refinement < MaxSteps; ++state.refinement)
  {
    Algorithm (status->statusPtr, &state, input, params);
    CHECKSTATUSPTR (status);

    integral[state.refinement] = state.integral;

    if (state.refinement >= Order)
    {
      InterpolatePar intpar;
      InterpolateOut intout;

      intpar.n = Order + 1;
      intpar.x = stepSize + state.refinement - Order;
      intpar.y = integral + state.refinement - Order;

      /* extrapolate to continuum limit (stepSize = 0) */
      PolynomialInterpolation (status->statusPtr, &intout, 0, &intpar);
      CHECKSTATUSPTR (status);

      if (fabs(intout.dy) < epsilon*fabs(intout.y))
      {
        *result = intout.y;
        DETATCHSTATUSPTR (status);
        RETURN (status);
      }
    }

    integral[state.refinement + 1] = state.integral;
    stepSize[state.refinement + 1] = refineFactor*stepSize[state.refinement];
  }

  ABORT (status, INTEGRATE_EMXIT, INTEGRATE_MSGEMXIT);
}
