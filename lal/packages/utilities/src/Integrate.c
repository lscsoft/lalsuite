/*----------------------------------------------------------------------- 
 * 
 * File Name: Integrate.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include "LALStdlib.h"
#include "Integrate.h"
#include "Interpolate.h"

NRCSID (INTEGRATEC, "$Id$");

typedef struct
tagSIntegralState
{
  REAL4 integral;    /* value of integral at current refinement level */
  INT4  refinement;  /* levels of refinement                          */
}
SIntegralState;

typedef struct
tagDIntegralState
{
  REAL8 integral;    /* value of integral at current refinement level */
  INT4  refinement;  /* levels of refinement                          */
}
DIntegralState;

static void
STrapezoid (
    LALStatus         *status,
    SIntegralState *output,
    SIntegrateIn   *input,
    void           *params
    )
{
  INITSTATUS (status, "STrapezoid", INTEGRATEC);
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
    output->integral = (input->xmax - input->xmin)*(y1 + y0)/2;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


static void
DTrapezoid (
    LALStatus         *status,
    DIntegralState *output,
    DIntegrateIn   *input,
    void           *params
    )
{
  INITSTATUS (status, "DTrapezoid", INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  if (output->refinement)
  {
    INT4  i;
    INT4  n   = 1 << (output->refinement - 1);
    REAL8 dx  = (input->xmax - input->xmin)/n;
    REAL8 x   = input->xmin + dx/2;
    REAL8 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL8 y;
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
    REAL8 y0;
    REAL8 y1;
    input->function (status->statusPtr, &y0, input->xmin, params);
    CHECKSTATUSPTR (status);
    input->function (status->statusPtr, &y1, input->xmax, params);
    CHECKSTATUSPTR (status);
    output->integral = (input->xmax - input->xmin)*(y1 + y0)/2;
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
SEqualsX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 1;
  return x;
}

static REAL4
SEqualsInvX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  REAL4 invx = 1/x;
  *jac = invx*invx;
  return invx;
}

static REAL4
SEqualsAPlusXSq (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 2*x;
  return a + x*x;
}

static REAL4
SEqualsBMinusXSq (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 2*x;
  return b - x*x;
}

static REAL4
SEqualsMinusLogX (REAL4 x, REAL4 a, REAL4 b, REAL4 *jac)
{
  *jac = 1/x;
  return -log(x);
}

static void
SMidpoint (
    LALStatus         *status,
    SIntegralState *output,
    SIntegrateIn   *input,
    void           *params
    )
{
  REAL4 (*ChangeOfVariables)(REAL4, REAL4, REAL4, REAL4 *);
  REAL4 a = input->xmin;
  REAL4 b = input->xmax;
  REAL4 xmax;
  REAL4 xmin;

  INITSTATUS (status, "SMidpoint", INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  switch (input->type)
  {
    case OpenInterval:
      ChangeOfVariables = SEqualsX;
      xmax = b;
      xmin = a;
      break;

    case SingularLowerLimit:
      ChangeOfVariables = SEqualsAPlusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case SingularUpperLimit:
      ChangeOfVariables = SEqualsBMinusXSq;
      xmax = sqrt(b - a);
      xmin = 0;
      break;

    case InfiniteDomainPow:
      ASSERT ((b > 0 && a > 0) || (b < 0 && a < 0), status,
              INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);
      ChangeOfVariables = SEqualsInvX;
      xmax = 1/a;
      xmin = 1/b;
      break;

    case InfiniteDomainExp:
      ChangeOfVariables = SEqualsMinusLogX;
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


static REAL8
DEqualsX (REAL8 x, REAL8 a, REAL8 b, REAL8 *jac)
{
  *jac = 1;
  return x;
}

static REAL8
DEqualsInvX (REAL8 x, REAL8 a, REAL8 b, REAL8 *jac)
{
  REAL8 invx = 1/x;
  *jac = invx*invx;
  return invx;
}

static REAL8
DEqualsAPlusXSq (REAL8 x, REAL8 a, REAL8 b, REAL8 *jac)
{
  *jac = 2*x;
  return a + x*x;
}

static REAL8
DEqualsBMinusXSq (REAL8 x, REAL8 a, REAL8 b, REAL8 *jac)
{
  *jac = 2*x;
  return b - x*x;
}

static REAL8
DEqualsMinusLogX (REAL8 x, REAL8 a, REAL8 b, REAL8 *jac)
{
  *jac = 1/x;
  return -log(x);
}

static void
DMidpoint (
    LALStatus         *status,
    DIntegralState *output,
    DIntegrateIn   *input,
    void           *params
    )
{
  REAL8 (*ChangeOfVariables)(REAL8, REAL8, REAL8, REAL8 *);
  REAL8 a = input->xmin;
  REAL8 b = input->xmax;
  REAL8 xmax;
  REAL8 xmin;

  INITSTATUS (status, "DMidpoint", INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  switch (input->type)
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
      ASSERT ((b > 0 && a > 0) || (b < 0 && a < 0), status,
              INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);
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
      ABORT (status, INTEGRATE_ETYPE, INTEGRATE_MSGETYPE);
  }

  if (output->refinement)
  {
    INT4  i;
    INT4  n   = ThreePow (output->refinement - 1);
    REAL8 dx1 = (xmax - xmin)/(3*n);
    REAL8 dx2 = 2*dx1;
    REAL8 x   = xmin + dx1/2;
    REAL8 sum = 0;
    for (i = 0; i < n; ++i)
    {
      REAL8 y;
      REAL8 jac;

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
    REAL8 x = (xmax + xmin)/2;
    REAL8 y;
    REAL8 jac;
    input->function (status->statusPtr, &y,
                     ChangeOfVariables (x, a, b, &jac), params);
    CHECKSTATUSPTR (status);
    output->integral = (xmax - xmin)*y*jac;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}




void
LALSRombergIntegrate (
    LALStatus       *status,
    REAL4        *result,
    SIntegrateIn *input,
    void         *params
    )
{
  const REAL4 epsilon = 1e-6;
  enum { MaxSteps     = 20 };
  enum { Order        = 4  };

  void (*Algorithm)(LALStatus *, SIntegralState *, SIntegrateIn *, void *);

  SIntegralState state;
  REAL4          integral[MaxSteps + 1];
  REAL4          stepSize[MaxSteps + 1];
  REAL4          refineFactor;

  INITSTATUS (status, "LALSRombergIntegrate", INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  ASSERT (result, status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->function,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->xmax > input->xmin, status,
          INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);

  switch (input->type)
  {
    case ClosedInterval:
      Algorithm    = STrapezoid;
      refineFactor = 0.25;
      break;

    case OpenInterval:
    case SingularLowerLimit:
    case SingularUpperLimit:
    case InfiniteDomainPow:
    case InfiniteDomainExp:
      Algorithm    = SMidpoint;
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
      SInterpolatePar intpar;
      SInterpolateOut intout;

      intpar.n = Order + 1;
      intpar.x = stepSize + state.refinement - Order;
      intpar.y = integral + state.refinement - Order;

      /* extrapolate to continuum limit (stepSize = 0) */
      LALSPolynomialInterpolation (status->statusPtr, &intout, 0, &intpar);
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


void
LALDRombergIntegrate (
    LALStatus       *status,
    REAL8        *result,
    DIntegrateIn *input,
    void         *params
    )
{
  const REAL8 epsilon = 1e-15;
  enum { MaxSteps     = 20 };
  enum { Order        = 4  };

  void (*Algorithm)(LALStatus *, DIntegralState *, DIntegrateIn *, void *);

  DIntegralState state;
  REAL8          integral[MaxSteps + 1];
  REAL8          stepSize[MaxSteps + 1];
  REAL8          refineFactor;

  INITSTATUS (status, "LALDRombergIntegrate", INTEGRATEC);
  ATTATCHSTATUSPTR (status);

  ASSERT (result, status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->function,  status, INTEGRATE_ENULL, INTEGRATE_MSGENULL);
  ASSERT (input->xmax > input->xmin, status,
          INTEGRATE_EIDOM, INTEGRATE_MSGEIDOM);

  switch (input->type)
  {
    case ClosedInterval:
      Algorithm    = DTrapezoid;
      refineFactor = 0.25;
      break;

    case OpenInterval:
    case SingularLowerLimit:
    case SingularUpperLimit:
    case InfiniteDomainPow:
    case InfiniteDomainExp:
      Algorithm    = DMidpoint;
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
      DInterpolatePar intpar;
      DInterpolateOut intout;

      intpar.n = Order + 1;
      intpar.x = stepSize + state.refinement - Order;
      intpar.y = integral + state.refinement - Order;

      /* extrapolate to continuum limit (stepSize = 0) */
      LALDPolynomialInterpolation (status->statusPtr, &intout, 0, &intpar);
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
