/*----------------------------------------------------------------------- 
 * 
 * File Name: FindRoot.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include "LALStdlib.h"
#include "LALConstants.h"
#include "FindRoot.h"

NRCSID (FINDROOTC, "$Id$");

void
LALSBracketRoot (
    LALStatus      *status,
    SFindRootIn *inout,
    void        *params
    )
{
  const REAL4 fac  = LAL_SQRT2;
  const INT4  imax = 64;

  INT4  i = 0;
  REAL4 y1;
  REAL4 y2;

  INITSTATUS (status, "LALSBracketRoot", FINDROOTC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (inout, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (inout->function, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  /* params can be NULL ... */

  ASSERT (inout->xmax != inout->xmin, status,
          FINDROOT_EIDOM, FINDROOT_MSGEIDOM);

  /* evaluate function at endpoints */

  inout->function (status->statusPtr, &y1, inout->xmin, params);
  CHECKSTATUSPTR (status);

  inout->function (status->statusPtr, &y2, inout->xmax, params);
  CHECKSTATUSPTR (status);

  while (1)
  {
    /* break out if root has been bracketed */
    if (y1*y2 < 0)
    {
      break;
    }

    /* increment iteration count */
    ASSERT (i < imax, status, FINDROOT_EMXIT, FINDROOT_MSGEMXIT);
    ++i;

    if (fabs(y1) < fabs(y2))
    {
      /* expand lower limit */
      inout->xmin += fac*(inout->xmin - inout->xmax);
      inout->function (status->statusPtr, &y1, inout->xmin, params);
      CHECKSTATUSPTR (status);
    }
    else
    {
      /* expand upper limit */
      inout->xmax += fac*(inout->xmax - inout->xmin);
      inout->function (status->statusPtr, &y2, inout->xmax, params);
      CHECKSTATUSPTR (status);
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALDBracketRoot (
    LALStatus      *status,
    DFindRootIn *inout,
    void        *params
    )
{
  const REAL8 fac  = LAL_SQRT2;
  const INT4  imax = 64;

  INT4  i = 0;
  REAL8 y1;
  REAL8 y2;

  INITSTATUS (status, "LALDBracketRoot", FINDROOTC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (inout, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (inout->function, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  /* params can be NULL ... */

  ASSERT (inout->xmax != inout->xmin, status,
          FINDROOT_EIDOM, FINDROOT_MSGEIDOM);

  /* evaluate function at endpoints */

  inout->function (status->statusPtr, &y1, inout->xmin, params);
  CHECKSTATUSPTR (status);

  inout->function (status->statusPtr, &y2, inout->xmax, params);
  CHECKSTATUSPTR (status);

  while (1)
  {
    /* break out if root has been bracketed */
    if (y1*y2 < 0)
    {
      break;
    }

    /* increment iteration count */
    ASSERT (i < imax, status, FINDROOT_EMXIT, FINDROOT_MSGEMXIT);
    ++i;

    if (fabs(y1) < fabs(y2))
    {
      /* expand lower limit */
      inout->xmin += fac*(inout->xmin - inout->xmax);
      inout->function (status->statusPtr, &y1, inout->xmin, params);
      CHECKSTATUSPTR (status);
    }
    else
    {
      /* expand upper limit */
      inout->xmax += fac*(inout->xmax - inout->xmin);
      inout->function (status->statusPtr, &y2, inout->xmax, params);
      CHECKSTATUSPTR (status);
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALSBisectionFindRoot (
    LALStatus      *status,
    REAL4       *root,
    SFindRootIn *input,
    void        *params
    )
{
  const INT4 imax = 40;

  INT4  i = 0;
  REAL4 y1;
  REAL4 y2;
  REAL4 x;
  REAL4 dx;

  INITSTATUS (status, "LALSBisectionFindRoot", FINDROOTC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (root, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (input, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (input->function, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  /* params can be NULL ... */

  /* evaluate function at endpoints */

  input->function (status->statusPtr, &y1, input->xmin, params);
  CHECKSTATUSPTR (status);

  input->function (status->statusPtr, &y2, input->xmax, params);
  CHECKSTATUSPTR (status);

  ASSERT (y1*y2 < 0, status, FINDROOT_EBRKT, FINDROOT_MSGEBRKT);

  if (y1 < 0)
  {
    /* start search at xmin and increase */
    x  = input->xmin;
    dx = input->xmax - input->xmin;
  }
  else
  {
    /* start search at xmax and decrease */
    x  = input->xmax;
    dx = input->xmin - input->xmax;
  }

  /* infinite loop to locate root */
  while (1)
  {
    REAL4 xmid;
    REAL4 ymid;

    /* increment iteration count */
    ASSERT (i < imax, status, FINDROOT_EMXIT, FINDROOT_MSGEMXIT);
    ++i;

    /* locate midpoint of domain */
    dx   /= 2;
    xmid  = x + dx;
    
    /* evaluate function at midpoint */
    input->function (status->statusPtr, &ymid, xmid, params);
    CHECKSTATUSPTR (status);

    if (ymid < 0)
    {
      /* function is in second half of domain */
      x = xmid;
    }
    else if (ymid == 0)
    {
      /* root has been found */
      *root = xmid;
      break;
    }

    if (fabs(dx) < input->xacc)
    {
      /* domain has shrunk to acceptably small size */
      *root = xmid;
      break;
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALDBisectionFindRoot (
    LALStatus      *status,
    REAL8       *root,
    DFindRootIn *input,
    void        *params
    )
{
  const INT4 imax = 80;

  INT4  i = 0;
  REAL8 y1;
  REAL8 y2;
  REAL8 x;
  REAL8 dx;

  INITSTATUS (status, "LALDBisectionFindRoot", FINDROOTC);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (root, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (input, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  ASSERT (input->function, status, FINDROOT_ENULL, FINDROOT_MSGENULL);
  /* params can be NULL ... */

  /* evaluate function at endpoints */

  input->function (status->statusPtr, &y1, input->xmin, params);
  CHECKSTATUSPTR (status);

  input->function (status->statusPtr, &y2, input->xmax, params);
  CHECKSTATUSPTR (status);

  ASSERT (y1*y2 < 0, status, FINDROOT_EBRKT, FINDROOT_MSGEBRKT);

  if (y1 < 0)
  {
    /* start search at xmin and increase */
    x  = input->xmin;
    dx = input->xmax - input->xmin;
  }
  else
  {
    /* start search at xmax and decrease */
    x  = input->xmax;
    dx = input->xmin - input->xmax;
  }

  /* infinite loop to locate root */
  while (1)
  {
    REAL8 xmid;
    REAL8 ymid;

    /* increment iteration count */
    ASSERT (i < imax, status, FINDROOT_EMXIT, FINDROOT_MSGEMXIT);
    ++i;

    /* locate midpoint of domain */
    dx   /= 2;
    xmid  = x + dx;
    
    /* evaluate function at midpoint */
    input->function (status->statusPtr, &ymid, xmid, params);
    CHECKSTATUSPTR (status);

    if (ymid < 0)
    {
      /* function is in second half of domain */
      x = xmid;
    }
    else if (ymid == 0)
    {
      /* root has been found */
      *root = xmid;
      break;
    }

    if (fabs(dx) < input->xacc)
    {
      /* domain has shrunk to acceptably small size */
      *root = xmid;
      break;
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

