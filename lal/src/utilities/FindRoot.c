/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
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

/* ---------- see FindRoot.h for doxygen documentation ---------- */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindRoot.h>

/** \see See \ref FindRoot_h for documentation */
void
LALSBracketRoot (
    LALStatus      *status,
    SFindRootIn *inout,
    void        *params
    )
{
  const REAL4 fac  = LAL_SQRT2;

  INT4  i = 0;
  REAL4 y_1;
  REAL4 y_2;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (inout, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (inout->function, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  /* params can be NULL ... */

  ASSERT (inout->xmax != inout->xmin, status,
          FINDROOTH_EIDOM, FINDROOTH_MSGEIDOM);

  /* evaluate function at endpoints */

  inout->function (status->statusPtr, &y_1, inout->xmin, params);
  CHECKSTATUSPTR (status);

  inout->function (status->statusPtr, &y_2, inout->xmax, params);
  CHECKSTATUSPTR (status);

  while (1)
  {
    /* break out if root has been bracketed */
    if (y_1*y_2 < 0)
    {
      break;
    }

    /* increment iteration count */
#ifndef LAL_NDEBUG
    INT4 imax = 64;
    ASSERT (i < imax, status, FINDROOTH_EMXIT, FINDROOTH_MSGEMXIT);
#endif
    ++i;

    if (fabs(y_1) < fabs(y_2))
    {
      /* expand lower limit */
      inout->xmin += fac*(inout->xmin - inout->xmax);
      inout->function (status->statusPtr, &y_1, inout->xmin, params);
      CHECKSTATUSPTR (status);
    }
    else
    {
      /* expand upper limit */
      inout->xmax += fac*(inout->xmax - inout->xmin);
      inout->function (status->statusPtr, &y_2, inout->xmax, params);
      CHECKSTATUSPTR (status);
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/** \see See \ref FindRoot_h for documentation */
int
XLALDBracketRoot(
    REAL8 (*y)(REAL8, void *),
    REAL8 *xmin,
    REAL8 *xmax,
    void *params
)
{
  const INT4 imax = 64;
  INT4 i;
  REAL8 y_1;
  REAL8 y_2;

  /* put xmin and xmax in the correct order (use y_1 as temporary storage) */
  if(*xmin > *xmax)
  {
    y_1 = *xmin;
    *xmin = *xmax;
    *xmax = y_1;
  }

  /* evaluate function at endpoints */
  y_1 = y(*xmin, params);
  y_2 = y(*xmax, params);

  /* loop until iteration limit exceeded or root bracketed */
  for(i = 0; y_1 * y_2 >= 0.0; i++)
  {
    if(XLALIsREAL8FailNaN(y_1) || XLALIsREAL8FailNaN(y_2))
      XLAL_ERROR(XLAL_EFUNC);
    if(i >= imax)
      XLAL_ERROR(XLAL_EMAXITER);
    if(fabs(y_1) < fabs(y_2))
    {
      /* expand lower limit */
      *xmin -= LAL_SQRT2 * (*xmax - *xmin);
      y_1 = y(*xmin, params);
    }
    else
    {
      /* expand upper limit */
      *xmax += LAL_SQRT2 * (*xmax - *xmin);
      y_2 = y(*xmax, params);
    }
  }

  return(0);
}

/** \see See \ref FindRoot_h for documentation */
void
LALDBracketRoot (
    LALStatus      *status,
    DFindRootIn *inout,
    void        *params
    )
{
  const REAL8 fac  = LAL_SQRT2;

  INT4  i = 0;
  REAL8 y_1;
  REAL8 y_2;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (inout, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (inout->function, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  /* params can be NULL ... */

  ASSERT (inout->xmax != inout->xmin, status,
          FINDROOTH_EIDOM, FINDROOTH_MSGEIDOM);

  /* evaluate function at endpoints */

  inout->function (status->statusPtr, &y_1, inout->xmin, params);
  CHECKSTATUSPTR (status);

  inout->function (status->statusPtr, &y_2, inout->xmax, params);
  CHECKSTATUSPTR (status);

  while (1)
  {
    /* break out if root has been bracketed */
    if (y_1*y_2 < 0)
    {
      break;
    }

    /* increment iteration count */
#ifndef LAL_NDEBUG
    INT4 imax = 64;
    ASSERT (i < imax, status, FINDROOTH_EMXIT, FINDROOTH_MSGEMXIT);
#endif
    ++i;

    if (fabs(y_1) < fabs(y_2))
    {
      /* expand lower limit */
      inout->xmin += fac*(inout->xmin - inout->xmax);
      inout->function (status->statusPtr, &y_1, inout->xmin, params);
      CHECKSTATUSPTR (status);
    }
    else
    {
      /* expand upper limit */
      inout->xmax += fac*(inout->xmax - inout->xmin);
      inout->function (status->statusPtr, &y_2, inout->xmax, params);
      CHECKSTATUSPTR (status);
    }

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/** \see See \ref FindRoot_h for documentation */
void
LALSBisectionFindRoot (
    LALStatus      *status,
    REAL4       *root,
    SFindRootIn *input,
    void        *params
    )
{

  INT4  i = 0;
  REAL4 y_1;
  REAL4 y_2;
  REAL4 x;
  REAL4 dx;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (root, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (input, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (input->function, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  /* params can be NULL ... */

  /* evaluate function at endpoints */

  input->function (status->statusPtr, &y_1, input->xmin, params);
  CHECKSTATUSPTR (status);

  input->function (status->statusPtr, &y_2, input->xmax, params);
  CHECKSTATUSPTR (status);

  ASSERT (y_1*y_2 < 0, status, FINDROOTH_EBRKT, FINDROOTH_MSGEBRKT);

  if (y_1 < 0)
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
#ifndef LAL_NDEBUG
    INT4 imax = 40;
    ASSERT (i < imax, status, FINDROOTH_EMXIT, FINDROOTH_MSGEMXIT);
#endif
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


/** \see See \ref FindRoot_h for documentation */
REAL8
XLALDBisectionFindRoot (
    REAL8 (*y)(REAL8, void *),
    REAL8 xmin,
    REAL8 xmax,
    REAL8 xacc,
    void *params
)
{
  const INT4 imax = 80;
  INT4  i;
  REAL8 y_1;
  REAL8 y_2;
  REAL8 xmid;
  REAL8 ymid;

  /* check arguments */
  if(xacc < 0.0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  /* put xmin and xmax in the correct order, using y_1 as temporary storage */
  if(xmin > xmax) {
    y_1 = xmin;
    xmin = xmax;
    xmax = y_1;
  }

  /* evaluate function at endpoints */
  y_1 = y(xmin, params);
  y_2 = y(xmax, params);
  if(XLALIsREAL8FailNaN(y_1) || XLALIsREAL8FailNaN(y_2))
    XLAL_ERROR_REAL8(XLAL_EFUNC);

  /* loop until root found within requested accuracy or iteration limit
   * exceeded */
  for(i = 0; (xmax - xmin) > xacc; i++)
  {
    if(i >= imax)
      XLAL_ERROR_REAL8(XLAL_EMAXITER);

    /* evaluate function at midpoint */
    xmid = (xmin + xmax) / 2.0;
    ymid = y(xmid, params);
    if(XLALIsREAL8FailNaN(ymid))
      XLAL_ERROR_REAL8(XLAL_EFUNC);

    /* did we get lucky? */
    if(ymid == 0.0)
      break;

    if(y_1 * ymid < 0.0)
    {
      /* root is in first half of domain */
      xmax = xmid;
      y_2 = ymid;
    }
    else if(y_2 * ymid < 0.0)
    {
      /* root is in second half of domain */
      xmin = xmid;
      y_1 = ymid;
    }
    else
    {
      /* something's gone wrong */
      XLAL_ERROR_REAL8(XLAL_EFAILED);
    }
  }

  return((xmin + xmax) / 2.0);
}

/** \see See \ref FindRoot_h for documentation */
void
LALDBisectionFindRoot (
    LALStatus      *status,
    REAL8       *root,
    DFindRootIn *input,
    void        *params
    )
{

  INT4  i = 0;
  REAL8 y_1;
  REAL8 y_2;
  REAL8 x;
  REAL8 dx;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* check that arguments are reasonable */
  ASSERT (root, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (input, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  ASSERT (input->function, status, FINDROOTH_ENULL, FINDROOTH_MSGENULL);
  /* params can be NULL ... */

  /* evaluate function at endpoints */

  input->function (status->statusPtr, &y_1, input->xmin, params);
  CHECKSTATUSPTR (status);

  input->function (status->statusPtr, &y_2, input->xmax, params);
  CHECKSTATUSPTR (status);

  ASSERT (y_1*y_2 < 0, status, FINDROOTH_EBRKT, FINDROOTH_MSGEBRKT);

  if (y_1 < 0)
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
#ifndef LAL_NDEBUG
    INT4 imax = 80;
    ASSERT (i < imax, status, FINDROOTH_EMXIT, FINDROOTH_MSGEMXIT);
#endif
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
