/*  <lalVerbatim file="LALInspiralTiming2CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTiming2.c}}

Module used in solving the timing and phasing functions in quadrature for the
{\tt Approximant TaylorT2}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralTiming2CP}
\index{\verb&LALInspiralTiming2()&}

\subsubsection*{Description}

Given $t$ and $v$ this module computes the quantity
\begin{equation}
{\tt tofv} = t - t_C - t_N(v) \sum t_k v^k,
\end{equation}
where the coefficients $t_k$ and the Newtonian value $t_N$ are all defined
in Table~\ref{table:flux}.

\subsubsection*{Algorithm}
None


\subsubsection*{Uses}
None

\subsubsection*{Notes}
None


\vfill{\footnotesize\input{LALInspiralTiming2CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALTIMING2C, "$Id$");

/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_0PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v8 = pow(v,8.);

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_2PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v8 = v2*v2*v2*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2);


  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_3PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v8 = v3*v3*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_4PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_5PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_6PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v8 = v6*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(v)) * v6);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void 
LALInspiralTiming2_7PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   ) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v7, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM*f, oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  v8 = v7*v;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(v)) * v6
        + toffIn->t7 * v7);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
