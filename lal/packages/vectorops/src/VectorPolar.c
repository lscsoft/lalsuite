/**** <lalVerbatim file="VectorPolarCV">
 * Author: T. D. Creighton, A. M. Sintes
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{VectorPolar.c}}
 *
 * Convert complex vector components from rectangular coordinates to polar
 * coordinates.
 *
 * \subsubsection*{Prototypes}
 * \input{VectorPolarCP}
 * \idx{LALCVectorAbs()}
 * \idx{LALCVectorAngle()}
 * \idx{LALUnwrapREAL4Angle()}
 * \idx{LALZVectorAbs()}
 * \idx{LALZVectorAngle()}
 * \idx{LALUnwrapREAL8Angle()}
 * 
 * \subsubsection*{Description}
 *
 * Let \texttt{u} be an object of type \texttt{COMPLEX8Vector}, and let
 * \texttt{a} and \texttt{b} be objects of type \texttt{REAL4Vector}.
 * 
 * The \verb:LALCVectorAbs( &status, &a, &u ): function computes
 * the magnitude of a complex vector \texttt{u}.
 * $\mbox{\texttt{a.data[i]}}=\mbox{\texttt{sqrt}} (
 * \mbox{\texttt{u.data[i].re}}^2 + \mbox{\texttt{v.data[i].im}}^2 ) $.
 *
 * The \verb:LALCVectorAngle( &status, &a, &u ): function computes
 * the phase angle of a complex vector \texttt{u}
 * in the interval $[-\pi, \pi]$ radians.\\
 * $\mbox{\texttt{a.data[i]}}=\mbox{\texttt{atan2}} (
 * \mbox{\texttt{u.data[i].im}}, \mbox{\texttt{v.data[i].re}} ) $.
 *
 * The \verb:LALUnwrapREAL4Angle( &status, &a, &b ): function
 * corrects the radian phase angles of a real vector  \texttt{b}
 * by adding multiples of
 * $\pm\pi$ when the absolute jumps between consecutive
 * angle elements are greater than $\pi$ radians.
 * This function detects branch cut crossings, but it can be 
 * fooled by sparse, rapidly changing phase values.
 * 
 * The double-precision functions are similar.
 *  
 * \subsubsection*{Algorithm}
 *
 * 
 * The algorithm for LALUnwrapREAL4Angle and LALUnwrapREAL8Angle
 * (Inspired from the MATLAP function unwrap):
 * \begin{verbatim}
 * 
 *   a = in->data;
 *   b = out->data;
 *   n = out->length;
 *   
 *   cumsum = 0.0;
 *   phaseI = *a;
 *   *b = phaseI;
 *   --n;
 * 
 *   while (n-- > 0)
 *   {
 *     ++a;
 *     ++b;
 *     phaseII = *a;
 *     diffph = phaseII - phaseI;
 *     phaseI = phaseII;
 *     
 *     cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
 *     
 *     *b= phaseII + cumsum;
 *   }
 * 
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * For the LALUnwrapREAL4Angle and LALUnwrapREAL8Angle functions, \texttt{a},
 * and \texttt{b} should  not point to the same memory location (\texttt{a !=
 * b}).
 * 
 * \vfill{\footnotesize\input{VectorPolarCV}}
 * 
 **** </lalLaTeX> */ 

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

NRCSID (VECTORPOLARC, "$Id$");
 
/* <lalVerbatim file="VectorPolarCP"> */
void
LALCVectorAbs(
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{ /* </lalVerbatim> */
  COMPLEX8 *a;
  REAL4    *b;
  INT4      n;

  INITSTATUS (status, "LALCVectorAbs", VECTORPOLARC);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  a = in->data;
  b = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 ar = a->re;
    REAL4 ai = a->im;

    if (fabs(ar) > fabs(ai))
    {
      REAL4 ratio = ai/ar;
      *b=fabs(ar)*sqrt(1.0 + ratio*ratio);
    }
    else
    {
      if (fabs(ai) < LAL_REAL4_MIN)
      {
         *b=0.0;   /* to avoid NaN */
      }
      else
      {
         REAL4 ratio = ar/ai;
         *b=fabs(ai)*sqrt(1.0 + ratio*ratio);
      }
    }

    ++a;
    ++b;
  }

  RETURN (status);
}


/* <lalVerbatim file="VectorPolarCP"> */
void
LALZVectorAbs(
    LALStatus             *status,
    REAL8Vector           *out,
    const COMPLEX16Vector *in
    )
{ /* </lalVerbatim> */
  COMPLEX16 *a;
  REAL8    *b;
  INT4      n;

  INITSTATUS (status, "LALZVectorAbs", VECTORPOLARC);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  a = in->data;
  b = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;

    if (fabs(ar) > fabs(ai))
    {
      REAL8 ratio = ai/ar;
      *b=fabs(ar)*sqrt(1.0 + ratio*ratio);
    }
    else
    {
      if (fabs(ai) < LAL_REAL8_MIN)
      {
         *b=0.0;   /* to avoid NaN */
      }
      else
      {
         REAL8 ratio = ar/ai;
         *b=fabs(ai)*sqrt(1.0 + ratio*ratio);
      }
    }

    ++a;
    ++b;
  }

  RETURN (status);
}


/* <lalVerbatim file="VectorPolarCP"> */
void
LALCVectorAngle (
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{ /* </lalVerbatim> */
  COMPLEX8 *a;
  REAL4    *b;
  INT4      n;

  INITSTATUS (status, "LALCVectorAngle", VECTORPOLARC);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  a = in->data;
  b = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL4 ar = a->re;
    REAL4 ai = a->im;
    
    if ((fabs(ar)<LAL_REAL4_MIN) && (fabs(ai)<LAL_REAL4_MIN))
    {
      *b=0.0;   /* to avoid NaN when ar=ai=0.0 */
    }
    else
    {
      *b= atan2(ai,ar);
    }
    

    ++a;
    ++b;
  }

  RETURN (status);
}


/* <lalVerbatim file="VectorPolarCP"> */
void
LALZVectorAngle (
    LALStatus              *status,
    REAL8Vector            *out,
    const COMPLEX16Vector  *in
    )
{ /* </lalVerbatim> */
  COMPLEX16  *a;
  REAL8      *b;
  INT4        n;

  INITSTATUS (status, "LALZVectorAngle", VECTORPOLARC);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  a = in->data;
  b = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;
    
    
    if ((fabs(ar)<LAL_REAL8_MIN) && (fabs(ai)<LAL_REAL8_MIN))
    {
      *b=0.0;   /* to avoid NaN when ar=ai=0.0 */
    }
    else
    {
      *b= atan2(ai,ar);
    }


    ++a;
    ++b;
  }

  RETURN (status);
}


/* <lalVerbatim file="VectorPolarCP"> */
void
LALUnwrapREAL4Angle (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in
    )
{ /* </lalVerbatim> */
  REAL4    *a;
  REAL4    *b;
  INT4      n;
  REAL4     cumsum;
  REAL4     phaseI;
  REAL4     phaseII;
  REAL4     diffph;
  
  INITSTATUS (status, "LALUnwrapREAL4Angle", VECTORPOLARC );
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
          
  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPSH_ESAME, VECTOROPSH_MSGESAME);


  a = in->data;
  b = out->data;
  n = out->length;
  
  cumsum = 0.0;
  phaseI = *a;
  *b = phaseI;
  --n;

  while (n-- > 0)
  {
    ++a;
    ++b;
    phaseII = *a;
    diffph = phaseII - phaseI;
    phaseI = phaseII;
    
    cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
    
    *b= phaseII + cumsum;
  }

  RETURN (status);
}


/* <lalVerbatim file="VectorPolarCP"> */
void
LALUnwrapREAL8Angle (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in
    )
{ /* </lalVerbatim> */
  REAL8    *a;
  REAL8    *b;
  INT4      n;
  REAL8     cumsum;
  REAL8     phaseI;
  REAL8     phaseII;
  REAL8     diffph;
  
  INITSTATUS (status, "LALUnwrapREAL8Angle", VECTORPOLARC );
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
          
  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPSH_ESAME, VECTOROPSH_MSGESAME);


  a = in->data;
  b = out->data;
  n = out->length;
  
  cumsum = 0.0;
  phaseI = *a;
  *b = phaseI;
  --n;

  while (n-- > 0)
  {
    ++a;
    ++b;
    phaseII = *a;
    diffph = phaseII - phaseI;
    phaseI = phaseII;
    
    cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
    
    *b= phaseII + cumsum;
  }

  RETURN (status);
}

