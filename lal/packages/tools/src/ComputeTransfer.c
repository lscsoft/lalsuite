/********************************* <lalVerbatim file="ComputeTransferCV">
Author: Patrick Brady
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{ComputeTransfer.c}}
\label{ss:ComputeTransfer.c}

Computes the transfer function from zero-pole-gain representation.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComputeTransferCP}
\idx{LALComputeTransfer()}

\subsubsection*{Description}

A transfer function can either be specified as a list of coefficients or a 
list of poles and zeros. The function \verb@LALComputeTransfer()@ computes the 
frequency representation of the transfer function \verb+calrec->transfer+
described by the zeroes,  poles and gain in \verb@*calrec@.   The memory for 
the frequency series should be allocated before calling this routine which uses
\verb+calrec->transfer->deltaF+ and \verb+calrec->transfer->data->npoints+.

\subsubsection*{Algorithm}

The transfer function is deduced from the poles and zeros as follows:
\begin{equation}
T(f) = c_{\mathrm{gain}} 
{\prod_i \textrm{zero}(f,z_i)}{ \prod_{i} \textrm{pole}(f,p_i)}
\end{equation}
where 
\begin{equation}
\textrm{zero}(f,z) = \left\{ \begin{array}{ll}
(i f / z) + 1 & \textrm{ when } z\neq 0 \\
i f & \textrm{ when } z = 0
\end{array}
\right.
\end{equation}
and 
\begin{equation}
\textrm{pole}(f,p) = \left\{ \begin{array}{ll}
\frac{1}{(i f / p) + 1} & \textrm{ when } p \neq 0 \\
\frac{1}{i f} & \textrm{ when } p = 0
\end{array}
\right.
\end{equation}
For both the transfer function and the pole-zero notation the units for
frequency is Hz rather than rad/s (angular frequency).  In particular, poles
and zeros are specified by their location in frequency space

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}
The DC component of \verb+calrec->transfer+ is always filled with $1 + i 0$.  
In most cases,  this should be irrelevant for gravitational wave data analysis,  
but care should be taken if DC is relevant when this function is used.  

\vfill{\footnotesize\input{ComputeTransferCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>

NRCSID( COMPUTETRANSFERC, "$Id$" );

static void product(COMPLEX8 *c,COMPLEX8 *a, COMPLEX8 *b) {
  
  c->re = a->re * b->re - a->im * b->im;
  c->im = a->re * b->im + a->im * b->re;
  
  return;
}

static void ratio(COMPLEX8 *c,COMPLEX8 *a, COMPLEX8 *b) {
  REAL4 norm;

  norm = b->re * b->re + b->im * b->im;
  
  c->re = (a->re * b->re + a->im * b->im)/norm;
  c->im = (- a->re * b->im + a->im * b->re)/norm;
  
  return;
}

/* <lalVerbatim file="ComputeTransferCP"> */
void
LALComputeTransfer( LALStatus                 *stat,
                    CalibrationRecord         *calrec
                    )
/* </lalVerbatim> */
{ 
  UINT4         i, j;                    /* indexes               */
  UINT4         jmin;                    /* used to handle DC     */
  COMPLEX8      dummy, dummyC, factor;   /* dummy complex numbers */
  REAL4         f,df;                    /* freq and interval     */
  REAL8         norm;

  INITSTATUS( stat, "LALComputeTransfer", COMPUTETRANSFERC );
  ATTATCHSTATUSPTR (stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT( calrec, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->zeros, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->poles, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );  
  ASSERT( calrec->transfer, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );  

  /* initialise everything */ 
  dummyC.re = factor.re = 1.0;
  dummyC.im = factor.im = 0.0;
  df = calrec->transfer->deltaF;
  jmin = 0;

  /* compute the normalization constant */
  norm = calrec->gain;
  
  for ( i=0 ; i < calrec->poles->length ; i++)
    if ( calrec->poles->data[i] != 0 )
    {
      norm *= calrec->poles->data[i];
    }

  for ( i=0 ; i < calrec->zeros->length ; i++)
    if ( calrec->zeros->data[i] != 0 )
    {
      norm /= calrec->zeros->data[i];
    }

  /* Handle DC if necessary */
  if ( calrec->transfer->f0 == 0.0 )
  {
    calrec->transfer->data->data[0].re = 1.0;
    calrec->transfer->data->data[0].im = 0.0;
    jmin = 1;
  }

  /* loop over frequency in the output */
  for ( j=jmin ; j<calrec->transfer->data->length ; j++)
  {
    /* the frequency */
    f = calrec->transfer->f0 + (REAL4) j * df;
    dummyC.re=1.0;
    dummyC.im=0.0;

    /* loop over zeroes */
    for (i = 0 ; i < calrec->zeros->length ; i++)
    {
      factor.re = (calrec->zeros->data[i]);
      factor.im = f;
      product( &dummy, &dummyC, &factor);
      dummyC.re=dummy.re;
      dummyC.im=dummy.im;
    }

    /* loop over poles */
    for (i = 0 ; i < calrec->poles->length ; i++)
    {
      factor.re = (calrec->poles->data[i]);
      factor.im = f;
      ratio( &dummy, &dummyC, &factor);
      dummyC.re=dummy.re;
      dummyC.im=dummy.im;
    }

    /* fill the frequency series */
    calrec->transfer->data->data[j].re = norm * dummyC.re;
    calrec->transfer->data->data[j].im = norm * dummyC.im;
  }

 
  /* we're out of here */
  DETATCHSTATUSPTR (stat);
  RETURN( stat );
}
