/*----------------------------------------------------------------------- 
 * 
 * File Name: HarmonicFinder.c
 * 
 * Author: Sintes, A. M.
 * 
 * Revision: $Id: 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 *  HarmonicFinder.c
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 *   Given certain harmonic indices (k)  finds the frequency interval
 *   locations (in bins) of the interference (around k*fLine)
 * 
 * DIAGNOSTICS 
 *
 * CALLS
 * 
 * NOTES
 * 
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="HarmonicFinderCV">
Author: Sintes, A. M. 
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{HarmonicFinder.c}}
\label{ss:HarmonicFinder.c}
 Given certain harmonic indices $\{ k\} $  finds the frequency interval
   location (in bins) of the interference (around $k\times f_0$).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{HarmonicFinderD}
\idx{LALHarmonicFinder()}

\subsubsection*{Description}
This routine determines the lower and upper frequency limit (in bins)
of each harmonic line considered, $(\nu_{ik}, \nu_{fk})$, from the power
spectrum.

The harmonic indices  are given as an input \verb@*in1@.
\begin{description}
\item[\texttt{in1->length}] Number of harmonics.
\item[\texttt{in1->data}]   List of harmonics to consider, e.g., 
       $\{ k\} $ =  $\{ 3, 5, 9, 11 \ldots\}$
\end{description}

The power spectrum, $\vert\tilde x(\nu)\vert^2$, together with the approximate
frequency $f_0$ (in Hz) of the interference fundamental harmonic and the 
frequency
resolution are also given as an input  \verb@*in2@.
\begin{description}
\item[\texttt{in2->length}] The number of elements in \texttt{in2->data}.
\item[\texttt{in2->data}]   The spectrum,  $\vert\tilde x(\nu)\vert^2$.
\item[\texttt{in2->deltaF}] The $\Delta$F offset between samples (in Hz).
\item[\texttt{in2->fLine}]  The interference fundamental frequency $f_0$ 
              (in Hz), e.g., 60 Hz.
\end{description}

The  output  \verb@*out@ is a vector whose length is 
\verb@out->length@ = $3\times$~\verb@in1->length@,
 and contains  for each considered harmonic, in the following order,
 its index $k$ and the bin location of $\nu_{ik}$ and  $\nu_{fk}$.
\begin{description}
\item[\texttt{out->length}] The number of elements in \texttt{out->data}.
\item[\texttt{out->data}]    $\{ k,\nu_{ik}, \nu_{fk} \} $,
      e.g.,  $\{$3, 9868, 9894, 5, 16449, 16487, 9, 29607, 29675$\ldots\}$.
\end{description}




\subsubsection*{Algorithm}

It looks for the location of interference harmonics assuming that
the fundamental harmonic is located somewhere in the interval
\verb@in2->fLine@$ - 0.7$ Hz and \verb@in2->fLine@$+ 0.7$ Hz.
First, the power spectrum is smoothed by averaging neighboring bins.
Then, the corresponding frequency  intervals
of the harmonics considered are reduced by finding the central
bin position of the lines and their standard deviation. This is done
using the smooth power spectrum as a probability density distribution.
The limits of the lines are set initially  at 1 or 2 sigma from the
central bin location and, later, they are moved until they hit a local
 minimum in a selected interval.
See the code for details.

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVector()
LALSDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{HarmonicFinderCV}}

</lalLaTeX> */



#include <math.h>
#include <lal/LALConstants.h>
#include <lal/CLR.h>

#define MIN(A, B)       ((A) < (B) ? (A) : (B))
#define MAX(A, B)       ((A) > (B) ? (A) : (B))

#define log2( x )       ( log( x ) / LAL_LN2 )

NRCSID (HARMONICFINDERC, "$Id$");


/* <lalVerbatim file="HarmonicFinderD"> */
void LALHarmonicFinder (LALStatus  *status,
         INT4Vector         *out,   /* harmonic index and location, size 3*l */
         REAL4FVectorCLR    *in2,   /* |x(f)|^2, data + information */
         INT4Vector         *in1)   /* the harmonic index, size l */
{ /* </lalVerbatim> */

  INT4    n,l;
  INT4    k,binini,binfin;
  INT4    nInf1s,nInf2s,nSup1s,nSup2s;
  INT4    nBins;
  INT4    *kv;
  INT4    *kf;

  INT4    i,j;

  REAL4   devF,fL,myfmax;
  REAL4   Tobs;
  
  REAL8   sumpx,mean1,std1,mn2,sn2;
  REAL4   llindar,cc,invk,norma;
  
  REAL4   sInf1,sInf2,sSup1,sSup2;

  REAL4        *px;
  REAL4Vector  *pxs  = NULL;
    
/* --------------------------------------------- */

  INITSTATUS (status, "LALHarmonicFinder", HARMONICFINDERC);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (out, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in2, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure the data pointers are not NULL: */ 
  ASSERT (out->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in2->data, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure that the size is correct:  */
  ASSERT (out->length > 0, status, CLRH_ESIZE, CLRH_MSGESIZE);
  ASSERT (in2->length > 2, status, CLRH_ESIZE, CLRH_MSGESIZE);

  /*   Make sure that the lengths are correct (size mismatch): */  
  ASSERT (3 *in1->length == out->length, status, CLRH_ESZMM, CLRH_MSGESZMM);


   /*  Make sure F_max > fLine  */ 
  ASSERT (fabs(1.01*  in2->fLine) <  fabs(in2->deltaF *(in2->length - 1) ), 
               status, CLRH_EFREQ, CLRH_MSGEFREQ); 

  /*   Make sure deltaF!= 0 */
  ASSERT (fabs(in2->deltaF) != 0, status, CLRH_EFREQ, CLRH_MSGEFREQ); 

  /* -------------------------------------------   */
 
  devF = 0.7;
  fL = fabs( in2->fLine);
  
  px = in2->data;
  n  = in2->length;

  l  = in1->length;
  kv = in1->data;
  kf = out->data;

  Tobs = fabs(1.0 / in2->deltaF );
  myfmax = fabs(in2->deltaF) *(n - 1.0);
 
  /* create extra vectors */
  TRY(LALSCreateVector(status->statusPtr, &pxs,  n), status);
 
 /* -------------------------------------------   */
  /* smoothing the spectrum by averaging neighboring bins */
  /* ignoring border effects */

  nBins = 1;
  if (Tobs > 1.0 )
    { nBins = ceil ( log2(Tobs)) ;  }
  
  /* a normalization to control overflow */
  norma = 1.0/( (2.0*n -2.0)*(2.0*nBins +1.0)  );
  
  /* Border effects */
  sumpx = 0.0;
  for ( i=0; i<= 2* nBins; ++i)
    { sumpx += px[i];}

  for (i=0; i<= nBins; ++i )
    { pxs->data[i] = sumpx*norma;  }
  
  /* intermediate values */
  for (i=nBins+1; i< n-nBins; ++i)
    { pxs->data[i] = pxs->data[i-1] + (px[i+nBins] - px[i-nBins-1])*norma; }
  
  /* the other border */
  for (i=n-nBins; i<n; ++i)
    { pxs->data[i] = pxs->data[i-1];}
  /* -------------------------------------------   */
  
  /* the first line */
  
  /* k = fabs(*kv);  */
  
  k= kv[0];
  if ( k < 0 ) { k = -k; }
  
  binini = floor( k*( fL - devF)*Tobs );
  binfin = ceil(  k*( fL + devF)*Tobs );
  
  ASSERT ( k != 0, status, CLRH_EFREQ, CLRH_MSGEFREQ); 
  ASSERT ( binfin <  n, status, CLRH_EFREQ, CLRH_MSGEFREQ); 
  ASSERT ( fL > devF, status, CLRH_EFREQ, CLRH_MSGEFREQ); 

  sumpx = 0.0;
  mean1 = 0.0;
  std1  = 0.0;
  for (i=binini; i<= binfin; ++i) {
    sumpx += pxs->data[i]; /* sum of the distribution, normalize it */
    mean1 += i*pxs->data[i]; 
  }
  mean1 = mean1/sumpx;   /* central point of the first line considered */
  
  /* separated loops due to precision problems */
  for (i=binini; i<= binfin; ++i) {
    std1  += (i-mean1)*(i-mean1)*pxs->data[i]; 
  }
  std1 = sqrt( std1/sumpx );      /* std of the line */
  
  invk = 1.0/k;
  
  /* -------------------------------------------   */
  
  /* for all the  lines */
  llindar = 1.1;
  cc = 4.0;
  
  for (j=0; j<l; ++j) {
    /* k = fabs( kv[j] ); */ 
    
    k = kv[j];
    if ( k < 0 ) { k = -k; }
    binini = floor( k*(mean1 - std1*cc)*invk );
    binfin = ceil ( k*(mean1 + std1*cc)*invk );
    ASSERT ( binfin <  n, status, CLRH_EFREQ, CLRH_MSGEFREQ); 

    kf[3*j] = k;
    
    sumpx = 0.0;
    mn2 = 0.0;
    sn2 = 0.0;
    
    for (i=binini; i<= binfin; ++i) {
      sumpx += pxs->data[i]; /* sum of the distribution, normalize it */
      mn2 += i*pxs->data[i]; 
    }
    mn2 = mn2/sumpx;   /* central point of the first line considered */
    
    
    /* separated loops due to precision problems */
    for (i=binini; i<= binfin; ++i) {
      sn2 += (i- mn2)*(i- mn2)*pxs->data[i]; 
    }
    
    sn2 = sqrt( sn2/sumpx );      /* std of the line */
    
    /* initial settings for the limits */
    nInf1s= MAX(binini, floor(mn2-  sn2 ) );
    nInf2s= MAX(binini, floor(mn2-2*sn2 ) );
    
    nSup1s= MIN(binfin, ceil(mn2+  sn2 ) );
    nSup2s= MIN(binfin, ceil(mn2+2*sn2 ) );
    
    /* shifting from 1 sigma to 2 sigma if there is a big difference */
    /* and setting closer intervals where there is no big difference */
    sInf1 = 0.0;
    sInf2 = 0.0;
    sSup1 = 0.0;
    sSup2 = 0.0;
   
    for (i=0; i<= nBins; ++i) {
      sInf1 += pxs->data[nInf1s -i];
      sInf2 += pxs->data[nInf2s -i];
      sSup1 += pxs->data[nSup1s +i];
      sSup2 += pxs->data[nSup2s +i];
    }  
    
    if (sInf1/sInf2 > llindar ) 
      {  nInf1s = nInf2s;  }
    else
      {  binini = nInf2s; } 
    
    
    if (sSup1/sSup2 > llindar)
      {  nSup1s = nSup2s; }
    else
      {  binfin = nSup2s; }
    
    /* move to the next local minimum in the selected interval*/
    
    i=nInf1s;
    while ( i>binini ) {
      if ( pxs->data[i] < pxs->data[i-1])
	break;
      --i;
    }
    nInf1s = i;
    kf[3*j+1] = nInf1s;
    
    i=nSup1s;
    while ( i<binfin ) {
      if ( pxs->data[i] < pxs->data[i+1])
	break;
      ++i;
    }
    nSup1s = i;
    kf[3*j+2] = nSup1s;
    
  } /* closes loop of all lines */
  
  /* -------------------------------------------   */
  
  /* Destroy Vectors  */
  TRY(LALSDestroyVector(status->statusPtr, &pxs), status);
  
  /* -------------------------------------------   */
  
  DETATCHSTATUSPTR (status);
  
  /* normal exit */
  RETURN (status);
}

#undef MIN
#undef MAX
