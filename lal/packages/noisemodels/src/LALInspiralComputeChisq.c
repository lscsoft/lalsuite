/*  <lalVerbatim file="LALInspiralComputeChisqCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralComputeChisq.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralComputeChisqCP}
\index{\verb&LALInspiralComputeChisq()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALREAL4VectorFFT
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralComputeChisqCV}}
</lalLaTeX>  */

#include <lal/AVFactories.h>
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALCOMPUTECHISQ, "$Id$");

/*  <lalVerbatim file="LALInspiralComputeChisqCP"> */
void
LALInspiralComputeChisq
   (
   LALStatus     *status,
   REAL4         *chisq,
   InspiralChisqDataVec  *input,
   InspiralChisqParams   *params
   )
{  /*  </lalVerbatim>  */

  REAL8 c2, flso, df, duration, totalModelSNR, SNRPerBin, binSNR, totalSNR, diffSNR, mSevenBy3;
  INT4 n, k_0, k_N, /*ki,*/ count, k, k_start, k_end, i;
  UINT4 binIndexLength;
  INT4Vector *binIndexes = NULL;
  REAL8Vector *binSNRs = NULL;

  INITSTATUS (status, "LALInspiralComputeChisq", LALINSPIRALCOMPUTECHISQ);
  ATTATCHSTATUSPTR(status);

/* Need to calculate the total SNR from the stationary phase approximation */



  binIndexLength = params->nBins + 1;
  LALI4CreateVector(status->statusPtr, &binIndexes, binIndexLength);
  LALDCreateVector(status->statusPtr, &binSNRs, binIndexLength-1);
  
  flso = 1.0/(pow(6.0,3.0/2.0)*LAL_PI*params->totalMass);

  duration = 2.L * (input->psd->length-1)*params->deltaT;
  df = 1.L/duration;

  k_0 = floor(params->fLower/df);  /* The lower discrete frequency in summation */
  k_N = floor(flso/df);            /* The upper discrete frequency in summation */

  n = input->SNRIntegrand->length;

  /*
  fprintf(stderr, "%d %d %d\n", n, k_0, k_N);
   */
  
  totalModelSNR = 0.0;
  mSevenBy3 = -7.L/3.L;
  for(k=k_0; k<=k_N; k++)
  {
     if (input->psd->data[k]) totalModelSNR += pow( (float)k*df, mSevenBy3) / input->psd->data[k];
  }

  c2 = SNRPerBin = totalModelSNR/(double)params->nBins;
  /*
  fprintf(stderr, "SNRPerBin=%e\n", SNRPerBin);
  */

  count = 0;
  binIndexes->data[count] = k_0;
  k_start = k_0;

  /* At the end of this loop we should have the INT4Vector binIndexes filled up
   * with the indexes k corresponding to the boundaries between the chi-squared bins
   */

  k = k_start;
  while(count < params->nBins-1)
  {
     binSNR = 0.0;
     for(k=k_start; k<=k_N; k++)
     {
	double dSNR=0., ePlus=0., eMinus=0.;
        if (input->psd->data[k]) binSNR += (dSNR=pow( (float)k*df, mSevenBy3) / input->psd->data[k]);

           if (binSNR > SNRPerBin)
	   {  
		   count++;
		   /* The following is added so that we are as close to the SNRPerBin
		    * as possible; bug noted by S. Babak on April 22, 2003. Also, instead of
		    * using SNRPerBin we shall use the actual SNR in the bin
		    */
		   ePlus = binSNR - SNRPerBin;
		   eMinus = SNRPerBin + dSNR - binSNR;
		   if (ePlus < eMinus)
		   {
			   binIndexes->data[count] = k;
			   binSNRs->data[count-1] = binSNR;
			   k_start = k+1;
			   /*
			   fprintf(stderr, "I%d %d %e %e %e\n", count, k, binSNRs->data[count-1], ePlus, eMinus);
			    */
		   }
		   else
		   {
			   binIndexes->data[count] = k-1;
			   binSNRs->data[count-1] = binSNR-dSNR;
			   k_start = k;
			   /*
			   fprintf(stderr, "E%d %d %e %e %e\n", count, k, binSNRs->data[count-1], ePlus, eMinus);
			   */
		   }
		   break;
	   }
     }
  }
  count++;
  binSNR = 0.;
  for(k=k_start; k<=k_N; k++)
  {
	  if (input->psd->data[k]) binSNR += pow( (float)k*df, mSevenBy3) / input->psd->data[k];
  }
			   
  binIndexes->data[count] = k_N;
  binSNRs->data[count-1] = binSNR;
  /*
  fprintf(stderr, "E%d %d %e\n", count, k, binSNRs->data[count-1]);
  */

  if (count != params->nBins) 
  {
	  fprintf(stderr, "count=%d not same as binIndexLength=%d\n", count, binIndexLength);
  }


  /*
  for(k=1; k < input->psd->length; k++)
  {
     if (input->psd->data[k]) 
	     fprintf(stdout, "%e %e \n", (float)k*df, (pow((float)k * df, mSevenBy3) / input->psd->data[k]));
  }
     fprintf(stdout, "&\n");
   for(i=0; i < binIndexLength; i++)
   {
      fprintf(stdout, "%e %d\n", binIndexes->data[i]*df, i);
   } 
  
   */

  /* Now get the total SNR from the data */

  totalSNR = 0.0;
  for(k=k_0; k<=k_N; k++)
  {
     totalSNR += input->SNRIntegrand->data[k] + input->SNRIntegrand->data[n-k];
  }
  /* Now count the negative frequencies too */
  totalSNR *= 2.;

  /* the expected SNR per chisq bin in the frequency space */

  SNRPerBin = totalSNR / (double)params->nBins;
  c2 /= SNRPerBin;

  /* Now we need to sum up the SNR in each of the bins defined and stored in binIndexes->data[k] 
   * Then we take the difference between this figure and the SNRPerBin calculated immediately
   * above, and use this to compute the chi-squared statistic.
   */

   *chisq = 0.0;
   for(i=0; i<params->nBins; i++)
   {
      k_start = binIndexes->data[i];
      k_end = binIndexes->data[i+1];

      binSNR = 0.0;
      for(k=k_start; k<=k_end; k++)
      {
         binSNR += input->SNRIntegrand->data[k] + input->SNRIntegrand->data[n-k];
      }

      binSNR *= 2.;
      diffSNR = 1. - c2*binSNR/binSNRs->data[i];
      /* diffSNR = 1. - binSNR/SNRPerBin; */
      *chisq += diffSNR * diffSNR;
   }
   /*
   fprintf(stderr, "snr=%e chisq=%e\n", totalSNR, *chisq);
   */
      
  LALI4DestroyVector(status->statusPtr, &binIndexes);
  LALDDestroyVector(status->statusPtr, &binSNRs);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
