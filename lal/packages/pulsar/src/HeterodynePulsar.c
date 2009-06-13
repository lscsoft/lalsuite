/*
*  Copyright (C) 2007 Jolien Creighton, Matt Pitkin
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

/************************************ <lalVerbatim file="HeterodynePulsarCV">
Author: Dupuis, R. J.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{HeterodynePulsar.c}}

The routines in this module heterodyne, average and resample a time series for known pulsars.

\subsubsection*{Prototypes}
\input{HeterodynePulsarCP}
\idx{LALCoarseHeterodyne()}
\idx{LALFineHeterodyneToPulsar()}

\subsubsection*{Description}

\noindent The function \texttt{LALCoarseHeterodyne()}  ...
\newline

\noindent The function \texttt{LALFineHeterodyneToPulsar()} ...

\subsubsection*{Algorithm}

To be completed.

\subsubsection*{Uses}
\begin{verbatim}
LALSIIRFilter()
LALCCreateVector()
LALCDestroyVector()
LALBarycenterEarth()
LALBarycenter()
\end{verbatim}


\subsubsection*{Notes}

\vfill{\footnotesize\input{HeterodynePulsarCV}}

******************************************************* </lalLaTeX> */

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/HeterodynePulsar.h>
#include <lal/BinaryPulsarTiming.h>

#include <lal/LALConstants.h>
/******* DEFINE RCS ID STRING ************/
NRCSID( HETERODYNEPULSARC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

/******* DECLARE LOCAL (static) FUNCTIONS ************/
/* (definitions can go here or at the end of the file) */

/******* DEFINE GLOBAL FUNCTIONS ************/

/* <lalVerbatim file="HeterodynePulsarCP"> */
void
LALCoarseHeterodyne       ( LALStatus                      *status,
		  	    CoarseHeterodyneOutput   	   *output,
		   	    CoarseHeterodyneInput	   *input,
		   	    CoarseHeterodyneParams	   *params )
/* </lalVerbatim> */
{
  REAL8   		t0;
  REAL8			f0;
  UINT4 		i,j,k,n;
  UINT4			nbox, niir;
  REAL8 		t;
  COMPLEX8Vector    	*Vh;
  COMPLEX8Vector        *Vbox;
  COMPLEX16 		sum;
  REAL8 phase;


  INITSTATUS( status, "LALCoarseHeterodyne", HETERODYNEPULSARC );
  ATTATCHSTATUSPTR (status);

  /******* CHECK VALIDITY OF ARGUMENTS  ************/

 ASSERT(input != (CoarseHeterodyneInput *)NULL, status,
	HETERODYNEPULSARH_ENULLINPUT, HETERODYNEPULSARH_MSGENULLINPUT);

 ASSERT(output != (CoarseHeterodyneOutput *)NULL, status,
	HETERODYNEPULSARH_ENULLOUTPUT, HETERODYNEPULSARH_MSGENULLOUTPUT);

 ASSERT(params != (CoarseHeterodyneParams *)NULL, status,
	HETERODYNEPULSARH_ENULLPARAMS, HETERODYNEPULSARH_MSGENULLPARAMS);

 ASSERT((params->boxM > 0 && params->boxM  <= input->V.data->length), status,
	HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);

 ASSERT((params->iirM > 0 && params->iirM  <= input->V.data->length/params->boxM), status,
	HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);

 if (fmod(input->V.data->length, (params->boxM*params->iirM)) != 0.0)
 {
   ABORT(status, HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);
 }

  ASSERT(!(input->f0 < 0), status,
        HETERODYNEPULSARH_EINVALIDF0, HETERODYNEPULSARH_MSGEINVALIDF0);


  /******* EXTRACT INPUTS AND PARAMETERS ************/

   n = input->V.data->length;

   nbox = n / params->boxM;
   niir = nbox / params->iirM;

  /******* ALLOCATE MEMORY *************/

   Vh = NULL;
   LALCCreateVector( status->statusPtr, &Vh, n);

   Vbox = NULL;
   LALCCreateVector( status->statusPtr, &Vbox, nbox);

  /******* DO ANALYSIS ************/

  f0 = (REAL8)input->f0;
  t0 = (REAL8)input->V.epoch.gpsSeconds + (REAL8)input->V.epoch.gpsNanoSeconds*1.e-9;

  /* calculate initial phase */
  if (f0 == 0)
  {
    phase = 0.0;
  }
  else
  {
   phase = 2.0*LAL_PI*t0*f0;
   phase = fmod(phase, 2.0*LAL_PI);
  }

  /* heterodyne the data to f0 and apply 1st IIR filter */
  for (i = 0; i < n; i++)
  {
    t = (REAL8)i * input->V.deltaT;
    Vh->data[i].re = input->V.data->data[i]* cos(-2.0*(REAL8)LAL_PI*f0*t - phase);
    Vh->data[i].re = LALSIIRFilter( Vh->data[i].re, params->iirFilter1Re );

    Vh->data[i].im = input->V.data->data[i]* sin(-2.0*(REAL8)LAL_PI*f0*t - phase);
    Vh->data[i].im = LALSIIRFilter( Vh->data[i].im, params->iirFilter1Im );
  }

   /* average (apply boxcar) and resample */

  k=0;
  for (i = 0; i < n; i += params->boxM)
  {
    sum.re = 0.0;
    sum.im = 0.0;

    for (j = i; j < i + params->boxM; j++)
    {
      sum.re += (REAL8)Vh->data[j].re;
      sum.im += (REAL8)Vh->data[j].im;
    }

    Vbox->data[k].re = (REAL4)sum.re / (REAL4) params->boxM;
    Vbox->data[k].im = (REAL4)sum.im / (REAL4) params->boxM;

    k++;
  }

  /* apply second iir filter */

  for (i=0; i<nbox;i++)
  {
    Vbox->data[i].re = LALSIIRFilter( Vbox->data[i].re, params->iirFilter2Re );
    Vbox->data[i].im = LALSIIRFilter( Vbox->data[i].im, params->iirFilter2Im );
  }

  /* resample without boxcar filter for sharper edge */
  k=0;
  for (i=0; i < nbox; i += params->iirM)
  {
    output->Vh.data->data[k].re = Vbox->data[i].re;
    output->Vh.data->data[k].im = Vbox->data[i].im;
    k++;
  }

  /* calculate average */
  if (params->stats == 1 || params->stats == 2)
  {
    sum.re = 0.0;
    sum.im = 0.0;

     for (i=0;i<niir;i++)
     {
       sum.re += (REAL8)output->Vh.data->data[i].re;
       sum.im += (REAL8)output->Vh.data->data[i].im;
     }
     output->avg.re = sum.re / (REAL8) niir;
     output->avg.im = sum.im / (REAL8) niir;
   }
   else
   {
     output->avg.re = 0.0;
     output->avg.im = 0.0;
   }

 /* calculate stdev */
  if (params->stats == 1  || params->stats == 2)
  {
    sum.re = 0.0;
    sum.im = 0.0;

    for (i=0;i<niir;i++)
    {
      sum.re += ((REAL8)output->Vh.data->data[i].re - output->avg.re)*((REAL8)output->Vh.data->data[i].re - output->avg.re);
      sum.im += ((REAL8)output->Vh.data->data[i].im - output->avg.im)*((REAL8)output->Vh.data->data[i].im - output->avg.im);
    }

    output->varh.re = sum.re / (REAL8) niir;
    output->varh.im = sum.im / (REAL8) niir;
  }
  else
  {
    output->varh.re = 0.0;
    output->varh.im = 0.0;
  }

   /* calculate skewness */
   if (params->stats == 2)
   {
     sum.re = 0.0;
     sum.im = 0.0;

     for (i=0;i<niir;i++)
      {
	sum.re += ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re)
	* ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re)
	* ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re);

        sum.im += ((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im)
	*((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im)
	*((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im);

      }

      output->skew.re = sum.re /
        (pow((REAL8)output->varh.re, 1.5)*((REAL8)niir - 1.0));
      output->skew.im = sum.im /
        (pow((REAL8)output->varh.im, 1.5)*((REAL8)niir - 1.0));
   }
   else
   {
     output->skew.re = 0.0;
     output->skew.im = 0.0;
   }

   /* calculate excess kurtosis */
   if (params->stats == 2)
   {
     sum.re = 0.0;
     sum.im = 0.0;
     for (i=0;i<niir;i++)
     {

       sum.re += ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re)
               * ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re)
	       * ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re)
	       * ((REAL8)output->Vh.data->data[i].re - (REAL8)output->avg.re);

       sum.im += ((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im)
	       * ((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im)
	       * ((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im)
	       * ((REAL8)output->Vh.data->data[i].im - (REAL8)output->avg.im);
     }

     output->kurt.re = sum.re*pow((REAL8)niir, 4.0)
     / (pow((REAL8)output->varh.re, 2.0) * pow(((REAL8)niir - 1.0), 5.0)) - 3.0;

     output->kurt.im = sum.im*pow((REAL8)niir, 4.0)
     / (pow((REAL8)output->varh.im, 2.0) * pow(((REAL8)niir - 1.0), 5.0)) - 3.0;
   }
   else
   {
     output->kurt.re = 0.0;
     output->kurt.im = 0.0;
   }


   /* calculate the first covariance term assuming mean = 0 */
   if (params->stats == 2)
   {
     sum.re = 0.0;
     sum.im = 0.0;

     for (i=0;i<niir;i++)
     {
        if (i+1 < niir)
	{
          sum.re += (REAL8)output->Vh.data->data[i].re*(REAL8)output->Vh.data->data[i+1].re;
          sum.im += (REAL8)output->Vh.data->data[i].im*(REAL8)output->Vh.data->data[i+1].im;
	}
	else /* just use first point as niir+1 */
	{
	  sum.re += (REAL8)output->Vh.data->data[i].re*(REAL8)output->Vh.data->data[0].re;
          sum.im += (REAL8)output->Vh.data->data[i].im*(REAL8)output->Vh.data->data[0].im;
	}
     }

     output->covar.re = sum.re / (REAL8) niir;
     output->covar.im = sum.im / (REAL8) niir;
   }
   else
   {
     output->covar.re = 0.0;
     output->covar.im = 0.0;
   }

  output->Vh.data->length = niir;

  output->Vh.epoch = input->V.epoch;

  output->Vh.deltaT = input->V.deltaT * params->boxM * params->iirM;

  output->Vh.f0 = f0;

  /****** CLEAN UP *******/

  LALCDestroyVector(status->statusPtr, &Vh);
  LALCDestroyVector(status->statusPtr, &Vbox);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*****************************************************************************/

/* <lalVerbatim file="HeterodynePulsarCP"> */
void
LALFineHeterodyneToPulsar ( LALStatus                      *status,
		  	    FineHeterodyneOutput	   *output,
		   	    FineHeterodyneInput 	   *input,
		   	    FineHeterodyneParams	   *params )
/* </lalVerbatim> */
{
  UINT4		n;
  UINT4		i;
  UINT4		npoints;
  REAL4Vector   *phase;
  REAL4		f0;
  REAL4 	f1;
  REAL4 	f2;
  REAL8 	t,T0;
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime  emit;
  REAL4		deltaT;
  COMPLEX16Vector   	*Xh;
  UINT4	k,j;
  REAL8 tdt,ph, ph0;
  COMPLEX16 sum, sumvar;
  /* REAL8 binaryDeltaT; */

  BinaryPulsarParams orbit;
  /*BinaryPulsarTiming outputBinary;*/

  INITSTATUS( status, "LALFineHeterodyneToPulsar", HETERODYNEPULSARC );
  ATTATCHSTATUSPTR (status);

 /******* CHECK VALIDITY OF ARGUMENTS  ************/

  ASSERT(input != (FineHeterodyneInput *)NULL, status,
	HETERODYNEPULSARH_ENULLINPUT, HETERODYNEPULSARH_MSGENULLINPUT);

  ASSERT(output != (FineHeterodyneOutput *)NULL, status,
	HETERODYNEPULSARH_ENULLOUTPUT, HETERODYNEPULSARH_MSGENULLOUTPUT);

  ASSERT(params != (FineHeterodyneParams *)NULL, status,
	HETERODYNEPULSARH_ENULLPARAMS, HETERODYNEPULSARH_MSGENULLPARAMS);

  ASSERT((params->M > 0 && params->M <= input->Vh.data->length), status,
	HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);

  ASSERT(input->Vh.data->length == input->varh.data->length, status,
	 HETERODYNEPULSARH_ELENGTH, HETERODYNEPULSARH_MSGELENGTH);

 if (fmod(input->Vh.data->length, params->M) != 0.0)
 {
   ABORT(status, HETERODYNEPULSARH_ERFACTOR, HETERODYNEPULSARH_MSGERFACTOR);
 }

 /* since binary model not yet implemented, check to make sure it is an isolated pulsar */
  ASSERT(input->model == 0, status, HETERODYNEPULSARH_EBINARY, HETERODYNEPULSARH_MSGEBINARY);


  n = input->Vh.data->length;
  npoints = n / params->M;

 /******* EXTRACT INPUTS AND PARAMETERS ************/

  deltaT = input->Vh.deltaT;

  f0 = input->f0;
  f1 = input->f1;
  f2 = input->f2;

  /* JC: THIS BETTER NOT DO ANYTHING ... BECAUSE IT IS A SYNTAX ERROR! */
  /* orbit.model = input->model; */

 /* if (orbit.model != 0)
  {
    orbit.Pb = input->Pb;
    orbit.e = input->e;
    orbit.w = input->w;
    orbit.x = input->x;
    orbit.T0 = input->T0;
  }
  */
  baryinput.site.location[0] = params->detector.location[0]/LAL_C_SI;
  baryinput.site.location[1] = params->detector.location[1]/LAL_C_SI;
  baryinput.site.location[2] = params->detector.location[2]/LAL_C_SI;

 /* corrections for proper motion */
  baryinput.tgps = input->Vh.epoch;
  baryinput.alpha = input->source.longitude +((REAL8)input->Vh.epoch.gpsSeconds - input->posEpochGPS)*input->pmRA/LAL_YRSID_SI;
  baryinput.delta = input->source.latitude + ((REAL8)input->Vh.epoch.gpsSeconds - input->posEpochGPS)*input->pmDEC/LAL_YRSID_SI;

  /******* ALLOCATE MEMORY *************/

  phase = NULL;
  LALCreateVector(status->statusPtr, &phase, n);

  Xh = NULL;
  LALZCreateVector(status->statusPtr, &Xh, n);

 /******* DO ANALYSIS ************/

 /* calculate instantaneous phase */
  if (orbit.model == 0) /* isolated pulsar */
  {
    for (i=0;i<n;i++)
    {
      t = (REAL8)input->Vh.epoch.gpsSeconds + (REAL8)input->Vh.epoch.gpsNanoSeconds*1e-9 + (REAL8)i*deltaT;
      baryinput.tgps.gpsSeconds = (INT4)floor(t);
      baryinput.tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
      LALBarycenterEarth(status->statusPtr, &earth, &baryinput.tgps, params->edat);
      LALBarycenter(status->statusPtr, &emit, &baryinput, &earth);

      T0 =input->fEpochGPS;

      tdt = t + emit.deltaT - T0;

      ph0 = input->Vh.f0*t;

      ph = f0*tdt + 0.5*f1*tdt*tdt + f2*tdt*tdt*tdt/6.0;

      phase->data[i] =  2.0*LAL_PI*fmod((ph-ph0), 1.0);
    }
  }
  else /* binary pulsar */
  {
  /*** NOT COMPLETED ***
      for (i=0;i<n;i++)
    {
      t = (REAL8)input->Vh.epoch.gpsSeconds + (REAL8)input->Vh.epoch.gpsNanoSeconds*1e-9 + (REAL8)i*deltaT;
      baryinput.tgps.gpsSeconds = (INT4)floor(t);
      baryinput.tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
      LALBarycenterEarth(status->statusPtr, &earth, &baryinput.tgps, params->edat);
      LALBarycenter(status->statusPtr, &emit, &baryinput, &earth);

      T0 =input->fEpochGPS;
     LALBinaryPulsarDeltaT(status->statusPtr, &outputBinary, (t+emit.deltaT),&orbit);
     tdt = t + emit.deltaT +outputBinary.deltaT  - T0;


      ph0 = input->Vh.f0*t;

      ph = f0*tdt + 0.5*f1*tdt*tdt + f2*tdt*tdt*tdt/6.0;

      phase->data[i] =  2.0*LAL_PI*fmod((ph-ph0), 1.0);
    }
  *****/

  }


  /* do heterodyning, iir filtering, and resampling */
 for (i = 0; i < n; i++)
 {
  Xh->data[i].re = (REAL8)input->Vh.data->data[i].re* cos(-phase->data[i])
		    - (REAL8)input->Vh.data->data[i].im*sin(-phase->data[i]);

  if (params->iirFlag == 1)
    Xh->data[i].re = LALSIIRFilter( Xh->data[i].re, params->iirFilterRe );

  Xh->data[i].im = (REAL8)input->Vh.data->data[i].re * sin(-phase->data[i])
		    + (REAL8)input->Vh.data->data[i].im * cos(-phase->data[i]);

  if (params->iirFlag == 1)
    Xh->data[i].im = LALSIIRFilter( Xh->data[i].im, params->iirFilterIm );
 }

 k=0;
 for (i = 0; i < n; i+=params->M)
 {
   sum.re = 0.0;
   sum.im = 0.0;

   for (j=i;j<i+params->M;j++)
   {
     sum.re += Xh->data[j].re;
     sum.im += Xh->data[j].im;
   }

   output->B.data->data[k].re = sum.re / (REAL4)params->M;
   output->B.data->data[k].im = sum.im / (REAL4)params->M;

   sumvar.re = 0.0;
   sumvar.im = 0.0;

   for (j=i;j<i+params->M;j++)
   {
     sumvar.re += 1.0/ input->varh.data->data[i].re;
     sumvar.im += 1.0 / input->varh.data->data[i].im ;
   }
   output->var.data->data[k].re = 1.0 / (REAL4)sumvar.re;
   output->var.data->data[k].im = 1.0 / (REAL4)sumvar.im;

   k++;
 }

 output->B.data->length = n/params->M;
 output->var.data->length = n/params->M;

 output->B.epoch = input->Vh.epoch;
 output->var.epoch = input->Vh.epoch;

 output->B.deltaT = input->Vh.deltaT * params->M;
 output->var.deltaT = input->Vh.deltaT * params->M;

 output->phase = phase->data[0];

  /****** CLEAN UP *******/

  LALDestroyVector(status->statusPtr, &phase);
  LALZDestroyVector(status->statusPtr, &Xh);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*****************************************************************************/
