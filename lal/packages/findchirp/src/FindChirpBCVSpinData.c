/*
*  Copyright (C) 2007 Duncan Brown, Gareth Jones, Patrick Brady
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCVSpinData.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVSpinDataCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinData.c}}
\label{ss:FindChirpBCVSpinData.c}

Provides functions to condition data prior to filtering with spinning BCV
detection templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpBCVDataCP}
\idx{LALFindChirpBCVData()}

The function \texttt{LALFindChirpBCVSpinData()} constions the data
as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpBCVSpinDataCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCVSpin.h>


NRCSID (FINDCHIRPBCVSPINDATAC, "$Id$");

/* <lalVerbatim file="FindChirpBCVSpinDataCP"> */
void
LALFindChirpBCVSpinData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, k;
  UINT4                 cut;

  REAL4                *w;
  REAL4                *amp;
  REAL4                *ampBCV;
  COMPLEX8             *wtilde;
  REAL4		       *tmpltPower;
  REAL4		       *tmpltPowerBCV;

  REAL4Vector          *dataVec;
  REAL4                *spec = NULL;
  COMPLEX8             *resp = NULL;

  COMPLEX8             *outputData;

  UINT4                *chisqBin    = NULL;
  UINT4                *chisqBinBCV = NULL;
  UINT4                 numChisqBins;

  REAL8                 *ampBCVSpin2;

  FindChirpSegment     *fcSeg   = NULL;
  DataSegment          *dataSeg = NULL;

  /* REMOVE THIS */
  /*FILE                 *fpDataIn     = NULL;
  char                  filename[10];
  char                  suffix[10];*/
  /* REMOVE THIS */


  /*declaration*/
  INITSTATUS( status, "LALFindChirpBCVSpinData", FINDCHIRPBCVSPINDATAC );
  ATTATCHSTATUSPTR( status );

  /* check that the output exists */
  ASSERT( fcSegVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec" );
  ASSERT( fcSegVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data" );
  ASSERT( fcSegVec->data->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->data" );
  ASSERT( fcSegVec->data->data->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->data->data" );

  /* check that the parameter structure exists */
  ASSERT( params, status,
	  FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
	  ": params" );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCVSpin, status,
      FINDCHIRPBCVSPINH_EMAPX, FINDCHIRPBCVSPINH_MSGEMAPX );

  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVecBCVSpin1, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVecBCVSpin1->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVecBCVSpin2, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVecBCVSpin2->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->wVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->wVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->invPlan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPBCVSPINH_EFLOW, FINDCHIRPBCVSPINH_MSGEFLOW );
  ASSERT( params->dynRange > 0, status,
      FINDCHIRPBCVSPINH_EDYNR, FINDCHIRPBCVSPINH_MSGEDYNR );

  /* check that the input exists */
  ASSERT( dataSegVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec" );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data" );
  ASSERT( dataSegVec->data->chan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data->chan" );
  ASSERT( dataSegVec->data->chan->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data->chan->data" );

  /*
   *
   * set up local segment independent pointers
   * will I need the last two? maybe for chisq
   *
   */


  w             = params->wVec->data;
  amp           = params->ampVec->data;
  ampBCV        = params->ampVecBCV->data;
  wtilde        = params->wtildeVec->data;
  tmpltPower    = params->tmpltPowerVec->data;
  tmpltPowerBCV = params->tmpltPowerVecBCV->data;
  ampBCVSpin2   = params->ampVecBCVSpin2->data;

  /*
   *
   * loop over data segments
   *
   */

  for ( i = 0; i < dataSegVec->length; ++i )
  {

	 /* REMOVE THIS */
	 /*snprintf (suffix, "%d", i);
         strcpy (filename, "dataSegment.");
         strncat (filename, suffix, 5);

         {
            fpDataIn     = fopen (filename,"w");
         }*/
         /* REMOVE THIS */

    /*
     	 *
     	 * set up segment dependent pointers
     	 *
     	 */

    	dataSeg      = &(dataSegVec->data[i]);
    	fcSeg        = &(fcSegVec->data[i]);

    	dataVec      = dataSeg->chan->data;
    	spec         = dataSeg->spec->data->data;
    	resp         = dataSeg->resp->data->data;

    	outputData   = fcSeg->data->data->data;


    	if ( fcSeg->chisqBinVec->length )
    	{
      		chisqBin     = fcSeg->chisqBinVec->data;
      		chisqBinBCV  = fcSeg->chisqBinVecBCV->data;
      		numChisqBins = fcSeg->chisqBinVec->length - 1;
    	}
    	else
    	{
      		numChisqBins = 0;
    	}

    	ASSERT( params->wtildeVec->length == fcSeg->data->data->length, status,
        FINDCHIRPBCVSPINH_EMISM, FINDCHIRPBCVSPINH_MSGEMISM );

	/* store the waveform approximant in the data segment */
    	fcSeg->approximant = BCVSpin;


       /* REMOVE THIS */
       /*{
                fprintf (stdout, "Writing input data (time domain) to file %s\n", filename );

                for ( k = 0; k < dataVec->length; ++k )
                {
                     fprintf (fpDataIn, "%d\t%e\n",  k, dataVec->data[k]);
                }
       }*/
       /* REMOVE THIS */

    	/*
     	 *
     	 * compute htilde and store in fcSeg
     	 *
     	 */

    	LALForwardRealFFT( status->statusPtr, fcSeg->data->data,
        	dataVec, params->fwdPlan );
    	CHECKSTATUSPTR( status );

  	/* compute strain */
  	for ( k = 0; k < fcSeg->data->data->length; ++k )
  	{
    		REAL4 p = outputData[k].re;
   	 	REAL4 q = outputData[k].im;
   	 	REAL4 x = resp[k].re * params->dynRange;
    		REAL4 y = resp[k].im * params->dynRange;

    		outputData[k].re =  (p*x) - (q*y);
    		outputData[k].im =  (p*y) + (q*x);
  	}


     /* set output frequency series parameters */
        strncpy( fcSeg->data->name, dataSeg->chan->name, LALNameLength );

        fcSeg->data->epoch.gpsSeconds     = dataSeg->chan->epoch.gpsSeconds;
        fcSeg->data->epoch.gpsNanoSeconds = dataSeg->chan->epoch.gpsNanoSeconds;

        fcSeg->data->f0     = dataSeg->chan->f0;
        fcSeg->data->deltaF = 1.0 /
 	       ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;

	fcSeg->deltaT       = dataSeg->chan->deltaT;
        fcSeg->number       = dataSeg->number;
        fcSeg->analyzeSegment = dataSeg->analyzeSegment;

        /* store low frequency cutoff and invSpecTrunc in segment */
        fcSeg->fLow         = params->fLow;
        fcSeg->invSpecTrunc = params->invSpecTrunc;

  } /* end of loop over data segments */


    	/*
     	 *
     	 * compute inverse power spectrum
      	 *
     	 */

    	/* set low frequency cutoff inverse power spectrum */
    	cut = params->fLow / dataSeg->spec->deltaF > 1 ?
      	params->fLow / dataSeg->spec->deltaF : 1;

    	/* set inverse power spectrum to zero */
    	memset( wtilde, 0, params->wtildeVec->length * sizeof(COMPLEX8) );

    	/* compute inverse of S_v */
   	for ( k = cut; k < params->wtildeVec->length; ++k )
    	{
      		if ( spec[k] == 0 )
      		{
        	ABORT( status, FINDCHIRPBCVSPINH_EDIVZ,
		FINDCHIRPBCVSPINH_MSGEDIVZ );
      		}

		wtilde[k].re = 1.0 / spec[k];
    	}

    	/*
     	 *
     	 * truncate inverse power spectrum in time domain if required
     	 *
     	 */

 	if ( params->invSpecTrunc )
 	{
      		/* compute square root of inverse power spectrum */
      		fprintf(stdout, "truncating wtilde! \n");

      		for ( k = cut; k < params->wtildeVec->length; ++k )
      		{
        	wtilde[k].re = sqrt( wtilde[k].re );
      		}

      		/* set nyquist and dc to zero */
      		wtilde[params->wtildeVec->length - 1].re = 0.0;
      		wtilde[0].re                             = 0.0;

      		/* transform to time domain */
      		LALReverseRealFFT( status->statusPtr,
			params->wVec,
			params->wtildeVec,
			params->invPlan );
      		CHECKSTATUSPTR (status);

      		/* truncate in time domain */
      		memset( w + params->invSpecTrunc/2, 0,
          	(params->wVec->length - params->invSpecTrunc) * sizeof(REAL4) );

      		/* transform to frequency domain */
      		LALForwardRealFFT( status->statusPtr, params->wtildeVec,
			params->wVec,
          		params->fwdPlan );
     	 	CHECKSTATUSPTR (status);

      		/* normalise fourier transform and square */
      		{
        		REAL4 norm = 1.0 / (REAL4) params->wVec->length;
        		for ( k = cut; k < params->wtildeVec->length; ++k )
        		{
          		wtilde[k].re *= norm;
          		wtilde[k].re *= wtilde[k].re;
          		wtilde[k].im = 0.0;
        		}
      		}

      		/* set nyquist and dc to zero */
      		wtilde[params->wtildeVec->length - 1].re = 0.0;
      		wtilde[0].re                             = 0.0;
 	}

    	/* set inverse power spectrum below cut to zero */
    	memset( wtilde, 0, cut * sizeof(COMPLEX8) );

    	/* convert from S_v to S_h */
    	for ( k = cut; k < params->wtildeVec->length; ++k )
    	{
      		REAL4 respRe = resp[k].re * params->dynRange;
      		REAL4 respIm = resp[k].im * params->dynRange;
      		REAL4 modsqResp = (respRe * respRe + respIm * respIm);
      		REAL4 invmodsqResp;

      		if ( modsqResp == 0 )
      		{
        		ABORT( status, FINDCHIRPBCVSPINH_EDIVZ,
				FINDCHIRPBCVSPINH_MSGEDIVZ );
      		}
      		invmodsqResp = 1.0 / modsqResp;

      		wtilde[k].re *= invmodsqResp;
    	}

	/* REMOVE THIS */
	/*fclose (fpDataIn);*/
        /* REMOVE THIS */


  DETATCHSTATUSPTR( status );
  RETURN( status );

}
