/*
*  Copyright (C) 2007 Jolien Creighton, Julien Sylvestre
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
 * File Name: TFClustersTest1.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TFClustersTest1
 *
 * SYNOPSIS
 * Implement a simple use of the algorithm on a white noise time series.
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


/* <lalLaTeX>
Sample usage of the routines in \texttt{TFClusters.h}.

\subsection*{Usage}
\texttt{TFClustersTest1}

\subsection*{Description}
Generates 128 seconds of white Gaussian noise at 2048 Hz, and apply the three levels of thresholding. Writes a list of clusters on stdout.

\subsection*{Exit codes}
Returns 0 on success, otherwise returns 1.

\subsection*{Uses}

\subsection*{Notes}

</lalLaTeX> */

#include "lal/LALRCSID.h"

NRCSID (MAIN, "$Id$");


#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/TFClusters.h>
#include <lal/TFCWavelet.h>
#include <lal/Random.h>
#include <string.h>

#define CHKST if(status.statusCode != 0) return -1

int lalDebugLevel = LALMSGLVL3;



int main(void) {

  static LALStatus status;
  REAL4TimeSeries tseries;
  CListDir dir;
  CList clist/*, list*/;
  TFCWParams twav;
  TFPlaneParams tspec;
  TFCSpectrogram spower;

  RandomParams *params = NULL;
  REAL4Vector *vect = NULL;

  REAL8 /*T,*/ P;
  UINT4 i, j, N;
  INT4 seed = 0;


  /* first generate a time series of white Gaussian noise of unit variance */
  N = 1024 * 128;

  LALCreateVector(&status, &vect, N);
  CHKST;

  LALCreateRandomParams(&status, &params, seed);
  CHKST;

  LALNormalDeviates(&status, vect, params);
  CHKST;

  tseries.epoch.gpsSeconds = 0;
  tseries.epoch.gpsNanoSeconds = 0;
  tseries.deltaT = 1.0/1024.0;
  tseries.f0 = 0.0;
  tseries.data = vect;

  /*
  for(i=1024*63;i<1024*63+16;i++) {
    vect->data[i] += (float)(i - 1024*63 + 1) / (float)1.0;
  }
  for(i=1024*63+16;i<1024*63+32;i++) {
    vect->data[i] += (float)(i - 1024*63) / (float)1.0;
  }
  */

  /*
  vect->data[1024*64+27] = 5.0;
  vect->data[1024*64+28] = 5.0;
  vect->data[1024*64+29] = 5.0;
  vect->data[1024*64+30] = 5.0;
  vect->data[1024*64+31] = 5.0;
  */

  twav.minScale = 4;
  twav.maxScale = 1024;
  tspec.timeBins = twav.timeBins = N / twav.minScale;
  tspec.freqBins = twav.freqBins = 9;
  tspec.deltaT = twav.deltaT = twav.minScale * tseries.deltaT;

  twav.wavelet = NULL;
  LALCreateVector( &status, &twav.wavelet, 4 );
  CHKST;
  memcpy( twav.wavelet->data, daub2,
      twav.wavelet->length * sizeof( *twav.wavelet->data ) );

  tspec.flow = 1.0 / twav.maxScale;

  spower.power = NULL;
  spower.params = NULL;

  LALComputeWaveletTFCSpectrogram(&status, &spower, &twav, &tseries);
  CHKST;


  /* Set thresholds */

  dir.freqBins = tspec.freqBins; /* number of frequency bins in spectrogram */
  dir.sigma = 5; /* threshold on cluster size */
  dir.minf = 0;
  dir.maxf = 1E30; /* max frequency to consider (Hz) */

  LALFillCListDir(&status, &dir, -log(1E-4)); /* allocate memory and set the threshold on power so that 1 every 10 pixel is black */
  CHKST;


  /* set thresholds on distance for different size couples */
  dir.d[0] = 0; /* 1,1 */
  dir.d[1] = 0; /* ... */
  dir.d[2] = 0;
  dir.d[3] = 0; /* 1,4 */
  dir.d[4] = 0; /* 2,2 */
  dir.d[5] = 0;
  dir.d[6] = 2; /* 2,4 */
  dir.d[7] = 3; /* 3,3 */
  dir.d[8] = 4; /* 3,4 */
  dir.d[9] = 4; /* 4,4 */

  dir.mdist = 4; /* no need to worry about things that are more than 4 units away from each other */

  /* run cluster threshold algorithm */
  LALInitCList(&status, &clist, &tspec); /* initialize clist */
  CHKST;

  LALGetClusters(&status, &clist, &spower, &dir); /* generate list of clusters */
  CHKST;

  LALFreeSpecgram(&status, &spower); /* spectrogram no longer useful */
  CHKST;



  LALFreeCListDir(&status, &dir);
  CHKST;




  /* display results to stdout */
  printf("Id\t\tSize\t\tPower\n");
  for(i=0; i<clist.nclusters; i++) {
    for(P=0, j=0; j<clist.sizes[i]; j++) P += clist.P[i][j];
    printf("%i\t\t%i\t\t%g\n",i,clist.sizes[i],P);
  }

  /* clean up */
  LALFreeCList(&status, &clist);
  CHKST;


  LALDestroyVector(&status, &vect);
  CHKST;

  LALDestroyVector( &status, &twav.wavelet );

  return 0;

}




