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

#include <stdio.h>

#include "lal/LALRCSID.h"

NRCSID (MAIN, "$Id$");


#include <lal/AVFactories.h>
#include <lal/TFClusters.h>
#include <lal/Random.h>

#define CHKST if(status.statusCode != 0) return -1

int lalDebugLevel = LALMSGLVL3 | LALNMEMDBG;


int main(void) {

  static LALStatus status;
  REAL4TimeSeries tseries;
  CListDir dir;
  CList clist;
  CList *list, * mainlist;
  TFPlaneParams tspec;
  TFCSpectrogram spower;

  RandomParams *params = NULL;
  REAL4Vector *vect = NULL;

  REAL8 T, P;
  UINT4 i, j, k, N, Nlists;
  INT4 seed = 0;


  /* first generate a time series of white Gaussian noise of unit variance */
  N = 1024 * 128;

  LALCreateVector(&status, &vect, N);
  CHKST;

  LALCreateRandomParams(&status, &params, seed);
  CHKST;

  tseries.epoch.gpsSeconds = 0;
  tseries.epoch.gpsNanoSeconds = 0;
  tseries.deltaT = 1.0/1024.0;
  tseries.f0 = 0.0;
  tseries.data = vect;

  /* Next compute a spectrogram for the time series */
  T = 1.0; /* this is the resolution in seconds of the spectrogram */

  LALPlainTFCSpectrogram(&status, &tspec, &tseries, T); /* this creates spectrogram parameters at the 'Heisenberg limit' from DC+1/T to the Nyquist frequency */
  CHKST;

  /* Set thresholds */

  dir.freqBins = tspec.freqBins; /* number of frequency bins in spectrogram */
  dir.sigma = 5; /* threshold on cluster size */
  dir.minf = 0;
  dir.maxf = 512; /* max frequency to consider (Hz) */

  LALFillCListDir(&status, &dir, -log(0.1)); /* allocate memory and set the threshold on power so that 1 every 10 pixel is black */
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


  mainlist = NULL;

  Nlists = 5;

  list = (CList *)LALMalloc(Nlists * sizeof(CList));

  for(k=0; k<Nlists; k++) { /* compute a number of lists */

  LALNormalDeviates(&status, vect, params);
  CHKST;

  tseries.data = vect;

  spower.power = NULL;
  spower.params = NULL;

  LALComputeTFCSpectrogram(&status, &spower, &tspec, &tseries);
  CHKST;


  /* run cluster threshold algorithm */
  LALInitCList(&status, &clist, &tspec); /* initialize clist */
  CHKST;

  LALGetClusters(&status, &clist, &spower, &dir); /* generate list of clusters */
  CHKST;

  LALFreeSpecgram(&status, &spower); /* spectrogram no longer useful */
  CHKST;



  /* run threshold on cluster total power */
  dir.alpha = 0.25; /* only 1/4 of all clusters from white noise will make it */

  LALInitCList(&status, list + k, &tspec); /* initialize list */
  CHKST;

  LALClustersPowerThreshold(&status, list + k, &clist, &dir); /* generate new list */
  CHKST;


  /* clean up a bit */
  LALFreeCList(&status, &clist);
  CHKST;

  /* merge */
  if(!mainlist) {mainlist = list + k;}
  else {
    LALMergeClusterLists(&status, mainlist, mainlist, list + k);
    CHKST;
  }
  }




  /* display results to stdout */
  printf("Id\t\tSize\t\tPower\n");
  for(i=0; i<mainlist->nclusters; i++) {
    for(P=0, j=0; j<mainlist->sizes[i]; j++) P += mainlist->P[i][j];
    printf("%i\t\t%i\t\t%g\n",i,mainlist->sizes[i],P);
  }

  /* clean up */
  LALFreeCListDir(&status, &dir);
  CHKST;

  for(k=0;k<Nlists;k++) {
    LALFreeCList(&status, list + k);
    CHKST;
  }
  LALFree(list);

  LALDestroyVector(&status, &vect);
  CHKST;

  return 0;

}




