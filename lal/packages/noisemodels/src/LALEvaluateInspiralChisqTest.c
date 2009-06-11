/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Anand Sengupta
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

#include <lal/LALStdlib.h>
#include <lal/LALNoiseModels.h>

NRCSID( LALEVALUATEINSPIRALCHISQTESTC, "$Id$" );

/* Local Functions used only in this file */
static REAL4 **matrix(long nrl, long nrh, long ncl, long nch);
static void  free_matrix(REAL4 **m, long nrl, long nrh, long ncl, long nch);


void LALEvaluateInspiralChisqTest ( LALStatus                *status,
        InspiralChisqOutput   *chisqOut,
        InspiralChisqInput    *chisqIn
        )

{
    INT4                      i, j, imax, imin, k;
    REAL8                     f, df, sum, powNorm, DeltaT;
    REAL8                     *freqBndry=NULL;
    REAL4Vector               output1, output2;
    REAL4                     **pBandCorrelation1 = NULL, **pBandCorrelation2 = NULL;
    REAL4                     *pBandRho1=NULL, *pBandRho2=NULL;
    InspiralWaveCorrelateIn   corrin;

    INITSTATUS (status, "LALEvaluateInspiralChisqTest", LALEVALUATEINSPIRALCHISQTESTC);
    ATTATCHSTATUSPTR(status);

    /* If fCutoff > flso, we need to modify imax
     */
    DeltaT = (REAL8) chisqIn->findEventsIn->signal.length /  chisqIn->findEventsIn->param.tSampling;
    df     = 1./DeltaT;

    imin   = (INT4) (chisqIn->findEventsIn->param.fLower/ df);
    /*imin = 1;
     */
    if (chisqIn->findEventsIn->param.fCutoff > chisqIn->flso)
          imax = (INT4) ceil(chisqIn->flso / df);
    else
          imax = (INT4) ceil(chisqIn->findEventsIn->param.fCutoff / df);

    /* What is the total powNorm till flso
     */
    f       = 0.0;
    powNorm = 0.0;
    for (i=imin; i<=imax; i++) {
        f = df *i;
        powNorm += pow (f, -7./3.) / chisqIn->findEventsIn->psd.data[i];
    }

    /* freqBndry stores the obvious - note that we need one xtra element to store the last
     */
    freqBndry = (REAL8 *) LALCalloc(1, sizeof(REAL8) * (chisqIn->chisqBins+1));

    freqBndry[0] = df * imin;
    k = 1;
    f = 0.0;
    sum = 0.0;
    for (i=imin; i<=imax; i++){
        f = (REAL8)(i) * df;
        sum += pow (f, -7./3.) / chisqIn->findEventsIn->psd.data[i];
        if (sum >= (powNorm/chisqIn->chisqBins)) {
            freqBndry[k] = f;
            k += 1;
            sum = 0.0;
        }
    }
    freqBndry[chisqIn->chisqBins] = df * imax;

    /* output1 and 2 will store the output of overlap with filter1 and filter2 respectively.
     */
    output1.length = output2.length = chisqIn->findEventsIn->signal.length;

    if (!(output1.data = (REAL4*) LALMalloc(sizeof(REAL4)*output1.length))) {
        ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
    }

    if (!(output2.data = (REAL4*) LALMalloc(sizeof(REAL4)*output2.length))) {
        LALFree(output1.data);
        output1.data = NULL;
        ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
    }

    /* Init the corrin structure except fCutoff and signal2
     */
    corrin.df            = df;
    corrin.psd           = chisqIn->findEventsIn->psd;
    corrin.revp          = chisqIn->findEventsIn->revp;
    corrin.signal1       = chisqIn->findEventsIn->signal;
    corrin.samplingRate  = chisqIn->findEventsIn->param.tSampling;

    /* The correlation time series for p bands will be stored in a 2 dimensional array
     */
    /* called pBandCorrelation[0][...] to pBandCorrelation[p-1][...]
       Note the fact that pBandCorrelation[0][...] stores the correlation from fmin to freqBndry[0],
       pBandCorrelation[1][...] stores it from fmin to freqBndry[1] and so on. Thus in order to extract
       what is the actual correlation in the band freqBndry[0] to freqBndry[1] you have to subtract
       pBandCorrelation[1][...] - pBandCorrelation[0][...] term by term.
       First index is chisqbin, Second index is time bin
     */

    pBandCorrelation1 = matrix (1, chisqIn->chisqBins, 0,   (long)(chisqIn->findEventsIn->signal.length-1));
    pBandCorrelation2 = matrix (1, chisqIn->chisqBins, 0,   (long)(chisqIn->findEventsIn->signal.length-1));

    /* Loop over the freqBndry and do the IFFTs
     */
    for (i=1; i<=chisqIn->chisqBins; i++) {
        corrin.fCutoff       = freqBndry[i];

        corrin.signal2 = chisqIn->filter1;
        LALInspiralWaveCorrelate(status->statusPtr, &output1, corrin);
        CHECKSTATUSPTR(status);

        corrin.signal2 = chisqIn->filter2;
        LALInspiralWaveCorrelate(status->statusPtr, &output2, corrin);
        CHECKSTATUSPTR(status);

        /* Store it in pBandCorrelation
         */
        for (j=0; j<chisqIn->findEventsIn->signal.length; j++){
            pBandCorrelation1[i][j] = output1.data[j];
            pBandCorrelation2[i][j] = output2.data[j];
        }
    }

    pBandRho1 = (REAL4 *) LALCalloc (1, sizeof(REAL4) * (chisqIn->chisqBins+1));
    pBandRho2 = (REAL4 *) LALCalloc (1, sizeof(REAL4) * (chisqIn->chisqBins+1));

    /* OK Now all the IFFTs are over. We now need to calculate the chisquares properly.
     */
    for (j=0; j<chisqIn->findEventsIn->signal.length; j++) {

        pBandRho1[1] = pBandCorrelation1[1][j];
        pBandRho2[1] = pBandCorrelation2[1][j];

        for (i=2; i<=chisqIn->chisqBins; i++) {
            pBandRho1[i] = pBandCorrelation1[i][j] - pBandCorrelation1[i-1][j];
            pBandRho2[i] = pBandCorrelation2[i][j] - pBandCorrelation2[i-1][j];
        }

        chisqOut->chisqZERO[j] = chisqOut->chisqPIbyTWO[j]= 0.0;
        for (i=1; i<=chisqIn->chisqBins; i++) {
            chisqOut->chisqZERO[j]    += pow ( ((REAL8)(pBandRho1[i]) - chisqIn->rho1.data[j]/(REAL8)(chisqIn->chisqBins)), 2.0 );
            chisqOut->chisqPIbyTWO[j] += pow ( ((REAL8)(pBandRho2[i]) - chisqIn->rho2.data[j]/(REAL8)(chisqIn->chisqBins)), 2.0 );
        }

        chisqOut->chisqZERO[j]    *= (REAL8)(chisqIn->chisqBins);
        chisqOut->chisqPIbyTWO[j] *= (REAL8)(chisqIn->chisqBins);
        chisqOut->chisq[j] = chisqOut->chisqZERO[j] + chisqOut->chisqPIbyTWO[j];

    }

    LALFree (freqBndry);
    LALFree (output1.data);
    LALFree (output2.data);
    free_matrix (pBandCorrelation1, 1, chisqIn->chisqBins, 0,   (long)(chisqIn->findEventsIn->signal.length-1));
    free_matrix (pBandCorrelation2, 1, chisqIn->chisqBins, 0,   (long)(chisqIn->findEventsIn->signal.length-1));
    LALFree (pBandRho1);
    LALFree (pBandRho2);

    /* Normal exit
     */
    DETATCHSTATUSPTR(status);
    RETURN(status);
}


/*-----------------------------------------------------------
  SUBROUTINES ADAPTED FROM NUMERICAL RECIPES
  Local functions to declare a 2 dimensional array and free it
  These subroutines have been taken from Numerical Recipes and
  modified to suit our purposes.
  ------------------------------------------------------------*/

#define NR_END 1
#define FREE_ARG char*

static REAL4 **matrix(long nrl, long nrh, long ncl, long nch)
    /* allocate a REAL4 matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    REAL4 **m;

    /* allocate pointers to rows */
    m=(REAL4 **) LALMalloc((size_t)((nrow+NR_END)*sizeof(REAL4*)));
    if (!m) {
        fprintf (stderr, "allocation failure 1 in matrix()\n");
        abort ();
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl]=(REAL4 *) LALMalloc((size_t)((nrow*ncol+NR_END)*sizeof(REAL4)));
    if (!m[nrl]) {
        fprintf (stderr, "allocation failure 2 in REAL4 matrix()");
        abort ();
    }

    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    /* return pointer to array of pointers to rows */
    return m;
}


static void free_matrix(REAL4 **m, long nrl, long nrh, long ncl, long nch)
    /* free a REAL4 matrix allocated by matrix() */
{
    LALFree((FREE_ARG) (m[nrl]+ncl-NR_END));
    LALFree((FREE_ARG) (m+nrl-NR_END));
}


