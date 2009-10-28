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
#include <lal/RealFFT.h>
#include <lal/LALNoiseModels.h>

NRCSID (LALTRUNCATEINVSPECTRUMC, "$Id$");

void LALTruncateInvSpectrum(
        LALStatus                *status,
        REAL8Vector              *inputVec,
        InvSpecTruncationParams    *params
        )
{
    INT4              n, trunc_n;
    UINT4             i;
    REAL8             df;
    COMPLEX8Vector    *Hvec = NULL;
    REAL4Vector       *hvec = NULL;
    REAL4             *w, norm;

    INITSTATUS( status, "LALTruncateInvSpectrum", LALTRUNCATEINVSPECTRUMC);
    ATTATCHSTATUSPTR( status );

    ASSERT (inputVec,         status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (params,           status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (params->fwdp,     status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (params->revp,     status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
    ASSERT (params->df > 0.0, status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

    n       = 2 * (params->n - 1);
    df      = params->df;
    trunc_n = (INT4) params->psdTruncTime * (n*df);

    /* Hvec will contain the FFT i.e 1./sqrt(shf)
     */
    LALCCreateVector ( status->statusPtr, &Hvec, n/2 + 1 );
    CHECKSTATUSPTR( status );

    /* hvec will contain time domain IFFT
     */
    LALSCreateVector ( status->statusPtr, &hvec, n );
    CHECKSTATUSPTR( status );

    /* w should point at hvec->data
     */
    w = hvec->data;

    /* Store sqrt (inv spectrum)  (f-domain) in Hvec
     */
    for (i=0; i<Hvec->length; i++){
        if (inputVec->data[i] > 0.0) {
            Hvec->data[i].re = (REAL4) (1./sqrt(inputVec->data[i]));
            Hvec->data[i].im = 0.0;
        } else{
            ABORT( status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE );
        }
    }

    /* ---- If we are debugging then print the inv sqrt spectrum
       as it was before (input) --- */
    if (params->ifDebug) {
        FILE  *out=NULL;
        /*REAL8 dt;*/
        /* We should print the frequency domain sqrt (inv spectrum) */
        out = fopen("FB_sqrt_invSpectrum_before.f","w");
        for (i=0; i<Hvec->length; i++){
            fprintf (out, "%e %e\n",
                    params->df*i,
                    Hvec->data[i].re);
        }
        fclose(out);
    }

    /* Set Nyquist and zero frequency components of the input spectrum to zero
     */
    Hvec->data[Hvec->length - 1].re = 0.0;
    Hvec->data[0].re                = 0.0;

    /* Inverse Fourier Transform to time domain
     */
    if (XLALREAL4ReverseFFT(hvec, Hvec, params->revp) != 0)
      ABORTXLAL(status);

    /* ---- If debugging then print the spectrum in time domain
       after its inverse Fourier transform --- */
    if (params->ifDebug) {
        FILE  *out=NULL;
        REAL8 dt;
        dt = 1.0 / (params->df * n);
        /* We should print the time domain sqrt(inv spectrum) */
        out = fopen("FB_sqrt_invSpectrum.t","w");
        for (i=0; i<hvec->length; i++){
            fprintf (out, "%e %e\n",
                    dt*i,
                    hvec->data[i]);
        }
        fclose(out);
    }

    /* truncate in time domain
     */
    memset( w + trunc_n/2, 0, (hvec->length - trunc_n) * sizeof(REAL4) );

    /* ---- If debugging print the time domain spectrum after truncation
       in the time domain --- */
    if (params->ifDebug) {
        FILE  *out=NULL;
        REAL8 dt;
        dt = 1. / (params->df * n);
        /* We should print the time domain sqrt (inv spectrum) after padding */
        out = fopen("FB_sqrt_invSpectrum_pad.t","w");
        for (i=0; i<hvec->length; i++){
            fprintf (out, "%e %e\n",
                    dt*i,
                    hvec->data[i]);
        }
        fclose(out);
    }

    /* transform to frequency domain
     */
    if (XLALREAL4ForwardFFT(Hvec, hvec, params->fwdp) != 0)
      ABORTXLAL(status);

    /* normalise fourier transform and square
     */
    norm = 1.0 / (REAL4) (n);
    for ( i = 0; i < Hvec->length; i++ )
    {
        Hvec->data[i].re *= norm;
        Hvec->data[i].re *= Hvec->data[i].re;
        Hvec->data[i].im  = 0.0;
    }

    /* ---- If debugging print the inv sqrt spectrum in Fourier domain
       after time domain truncation ---- */
    if (params->ifDebug) {
        FILE  *out=NULL;
        /* We should print the frequency domain sqrt (inv spectrum) after inv truncation */
        out = fopen("FB_sqrt_invSpectrum_after.f","w");
        for (i=0; i<Hvec->length; i++){
            fprintf (out, "%e %e\n",
                    params->df*i,
                    sqrt(Hvec->data[i].re));
        }
        fclose(out);
    }

    /* populate the input vector structure back with truncated psd
     *
     */
    for (i=0; i<Hvec->length; i++){
        if (Hvec->data[i].re > 0.0)
              inputVec->data[i] = (REAL8) (1./Hvec->data[i].re);
    }

    /* Clear work space
     */
    LALCDestroyVector (status->statusPtr, &Hvec);
    CHECKSTATUSPTR(status);

    LALSDestroyVector (status->statusPtr, &hvec);
    CHECKSTATUSPTR( status );

    /* normal exit
     */
    DETATCHSTATUSPTR( status );
    RETURN( status );
}
