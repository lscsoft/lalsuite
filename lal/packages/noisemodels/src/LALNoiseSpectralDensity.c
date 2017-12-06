/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

#include <lal/LALNoiseModels.h>

/**
 * \author Sathyaprakash, B. S.
 * \ingroup LALNoiseModels_h
 *
 * \brief This function generates an array of size specified in the vector \c psd, that is <tt>psd.length</tt>.
 * The inputs are
 * <ol>
 * <li> The lenght of the psd array required: this must be
 * given as a non-zero positive integer by setting the \c length
 * of the \c psd vector to the desired value;</li>
 * <li> Frequency resolution \c df in Hz.</li>
 * <li> Pointer to a function that should be used in generating the
 * power spectral density values in units of Hz\f$^{-1}.\f$ This function
 * must necessarily be of the type:
 * <tt> void  (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f).</tt></li>
 * </ol>
 * Presently, there are four such functions in the \c noisemodels
 * package. These are <tt>LALGEOPsd, LALLIGOIPsd, LALTAMAPsd, LALVIRGOPsd.</tt>
 * These four packages return a scaled PSD while this module returns the
 * correctly scaled version. It is assumed that new PSD modules return
 * unscaled versions. (Note, however, that it might be better to use the
 * scaled versions of the PSD when computing the metric on the signal
 * manifold; this is because computing the metric involves calculation of
 * many moments of the noise PSD and one might encounter round-off errors
 * if un-scaled version of PSD is used; I have not checked this to be
 * the case but suspect that there might be some problems.)
 *
 */
void
LALNoiseSpectralDensity
   (
   LALStatus    *status,
   REAL8Vector  *psd,
   void         (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f),
   REAL8        df
   )
{

    REAL8 f, shf, fs, s0;
    INT4 i, n;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (psd->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (NoisePsd, status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (df > 0., status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (psd->length > 1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   /*
   switch (NoisePsd)
   {
	case LALGEOPsd:
           s0 = 1.e-46;
           fs = 10.;
	   break;
        case LALLIGOIPsd:
           s0 = 9.0e-46;
           fs = 10.;
	   break;
        case LALVIRGOPsd:
           s0 = 3.24e-46;
           fs = 10.;
	   break;
        case LALTAMAPsd:
           s0 = 75.e-46;
           fs = 75.;
	   break;
	default:
           s0 = 1.;
           fs = 0.;
	   break;
   }
   */
   if (NoisePsd == LALGEOPsd) {
           s0 = 1.e-46;
           fs = 30.;
   } else if(NoisePsd == LALLIGOIPsd) {
           s0 = 9.0e-46;
           fs = 10.;
   } else if(NoisePsd == LALTAMAPsd) {
           s0 = 75.e-46;
           fs = 75.;
   } else if(NoisePsd == LALVIRGOPsd) {
           s0 = 1;
           fs = 10.;
   } else if(NoisePsd == LALAdvLIGOPsd) {
           s0 = 1e-49;
           fs = 10.;
   } else {
           s0 = 1.;
           fs = 0.;
   }

   n = (INT4) psd->length-1;

   /*
    * Set DC and Nyquist components to zero
    */

   psd->data[0] = 0.;

   for (i=0; i<=n; i++)
   {
      f = i*df;

      if (f>=fs)
      {
         (*NoisePsd)(status->statusPtr, &shf, f);
         psd->data[i] = s0 * shf;

      } else {

         psd->data[i] = 0.;
      }
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
