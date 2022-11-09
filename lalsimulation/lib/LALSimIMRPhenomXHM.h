/*
* Copyright (C) 2019 Marta Colleoni, Cecilio García Quirós
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef _LALSIM_IMR_PHENOMXHM_H       //ifndef is for conditional compilation - determines it doesnt exist before including in compilation
#define _LALSIM_IMR_PHENOMXHM_H       //creates definition  


#ifdef __cplusplus        //ifdef is for conditional compilation - determines it exists before including in compilation
extern "C" {              //Ensures names arent mangled within C++ compiler so names arent changed 
#endif                    //just closes ifdef statement 

#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>
   
//These are all other files contained within the lalsuite repository which are just being imported


// /* Returns the Fourier domain strain of just one negative mode: h_l-m. This quantity is zero for negative frequencies.
// However the positive mode h_lm is zero for positive frequencies and for the negative frequencies is equal to (-1)^l h*_l-m(-f).
// This is a wrapper function that use XLALSimIMRPhenomXASGenerateFD for the 22 mode and XLALSimIMRPhenomXHMOneMode for the higher modes. */


#include "LALSimIMRPhenomXHM_structs.h"      

int IMRPhenomXHMGenerateFDOneMode(                                                //int gives the variable an integer type
   COMPLEX16FrequencySeries **htildelm,  /**< [out] hlm for one mode **/          //** is pointer to pointer - so stores the address of another pointer variable
   const REAL8Sequence *freqs_In,        /**< fmin, fmax [Hz] **/                 //const specifies that th variable is constant so tells compiler to prevent programmer from modifying 
   IMRPhenomXWaveformStruct *pWF,        /**< structure of the 22 mode **/        //* is a pointer and will store addresses
   UINT4 ell,                            /**< first index of the mode **/         //comma evaluates first thing, discards then evaluates second and gives result 
   UINT4 emm,                            /**< second index of the mode **/
   LALDict *lalParams                    /**< extra params **/
);

/* Compute the frequency array and initialize htildelm to the corresponding length. */
int SetupWFArrays(
  REAL8Sequence **freqs,                /**< [out] frequency grid [Hz] */
  COMPLEX16FrequencySeries **htildelm,  /**< [out] Frequency domain hlm GW strain */
  const REAL8Sequence *freqs_In,        /**< fmin, fmax [Hz] */
  IMRPhenomXWaveformStruct *pWF,        /**< Waveform structure with parameters */
  LIGOTimeGPS ligotimegps_zero          /**< = {0,0} */
);

#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMXHM_H */
