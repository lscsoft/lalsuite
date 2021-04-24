/*
  *  Copyright (C) 2019 Andrew Matas
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRv4ROM_NSBHAmplitudeCorrection.h"

/**
 * Tanh window function 
 * w^{+/-}_{f0,sigma}(f) = 0.5 * (1 +/- tanh(4*(f-f0)/sigma))
 * Auxiliary function used by SEOBNRv4_ROM_NRTidalv2_NSBH.
 */
REAL8 TanhWindow(
       const REAL8 f, /**< frequency at which to evaluate window function */ 
       const int sign, /**< sign (+1 or -1), used to determine whether window is "on" or "off"*/
       const REAL8 f0, /**< central frequency of window function */
       const REAL8 sigma) /**< width of window function */
{
  REAL8 window;
  REAL8 x;
  x=4*(f-f0)/sigma;
  window =0.5 * (1+sign*tanh(x));
  return window;
}


/**
 * Returns ringdown frequency in units of the total mass using 
 * fits from gr-qc/0512160/. Auxiliary function used by SEOBNRv4_ROM_NRTidalv2_NSBH.
 */ 
REAL8 CalcRDFrequency(
      REAL8 MF, /**< Final mass (in solar masses) */ 
      REAL8 chiF, /**< Final spin (dimensionless */
      REAL8 Mtot) /**< Iniital toal mass (in solar masses */
{
    REAL8 k1=1.5251;
    REAL8 k2=-1.1568;
    REAL8 k3=0.1292;
    REAL8 fRD=1/(2*LAL_PI) * (k1+k2*pow(1-chiF,k3)) * Mtot/MF;
    return fRD;
}

/**
 * @addtogroup LALSimIMRTIDAL_c
 * 
 * @{
 *
 * @name SEOBNRv4_ROM_NRTidalv2_NSBH
 *
 * @author Andrew Matas
 *
 * @brief C code for SEOBNRv4_ROM_NRTidalv2_NSBH model. A technical note deescribing the model can be found at https://dcc.ligo.org/LIGO-T1900723.
 *
 * SEOBNRv4_ROM_NRTidalv2_NSBH is a frequency domain model that applies amplitude corrections due to tidal disruption to the SEOBNRv4ROM model. It is based on the SEOBNRv4_ROM_NRTidalv2.
 *
 * @note Parameter ranges:
 *  * 1 <= q = m1/m2 <= 100 
 *  * m2 <= 3 Msun
 *  * lambda1 = 0
 *  * 0 <= lambda2 <= 5000
 *  mi = mass of object i, lambdai = tidal parameter of object i (i={1,2}). i=1 is the black hole, i=2 is the neutron star. The model was compared with NR simulations with BH spin magnitudes up to 0.9.
 *
 *  @note A warning is issued when
 *  * chi2 is not 0 (model was fit to NR simulations with chi2=0, but checked against existing simulations with chi2=-0.2)
 *  * m1 < 1 Msun 
 *
 *  chi2 = neutron star spin
 *
  * @review SEOBNRv4_ROM_NRTidalv2_NSBH review by Frank Ohme, Tim Dietrich, Shrobana Ghosh,
  * Andrew Matas, Jonathan Thompson, Edward Fauchon-Jones. The review concluded
  * on 3 February 2020. The review documentation, resources, and final git hash
  * can be found at https://git.ligo.org/waveforms/reviews/nsbh-models/wikis/home.
  *
 * @{
 */

 /**
  * Compute amplitude correction to SEOBNRv4_ROM_NRTidalv2 in LAL format at specified 
  * frequencies for the SEOBNRv4_ROM_NRTidalv2_NSBH model, incorporating tidal disruption
  * effects. 
  *
  * XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries returns a frequency series with
  * real numbers between 0 and 1, which corrects the amplitude of SEOBNRv4_ROM_NRTidalv2.
  *
  */
int XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries(
    const REAL8Sequence *amp_tidal, /**< [out] tidal amplitude frequency series */
    const REAL8Sequence *fHz, /**< list of input Gravitational wave Frequency in Hz to evaluate */
    REAL8 m1_SI, /**< Mass of companion 1 (kg) */
    REAL8 m2_SI, /**< Mass of companion 2 (kg) */
    REAL8 chi1, /**< Spin of black hole */
    REAL8 lambda2 /**< (tidal deformability of mass 2) / m2^5 (dimensionless) */
) {

   // Change masses to natural units
    REAL8 MBH = m1_SI / LAL_MSUN_SI;
    REAL8 MNS = m2_SI / LAL_MSUN_SI;

    // Derive useful quantities used later

    // Properties of progenitor
    REAL8 Mtot=MBH+MNS;
    REAL8 MtotSI=m1_SI+m2_SI; // Total mass
    REAL8 fmtotSI=pow(LAL_C_SI,3) / ( LAL_G_SI * MtotSI ); // frequency for total mass
    REAL8 nu = (MBH*MNS)/pow(MBH+MNS,2); // symmetric mass rtio
    REAL8 Q=MBH/MNS;
    REAL8 CNS=XLALSimNSBH_compactness_from_lambda(lambda2);
    REAL8 RNS=MNS/CNS; // NS radius

    // Properties of remnant
    REAL8 MF=XLALBHNS_mass_aligned(MBH, MNS, chi1, lambda2); // final mass of remnant
    REAL8 chiF=XLALBHNS_spin_aligned(MBH, MNS, chi1, lambda2); // final spin of remnant
    REAL8 fRD=CalcRDFrequency(MF,chiF,Mtot); // ringdown frequency of remant in units of total mass


    // Mass shedding
    REAL8 xiTide=XLALSimNSBH_xi_tide(Q,chi1,MBH/RNS); // Relativistic correction to mass shedding
    REAL8 ZERO = 1e-15; // to prevent tidal frequency from going to infinity in BH limit
    REAL8 RTidal=xiTide*RNS*(1-2*CNS) + ZERO;
    REAL8 fTidal=XLALSimNSBH_fGWinKerr(RTidal,MF,chiF) * Mtot;
    REAL8 MbTorus=XLALSimNSBH_torus_mass_fit(Q, chi1, CNS);

    fTidal=fabs(fTidal);

    // Construct correction functions

    // Compute non-disrupitve correction 
    // parameters determined by fit to NR data
    REAL8 x1=-0.09236597801342522; 
    REAL8 x2=-0.1773927624795226; 
    REAL8 xND_C=-0.4865330927898738;
    REAL8 xND_chi=-0.03143937714260868; 
    REAL8 xNDprime_C=0.4933764101669873;
    REAL8 xNDprime_chi=0.05691547067814197;
    REAL8 d1=0.01871545791809104;
    REAL8 d2=0.771909557448921;
    REAL8 dND=0.022500562246265655;

    REAL8 fratio=(fTidal-fRD)/fRD;
    REAL8 xND=pow(fratio,2)+xND_C*CNS+xND_chi*chi1;
    REAL8 xNDprime=pow(fratio,2)+xNDprime_C*CNS+xNDprime_chi*chi1;
    REAL8 eTide=TanhWindow(xND,+1,x1,d1);
    REAL8 sigma_tide=2*TanhWindow(xNDprime,-1,x2,d2);

    REAL8 f0_ND=fRD;
    REAL8 sigma_ND=dND+sigma_tide;

    // Compute and apply disrupitve correction 
    // parameters determined by fit to NR data
    REAL8 a1=1.2728043573489636;
    REAL8 a2=0.1853261083544252;
    REAL8 b1=-1.6873457237092873;
    REAL8 b2=-0.25347578534406;
    REAL8 xD_C=0.8496732940251721;
    REAL8 xD_nu=0.3022694700157108;
    REAL8 xDprime_C=-0.9904717980366731;
    REAL8 xDprime_nu=1.1227719410457802;
    REAL8 xD_chi=-0.16594256718148745;
    REAL8 xDprime_chi1=0.002986871614045452;
    REAL8 xDprime_chi2=-0.07136411471590108;
    REAL8 xDprime_chi3=-0.11261503453409044;


    REAL8 xD=MbTorus + xD_C*CNS + xD_nu*pow(nu,0.5) + xD_chi*chi1;
    REAL8 xDprime=MbTorus + xDprime_C*CNS + xDprime_nu*pow(nu,0.5) 
                   + xDprime_chi1*chi1 + xDprime_chi2*pow(chi1,2) 
                   + xDprime_chi3*pow(chi1,3);

    REAL8 eins=a1+b1*xD;
    REAL8 sigma_tide2=a2+b2*xDprime;

    REAL8 f0_D=eins*fTidal;
    REAL8 sigma_D=sigma_tide2;

    
    REAL8 f0=0;
    REAL8 sigma=0;
    REAL8 eRD=0;
    if (fTidal>=fRD && MbTorus==0){ // Case 1 (non-disruptive)
        f0=f0_ND;
        sigma=sigma_ND;
        eRD=eTide;
    } 
    else if (fTidal<fRD && MbTorus>0){ // Case 2 (disruptive)
        f0=f0_D;
        sigma=sigma_D;
        eRD=0;
    }
    else if (fTidal<fRD && MbTorus==0){ // Case 3 (mildly disruptive without remnant torus)
        f0 = (1-pow(Q,-1))*f0_ND + pow(Q,-1)*f0_D;
        sigma=0.5*(sigma_ND+sigma_D);
        eRD=0;
    }
    else if (fTidal>=fRD && MbTorus>0){ // Case 4 (mildly disruptive with remnant torus)
        f0 = eins*fRD; 
        sigma=sigma_ND;
        eRD=eTide;
    }
   

   // Apply window corrections
   for (UINT8 ii=0;ii<(*fHz).length;ii++){
       (*amp_tidal).data[ii] = TanhWindow((*fHz).data[ii]/fmtotSI,-1,f0,sigma)
                                + eRD*TanhWindow((*fHz).data[ii]/fmtotSI,+1,f0,sigma);
   }
   
   return XLAL_SUCCESS;
}

/** @} */
/** @} */

