/*
*  Copyright (C) 2007 David Churches, Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer, Duncan Brown, Riccardo Sturani, Laszlo Vereb, Drew Keppel
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

/**
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Module to set the pointers to the required energy and flux functions.
 * Normally, a user is not required to call this function to generate a waveform.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralChooseModel()</tt>
 * <ul>
 * <li> \c f: Output containing the pointers to the appropriate
 * energy, flux, frequency, timing and phasing functions.</li>
 * <li> \c ak: Output containing the PN expnasion coefficients.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * ### Description ###
 *
 * This module gives the post-Newtonian expansions and/or P-approximants
 * to the energy, its derivative and gravitational-wave flux functions. More
 * specifically, the <tt>static REAL8</tt> functions below give Taylor expansions
 * of \f$dE/dv,\f$ and \f${\cal F}(v),\f$ P-approximants of \f$e(v),\f$ \f$dE/dv\f$
 * (derived from \f$e(v)\f$) and \f${\cal F}(v).\f$
 *
 * \c LALInspiralChooseModel
 * is used to set pointers to the required energy and flux functions
 * \f$E^{\prime}_T(v),\f$ \f$\mathcal{F}_T(v),\f$ \f$E^{\prime}_P(v)\f$ and \f$\mathcal{F}_P(v),\f$
 * in <tt>expnFunc,</tt> as also the GW phasing and frequency fucntions used in
 * the various approximants to generate the waveform.
 * More specifically pointers are set to the following functions in the structure
 * \c expnFunc:
 * <ul>
 * <li> <tt>EnergyFunction *dEnergy</tt>
 * </li><li> <tt>FluxFunction *flux</tt>
 * </li><li> <tt>InspiralTiming2 *timing2</tt>
 * </li><li> <tt>InspiralPhasing2 *phasing2</tt>
 * </li><li> <tt>InspiralPhasing3 *phasing3</tt>
 * </li><li> <tt>InspiralFrequency3 *frequency3</tt></li>
 * </ul>
 * \c LALInspiralChooseModel also outputs in \c ak the
 * last stable orbit (LSO) velocity \f$v_\textrm{LSO}\f$ (as <tt>ak->vn</tt>)
 * defined by the equation \f$E'(v_\textrm{LSO})=0,\f$
 * the values of the GW frequency \f$f_\textrm{LSO}=v_\textrm{LSO}^3/(\pi m)\f$
 * (as <tt>ak->fn</tt>) and time (as <tt>ak->tn</tt>) elapsed from <tt>params->fLower</tt>
 * to smaller of \c fCutOff and <tt>ak->fn</tt> by evaluating the integral
 * \f{equation}{
 * t_n = t_{0} - m \int^{v_n}_{v_0} \frac{E^{\prime}(v)}{\mathcal{F}(v)} \, dv\,,
 * \f}
 * where \f$t_{0}\f$ (usually equal to zero) is the user specified starting
 * time for the waveform when the wave frequency reaches <tt>params->fLower</tt>
 * and \f$v_{0}= (\pi m f)^{1/3}\f$ (with \f$f=<tt>params->fLower</tt>\f$) is the  velocity
 * at time \f$t_{0}.\f$  Note that \f$E'(v)\f$ and \f${\cal F}(v)\f$ are defined in
 * <tt>f->dEnergy</tt> and <tt>f->flux.</tt>
 *
 * ### Algorithm ###
 *
 * Numerical integration is used to compute <tt>ak->tn.</tt>
 *
 * ### Uses ###
 *
 * LALInspiralTofV
 *
 * ### Notes ###
 *
 * <ul>
 * <li> See Damour, Iyer and Sathyaprakash, PRD 57, 885, 1998 for further details.
 * Damour, Iyer and Sathyaprakash, PRD 63, 044023, 2001 is a resource paper that
 * summarizes how to generate waveforms in different approximations to the dynamics
 * of a compact binary under radiation reaction.</li>
 * <li> The Pade Approximant for the 1PN expansion is undefined as also
 * EOB at orders less than 2PN. BCV is independent of the PN order.
 * Spinning waveforms are only defined at the highest PN order.</li>
 * </ul>
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

static REAL8 dEt0(REAL8 v, expnCoeffs *ak)
{
   REAL8 dEnergy;
   dEnergy = ak->dETaN * v;
   return (dEnergy);
}


static REAL8 dEt2(REAL8 v, expnCoeffs *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x);
   return (dEnergy);
}


static REAL8 dEt4(REAL8 v, expnCoeffs *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x);
   return (dEnergy);
}



static REAL8 dEt6(REAL8 v, expnCoeffs *ak)
{
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x + ak->dETa3*x*x*x);
   return (dEnergy);
}






static REAL8 Ft0(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10;
   return (flux);
}


static REAL8 Ft2(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2);
   return (flux);
}


static REAL8 Ft3(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v);
   return (flux);
}


static REAL8 Ft4(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4);
   return (flux);
}


static REAL8 Ft5(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v);
   return (flux);
}


static REAL8 Ft6(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6);
   return (flux);
}


static REAL8 Ft7(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6 + ak->FTa7*v6*v);
   return (flux);
}


/*
static REAL8 ep0(REAL8 v, expnCoeffs *ak)
{
   REAL8 x, energy;
   ak = NULL;
   x = v*v;
   energy = -x;
   return (energy);
}
*/


static REAL8 ep2(REAL8 v, expnCoeffs *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1 * x);
   return (energy);
}


static REAL8 ep4(REAL8 v, expnCoeffs *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x));
   return (energy);
}


static REAL8 ep6(REAL8 v, expnCoeffs *ak)
{
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x /(1. + ak->ePa3*x)));
   return (energy);
}

/*
static REAL8 dEp0(REAL8 v, expnCoeffs *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, y;
   energy = ep0(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1;
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}
*/


static REAL8 dEp2(REAL8 v, expnCoeffs *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep2(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1. / ((1. + ak->ePa1*x)*(1. + ak->ePa1*x));
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


static REAL8 dEp4(REAL8 v, expnCoeffs *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = ep4(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;

   denom = 1. + (ak->ePa1 + ak->ePa2) * x;
   denom = denom * denom;
   denergy = (1. + 2.*ak->ePa2*x + ((ak->ePa1 + ak->ePa2) * ak->ePa2 * x*x))/denom;
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}



static REAL8 dEp6(REAL8 v, expnCoeffs *ak)
{
   REAL8 energy, denergy, Energy, dEnergy, denom, x, y;
   x = v*v;
   energy = ep6(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   
   denom = 1. + (ak->ePa1 + ak->ePa2 + ak->ePa3) * x + ak->ePa1*ak->ePa3*x*x;
   denom = denom * denom;

   denergy = (1. + 2.*(ak->ePa2+ak->ePa3)*x + (ak->ePa1*ak->ePa2
           + ak->ePa2*ak->ePa2 + 2.* ak->ePa2*ak->ePa3
           + ak->ePa3*ak->ePa3) * x*x)
           /denom;
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}



/*
static REAL8 Fp0(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10;
   return (flux);
}
*/

/*
static REAL8 Fp1(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v) * (1.-v/ak->vpoleP4));
   return (flux);
}
*/

/*
static REAL8 Fp2(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v / (1.+ak->fPa2*v)) * (1.-v/ak->vpoleP4));
   return (flux);
}
*/


static REAL8 Fp3(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v)))
	* (1.-v/ak->vpoleP4));
   return (flux);
}


static REAL8 Fp4(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v)))) * (1.-v/ak->vpoleP4));
   return (flux);
}


static REAL8 Fp5(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v / (1.+ak->fPa5*v))))) * (1.-v/ak->vpoleP4));
   return (flux);
}


static REAL8 Fp6(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v))))))
        * (1.-v/ak->vpoleP6));
   /* */
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}


static REAL8 Fp7(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v)))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}

/* Flux for the EOBNRv2 model */
static REAL8 Fp8PP(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10, l6, l8;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   l6 = ak->FTl6;
   l8 = ak->FTl8 - ak->FTa2*ak->FTl6;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v
        / (1.+ak->fPa8*v))))))))
        * (1.-v/ak->vpolePP));
   flux *= (1.+  log(v/ak->vlsoPP) * (l6*v6 + l8*v8) ) ;
   return (flux);
}

/*
static REAL8 Fp8(REAL8 v, expnCoeffs *ak)
{
   REAL8 flux,v2,v4,v6,v8,v10, l6, l8;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   l6 = ak->FTl6;
   l8 = ak->FTl8 - ak->FTa2*ak->FTl6;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v
	/ (1.+ak->fPa8*v))))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * (l6*v6 + l8*v8) ) ;
   return (flux);
}
*/


void
LALInspiralChooseModel(
   LALStatus        *status,
   expnFunc         *f,
   expnCoeffs       *ak,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralChooseModel");

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if (XLALInspiralChooseModel(f, ak, params))
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

int
XLALInspiralChooseModel(
   expnFunc         *f,
   expnCoeffs       *ak,
   InspiralTemplate *params
   )
{
   REAL8 vn, vlso;
   TofVIn in1;
   REAL8 tofv;
   void *in2;

   if (f == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (ak == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (params == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (params->order == LAL_PNORDER_HALF || (INT4)params->order < 0
      || (INT4)params->order > 8)
   {
      XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", __func__,
         ((INT4)params->order)/2, ((INT4)params->order)%2?".5":"");
      XLAL_ERROR(XLAL_EINVAL);
   }

   vlso = 0;

   switch (params->order)
   {
      case LAL_PNORDER_NEWTONIAN:
      switch (params->approximant)
      {
         case AmpCorPPN:
         case Eccentricity:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
            ak->vn = ak->vlso = vlso = ak->vlsoT0;
            f->dEnergy = dEt0;
            f->flux = Ft0;
            f->phasing2 = &XLALInspiralPhasing2_0PN;
            f->timing2 = &XLALInspiralTiming2_0PN;
            f->phasing3 = &XLALInspiralPhasing3_0PN;
            f->frequency3 = &XLALInspiralFrequency3_0PN;
            break;
         case PadeT1:
         case PadeF1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_HALF:
        XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
        break;
      case LAL_PNORDER_ONE:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:

            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft2;
            f->phasing2 = &XLALInspiralPhasing2_2PN;
            f->timing2 = &XLALInspiralTiming2_2PN;
            f->phasing3 = &XLALInspiralPhasing3_2PN;
            f->frequency3 = &XLALInspiralFrequency3_2PN;
            break;
         case PadeT1:
         case PadeF1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_ONE_POINT_FIVE:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft3;
            f->phasing3 = &XLALInspiralPhasing3_3PN;
            f->frequency3 = &XLALInspiralFrequency3_3PN;
            f->phasing2 = &XLALInspiralPhasing2_3PN;
            f->timing2 = &XLALInspiralTiming2_3PN;
            break;
         case PadeT1:
            ak->vn = ak->vlso = vlso = ak->vlsoP0;
            f->dEnergy = dEp2;
            f->flux = Fp3;
            break;
         case PadeF1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_TWO:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
/*
   The value vlsoT4 is too large and doesn't work sometimes;
   so we use vlsoT2.
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt4;
            f->flux = Ft4;
            f->phasing2 = &XLALInspiralPhasing2_4PN;
            f->timing2 = &XLALInspiralTiming2_4PN;
            f->phasing3 = &XLALInspiralPhasing3_4PN;
            f->frequency3 = &XLALInspiralFrequency3_4PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp4;
            break;
         case PadeF1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_TWO_POINT_FIVE:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
/*
   The value vlsoT4 is too large and doesn't work with 2.5 PN
   Taylor approximant; so we use vlsoT2.
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt4;
            f->flux = Ft5;
            f->phasing2 = &XLALInspiralPhasing2_5PN;
            f->timing2 = &XLALInspiralTiming2_5PN;
            f->phasing3 = &XLALInspiralPhasing3_5PN;
            f->frequency3 = &XLALInspiralFrequency3_5PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp5;
            break;
         case PadeF1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_THREE:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
	     case SpinQuadTaylor:
/*
   vlsoT6 is as yet undetermined and vlsoT4 is too large in
   certain cases (TaylorT2 crashes for (1.4,10)); using vlsoT2;
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt6;
            f->flux = Ft6;
            f->phasing2 = &XLALInspiralPhasing2_6PN;
            f->timing2 = &XLALInspiralTiming2_6PN;
            f->phasing3 = &XLALInspiralPhasing3_6PN;
            f->frequency3 = &XLALInspiralFrequency3_6PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp6;
            break;
         case PadeF1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
      }
      break;
      case LAL_PNORDER_THREE_POINT_FIVE:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
		 case SpinQuadTaylor:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt6;
            f->flux = Ft7;
            f->phasing2 = &XLALInspiralPhasing2_7PN;
            f->timing2 = &XLALInspiralTiming2_7PN;
            f->phasing3 = &XLALInspiralPhasing3_7PN;
            f->frequency3 = &XLALInspiralFrequency3_7PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case EOBNRv2:
         case EOBNRv2HM:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp7;
            break;
         case PadeF1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
      }
      break;
      case LAL_PNORDER_PSEUDO_FOUR:
      switch (params->approximant)
      {
         case Eccentricity:
            XLALPrintError("XLAL Error - %s: The PN order requested is not implemented for this approximant\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         case EOBNRv2:
         case EOBNRv2HM:
            ak->vn = ak->vlso = vlso = ak->vlsoPP;
            f->dEnergy = dEp6;
            f->flux = Fp8PP;
            break;
         case EOB:
         case EOBNR:
         case IMRPhenomA:
         case IMRPhenomB:
         case IMRPhenomFA:
         case IMRPhenomFB:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp7;
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case TaylorF2RedSpin:
         case SpinTaylorFrameless:
         case SpinTaylorT3:
         case SpinTaylor:
         case PhenSpinTaylorRD:
         case SpinQuadTaylor:
         case PadeT1:
         case PadeF1:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
         default:
            XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
      }
      break;
      default:
         XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
         XLAL_ERROR(XLAL_EINVAL);
   }

   switch (params->approximant)
   {
      case AmpCorPPN:
      case TaylorT1:
      case TaylorT2:
      case TaylorT3:
      case TaylorF1:
      case EOB:
      case EOBNR:
      case PadeT1:
      case PadeF1:
      case TaylorF2:
      case TaylorF2RedSpin:
      case SpinTaylorFrameless:
      case SpinTaylorT3:
      case SpinTaylor:
      case PhenSpinTaylorRD:
      case SpinQuadTaylor:
      case TaylorEt:
      case TaylorT4:
      case TaylorN:
         ak->flso = vlso * vlso * vlso /(LAL_PI * ak->totalmass);

         if (ak->fn)
         {
            vn = cbrt(LAL_PI * ak->totalmass * ak->fn);
            ak->vn = (vn < vlso) ? vn :  vlso;
         }

         in1.t=0.0;
         in1.v0=ak->v0;
         in1.t0=ak->t0;
         in1.vlso=ak->vlso;
         in1.totalmass = ak->totalmass;
         in1.dEnergy = f->dEnergy;
         in1.flux = f->flux;
         in1.coeffs = ak;

         in2 = (void *) &in1;

         tofv = XLALInspiralTofV(ak->vn, in2);
         if (XLAL_IS_REAL8_FAIL_NAN(tofv))
            XLAL_ERROR(XLAL_EFUNC);

         ak->tn = -tofv - ak->samplinginterval;
         params->fCutoff = ak->fn = pow(ak->vn, 3.)/(LAL_PI * ak->totalmass);
         /*
         for (v=0; v<ak->vn; v+=0.001)
         {
            FtN = Ft0(v,ak);
            printf("%e %e %e %e %e %e %e\n", v,
            Ft2(v,ak)/FtN, Ft3(v,ak)/FtN, Ft4(v,ak)/FtN, Ft5(v,ak)/FtN,
            Ft6(v,ak)/FtN, Ft7(v,ak)/FtN);
         }
         exit(0);
         */
         break;
      case BCV:
      case BCVSpin:
         ak->tn = 100.;
         break;
      case IMRPhenomA:
      case IMRPhenomB:
      case IMRPhenomFA:
      case IMRPhenomFB:
      case EOBNRv2:
      case EOBNRv2HM:
         ak->tn = 5.*ak->totalmass/(256.*ak->eta*pow(ak->v0,8.)) + 1000.*ak->totalmass;
         break;
      case Eccentricity:
         /* The eccentric waveforms contain harmonic, so similarly to amplitude corrected waveforms
          * the duration are longer than non eccentric waveform and starts at 2fl/3*/
         ak->tn = 5.*ak->totalmass/256./ak->eta/pow(LAL_PI*ak->totalmass*params->fLower/3.*2.,8./3.);
         ak->flso = vlso * vlso * vlso /(LAL_PI * ak->totalmass);
         break;
      default:
         XLALPrintError("XLAL Error - %s: Unknown case in PN approximant switch\n", __func__);
         XLAL_ERROR(XLAL_EINVAL);
   }

  return 0;
}

