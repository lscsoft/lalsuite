/*
*  Copyright (C) 2018 Riccardo Sturani
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

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>

#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralWaveformParams.h>
#include <lal/LALSimSphHarmMode.h>

int main(void){

  /* This codes checks that the SpinTaylor approximant constructed on
   *
   */

  REAL8TimeSeries *V=NULL;
  REAL8TimeSeries *Phi=NULL;
  REAL8TimeSeries *S1x=NULL;
  REAL8TimeSeries *S1y=NULL;
  REAL8TimeSeries *S1z=NULL;
  REAL8TimeSeries *S2x=NULL;
  REAL8TimeSeries *S2y=NULL;
  REAL8TimeSeries *S2z=NULL;
  REAL8TimeSeries *LNhx=NULL;
  REAL8TimeSeries *LNhy=NULL;
  REAL8TimeSeries *LNhz=NULL;
  REAL8TimeSeries *E1x=NULL;
  REAL8TimeSeries *E1y=NULL;
  REAL8TimeSeries *E1z=NULL;

  REAL8 f_low = 50.;
  REAL8 f_ref=50.;
  REAL8 sample_rate = 8192.;
  REAL8 dT=1./sample_rate;
  REAL8 m1=21.;
  REAL8 m2=11.;
  REAL8 m1_SI=m1*LAL_MSUN_SI;
  REAL8 m2_SI=m2*LAL_MSUN_SI;
  REAL8 dist_SI=1.e6*LAL_PC_SI;
  REAL8 phi_ref=0.4;
  REAL8 s1x=0.3;
  REAL8 s1y=0.2;
  REAL8 s1z=0.5;
  REAL8 s2x=0.;
  REAL8 s2y=-0.4;
  REAL8 s2z=0.;
  REAL8 incl=0.1;
  REAL8 lam1=1.e3;
  REAL8 lam2=5.e2;
  int ampO=2;//Test will fail if amplitude order >2, as it will involve l=5 modes which are not yet coded
  const UINT4 LMAX=4;

  LALDict *params=XLALCreateDict();
  int errCode=XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(params,ampO);
  errCode   +=XLALSimInspiralWaveformParamsInsertTidalLambda1(params,lam1);
  errCode   +=XLALSimInspiralWaveformParamsInsertTidalLambda2(params,lam2);

  REAL8TimeSeries *hp_std=NULL,*hc_std=NULL;
  Approximant apprx=XLALGetApproximantFromString("SpinTaylorT4");
  errCode +=XLALSimInspiralChooseTDWaveform(&hp_std,&hc_std,m1_SI,m2_SI,s1x,s1y,s1z,s2x,s2y,s2z,dist_SI,incl,phi_ref,0.,0.,0.,dT,f_low,f_ref,params,apprx);

  REAL8 inclination,spin1x,spin1y,spin1z,spin2x,spin2y,spin2z;
  errCode +=XLALSimInspiralInitialConditionsPrecessingApproxs(&inclination,&spin1x,&spin1y,&spin1z,&spin2x,&spin2y,&spin2z,incl,s1x,s1y,s1z,s2x,s2y,s2z,m1,m2,f_ref,phi_ref,XLALSimInspiralWaveformParamsLookupFrameAxis(params));

  REAL8 lnhx=sin(inclination);
  REAL8 lnhy=0.;
  REAL8 lnhz=cos(inclination);
  REAL8 e1x=0.;
  REAL8 e1y=1.;
  REAL8 e1z=0.;
  UINT4 idx=0;

  REAL8TimeSeries *hp_drv=NULL,*hc_drv=NULL;
  errCode+=XLALSimInspiralSpinTaylorDriver(&hp_drv,&hc_drv,&V,&Phi,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,&LNhx,&LNhy,&LNhz,&E1x,&E1y,&E1z, phi_ref, dT, m1_SI, m2_SI, f_low, f_ref, dist_SI, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, lnhx, lnhy, lnhz, e1x, e1y, e1z, params, apprx);
  // Additional rotation by polarization angle Pi/2 as done in ChooseWf
  for (idx=0;idx<hp_drv->data->length;idx++) {
    hp_drv->data->data[idx]*=-1.;
    hc_drv->data->data[idx]*=-1.;
  }

  UINT4 l;
  INT4 m;
  LALValue *modearray=XLALSimInspiralCreateModeArray();
  for (l=2; l<=LMAX; l++)
    XLALSimInspiralModeArrayActivateAllModesAtL(modearray, l);

  SphHarmTimeSeries *Hlms_orb=NULL;
  errCode+=XLALSimInspiralSpinTaylorHlmModesFromOrbit(&Hlms_orb,V,Phi,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,m1_SI,m2_SI,dist_SI,ampO,modearray);
  XLALDestroyREAL8TimeSeries(V);
  XLALDestroyREAL8TimeSeries(Phi);
  XLALDestroyREAL8TimeSeries(S1x);
  XLALDestroyREAL8TimeSeries(S1y);
  XLALDestroyREAL8TimeSeries(S1z);
  XLALDestroyREAL8TimeSeries(S2x);
  XLALDestroyREAL8TimeSeries(S2y);
  XLALDestroyREAL8TimeSeries(S2z);
  XLALDestroyREAL8TimeSeries(LNhx);
  XLALDestroyREAL8TimeSeries(LNhy);
  XLALDestroyREAL8TimeSeries(LNhz);
  XLALDestroyREAL8TimeSeries(E1x);
  XLALDestroyREAL8TimeSeries(E1y);
  XLALDestroyREAL8TimeSeries(E1z);

  SphHarmTimeSeries *Hlms_modes=XLALSimInspiralChooseTDModes(0.,dT, m1_SI, m2_SI, s1x, s1y, s1z, s2x, s2y, s2z, f_low, f_ref, dist_SI, params, LMAX, apprx);

  UINT4 lmax=XLALSphHarmTimeSeriesGetMaxL(Hlms_orb);
  COMPLEX16TimeSeries *hlm_tmp=XLALSphHarmTimeSeriesGetMode(Hlms_orb,2,2);
  REAL8TimeSeries *hp_orb=XLALCreateREAL8TimeSeries("H+ from orbit", &hlm_tmp->epoch, 0., dT, &lalStrainUnit, hlm_tmp->data->length);
  REAL8TimeSeries *hc_orb=XLALCreateREAL8TimeSeries("Hx from orbit", &hlm_tmp->epoch, 0., dT, &lalStrainUnit, hlm_tmp->data->length);
  hlm_tmp=XLALSphHarmTimeSeriesGetMode(Hlms_modes,2,2);
  REAL8TimeSeries *hp_modes=XLALCreateREAL8TimeSeries("H+ from mods", &hlm_tmp->epoch, 0., dT, &lalStrainUnit, hlm_tmp->data->length);
  REAL8TimeSeries *hc_modes=XLALCreateREAL8TimeSeries("Hx from modes ",&hlm_tmp->epoch, 0., dT, &lalStrainUnit, hlm_tmp->data->length);
  memset(hp_orb->data->data, 0, sizeof(REAL8)*hc_orb->data->length);
  memset(hc_orb->data->data, 0, sizeof(REAL8)*hc_orb->data->length);
  memset(hp_modes->data->data, 0, sizeof(REAL8)*hp_modes->data->length);
  memset(hc_modes->data->data, 0, sizeof(REAL8)*hc_modes->data->length);

  for (l=2; l<=lmax; l++) {
    for (m=-l; m<=(INT4)l; m++) {
      errCode += XLALSimAddMode(hp_orb,  hc_orb,  XLALSphHarmTimeSeriesGetMode(Hlms_orb,l,m), 0., 0., l, m, 0);
      errCode += XLALSimAddMode(hp_modes, hc_modes, XLALSphHarmTimeSeriesGetMode(Hlms_modes,l,m), incl, LAL_PI/2.-phi_ref, l, m, 0);
    }
  }

  REAL8 tmp;
  INT4 ret=0;
  // We take out of the comparison the last samples
  UINT4 minlen=hp_orb->data->length;
  if (minlen>hp_modes->data->length)
     minlen=hp_modes->data->length;
  for (idx=1; idx<minlen; idx++) {
    tmp=sqrt((hp_std->data->data[idx]*hp_std->data->data[idx])+(hc_std->data->data[idx]*hc_std->data->data[idx]));
    if (fabs(hp_std->data->data[idx]-hp_drv->data->data[idx])>0.01*tmp)
      ret+=1;
    if (fabs(hp_std->data->data[idx]-hp_orb->data->data[idx])>0.01*tmp)
      ret+=1;
    if ( fabs(hp_std->data->data[idx]-hp_modes->data->data[idx])>0.01*tmp)
      ret+=1;
    if (fabs(hc_std->data->data[idx]-hc_drv->data->data[idx])>0.01*tmp)
      ret+=1;
    if (fabs(hc_std->data->data[idx]-hc_orb->data->data[idx])>0.01*tmp)
      ret+=1;
    if (fabs(hc_std->data->data[idx]-hc_modes->data->data[idx])>0.01*tmp)
      ret+=1;
  }

  if ( (ret == 0) && (errCode == 0) ) {
    fprintf(stdout, "\n Mode vs. polarization test passed.\n");
  }
  else {
    fprintf(stderr, "\nFAILURE: %u %u Mode vs. polarization test failed.\n", ret,errCode);
  }

 return ret;

}
