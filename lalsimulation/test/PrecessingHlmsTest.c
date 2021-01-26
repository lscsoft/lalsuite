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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

//#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralPrecess.h>
#include <lal/LALSimSphHarmMode.h>

/*
#include <lal/LALSimInspiral.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiralWaveformFlags.h>
#include <lal/Units.h>
#include <lal/XLALError.h>*/

#define EPSILON 1.e-8

#define ROTATEZ(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) - vy*sin(angle);\
	tmp2 = vx*sin(angle) + vy*cos(angle);\
	vx = tmp1;\
	vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) + vz*sin(angle);\
	tmp2 = - vx*sin(angle) + vz*cos(angle);\
	vx = tmp1;\
	vz = tmp2

static int c_compare(COMPLEX16 val1,
		     COMPLEX16 val2)
{
  if ( (fabs(creal(val1) - creal(val2)) > EPSILON) || (fabs(cimag(val1) - cimag(val2)) > EPSILON) )
    return 1;
  else
    return 0;
}

static int compare(REAL8 val1,
		   REAL8 val2)
{
  if ( ( (fabs(val1 - val2)) > EPSILON) > EPSILON )
    return 1;
  else
    return 0;
}

int main (int argc, char **argv)
{
    /* Ignore unused parameters. */
    (void)argc;
    (void)argv;

    const UINT4 Ntest=10;
    const REAL8 NtestR=(REAL8) Ntest;
    int errCode=0;
    int ret=0;
    REAL8 idxr;
    REAL8 alpha,iota,psi,iota2,psi2,v;
    REAL8 tmp1,tmp2,tmpx,tmpy,tmpz;

    const INT4 ampO  = -1;
    const REAL8 m1   = 10.*LAL_MSUN_SI;
    const REAL8 m2   =  5.*LAL_MSUN_SI;
    const REAL8 eta  = m1/(m1+m2)*m2/(m1+m2);
    const REAL8 dist = LAL_PC_SI*1.e6;
    const REAL8 norm = eta*(m1+m2)*LAL_G_SI/pow(LAL_C_SI,2.) / dist;
    const REAL8 chi1 = 0.8;
    const REAL8 chi2 = 0.3;
    const REAL8 dt   = 1.;
    const LIGOTimeGPS time0=LIGOTIMEGPSZERO;

    REAL8TimeSeries *zts     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *vts     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *psits   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *mpsits  = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *psi0ts  = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhx0   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhy0   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhz0   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhx    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhy    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *LNhz    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1x0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1y0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1z0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1x     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1y     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *e1z     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *iotats  = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *alphats = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *iota2ts = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *miota2ts= XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *psi2ts  = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *mpsi2ts = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *miotats = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *malphats= XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1x0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1y0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1z0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2x0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2y0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2z0    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1x     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1y     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S1z     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2x     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2y     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *S2z     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalDimensionlessUnit,Ntest);
    REAL8TimeSeries *hpF     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hcF     = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hpFs    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hcFs    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hpFc    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hcFc    = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hpFc2   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hcFc2   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hpFc3   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);
    REAL8TimeSeries *hcFc3   = XLALCreateREAL8TimeSeries("",&time0,0.,dt,&lalStrainUnit,Ntest);

    SphHarmTimeSeries *hStart, *hFin, *hFin2, *hFinRotm1, *hFinCheck, *hFinCheck2;
    hStart=NULL;
    hFin=NULL;
    hFin2=NULL;
    hFinRotm1=NULL;
    hFinCheck=NULL;
    hFinCheck2=NULL;
    COMPLEX16TimeSeries *hStartm,*hFinm,*hFinRotm1m,*hFin2m,*hFinChkm,*hFinChk2m;

    for (INT4 l=2;l<=4;l++) {

      for (UINT4 idx=0;idx<Ntest;idx++) {

	idxr=((REAL8)idx)/((REAL8) Ntest);
	alpha = idxr*2.*LAL_PI/NtestR;
	iota  = idxr*LAL_PI/NtestR;
	iota2 = idxr*LAL_PI/NtestR*0.4;
	psi   = (NtestR-idxr-1.)*2.*LAL_PI/NtestR;
	psi2  = (NtestR-idxr-1.)*2.*LAL_PI/NtestR*0.4;
	v     = (idxr+1.)/NtestR;

	zts->data->data[idx]=0.;
	vts->data->data[idx]=v;
	psits->data->data[idx]=psi;
	psi2ts->data->data[idx]=psi2;
	mpsits->data->data[idx]=-psi;
	mpsi2ts->data->data[idx]=-psi2;
	psi0ts->data->data[idx]=0.;
	iotats->data->data[idx]=iota;
	iota2ts->data->data[idx]=iota2;
	alphats->data->data[idx]=alpha;
	miotats->data->data[idx]=-iota;
	miota2ts->data->data[idx]=-iota2;
	malphats->data->data[idx]=-alpha;
	LNhx0->data->data[idx]=0.;
	LNhy0->data->data[idx]=0.;
	LNhz0->data->data[idx]=1.;
	e1x0->data->data[idx]=0.;
	e1y0->data->data[idx]=1.;
	e1z0->data->data[idx]=0.;
	LNhx->data->data[idx]=sin(iota)*cos(alpha);
	LNhy->data->data[idx]=sin(iota)*sin(alpha);
	LNhz->data->data[idx]=cos(iota);
	e1x->data->data[idx]=-sin(alpha);
	e1y->data->data[idx]=cos(alpha);
	e1z->data->data[idx]=0.;
	hpF->data->data[idx]=0.;
	hcF->data->data[idx]=0.;
	hpFs->data->data[idx]=0.;
	hcFs->data->data[idx]=0.;
	hpFc->data->data[idx]=0.;
	hcFc->data->data[idx]=0.;
	hpFc2->data->data[idx]=0.;
	hcFc2->data->data[idx]=0.;
	hpFc3->data->data[idx]=0.;
	hcFc3->data->data[idx]=0.;

	S1x0->data->data[idx]=0.;
	S1y0->data->data[idx]=0.;
	S1z0->data->data[idx]=chi1;
	S2x0->data->data[idx]=0.;
	S2y0->data->data[idx]=0.;
	S2z0->data->data[idx]=chi2;

	tmpx=S1x0->data->data[idx];
	tmpy=S1y0->data->data[idx];
	tmpz=S1z0->data->data[idx];
	ROTATEZ(psi,  tmpx,tmpy,tmpz);
	ROTATEY(iota, tmpx,tmpy,tmpz);
	ROTATEZ(alpha,tmpx,tmpy,tmpz);
	S1x->data->data[idx]=tmpx;
	S1y->data->data[idx]=tmpy;
	S1z->data->data[idx]=tmpz;

	tmpx=S2x0->data->data[idx];
	tmpy=S2y0->data->data[idx];
	tmpz=S2z0->data->data[idx];
	ROTATEZ(psi,  tmpx,tmpy,tmpz);
	ROTATEY(iota, tmpx,tmpy,tmpz);
	ROTATEZ(alpha,tmpx,tmpy,tmpz);
	S2x->data->data[idx]=tmpx;
	S2y->data->data[idx]=tmpy;
	S2z->data->data[idx]=tmpz;

      }

      //We create the mode with trivial and non-trivial argument
      if (l==2) {
	errCode += XLALSimInspiralSpinPNMode2m(&hStart,vts,psi0ts,LNhx0,LNhy0,LNhz0,e1x0,e1y0,e1z0,S1x0,S1y0,S1z0,S2x0,S2y0,S2z0,m1,m2,dist,ampO);
	errCode += XLALSimInspiralSpinPNMode2m(&hFin,vts,psits,LNhx,LNhy,LNhz,e1x,e1y,e1z,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,dist,ampO);
      }
      else if (l==3) {
        errCode += XLALSimInspiralSpinPNMode3m(&hStart,vts,psi0ts,LNhx0,LNhy0,LNhz0,e1x0,e1y0,e1z0,S1x0,S1y0,S1z0,S2x0,S2y0,S2z0,m1,m2,dist,ampO);
	errCode += XLALSimInspiralSpinPNMode3m(&hFin,vts,psits,LNhx,LNhy,LNhz,e1x,e1y,e1z,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,dist,ampO);
      }
      else if (l==4) {
	errCode += XLALSimInspiralSpinPNMode4m(&hStart,vts,psi0ts,LNhx0,LNhy0,LNhz0,e1x0,e1y0,e1z0,S1x0,S1y0,S1z0,S2x0,S2y0,S2z0,m1,m2,dist,ampO);
	errCode += XLALSimInspiralSpinPNMode4m(&hFin,vts,psits,LNhx,LNhy,LNhz,e1x,e1y,e1z,S1x,S1y,S1z,S2x,S2y,S2z,m1,m2,dist,ampO);
      }

      //and we check that a rotation takes one into the other and back
      errCode  = XLALSimInspiralPrecessionRotateModesOut(&hFin2, hStart, psits, iotats, alphats);
      errCode += XLALSimInspiralPrecessionRotateModesOut(&hFinRotm1, hFin, malphats, miotats, mpsits);
      //Two more checks with additional rotations
      errCode += XLALSimInspiralPrecessionRotateModesOut(&hFinCheck,   hFin, zts, iota2ts, psi2ts);
      errCode += XLALSimInspiralPrecessionRotateModesOut(&hFinCheck2,  hFin, psi2ts, iota2ts, zts);

      for(INT4 m=-l;m<=l;m++) {
	hStartm    = XLALSphHarmTimeSeriesGetMode(hStart, l, m );
	hFinRotm1m = XLALSphHarmTimeSeriesGetMode(hFinRotm1, l, m );
	hFinm      = XLALSphHarmTimeSeriesGetMode(hFin, l, m );
	hFin2m     = XLALSphHarmTimeSeriesGetMode(hFin2, l, m );
	hFinChkm   = XLALSphHarmTimeSeriesGetMode(hFinCheck, l, m );
	hFinChk2m  = XLALSphHarmTimeSeriesGetMode(hFinCheck2, l, m );

	errCode += XLALSimAddMode(hpF, hcF, hFinm, 0., 0., l, m, 0);
	errCode += XLALSimAddModeAngleTimeSeries(hpFs, hcFs, hStartm, miotats, mpsits, l, m, 0);
	errCode += XLALSimAddModeAngleTimeSeries(hpFc, hcFc, hFinChkm, iota2ts, psi2ts, l, m, 0);
	errCode += XLALSimAddMode(hpFc2, hcFc2, hFinChk2m, 0., 0., l, m, 0);
	errCode += XLALSimAddModeAngleTimeSeries(hpFc3, hcFc3, hFinm, miota2ts , mpsi2ts, l, m, 0);

	for(UINT4 idx=0;idx<Ntest;idx++) {
	  ret+=c_compare(hStartm->data->data[idx]/norm,hFinRotm1m->data->data[idx]/norm);
	  ret+=c_compare(hFinm->data->data[idx]/norm,hFin2m->data->data[idx]/norm);
	}

	XLALDestroyCOMPLEX16TimeSeries(hStartm);
	hStartm=NULL;
	XLALDestroyCOMPLEX16TimeSeries(hFinRotm1m);
	hFinRotm1m=NULL;
	XLALDestroyCOMPLEX16TimeSeries(hFinm);
	hFinm=NULL;
	XLALDestroyCOMPLEX16TimeSeries(hFinChkm);
	hFinChkm=NULL;
	XLALDestroyCOMPLEX16TimeSeries(hFinChk2m);
	hFinChk2m=NULL;
      }

      hStart=NULL;
      hFin=NULL;
      hFin2=NULL;
      hFinRotm1=NULL;
      hFinCheck=NULL;
      hFinCheck2=NULL;

      for(UINT4 idx=0;idx<Ntest;idx++) {
	ret+=compare(hpF->data->data[idx]/norm,(cos(2.*alphats->data->data[idx])*hpFs->data->data[idx]-sin(2.*alphats->data->data[idx])*hcFs->data->data[idx])/norm);
	ret+=compare(hcF->data->data[idx]/norm,(cos(2.*alphats->data->data[idx])*hcFs->data->data[idx]+sin(2.*alphats->data->data[idx])*hpFs->data->data[idx])/norm);
	ret+=compare(hpF->data->data[idx]/norm,hpFc->data->data[idx]/norm);
	ret+=compare(hcF->data->data[idx]/norm,hcFc->data->data[idx]/norm);
	ret+=compare(hpFc2->data->data[idx]/norm,hpFc3->data->data[idx]/norm);
	ret+=compare(hcFc2->data->data[idx]/norm,hcFc3->data->data[idx]/norm);
      }

    }

    if ( (ret == 0) && (errCode == 0) ) {
      fprintf(stdout, "\n Precessing modes test passed.\n");
    }
    else {
      fprintf(stderr, "\nFAILURE: %u Precessing modes test failed.\n", ret+errCode);
    }

    return ret + errCode ;
}
