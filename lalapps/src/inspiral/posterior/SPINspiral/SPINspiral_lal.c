/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_lal.c:          interfaces to LAL routines
   
   
   Copyright 2007, 2008, 2009 Christian Roever, Marc van der Sluys, Vivien Raymond, Ilya Mandel
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>


//#include <lal/LALConstants.h>


#include <lal/DetectorSite.h>

#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
//#include <lal/DetectorSite.h>

//////////////////////////////////////////
#include <SPINspiral.h>
//////////////////////////////////////////


#include <stdio.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>

#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/Units.h>
//#include <lal/TimeFreqFFT.h>
#include <lal/LALDetectors.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/VectorOps.h>
#include <lal/Random.h>
//#include <lal/LALNoiseModels.h>
#include <lal/XLALError.h>
#include <lal/LIGOLwXMLRead.h>
//#include <lal/LIGOLwXMLInspiralRead.h>


// For templateLALPhenSpinTaylorRD():
#include <lal/NRWaveInject.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Inject.h>




/**
 * \file SPINspiral_lal.c
 * \brief Contains interfaces to LAL routines
 */



// ****************************************************************************************************************************************************  
/**
 * \brief Compute a 12-parameter (one spin) LAL waveform
 *
 * \todo Merge this routine with LALHpHc12() ?
 *
 * Use the LAL <=3.5 PN spinning waveform, with 1 spinning object (12 parameters)
 */
// ****************************************************************************************************************************************************  
void templateLAL12(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
  int i=0;
  double samplerate=0.0,inversesamplerate=0.0;
  int length=0;
  
  samplerate = (double)ifo[ifonr]->samplerate;
  inversesamplerate = 1.0/samplerate;
  length     = ifo[ifonr]->samplesize;
  //printf("length = %d\n",length);
  
  // double hplusLAL[length+2];
  // double hcrossLAL[length+2];
  double *wave = (double*)calloc(length+2,sizeof(double));
  int lengthLAL = 0;
  
  
  /*LAL thingies needed. Have to be freed later*/
  
  static LALStatus    status;
  CoherentGW          waveform;
  SimInspiralTable    injParams;
  PPNParamStruc       ppnParams;
  
  memset( &status, 0, sizeof(LALStatus) );
  memset( &waveform, 0, sizeof(CoherentGW) );
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  ppnParams.deltaT   = inversesamplerate;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;
  
  
  
  // Compute h_+ and h_x
  //LALHpHc(&thewaveform, hplusLAL, hcrossLAL, &lengthLAL, length, par, ifo[ifonr], ifonr);
  LALHpHc12(&status, &waveform, &injParams, &ppnParams, &lengthLAL, par, ifo[ifonr], injectionWF, run);  //ifonr is for debugging purposes
  
  
  
  // Compute the detector response
  //double delay = LALFpFc(&thewaveform, wave, &lengthLAL, length, par, ifonr);
  /* double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); */ //Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x) - JOHN: Commented out as never used.
  /* delay = delay; */ //MvdS: remove 'declared but never referenced' warnings
  
  // printf("LALdelay = %10.10f\n", delay);
  
  LALfreedomSpin(&waveform);
  
  // printf("localtc = %f\n", localtc);
  
  //  localtc = ((par->par[2] - ifo[ifonr]->FTstart) - delay);
  
  //  int indexstart;
  //  indexstart = (int) (localtc*samplerate);// - (double)lengthLAL);   //MvdS: lengthLAL is used here, but seems never to be set in LALFpFc. Is it set in LALHpHc? Vivien: Yes it is but it is also available in thewaveform and wave
  //  if (indexstart<0) indexstart = 0;
  
  //printf("localtc2 = %f\n", localtc2);
  //printf("i1 = %d\n", i1);
  //printf("indexstart = %d\n" , indexstart);
  //printf("ifo[%d]->FTstart = %f\n", ifonr, ifo[ifonr]->FTstart);
  //printf("lengthLAL = %d\n", lengthLAL);
  //printf("length: %d,  lengthLAL: %d,  i1: %d,  i2: %d\n", length, lengthLAL,indexstart, indexstart+lengthLAL);
  
  for (i=0; i<length; ++i){
    //MvdS:  Reduce the number of if statements:
    //if(i<indexstart) ifo[ifonr]->FTin[i]   = 0.0;
    //if(i>=indexstart && i<(indexstart+lengthLAL)) ifo[ifonr]->FTin[i] = wave[i-indexstart];
    //if(i>=(indexstart+lengthLAL)) ifo[ifonr]->FTin[i]   = 0.0;
    //ifo[ifonr]->FTin[i] = 0.0;
    // if(i>=indexstart && i<(indexstart+lengthLAL)) ifo[ifonr]->FTin[i] = wave[i-indexstart];
    ifo[ifonr]->FTin[i] = wave[i];
  }
  
  free(wave);
  
} // End of templateLAL12()
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Compute waveform for a 12-parameter (one spin) LAL waveform
 *
 * \todo Merge this routine with templateLAL12() ?
 *
 * Compute h_+ and h_x from the parameters in par and interferometer information in ifo.
 * l is a pointer to get the lenght of the waveform computed, this length is also available as waveform->phi->data->length.
 */
// ****************************************************************************************************************************************************  
void LALHpHc12(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parSet *par, struct interferometer *ifo, int injectionWF, struct runPar run) {
  
  //  static LALStatus    mystatus;
  
  //  SimInspiralTable    injParams;
  //  PPNParamStruc       ppnParams;
  
  INT4        i;
  int                   lengthLAL;
  //REAL8       a1, a2, phi, shift;
  
  
  ////////////////////////////////////////////////////////////initialisation of memory/////////////////////////////////
  
  //  memset( &mystatus, 0, sizeof(LALStatus) );
  //  memset( waveform, 0, sizeof(CoherentGW));
  //  memset( &injParams, 0, sizeof(SimInspiralTable) );
  //  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  
  ////////////////////////////////////////////////////////////conversion between the parameter set of the Apostolatos waveform (parameter stored in par) to the parameter set used in LAL//////////////
  
  
  double pMc=0.0,pEta=0.0,pTc=0.0,pLogDl=0.0,pSpin1=0.0,pSpCosTh1=0.0,pRA=0.0,pLongi=0.0,pSinDec=0.0,pPhase=0.0,pSinThJ0=0.0,pPhiJ0=0.0,pSpPhi1=0.0,PNorder=0.0;
  
  if(injectionWF==1) {                                               // Then this is an injection waveform template
    pMc       = par->par[run.injRevID[61]];                                            // 61: Mc
    pEta      = par->par[run.injRevID[62]];                                            // 62: eta
    pTc       = par->par[run.injRevID[11]];                                            // 11: t_c
    pLogDl    = par->par[run.injRevID[22]];                                            // 22: log(d_L)
    pSpin1    = par->par[run.injRevID[71]];                                            // 71: a_spin1
    pSpCosTh1 = par->par[run.injRevID[72]];                                            // 72: cos(theta_spin1)
    pRA       = par->par[run.injRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.injRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pSinThJ0  = par->par[run.injRevID[53]];                                            // 53: sin(theta_J0)
    pPhiJ0    = par->par[run.injRevID[54]];                                            // 54: phi_J0
    pSpPhi1   = par->par[run.injRevID[73]];                                            // 73: phi_spin1    
    
    PNorder   = run.injectionPNorder;                                                  // Post-Newtonian order
  } else {                                                           // Then this is an MCMC waveform template
    pMc       = par->par[run.parRevID[61]];                                            // 61: Mc
    pEta      = par->par[run.parRevID[62]];                                            // 62: eta
    pTc       = par->par[run.parRevID[11]];                                            // 11: t_c
    pLogDl    = par->par[run.parRevID[22]];                                            // 22: log(d_L)  
    pSpin1    = par->par[run.parRevID[71]];                                            // 71: a_spin1           
    pSpCosTh1 = par->par[run.parRevID[72]];                                            // 72: cos(theta_spin1)
    pRA       = par->par[run.parRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)     
    pPhase    = par->par[run.parRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pSinThJ0  = par->par[run.parRevID[53]];                                            // 53: sin(theta_J0)
    pPhiJ0    = par->par[run.parRevID[54]];                                            // 54: phi_J0         
    pSpPhi1   = par->par[run.parRevID[73]];                                            // 73: phi_spin1    
    
    PNorder   = run.mcmcPNorder;                                                       // Post-Newtonian order
  }
  
  pLongi = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);    // RA -> 'lon'
  
  //printf(" LAL one-spin WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, sthJo: %f, phiJo: %f,  a1: %f, cth1: %f, phi1: %f\n",
  //injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pSinThJ0, pPhiJ0,  pSpin1, pSpCosTh1, pSpPhi1);
  
  double x;
  double m1=0.0,m2=0.0,M=0.0,mu=0.0;
  double cvec1[3],cvec2[3],cvec3[3];
  double tvec1[3],tvec2[3],tvec4[3],tvec6[3],tvec7[3];
  double n_L[3];
  double spin=0.0,samplerate=0.0,inversesamplerate=0.0;
  for(i=0;i<3;i++) {
    cvec1[i] = 0.0;
    cvec2[i] = 0.0;
    cvec3[i] = 0.0;
    tvec1[i] = 0.0;
    tvec2[i] = 0.0;
    //tvec3[i] = 0.0;
    tvec4[i] = 0.0;
    //tvec5[i] = 0.0;
    tvec6[i] = 0.0;
    tvec7[i] = 0.0;
    n_L[i] = 0.0;
  }
  
  
  double f_lower=ifo->lowCut;
  
  samplerate = (double)ifo->samplerate;
  inversesamplerate = 1.0/samplerate;
  
  
  
  double n_z[3] = {0.0,0.0,1.0};                                                                                        //North in global coordinates
  double normalvec[3];                                                                                    
  for(i=0;i<3;i++) normalvec[i] = ifo->normalvec[i];                                                             //Detector position normal vector = local zenith vector z'
  //  for(i=0;i<3;i++) rightvec[i] = ifo->rightvec[i];
  //  for(i=0;i<3;i++) leftvec[i] = ifo->leftvec[i];
  
  // for(i=0;i<3;i++) normalvec[i]=n_z[i];
  
  
  //double D_L = exp(pLogDl)*Mpcs;                                                                                    //Source luminosity distance, in seconds
  double coslati = sqrt(1.0-pSinDec*pSinDec);
  
  double n_N[3] = {cos(pLongi)*coslati,sin(pLongi)*coslati,pSinDec};                                                                            //n_N: Position unit vector = N^
  
  double sthJ0   = pSinThJ0;                                                                                        //n_J0: 'total' AM unit vector, J0^  (almost equal to the real J, see Eq.15)
  double cthJ0   = sqrt(1. - sthJ0*sthJ0);
  double n_J0[3] = { cos(pPhiJ0)*cthJ0 , sin(pPhiJ0)*cthJ0 , sthJ0 };
  
  par->NdJ = dotProduct(n_N,n_J0);                                                                                                                                                                              //Inclination of J_0; only for printing purposes, should be removed from this routine
  
  //Get individual masses from Mch and eta  CHECK: use McEta2masses()
  double root = sqrt(0.25-pEta);
  double fraction = (0.5-root) / (0.5+root);
  double inversefraction = 1.0/fraction;
  double Mc = pMc*M0;                                                                                               //Chirp mass in seconds
  x = exp(0.6*log(fraction));
  m1 = Mc * (pow(1.0+fraction,0.2) / x);
  m2 = Mc * (pow(1.0+inversefraction,0.2) * x);
  M = m1+m2;                                                                                                            // Eq.16a
  mu = m1*m2/M;                                                                                                         // Eq.16b
  spin = pSpin1*m1*m1;
  
  
  double cst5 = spin*sqrt(1.0-pSpCosTh1*pSpCosTh1);
  
  facVec(n_J0,-sthJ0,tvec1);      
  addVec(n_z,tvec1,tvec2);                                                                                                      //cvec1 = (n_z - J0^*cos(theta_J0))/sin(theta_J0)
  facVec(tvec2,1.0/cthJ0,cvec1);
  
  //Constant vector 2 for the construction of Eq.41e
  crossProduct(n_J0,n_z,tvec1);                                                                                               //cvec2 = (J0^ x z^) / sin(theta_J0)
  facVec(tvec1,1.0/cthJ0,cvec2);
  
  
  //Constant vector 3 for the construction of Eq.12
  facVec(n_N,-dotProduct(normalvec,n_N),tvec1);                                                                          //tvec1 = -N^(z^'.N^)
  addVec(normalvec,tvec1,cvec3);                                                                                         //cvec3 = z^' - N^(z^'.N^)
  
  
  double alpha=0.0;
  double omega_orb=0.0,l_L=0.0,Y=0.0,Gsq=0.0,g=0.0,slamL=0.0,clamL=0.0,LdotN=0.0;
  double cst4=0.0,x1=0.0,x2=0.0,x3=0.0;
  
  omega_orb=pi*f_lower;
  
  l_L = m1*m2*exp(-c3rd*log(omega_orb*M));
  
  
  Y = spin/l_L;                                                                                                    //Y = |S|/|L|, Eq.43
  Gsq = 1.0 + 2.0*pSpCosTh1*Y + Y*Y;                                                                                      //G^2, Eq.46
  g   = sqrt(Gsq);
  
  cst4 = l_L+pSpCosTh1*spin;
  x = mu*M;
  x1 = x*x*x;
  x = g*l_L;
  x2 = x*x*x;
  x3 = spin*spin*spin;
  alpha = pSpPhi1 - 5.0/(96.0*x1) * (1.0+0.75*m2/m1) * 
    (2.0*x2 - 3.0*pSpCosTh1*spin*cst4*g*l_L - 3.0*pSpCosTh1*x3*(1.0-pSpCosTh1*pSpCosTh1) * asinh(cst4/cst5));                                            //Eq.47
  
  slamL = cst5/(l_L*g);                                                                                            //sin(lambda_L), Eq.48a
  clamL = cst4/(l_L*g);                                                                                            //cos(lambda_L), Eq.48b
  
  //Construct Eq.41e
  facVec(n_J0,clamL,tvec1);                                                                                        //tvec1 = J0^*cos(lambda_L)
  facVec(cvec1,slamL*cos(alpha),tvec4);                                                                            //tvec4 = (n_z - J0^*cos(theta_J0))*sin(lambda_L)*cos(alpha)/sin(theta_J0)
  facVec(cvec2,slamL*sin(alpha),tvec6);                                                                            //tvec6 = (J0^ x z^) * sin(lambds_L)*sin(alpha)/sin(theta_J0)
  addVec(tvec1,tvec4,tvec7);                                                                                       //Construct Eq.41e
  addVec(tvec7,tvec6,n_L);
  //Eq.41e: n_L=L^
  
  LdotN = dotProduct(n_L,n_N);
  
  double r = pow(M/(omega_orb*omega_orb),1.0/3.0);
  double e = (16.0/5.0)*sqrt((M/r)*(M/r)*(M/r)) / ((1.0+(3.0/4.0)*m2/m1)*(1.0+2.0*pSpCosTh1*Y+Y*Y));
  double e2 = e*e;
  double e3 = e2*e;
  
  
  double n_S[3];
  
  n_S[0] = -(-n_J0[0] - n_J0[1]*n_L[2]*e + n_L[1]*n_J0[2]*e - n_L[0]*n_L[0]*n_J0[0]*e2 - n_L[0]*n_L[1]*n_J0[1]*e2 - n_L[0]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[0];
  n_S[1] = -(-n_J0[1] + n_J0[0]*n_L[2]*e - n_L[0]*n_J0[2]*e - n_L[0]*n_J0[0]*n_L[1]*e2 - n_L[1]*n_L[1]*n_J0[1]*e2 - n_L[1]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[1];
  n_S[2] = -(-n_J0[2] - n_J0[0]*n_L[1]*e + n_L[0]*n_J0[1]*e - n_L[0]*n_J0[0]*n_L[2]*e2 - n_L[1]*n_J0[1]*n_L[2]*e2 - n_L[2]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[2];
  
  n_S[0] = n_S[0]*g*l_L - l_L*n_L[0];
  n_S[1] = n_S[1]*g*l_L - l_L*n_L[1];
  n_S[2] = n_S[2]*g*l_L - l_L*n_L[2];
  
  double xloc[3],yloc[3],zloc[3];                        // coordinates in the global frame (in which N is defined) of the local vectors (e.g. z=N)                                                                          
  for(i=0;i<3;i++) zloc[i] = n_N[i];                                                             
  for(i=0;i<3;i++) yloc[i] = 0.0;
  for(i=0;i<3;i++) xloc[i] = 0.0;
  
  crossProduct(zloc,n_L,yloc);
  normalise(yloc);
  crossProduct(yloc,zloc,xloc);
  
  
  double n_Sloc[3];
  normalise(n_S);
  
  n_Sloc[0] = dotProduct(n_S,xloc);
  n_Sloc[1] = dotProduct(n_S,yloc);
  n_Sloc[2] = dotProduct(n_S,zloc);
  
  for(i=0;i<3;i++) n_S[i] = n_Sloc[i];
  
  
  n_S[0] = pSpin1*n_S[0];
  n_S[1] = pSpin1*n_S[1];
  n_S[2] = pSpin1*n_S[2];
  
  ////////////////////////////////////////////////////////////now we fill the injParam structure with the converted parameters//////////////
  
  injParams->mass1 = (float)(m1/M0);
  injParams->mass2 = (float)(m2/M0);
  
  injParams->f_final = (float)ifo->highCut;  //It seems injParams->f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams->f_lower = (float)f_lower;
  
  //Remember we're in the 12-par routine here
  char* waveformApproximant = (char*)calloc(128,sizeof(char));
  getWaveformApproximant("SpinTaylor",128,PNorder,waveformApproximant);  //Spinning
  //snprintf(waveformApproximant,128,"SpinTayloronePointFivePN"); //Set it manually
  //printf("\n  %s\n\n",waveformApproximant);
  
  snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",waveformApproximant);
  //snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s","SpinTayloronePointFivePN");
  //snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s","SpinTaylortwoPN");
  
  /* this is given in Mpc */    
  injParams->distance = (float)exp(pLogDl);//d_L;
  
  injParams->inclination = (float)acos(LdotN);
  
  injParams->spin1x = (float)n_S[0];
  injParams->spin1y = (float)n_S[1];
  injParams->spin1z = (float)n_S[2];
  
  injParams->spin2x = 0.0;
  injParams->spin2y = 0.0;
  injParams->spin2z = 0.0;
  
  // 4 parameters used after the computation of h+,x ********************//
  injParams->coa_phase = (float)pPhase;
  injParams->longitude = (float)pLongi;
  injParams->latitude = (float)asin(pSinDec);
  injParams->polarization = (float)pSpPhi1;    
  
  ppnParams->deltaT = inversesamplerate;//1.0 / 4096.0;
  
  
  //Print output:
  //  printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //     injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //     injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  /* --- now we can call the injection function --- */
  
  LALGenerateInspiral( status, waveform, injParams, ppnParams );
  if(status->statusCode) {
    fprintf(stderr, "\n\n   LALHpHc12:  ERROR generating waveform\n" );
    REPORTSTATUS(status);
    exit(1);
  }
  
  lengthLAL  = waveform->phi->data->length;
  *l = lengthLAL;
  
  //Print output:
  //printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //     injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //     injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  ///////////////////////////////////////////////////////at this point the structure waveform is still allocated in memory and will have to be freed. See LALfreedomSpin below//////////
  
  free(waveformApproximant);
  
  
} // End of LALHpHc12()
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Compute waveform for a 15-parameter (two spins) LAL waveform
 *
 * \todo Merge this routine with LALHpHc15() ? - Done: new routine now called templateLAL15()
 *
 * Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
 */
// ****************************************************************************************************************************************************  
void templateLAL15old(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
  int i=0;
  double samplerate=0.0,inversesamplerate=0.0;
  int length=0;
  
  samplerate = (double)ifo[ifonr]->samplerate;
  inversesamplerate = 1.0/samplerate;
  length     = ifo[ifonr]->samplesize;
  //printf("length = %d\n",length);
  
  // double hplusLAL[length+2];
  // double hcrossLAL[length+2];
  double *wave = (double*)calloc(length+2,sizeof(double));
  int lengthLAL = 0;
  
  
  /*LAL thingies needed. Have to be freed later*/
  
  static LALStatus    status;
  CoherentGW          waveform;
  SimInspiralTable    injParams;
  PPNParamStruc       ppnParams;
  
  memset( &status, 0, sizeof(LALStatus) );
  memset( &waveform, 0, sizeof(CoherentGW) );
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  ppnParams.deltaT   = inversesamplerate;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;
  
  
  
  // Compute h_+ and h_x
  //LALHpHc(&thewaveform, hplusLAL, hcrossLAL, &lengthLAL, length, par, ifo[ifonr], ifonr);
  LALHpHc15(&status, &waveform, &injParams, &ppnParams, &lengthLAL, par, ifo[ifonr], injectionWF, run);
  
  
  
  
  // Compute the detector response
  //double delay = LALFpFc(&thewaveform, wave, &lengthLAL, length, par, ifonr);
  /* double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); */
  //Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x)
  /* delay = delay; */ //MvdS: remove 'declared but never referenced' warnings
  
  // printf("LALdelay = %10.10f\n", delay);
  
  LALfreedomSpin(&waveform);
  
  // printf("localtc = %f\n", localtc);
  
  // localtc = ((par->par[2] - ifo[ifonr]->FTstart) - delay);  //Can't use par[2] anymore...
  
  //  int indexstart;
  //  indexstart = (int) (localtc*samplerate);// - (double)lengthLAL);   //MvdS: lengthLAL is used here, but seems never to be set in LALFpFc. Is it set in LALHpHc? Vivien: Yes it is but it is also available in thewaveform and wave
  //  if (indexstart<0) indexstart = 0;
  
  //printf("localtc2 = %f\n", localtc2);
  //printf("i1 = %d\n", i1);
  //printf("indexstart = %d\n" , indexstart);
  //printf("ifo[%d]->FTstart = %f\n", ifonr, ifo[ifonr]->FTstart);
  //printf("lengthLAL = %d\n", lengthLAL);
  //printf("length: %d,  lengthLAL: %d,  i1: %d,  i2: %d\n", length, lengthLAL,indexstart, indexstart+lengthLAL);
  
  for (i=0; i<length; ++i){
    //MvdS:  Reduce the number of if statements:
    //if(i<indexstart) ifo[ifonr]->FTin[i]   = 0.0;
    //if(i>=indexstart && i<(indexstart+lengthLAL)) ifo[ifonr]->FTin[i] = wave[i-indexstart];
    //if(i>=(indexstart+lengthLAL)) ifo[ifonr]->FTin[i]   = 0.0;
    //ifo[ifonr]->FTin[i] = 0.0;
    // if(i>=indexstart && i<(indexstart+lengthLAL)) ifo[ifonr]->FTin[i] = wave[i-indexstart];
    ifo[ifonr]->FTin[i] = wave[i];
  }
  
  free(wave);
  
} // End of templateLAL15old()
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Compute a waveform for a 15-parameter (two spins) LAL waveform
 *
 * \todo Merge this routine with templateLAL15() ?  - Done: new routine is templateLAL15(), old one templateLAL15old()
 *
 * Compute h_+ and h_x form the parameters in par and interferometer information in ifo.
 * l is a pointer to get the lenght of the waveform computed, this length is also available in waveform->phi->data->length.
 */
// ****************************************************************************************************************************************************  
void LALHpHc15(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parSet *par, 
               struct interferometer *ifo, int injectionWF, struct runPar run) 
{
  
  // static LALStatus    mystatus;
  
  // SimInspiralTable    injParams;
  // PPNParamStruc       ppnParams;
  
  // INT4        i;
  int                   lengthLAL;
  //REAL8       a1, a2, phi, shift;
  
  
  ////////////////////////////////////////////////////////////initialisation of memory/////////////////////////////////
  
  
  // memset( waveform, 0, sizeof(CoherentGW));
  // memset( &injParams, 0, sizeof(SimInspiralTable) );
  // memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  ////////////////////////////////////////////////////////////conversion/////////////////////////////////
  
  
  double f_lower=ifo->lowCut;
  double samplerate = (double)ifo->samplerate;
  double inversesamplerate = 1.0/samplerate;
  
  // Get the 15 waveform parameters from their array:
  double pMc=0.0,pEta=0.0,pTc=0.0,pLogDl=0.0,pRA=0.0,pLongi=0.0,pSinDec=0.0,pPhase=0.0,pCosI=0.0,pPsi=0.0;
  double pSpin1=0.0,pSpCosTh1=0.0,pSpPhi1=0.0,pSpin2=0.0,pSpCosTh2=0.0,pSpPhi2=0.0,PNorder=0.0;
  
  if(injectionWF==1) {                                               // Then this is an injection waveform template
    pTc       = par->par[run.injRevID[11]];                                            // 11: t_c
    pLogDl    = par->par[run.injRevID[22]];                                            // 22: log(d_L)
    pMc       = par->par[run.injRevID[61]];                                            // 61: Mc
    pEta      = par->par[run.injRevID[62]];                                            // 62: eta
    
    pRA       = par->par[run.injRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.injRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.injRevID[51]];                                            // 51: cos(inclination)
    pPsi      = par->par[run.injRevID[52]];                                            // 52: psi: polarisation angle
    
    pSpin1    = par->par[run.injRevID[71]];                                            // 71: a_spin1
    pSpCosTh1 = par->par[run.injRevID[72]];                                            // 72: cos(theta_spin1)
    pSpPhi1   = par->par[run.injRevID[73]];                                            // 73: phi_spin1    
    pSpin2    = par->par[run.injRevID[81]];                                            // 81: a_spin2
    pSpCosTh2 = par->par[run.injRevID[82]];                                            // 82: cos(theta_spin2)
    pSpPhi2   = par->par[run.injRevID[83]];                                            // 83: phi_spin2    
    
    PNorder   = run.injectionPNorder;                                                  // Post-Newtonian order
  } else {                                                           // Then this is an MCMC waveform template
    pTc       = par->par[run.parRevID[11]];                                            // 11: t_c
    pLogDl    = par->par[run.parRevID[22]];                                            // 22: log(d_L)
    pMc       = par->par[run.parRevID[61]];                                            // 61: Mc
    pEta      = par->par[run.parRevID[62]];                                            // 62: eta
    
    pRA       = par->par[run.parRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.parRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.parRevID[51]];                                            // 51: cos(inclination) of the binary
    pPsi      = par->par[run.parRevID[52]];                                            // 52: psi: polarisation angle of the binary
    
    pSpin1    = par->par[run.parRevID[71]];                                            // 71: a_spin1
    pSpCosTh1 = par->par[run.parRevID[72]];                                            // 72: cos(theta_spin1)
    pSpPhi1   = par->par[run.parRevID[73]];                                            // 73: phi_spin1    
    pSpin2    = par->par[run.parRevID[81]];                                            // 81: a_spin2
    pSpCosTh2 = par->par[run.parRevID[82]];                                            // 82: cos(theta_spin2)
    pSpPhi2   = par->par[run.parRevID[83]];                                            // 83: phi_spin2    
    
    PNorder   = run.mcmcPNorder;                                                       // Post-Newtonian order
  }
  
  pLongi = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);    // RA -> 'lon'
  
  //printf(" LAL two-spin WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, cos(i): %f, psi: %f,  a1: %f, cth1: %f, phi1: %f,  a2: %f, cth2: %f, phi2: %f\n",
  //injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pCosI, pPsi,  pSpin1, pSpCosTh1, pSpPhi1,  pSpin2, pSpCosTh2, pSpPhi2);
  
  
  // Get masses from Mch and eta:
  double m1,m2;
  McEta2masses(pMc,pEta,&m1,&m2);
  
  //printf(" Mc: %f,  M1: %f,  M2: %f,  Mtot: %f\n",pMc,m1,m2,m1+m2);
  
  
  ////////////////////////////////////////////////////////////now we fill the injParam structure with the parameters//////////////
  
  injParams->mass1 = (float)m1;
  injParams->mass2 = (float)m2;
  
  injParams->f_final = (float)ifo->highCut;  // It seems injParams->f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams->f_lower = (float)f_lower;
  
  //Remember we're in the 15-par routine here
  char* waveformApproximant = (char*)calloc(128,sizeof(char));
  getWaveformApproximant("SpinTaylor",128,PNorder,waveformApproximant);  //Spinning
  //snprintf(waveformApproximant,128,"%s","SpinTaylorthreePointFivePN"); //Set it manually
  //printf("\n  %s\n\n",waveformApproximant);
  
  snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",waveformApproximant);
  //snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s","SpinTayloronePointFivePN");
  //snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s","SpinTaylortwoPN");
  //snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s","SpinTaylorthreePointFivePN");
  
  // This is given in Mpc:
  injParams->distance = (float)exp(pLogDl); // d_L;
  
  injParams->inclination = (float)acos(pCosI);                      // Inclination of the binary
  
  double pSpsinth1 = sqrt(1.0 - pSpCosTh1*pSpCosTh1);
  injParams->spin1x = (float)(pSpin1 * pSpsinth1 * cos(pSpPhi1));
  injParams->spin1y = (float)(pSpin1 * pSpsinth1 * sin(pSpPhi1));
  injParams->spin1z = (float)(pSpin1 * pSpCosTh1);
  
  double pSpsinth2 = sqrt(1.0 - pSpCosTh2*pSpCosTh2);
  injParams->spin2x = (float)(pSpin2 * pSpsinth2 * cos(pSpPhi2));
  injParams->spin2y = (float)(pSpin2 * pSpsinth2 * sin(pSpPhi2));
  injParams->spin2z = (float)(pSpin2 * pSpCosTh2);
  
  // 4 parameters used after the computation of h+,x ********************//
  injParams->coa_phase = (float)pPhase;                                      // GW phase at coalescence
  injParams->longitude = (float)pLongi;                                      // 'Longitude'
  injParams->latitude = (float)asin(pSinDec);                                // Declination
  injParams->polarization = (float)pPsi;                                     // Polarisation angle
  
  REAL8 geocent_end_time = pTc;
  
  XLALGPSSetREAL8( &(injParams->geocent_end_time), geocent_end_time );
  
  ppnParams->deltaT = inversesamplerate;
  
  
  //Print output:
  //  printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //     injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //     injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  /* --- now we can call the injection function --- */
  
  LALGenerateInspiral( status, waveform, injParams, ppnParams );
  if(status->statusCode) {
    fprintf(stderr, "\n\n   LALHpHc15:  ERROR generating waveform\n" );
    REPORTSTATUS(status);
    exit(1);
  }
  // printf("ppnParams->tc = %f\n",ppnParams->tc);
  // LALInfo( status, ppnParams.termDescription );
  
  lengthLAL  = waveform->phi->data->length;
  *l = lengthLAL;
  
  //Print output:
  //printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //     injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //     injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  ///////////////////////////////////////////////////////at this point the structure waveform is still allocated in memory and will have to be freed. See LALfreedomSpin below//////////
  
  free(waveformApproximant);
  
  
} // End of LALHpHc15()
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Compute waveform for a 15-parameter (two spins) LAL waveform
 *
 * Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
 */
// ****************************************************************************************************************************************************  
void templateLAL15(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
  
  // Get the 15 waveform parameters from their array:
  double pMc=0.0,pEta=0.0,pTc=0.0,pDl=0.0,pRA=0.0,pLongi=0.0,pSinDec=0.0,pPhase=0.0,pCosI=0.0,pPsi=0.0;
  double pSpin1=0.0,pSpCosTh1=0.0,pSpPhi1=0.0,pSpin2=0.0,pSpCosTh2=0.0,pSpPhi2=0.0,PNorder=0.0;
  
  if(injectionWF==1) {                                               // Then this is an injection waveform template
	pTc       = par->par[run.injRevID[11]];                                            // 11: t_c
	pDl       = pow(par->par[run.injRevID[21]],c3rd);                                            // 21: d_L^3
	pMc       = pow(par->par[run.injRevID[65]],6);                                            // 65: Mc^(1/6)    
	pEta      = par->par[run.injRevID[62]];                                            // 62: eta
	  
    pRA       = par->par[run.injRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.injRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.injRevID[51]];                                            // 51: cos(inclination)
    pPsi      = par->par[run.injRevID[52]];                                            // 52: psi: polarisation angle
    
    pSpin1    = par->par[run.injRevID[71]];                                            // 71: a_spin1
    pSpCosTh1 = par->par[run.injRevID[72]];                                            // 72: cos(theta_spin1)
    pSpPhi1   = par->par[run.injRevID[73]];                                            // 73: phi_spin1    
    pSpin2    = par->par[run.injRevID[81]];                                            // 81: a_spin2
    pSpCosTh2 = par->par[run.injRevID[82]];                                            // 82: cos(theta_spin2)
    pSpPhi2   = par->par[run.injRevID[83]];                                            // 83: phi_spin2    
    
    PNorder   = run.injectionPNorder;                                                  // Post-Newtonian order
  } else {                                                           // Then this is an MCMC waveform template
	pTc       = par->par[run.parRevID[11]];                                            // 11: t_c
	pDl       = pow(par->par[run.parRevID[21]],c3rd);                                            // 21: d_L^3
	pMc       = pow(par->par[run.parRevID[65]],6);                                            // 65: Mc^(1/6)    
	pEta      = par->par[run.parRevID[62]];                                            // 62: eta
	  
    pRA       = par->par[run.parRevID[31]];                                            // 31: RA
    pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.parRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.parRevID[51]];                                            // 51: cos(inclination) of the binary
    pPsi      = par->par[run.parRevID[52]];                                            // 52: psi: polarisation angle of the binary
    
    pSpin1    = par->par[run.parRevID[71]];                                            // 71: a_spin1
    pSpCosTh1 = par->par[run.parRevID[72]];                                            // 72: cos(theta_spin1)
    pSpPhi1   = par->par[run.parRevID[73]];                                            // 73: phi_spin1    
    pSpin2    = par->par[run.parRevID[81]];                                            // 81: a_spin2
    pSpCosTh2 = par->par[run.parRevID[82]];                                            // 82: cos(theta_spin2)
    pSpPhi2   = par->par[run.parRevID[83]];                                            // 83: phi_spin2    
    
    PNorder   = run.mcmcPNorder;                                                       // Post-Newtonian order
  }
  
  pLongi = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);    // RA -> 'lon'
  
  // Get masses from Mch and eta:
  double m1,m2;
  McEta2masses(pMc,pEta,&m1,&m2);
  
  //printf(" LAL two-spin WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, cos(i): %f, psi: %f,  a1: %f, cth1: %f, phi1: %f,  a2: %f, cth2: %f, phi2: %f\n",
  //injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pCosI, pPsi,  pSpin1, pSpCosTh1, pSpPhi1,  pSpin2, pSpCosTh2, pSpPhi2);
  //printf(" Mc: %f,  M1: %f,  M2: %f,  Mtot: %f\n",pMc,m1,m2,m1+m2);
  
  
  // Cannot compute templates for Mtot >~ 146Mo (?).  Use 140Mo.
  int length = ifo[ifonr]->samplesize;
  int i=0;
  
  
  double samplerate = (double)ifo[ifonr]->samplerate;
  double inversesamplerate = 1.0/samplerate;
  double *wave = (double*)calloc(length+2,sizeof(double));
  
  // LAL thingies needed. Have to be freed later:
  static LALStatus    status;
  CoherentGW          waveform;
  SimInspiralTable    injParams;
  PPNParamStruc       ppnParams;
  
  memset( &status, 0, sizeof(LALStatus) );
  memset( &waveform, 0, sizeof(CoherentGW) );
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  ppnParams.deltaT   = inversesamplerate;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;
  
  
  ////////////////////////////////////////////////////////////now we fill the injParam structure with the parameters//////////////
  
  injParams.mass1 = (float)m1;
  injParams.mass2 = (float)m2;
  
  injParams.f_final = (float)ifo[ifonr]->highCut;  // It seems injParams.f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams.f_lower = (float)ifo[ifonr]->lowCut;
  
  // Remember we're in the 15-par routine here:
  char* waveformApproximant = (char*)calloc(128,sizeof(char));
  getWaveformApproximant("SpinTaylor",128,PNorder,waveformApproximant);  //Spinning
  //snprintf(waveformApproximant,128,"%s","SpinTaylorthreePointFivePN"); //Set it manually
  
  snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",waveformApproximant);
  
  // This is given in Mpc:
  injParams.distance = (float)pDl; // d_L;
  
  injParams.inclination = (float)acos(pCosI);                      // Inclination of the binary
  
  double pSpsinth1 = sqrt(1.0 - pSpCosTh1*pSpCosTh1);
  injParams.spin1x = (float)(pSpin1 * pSpsinth1 * cos(pSpPhi1));
  injParams.spin1y = (float)(pSpin1 * pSpsinth1 * sin(pSpPhi1));
  injParams.spin1z = (float)(pSpin1 * pSpCosTh1);
  
  double pSpsinth2 = sqrt(1.0 - pSpCosTh2*pSpCosTh2);
  injParams.spin2x = (float)(pSpin2 * pSpsinth2 * cos(pSpPhi2));
  injParams.spin2y = (float)(pSpin2 * pSpsinth2 * sin(pSpPhi2));
  injParams.spin2z = (float)(pSpin2 * pSpCosTh2);
  
  // 4 parameters used after the computation of h+,x ********************//
  injParams.coa_phase = (float)pPhase;                                      // GW phase at coalescence
  injParams.longitude = (float)pLongi;                                      // 'Longitude'
  injParams.latitude = (float)asin(pSinDec);                                // Declination
  injParams.polarization = (float)pPsi;                                     // Polarisation angle
  
  REAL8 geocent_end_time = pTc;
  
  XLALGPSSetREAL8( &(injParams.geocent_end_time), geocent_end_time );
  
  ppnParams.deltaT = inversesamplerate;
  
  
  //Print output:
  //  printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //     injParams.mass1,injParams.mass2,injParams.f_final,injParams.f_lower,injParams.distance,injParams.inclination,injParams.spin1x,injParams.spin1y,injParams.spin1z,
  //     injParams.spin2x,injParams.spin2y,injParams.spin2z,injParams.coa_phase,injParams.longitude,injParams.latitude,injParams.polarization);
  
  
  // Call the injection function:
  LALGenerateInspiral( &status, &waveform, &injParams, &ppnParams );
	if(status.statusCode) {
    if(injectionWF==1) {
      fprintf(stderr, "\n\n   templateLAL15():  ERROR generating injection waveform %s (too high mass?)\n   Aborting...\n\n", waveformApproximant);
      REPORTSTATUS(&status);
      exit(1);
    }
    for(i=0; i<length; ++i) ifo[ifonr]->FTin[i] = 0.0;
		
    //LALfreedomNoSpin(&waveform);  //Why does this give a seg.fault here, but not at the end of the routine?
    free(wave);
    free(waveformApproximant);
    return;
  }
  // LALInfo( status, ppnParams.termDescription );
  
  
  
  // Compute the detector response:
  /* double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); */ //Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x)
  /* delay = delay; */ //MvdS: remove 'declared but never referenced' warnings
  
  
  for (i=0; i<length; ++i) ifo[ifonr]->FTin[i] = wave[i];
  
  free(wave);
  free(waveformApproximant);
  LALfreedomSpin(&waveform);
  
} // End of templateLAL15()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Compute waveform for a 15-parameter (two spins) LAL PhenSpinTaylorRD waveform
 *
 * Use the LAL 3.5/2.5 PN spinning hybrid PhenSpinTaylorRD waveform, with 2 spinning objects (15 parameters). Phenomenological ring-down attached
 */
// ****************************************************************************************************************************************************  
void templateLALPhenSpinTaylorRD(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
	
	// Get the 15 waveform parameters from their array:
	double pMc=0.0,pEta=0.0,pTc=0.0,pLogDl=0.0,pRA=0.0,pLongi=0.0,pSinDec=0.0,pPhase=0.0,pCosI=0.0,pPsi=0.0;
	double pSpin1=0.0,pSpCosTh1=0.0,pSpPhi1=0.0,pSpin2=0.0,pSpCosTh2=0.0,pSpPhi2=0.0,PNorder=0.0;
	//REAL8 h1, h2, phi, shift;
	REAL4TimeSeries *htData = NULL;
	CHAR ifoname[LIGOMETA_IFO_MAX];
	INT8 waveformStartTime;
	INT4 dataLength, wfmLength;
	float *x1;
	
	if(injectionWF==1) {                                               // Then this is an injection waveform template
		pTc       = par->par[run.injRevID[11]];                                            // 11: t_c
		pLogDl    = par->par[run.injRevID[22]];                                            // 22: log(d_L)
		pMc       = par->par[run.injRevID[61]];                                            // 61: Mc
		pEta      = par->par[run.injRevID[62]];                                            // 62: eta
		
		pRA       = par->par[run.injRevID[31]];                                            // 31: RA
		pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
		pPhase    = par->par[run.injRevID[41]];                                            // 41: phi_c - GW phase at coalescence
		pCosI     = par->par[run.injRevID[51]];                                            // 51: cos(inclination)
		pPsi      = par->par[run.injRevID[52]];                                            // 52: psi: polarisation angle
		
		pSpin1    = par->par[run.injRevID[71]];                                            // 71: a_spin1
		pSpCosTh1 = par->par[run.injRevID[72]];                                            // 72: cos(theta_spin1)
		pSpPhi1   = par->par[run.injRevID[73]];                                            // 73: phi_spin1    
		pSpin2    = par->par[run.injRevID[81]];                                            // 81: a_spin2
		pSpCosTh2 = par->par[run.injRevID[82]];                                            // 82: cos(theta_spin2)
		pSpPhi2   = par->par[run.injRevID[83]];                                            // 83: phi_spin2    
		
		PNorder   = run.injectionPNorder;                                                  // Post-Newtonian order
	} else {                                                           // Then this is an MCMC waveform template
		pTc       = par->par[run.parRevID[11]];                                            // 11: t_c
		pLogDl    = par->par[run.parRevID[22]];                                            // 22: log(d_L)
		pMc       = par->par[run.parRevID[61]];                                            // 61: Mc
		pEta      = par->par[run.parRevID[62]];                                            // 62: eta
		
		pRA       = par->par[run.parRevID[31]];                                            // 31: RA
		pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)
		pPhase    = par->par[run.parRevID[41]];                                            // 41: phi_c - GW phase at coalescence
		pCosI     = par->par[run.parRevID[51]];                                            // 51: cos(inclination) of the binary
		pPsi      = par->par[run.parRevID[52]];                                            // 52: psi: polarisation angle of the binary
		
		pSpin1    = par->par[run.parRevID[71]];                                            // 71: a_spin1
		pSpCosTh1 = par->par[run.parRevID[72]];                                            // 72: cos(theta_spin1)
		pSpPhi1   = par->par[run.parRevID[73]];                                            // 73: phi_spin1    
		pSpin2    = par->par[run.parRevID[81]];                                            // 81: a_spin2
		pSpCosTh2 = par->par[run.parRevID[82]];                                            // 82: cos(theta_spin2)
		pSpPhi2   = par->par[run.parRevID[83]];                                            // 83: phi_spin2    
		
		PNorder   = run.mcmcPNorder;                                                       // Post-Newtonian order
	}
	
	pLongi = pRA; //XLALComputeDetAMResponse() takes RA. fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);    // RA -> 'lon'
	
	double temp_ra = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);
	printf("temp_ra=%f\n",temp_ra);
	
	// Get masses from Mch and eta:
	double m1,m2;
	McEta2masses(pMc,pEta,&m1,&m2);
	
	//printf(" LAL two-spin WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, cos(i): %f, psi: %f,  a1: %f, cth1: %f, phi1: %f,  a2: %f, cth2: %f, phi2: %f\n",
	//injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pCosI, pPsi,  pSpin1, pSpCosTh1, pSpPhi1,  pSpin2, pSpCosTh2, pSpPhi2);
	//printf(" Mc: %f,  M1: %f,  M2: %f,  Mtot: %f\n",pMc,m1,m2,m1+m2);
	
	
	// Cannot compute templates for Mtot >~ 146Mo (?).  Use 140Mo.
	int length = ifo[ifonr]->samplesize;
	int i=0;
	for ( i = 0; i < length; i++ ) {ifo[ifonr]->FTin[i] = 0.0;}
	
	REAL8 samplerate = (double)ifo[ifonr]->samplerate;
	double inversesamplerate = 1.0/samplerate;
	double *wave = (double*)calloc(length+2,sizeof(double));
	
	// LAL thingies needed. Have to be freed later:
	static LALStatus    status;
	CoherentGW          waveform;
	SimInspiralTable    injParams;
	PPNParamStruc       ppnParams;
	
	memset( &status, 0, sizeof(LALStatus) );
	memset( &waveform, 0, sizeof(CoherentGW) );
	memset( &injParams, 0, sizeof(SimInspiralTable) );
	memset( &ppnParams, 0, sizeof(PPNParamStruc) );
	
	ppnParams.deltaT   = inversesamplerate;
	ppnParams.lengthIn = 0;
	ppnParams.ppn      = NULL;
	
	
	////////////////////////////////////////////////////////////now we fill the injParam structure with the parameters//////////////
	
	injParams.mass1 = (float)m1;
	injParams.mass2 = (float)m2;
	
	injParams.f_final = (float)ifo[ifonr]->highCut;  // It seems injParams.f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
	injParams.f_lower = (float)ifo[ifonr]->lowCut;
	
	// Remember we're in the 15-par routine here:
	char* waveformApproximant = (char*)calloc(128,sizeof(char));
	getWaveformApproximant("PhenSpinTaylorRD",128,PNorder,waveformApproximant);  //Spinning
	//snprintf(waveformApproximant,128,"%s","SpinTaylorthreePointFivePN"); //Set it manually
	
	snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",waveformApproximant);
	
	// This is given in Mpc:
	injParams.distance = (float)exp(pLogDl); // d_L;
	
	injParams.inclination = (float)acos(pCosI);                      // Inclination of the binary
	
	double pSpsinth1 = sqrt(1.0 - pSpCosTh1*pSpCosTh1);
	injParams.spin1x = (float)(pSpin1 * pSpsinth1 * cos(pSpPhi1));
	injParams.spin1y = (float)(pSpin1 * pSpsinth1 * sin(pSpPhi1));
	injParams.spin1z = (float)(pSpin1 * pSpCosTh1);
	
	double pSpsinth2 = sqrt(1.0 - pSpCosTh2*pSpCosTh2);
	injParams.spin2x = (float)(pSpin2 * pSpsinth2 * cos(pSpPhi2));
	injParams.spin2y = (float)(pSpin2 * pSpsinth2 * sin(pSpPhi2));
	injParams.spin2z = (float)(pSpin2 * pSpCosTh2);
	
	// 4 parameters used after the computation of h+,x ********************//
	injParams.coa_phase = (float)pPhase;                                      // GW phase at coalescence
	ppnParams.phi = injParams.coa_phase;
	injParams.longitude = (float)pLongi;                                      // 'Longitude'
	injParams.latitude = (float)asin(pSinDec);                                // Declination
//	injParams.theta0       = injParams.longitude;         //?
//    injParams.phi0         = injParams.latitude;          //?
	injParams.polarization = (float)pPsi;                                     // Polarisation angle
	
	REAL8 geocent_end_time = pTc;
	
	XLALGPSSetREAL8( &(injParams.geocent_end_time), geocent_end_time );
	
	waveform.position.system=COORDINATESYSTEM_GEOGRAPHIC;
	//ppnParams.deltaT = inversesamplerate;
	
	
	//Print output:
	  //printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
	//  printf("  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f\n",
	//     injParams.mass1,injParams.mass2,injParams.f_final,injParams.f_lower,injParams.distance,injParams.inclination,injParams.spin1x,injParams.spin1y,injParams.spin1z,
	//     injParams.spin2x,injParams.spin2y,injParams.spin2z,injParams.coa_phase,injParams.longitude,injParams.latitude,injParams.polarization);
	
	//  printf("%8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f  %8.10f\n",
	//		 pMc,pEta,pTc,pLogDl,pRA,pSinDec,pCosI,pPhase,pPsi,pSpin1,pSpCosTh1,pSpPhi1,pSpin2,pSpCosTh2,pSpPhi2);	
	
	XLALSetErrno( 0 );
	// Call the injection function:
	LALGenerateInspiral( &status, &waveform, &injParams, &ppnParams );
	if(status.statusCode || waveform.h == NULL) {
		if(injectionWF==1) {
			fprintf(stderr, "\n\n   templateLALPhenSpinTaylorRD():  ERROR generating injection waveform %s (too high mass?)\n   Aborting...\n\n", waveformApproximant);
			REPORTSTATUS(&status);
			exit(1);
		}
		if(waveform.h == NULL){fprintf(stderr,"waveform.h == NULL\n");}
		REPORTSTATUS(&status);
		for(i=0; i<length; ++i) ifo[ifonr]->FTin[i] = 0.0;
    //printf("wave at zero\n");
		//LALfreedomPhenSpinTaylorRD(&waveform);
		//XLALDestroyREAL4Vector ( htData->data);
		//LALFree(htData);
		
		//LALfreedomNoSpin(&waveform);  //Why does this give a seg.fault here, but not at the end of the routine?
		free(wave);
		free(waveformApproximant);
		return;
	}
	// LALInfo( status, ppnParams.termDescription );
	
//	DetectorResponse detector;   // the detector in question 
	//memset( &detector, 0, sizeof( DetectorResponse ) );
	
	
	//LALDetector site;
//	detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
	
	//site = lalCachedDetectors[LALDetectorIndexLLODIFF];
	
	//detector.site = &site;
//	detector.transfer = NULL;
	
	switch ( ifonr )
	{
		case 0:
		//	*(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
			strcpy(ifoname, "H1");
			break;
		case 1:
		//	*(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
			strcpy(ifoname, "L1");
			break;
		case 2:
		//	*(detector.site) = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
			strcpy(ifoname, "V1");
			break;
		default:
		//	LALFree( detector.site );
		//	detector.site = NULL;
			fprintf(stderr, "\n\n   Unknown detector site ifonr=%d\n   Aborting...\n\n",ifonr);
			REPORTSTATUS(&status);
			exit(1);
			break;
	}
	
	waveformStartTime = XLALGPSToINT8NS( &(injParams.geocent_end_time) );
	waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

	XLALINT8NSToGPS( &(waveform.h->epoch), waveformStartTime );
//	memcpy( &(waveform.f->epoch), &(waveform.h->epoch),
//		   sizeof(LIGOTimeGPS) );
//	memcpy( &(waveform.phi->epoch), &(waveform.h->epoch),
//		   sizeof(LIGOTimeGPS) );
//	memcpy( &(waveform.a->epoch), &(waveform.a->epoch),
//		   sizeof(LIGOTimeGPS) );

	wfmLength = waveform.h->data->length;
	dataLength = 2*wfmLength;
	x1 = (float *) LALMalloc(sizeof(x1)*dataLength);
	/*
	 * We are using functions from NRWaveInject which do not take a vector
	 * with alternating h+ & hx but rather one that stores h+ in the first
	 * half and hx in the second half.
	 *
	 * We also must multiply the strain by the distance to be compatible
	 * with NRWaveInject.
	 *
	 */
	for( i = 0; i < dataLength; i++)
	{
        x1[i] = waveform.h->data->data[i]*injParams.distance;
	}
	for( i = 0; i < wfmLength; i++)
	{
		waveform.h->data->data[i] = x1[2*i];
		waveform.h->data->data[wfmLength+i] = x1[2*i+1];
	}
	
	LALFree(x1);
	
	waveform.h->data->vectorLength = wfmLength;
	
	INT4 sampleRate = (INT4)samplerate;
	
	// Compute the detector response:
	//htData = XLALCalculateNRStrain( waveform.h, &injParams, ifoname, sampleRate );
	
	
	
	LALDetector            det;
	double                 fplus;
	double                 fcross;
	double                 tDelay;
	//REAL4TimeSeries       *htData = NULL;
	UINT4                  vecLength, k;
	InterferometerNumber   ifoNumber = LAL_UNKNOWN_IFO;
	
	/* get the detector information */
	memset( &det, 0, sizeof(LALDetector) );
	ifoNumber = (InterferometerNumber)XLALIFONumber( ifoname );
	XLALReturnDetector( &det, ifoNumber );
	
	if( ifoNumber == LAL_UNKNOWN_IFO )
	{
		XLALPrintWarning( "Unknown ifo! Injecting signal overhead!\n" );
		fplus  = 1.0;
		fcross = 0.0;
		tDelay = 0.0;
	}
	else
	{
		/* compute detector response */
		XLALComputeDetAMResponse(&fplus, &fcross, (const REAL4 (*)[3])det.response, injParams.longitude,
								 injParams.latitude, injParams.polarization, injParams.end_time_gmst);
		
		/* calculate the time delay */
		tDelay = XLALTimeDelayFromEarthCenter( det.location, injParams.longitude,
											  injParams.latitude, &(injParams.geocent_end_time) );
	}
	
	/* create htData */
	htData = LALCalloc(1, sizeof(*htData));
	if (!htData)
	{
		printf("hData unassigned\n");
		//ABORTXLAL( &status );
		exit(1);
	}
	vecLength = waveform.h->data->vectorLength;
	htData->data = XLALCreateREAL4Vector( vecLength );
	if ( ! htData->data )
	{
		LALFree( htData );
		printf("hData->data unassigned\n");
		//ABORTXLAL( &status );
		exit(1);
	}
	
	/* store the htData */
	
	/* add timedelay to inj->geocent_end_time */
	{
		REAL8 tEnd;
		tEnd =  XLALGPSGetREAL8(&(injParams.geocent_end_time) );
		tEnd += tDelay;
		XLALGPSSetREAL8(&(htData->epoch), tEnd );
	}
	/* Using XLALGPSAdd is bad because that changes the GPS time! */
	/*  htData->epoch = *XLALGPSAdd( &(inj->geocent_end_time), tDelay ); */
	
	htData->deltaT = waveform.h->deltaT;
	htData->sampleUnits = lalADCCountUnit;
	
	for ( k = 0; k < vecLength; ++k )
	{
	  htData->data->data[k] = (float)(fplus * waveform.h->data->data[k]  +
				   fcross * waveform.h->data->data[vecLength + k]) / injParams.distance;
	}
	
	/*interpolate to given sample rate */
	//htData = XLALInterpolateNRWave( htData, sampleRate);


	
	
	
	
	
	REAL4TimeSeries *ret=NULL;
	REAL8 deltaTin, deltaTout, r, y_1, y_2;
	REAL8 tObs; /* duration of signal */
	UINT4 lo, numPoints;
	
	deltaTin = htData->deltaT;
	tObs = deltaTin * htData->data->length;
	
	/* length of output vector */
	numPoints = (UINT4) (sampleRate * tObs);
	
	/* allocate memory */
	ret = LALCalloc(1, sizeof(*ret));
	if (!ret)
	{
		printf("ret unassigned\n");
		//ABORTXLAL( &status );
		exit(1);
	}
	
	ret->data = XLALCreateREAL4Vector( numPoints );
	if (! ret->data)
	{
		printf("ret->data unassigned\n");
		//ABORTXLAL( &status );
		exit(1);
	}
	
	ret->deltaT = 1./sampleRate;
	deltaTout = ret->deltaT;
	
	/* copy stuff from in which should be the same */
	ret->epoch = htData->epoch;
	ret->f0 = htData->f0;
	ret->sampleUnits = htData->sampleUnits;
	strcpy(ret->name, htData->name);
	
	/* go over points of output vector and interpolate linearly
     using closest points of input */
	for (k = 0; k < numPoints; k++) {
		
		lo = (UINT4)( k*deltaTout / deltaTin);
		
		/* y_1 and y_2 are the input values at x1 and x2 */
		/* here we need to make sure that we don't exceed
		 bounds of input vector */
		if ( lo < htData->data->length - 1) {
			y_1 = htData->data->data[lo];
			y_2 = htData->data->data[lo+1];
			
			/* we want to calculate y_2*r + y_1*(1-r) where
			 r = (x-x1)/(x2-x1) */
			r = k*deltaTout / deltaTin - lo;
			
			ret->data->data[k] = (float)(y_2 * r + y_1 * (1 - r));
		}
		else {
			ret->data->data[k] = 0.0;
		}
	}
	
	if ( !ret )
	{
		printf("ret unassigned\n");
		//ABORTXLAL( &status );
		exit(1);
	}
	
	REAL8 offset;
	XLALFindNRCoalescenceTime( &offset, waveform.h);
	XLALGPSAdd( &(ret->epoch), -offset);
	
	REAL4TimeSeries signalvec;        // GW signal 
	
	memset( &signalvec, 0, sizeof(REAL4TimeSeries) );
	
	/* set the parameters for the signal time series */
	signalvec.deltaT = ret->deltaT;
	
	signalvec.sampleUnits = lalADCCountUnit;
	//signal.f0 = 0.0;
	signalvec.data = NULL;
	/* simulate the detectors response to the inspiral */
	LALSCreateVector( &status, &(signalvec.data), (UINT4)length );
	XLALGPSSetREAL8( &(signalvec.epoch), ifo[ifonr]->FTstart);	
	for ( i = 0; i < (int)signalvec.data->length; i++ ) {signalvec.data->data[i]=0.0;}
	LALSSInjectTimeSeries( &status, &signalvec, ret );
	
	for ( i = 0; i < (int)signalvec.data->length && i < length; i++ ) {ifo[ifonr]->FTin[i] = signalvec.data->data[i];
																//printf("signal.data->data[%d] = %10.10e\n",i,signal.data->data[i]);
																}
	
	free(wave);
	free(waveformApproximant);
	LALfreedomPhenSpinTaylorRD(&waveform);
	XLALDestroyREAL4Vector ( htData->data);
	LALFree(htData);
	XLALDestroyREAL4Vector ( ret->data);
	LALFree(ret);
	LALSDestroyVector( &status, &( signalvec.data ) );
//	if ( detector.site ) LALFree( detector.site );
//	LALFree( detector.transfer );

	
} // End of templateLALPhenSpinTaylorRD()
// ****************************************************************************************************************************************************  









// ****************************************************************************************************************************************************  
/**
 * \brief Compute waveform template for a LAL non-spinning inspiral waveform
 *
 * Uses GeneratePPN approximant.
 */
// ****************************************************************************************************************************************************  
void templateLALnonSpinning(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
  int i=0;
  
  // Get the 9 waveform parameters from their array:
  double pMc=0.0,pEta=0.0,pTc=0.0,pDl=0.0,pRA=0.0,pLongi=0.0,pSinDec=0.0,pPhase=0.0,pCosI=0.0,pPsi=0.0,PNorder=0.0;
  if(injectionWF==1) {                                               // Then this is an injection waveform template
    pTc       = par->par[run.injRevID[11]];                                            // 11: t_c
    pDl       = pow(par->par[run.injRevID[21]],c3rd);                                            // 21: d_L^3
    pMc       = pow(par->par[run.injRevID[65]],6);                                            // 65: Mc^(1/6)
    pEta      = par->par[run.injRevID[62]];                                            // 62: eta
    
    pRA       = par->par[run.injRevID[31]];                                            // 31: 'longi' := RA (?)
    pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.injRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.injRevID[51]];                                            // 51: cos(inclination)
    pPsi      = par->par[run.injRevID[52]];                                            // 52: psi: polarisation angle
    
    PNorder   = run.injectionPNorder;                                                  // Post-Newtonian order
  } else {                                                           // Then this is an MCMC waveform template
    pTc       = par->par[run.parRevID[11]];                                            // 11: t_c
    pDl       = pow(par->par[run.parRevID[21]],c3rd);                                            // 21: d_L^3
    pMc       = pow(par->par[run.parRevID[65]],6);                                            // 65: Mc^(1/6)
    pEta      = par->par[run.parRevID[62]];                                            // 62: eta
    
    pRA       = par->par[run.parRevID[31]];                                            // 31: longi := RA (?)
    pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)
    pPhase    = par->par[run.parRevID[41]];                                            // 41: phi_c - GW phase at coalescence
    pCosI     = par->par[run.parRevID[51]];                                            // 51: cos(inclination) of the binary
    pPsi      = par->par[run.parRevID[52]];                                            // 52: psi: polarisation angle of the binary
    
    PNorder   = run.mcmcPNorder;                                                       // Post-Newtonian order
  }
  
  pLongi = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);    // RA -> 'lon'
  //pLongi = pRA;                                          // For non-spinning LAL waveforms 'longi' = RA (?)  NO! - It's probably 'longitude'
  
  
  //printf(" LAL nS WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, cos(i): %f, psi: %f\n",
  //injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pCosI, pPsi);
  
  
  // Get masses from Mch and eta:
  double m1=0.0,m2=0.0;
  McEta2masses(pMc,pEta,&m1,&m2);
  
  //printf(" LAL nS WF pars:  injWF: %i, Mc: %f, eta: %f, M1: %f, M2: %f, Mtot: %f\n",injectionWF,pMc, pEta, m1, m2, m1+m2);
  
  int length = ifo[ifonr]->samplesize;
  
  
  // LAL structs needed. Have to be freed later
  static LALStatus    status;
  CoherentGW          waveform;  // i.e. output
  SimInspiralTable    injParams;  // Physical input parameters
  PPNParamStruc       ppnParams;  // 'non-physical' input parameters, e.g. f cuts, Delta-t, etc.
  
  memset( &status, 0, sizeof(LALStatus) );
  memset( &waveform, 0, sizeof(CoherentGW) );
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  double f_lower=ifo[ifonr]->lowCut;
  double samplerate = (double)ifo[ifonr]->samplerate;
  double inversesamplerate = 1.0/samplerate;
  double *wave = (double*)calloc(length+2,sizeof(double));
  
  
  
  
  // Store waveform family and pN order in injParams.waveform
  // Remember we're in the non-spinning LAL routine here
  char* waveformApproximant = (char*)calloc(128,sizeof(char));
  getWaveformApproximant("GeneratePPN",128,PNorder,waveformApproximant);  //Non-spinning
  //getWaveformApproximant("SpinTaylor",128,PNorder,waveformApproximant);  //Spinning
  //printf("\n  %s\n\n",waveformApproximant);
  
  snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",waveformApproximant);
  Approximant injapprox;
  injapprox = XLALGetApproximantFromString(injParams.waveform);
  if(injapprox!=GeneratePPN) fprintf(stderr,"\n *** Warning:  not using GeneratePPN approximant causes incoherent injections. Used: %s\t%s***\n",waveformApproximant,injParams.waveform);
  
  // Fill injParam with the waveform parameters:
  injParams.mass1 = (float)m1;
  injParams.mass2 = (float)m2;
  injParams.mchirp = (float)pMc;  // Get a seg.fault when setting both pairs?!?!?
  injParams.eta = (float)pEta;
  
  injParams.distance = (float)pDl;                                  // Distance in Mpc
  injParams.inclination = (float)acos(pCosI);                               // Inclination of the binary
  
  
  // 4 parameters used after the computation of h+,x ********************//
  injParams.coa_phase = (float)pPhase;                                      // GW phase at coalescence
  injParams.longitude = (float)pLongi;                                      // 'Longitude'  CHECK: is this actually RA?!?!?
  injParams.latitude = (float)asin(pSinDec);                                // Declination
  injParams.polarization = (float)pPsi;                                     // Polarisation angle
  
  injParams.f_final = (float)ifo[ifonr]->highCut;  // It seems injParams.f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams.f_lower = (float)f_lower;
  
  ppnParams.fStartIn = (float)f_lower;  //May be needed here as well...(?)
  ppnParams.deltaT   = inversesamplerate;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;
  
  
  
  
  REAL8 geocent_end_time = pTc;
  XLALGPSSetREAL8( &injParams.geocent_end_time, geocent_end_time );
  
  
  // Call the injection function; compute h_+ and h_x:
  LALGenerateInspiral(&status, &waveform, &injParams, &ppnParams );
  if(status.statusCode) {
    if(injectionWF==1) {
      fprintf(stderr, "\n\n   templateLALnonSpinning():  ERROR generating injection waveform %s (too high mass?)\n   Aborting...\n\n", waveformApproximant);
      REPORTSTATUS(&status);
      exit(1);
    }
    for(i=0; i<length; ++i) ifo[ifonr]->FTin[i] = 0.0;
    
    //LALfreedomNoSpin(&waveform);  //Why does this give a seg.fault here, but not at the end of the routine?
    free(wave);
    free(waveformApproximant);
    return;
  }
  
  
  // Compute the detector response  -  Is this done by LALGenerateInspiral for the non-spinning case?
  /* double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); */ //Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x)
  /* delay = delay; */ //MvdS: remove 'declared but never referenced' warnings
  
  
  
  for(i=0; i<length; ++i) ifo[ifonr]->FTin[i] = wave[i];
  
  LALfreedomNoSpin(&waveform);
  free(wave);
  free(waveformApproximant);
  
} // End of templateLALnonSpinning()
// ****************************************************************************************************************************************************  

















// ****************************************************************************************************************************************************  
/**
 * \brief Compute detector response for a given detector and given h_+,h_x
 *
 * Compute the detector response for a given detector (ifonr) and h_+,h_x. the structure waveform must already hold the computed values of h+,x (or just a1, a2, phi and shift as a function of time)
 */
// ****************************************************************************************************************************************************  
double LALFpFc(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, double *wave, int length, struct parSet *par, struct interferometer *ifo, int ifonr) 
{
  par->par[2] = par->par[2]; //MvdS: CHECK: remove 'never referenced' warning (or remove par from argument list)
  
  // static LALStatus stat;     // status structure
  
  // memset( &stat, 0, sizeof(LALStatus) );
  
  int i;
  
  DetectorResponse detector;   // the detector in question 
  //memset( &detector, 0, sizeof( DetectorResponse ) );
  
  
  //LALDetector site;
  detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  
  if(ifonr==0) *(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  if(ifonr==1) *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(ifonr==2) *(detector.site) = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  
  //site = lalCachedDetectors[LALDetectorIndexLLODIFF];
  
  //detector.site = &site;
  detector.transfer = NULL;
  
  /* set up units for the transfer function */
  /*    RAT4 negOne = { -1, 0 };
        LALUnit unit;
        LALUnitPair pair;
        pair.unitOne = &lalADCCountUnit;
        pair.unitTwo = &lalStrainUnit;
        LALUnitRaise( &stat, &unit, pair.unitTwo, &negOne );
        pair.unitTwo = &unit;
        LALUnitMultiply( &stat, &(detector.transfer->sampleUnits), &pair );*/
  
  
  //detector.transfer = NULL;
  
  
  /* invert the response function to get the transfer function */
  /*  LALCCreateVector( &stat, &( detector.transfer->data ), resp->data->length );
      
      LALCCreateVector( &stat, &unity, resp->data->length );
      for ( k = 0; k < resp->data->length; ++k ) 
      {
      unity->data[k].re = 1.0;
      unity->data[k].im = 0.0;
      }
      
      LALCCVectorDivide( &stat, detector.transfer->data, unity,
      resp->data );
      
      LALCDestroyVector( &stat, &unity );
  */
  
  
  
  
  REAL4TimeSeries signalvec;        // GW signal 
  
  memset( &signalvec, 0, sizeof(REAL4TimeSeries) );
  
  //REAL4TimeSeries chan;        // channel
  
  //memset( &chan, 0, sizeof(REAL4TimeSeries) );
  
  
  //signalvec.epoch.gpsSeconds = (INT4)par->par[2];  // Can't use par[i] anymore...
  //signalvec.epoch.gpsNanoSeconds = (INT4)(100000000.0*(par->par[2] - (double)signalvec.epoch.gpsSeconds));  // Can't use par[i] anymore...
  
  //waveform->f->epoch = waveform->phi->epoch = waveform->a->epoch = signalvec.epoch; 
  INT8 waveformStartTime;
  
  //  waveformStartTime = par->par[2];  // Can't use par[i] anymore...
  
  /* set the start times for injection */
  //LALINT8toGPS( &stat, &(waveform->a->epoch), &waveformStartTime );
  //  LALFloatToGPS( &stat, &(waveform->a->epoch), &waveformStartTime);
  
  
  waveformStartTime = XLALGPSToINT8NS( &(injParams->geocent_end_time) );
  
  
  waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams->tc );
  
  XLALINT8NSToGPS(&(waveform->a->epoch), waveformStartTime );
  
  //memcpy( &(waveform->f->epoch), &(waveform->a->epoch), sizeof(LIGOTimeGPS) );
  //memcpy( &(waveform->phi->epoch), &(waveform->a->epoch), sizeof(LIGOTimeGPS) );

  waveform->f->epoch = waveform->a->epoch;
  waveform->phi->epoch = waveform->a->epoch;

  if ( waveform->shift )
    {
      waveform->shift->epoch = waveform->a->epoch;
    }
  
  
	/* set the start time of the signal vector to the start time of ifo[ifonr]->FTstart */
  
  
  
  
  
  /* set the parameters for the signal time series */
  signalvec.deltaT = waveform->phi->deltaT;
  
  signalvec.sampleUnits = lalADCCountUnit;
  //signalvec.f0 = 0.0;
  signalvec.data = NULL;
  /* simulate the detectors response to the inspiral */
  LALSCreateVector( status, &(signalvec.data), (UINT4)length );
  XLALGPSSetREAL8( &(signalvec.epoch), ifo->FTstart);
  
  
  /* set the parameters for the signal time series */
  //   chan.deltaT = waveform->phi->deltaT;
  
  //   chan.sampleUnits = lalADCCountUnit;
  // chan.f0 = 0.0;
  // chan.data = NULL;
  /* simulate the detectors response to the inspiral */
  //  LALSCreateVector( status, &(chan.data), (UINT4)length );
  //  LALFloatToGPS( status, &(chan.epoch), &(ifo->FTstart));
  
  waveform->position.system=COORDINATESYSTEM_GEOGRAPHIC;
  
  LALSimulateCoherentGW( status, &signalvec, waveform, &detector );//////////////////this is were F+,x are being computed.
  
  //LALFloatToGPS( status, &(chan.epoch), &(ifo->FTstart));
  
  
  
  
  
  
  
  // signalvec.deltaT = waveform->phi->deltaT;
  // signalvec.f0 = 0.0;
  // signalvec.data = NULL;
  
  //      LALSSInjectTimeSeries(status, &chan, &signalvec );
  
  for ( i = 0; i < (int)signalvec.data->length && i < length; i++ ){
    
    //printf("%d\t%10.10e\n", i, chan.data->data[i]);
    
    wave[i] = signalvec.data->data[i]; // wave is my array of doubles to send back the waveform to the rest of SPINspiral.
	  //printf("%d\t%10.10e\n", i, signalvec.data->data[i]);
  }
  //printf("1000\t%10.10e\n",wave[1000]);
  /*********TIME DELAY***********/
  
  
  REAL8            delay=0.0;
  //  DetTimeAndASource     det1_and_source;
  //  LALPlaceAndGPS        det1_and_gps;
  
  //  det1_and_gps.p_detector = detector.site;
  // det1_and_gps.p_gps      = &(signalvec.epoch);
  
  // det1_and_source.p_det_and_time = &det1_and_gps;
  //  det1_and_source.p_source       = &(waveform->position);
  
  // LALTimeDelayFromEarthCenter(status, &delay, &det1_and_source);
  
  LALSDestroyVector( status, &( signalvec.data ) );
  // LALSDestroyVector( status, &( chan.data ) );
  
  // if(waveform->position.system==COORDINATESYSTEM_EQUATORIAL) printf("youpi\n");
  // printf("position=%f,%f\n",waveform->position.longitude,waveform->position.latitude);
  
  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );
  
  return delay;
  
} // End of LALFpFc()
// ****************************************************************************************************************************************************  


























// ****************************************************************************************************************************************************  
/**
 * \brief Compose the waveform approximant from the family name and pN order
 */
// ****************************************************************************************************************************************************  
void getWaveformApproximant(const char* familyName, int length, double PNorder, char* waveformApproximant) {
  int PNorderTimesTwo = (int)rint(PNorder*2.0);
  switch(PNorderTimesTwo) {
  case 2:
    snprintf(waveformApproximant,length,"%s%s",familyName,"onePN");
    break;
    
  case 3:
    snprintf(waveformApproximant,length,"%s%s",familyName,"onePointFivePN");
    break;
    
  case 4:
    snprintf(waveformApproximant,length,"%s%s",familyName,"twoPN");
    break;
    
  case 5:
    snprintf(waveformApproximant,length,"%s%s",familyName,"twoPointFivePN");
    break;
    
  case 6:
    snprintf(waveformApproximant,length,"%s%s",familyName,"threePN");
    break;
    
  case 7:
    snprintf(waveformApproximant,length,"%s%s",familyName,"threePointFivePN");
    break;
    
  case 8:
    snprintf(waveformApproximant,length,"%s%s",familyName,"fourPN");
    break;
    
  default:
    fprintf(stderr, "\n\n   ERROR:  pN order%4.1f is not (yet) supported, aborting.\n\n\n",PNorder);
    exit(1);
  }
}
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Free spinning LAL variables
 */
// ****************************************************************************************************************************************************  
void LALfreedomSpin(CoherentGW *waveform) {
  // Free LAL stuff  
  static LALStatus stat;     // status structure
  
  memset( &stat, 0, sizeof(LALStatus) );
  
  LALSDestroyVectorSequence(&stat, &( waveform->a->data ));
  LALSDestroyVector(&stat, &( waveform->f->data ));
  LALDDestroyVector(&stat, &( waveform->phi->data ));
  LALSDestroyVector(&stat, &( waveform->shift->data ));
  
  LALFree( waveform->a ); waveform->a = NULL;
  LALFree( waveform->f ); waveform->f = NULL;
  LALFree( waveform->phi) ; waveform->phi = NULL;
  LALFree( waveform->shift ); waveform->shift = NULL;
  
} // End of LALfreedomSpin
// ****************************************************************************************************************************************************  

// ****************************************************************************************************************************************************  
/**
 * \brief Free spinning LAL variables
 */
// ****************************************************************************************************************************************************  
void LALfreedomPhenSpinTaylorRD(CoherentGW *waveform) {
	// Free LAL stuff  
	static LALStatus stat;     // status structure
	
	memset( &stat, 0, sizeof(LALStatus) );
	
	LALSDestroyVectorSequence(&stat, &( waveform->h->data ));
	LALSDestroyVector(&stat, &( waveform->f->data ));
	LALDDestroyVector(&stat, &( waveform->phi->data ));
	LALSDestroyVector(&stat, &( waveform->shift->data ));
	
	LALFree( waveform->h ); waveform->h = NULL;
	LALFree( waveform->f ); waveform->f = NULL;
	LALFree( waveform->phi) ; waveform->phi = NULL;
	LALFree( waveform->shift ); waveform->shift = NULL;
	
} // End of LALfreedomSpin
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Free non-spinning LAL variables
 */
// ****************************************************************************************************************************************************  
void LALfreedomNoSpin(CoherentGW *waveform) {
  // Free LAL stuff  
  static LALStatus stat;     // status structure
  
  memset( &stat, 0, sizeof(LALStatus) );
  
  LALSDestroyVectorSequence(&stat, &( waveform->a->data ));
  LALSDestroyVector(&stat, &( waveform->f->data ));
  LALDDestroyVector(&stat, &( waveform->phi->data ));
  //LALSDestroyVector(&stat, &( waveform->shift->data ));
  
  LALFree( waveform->a );
  LALFree( waveform->f ); 
  LALFree( waveform->phi) ;
  LALFree( waveform->shift );
  
} // End of LALfreedomNoSpin
// ****************************************************************************************************************************************************  



