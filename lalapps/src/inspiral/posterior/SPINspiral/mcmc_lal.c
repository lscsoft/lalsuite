/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   mcmc_lal.c:                interfaces to LAL routines
   
   
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


#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>

//#include <lal/LALConstants.h>

#include <mcmc_lal.h>

#include <lal/DetectorSite.h>

#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
//#include <lal/DetectorSite.h>

//////////////////////////////////////////
#include <mcmc.h>
//////////////////////////////////////////








void templateLAL12(struct parset *par, struct interferometer *ifo[], int ifonr)
//Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
{
  int i=0;
  double localtc=0.0,samplerate=0.0,inversesamplerate=0.0;
  int length=0;
  
  samplerate = (double)ifo[ifonr]->samplerate;
  inversesamplerate = 1.0/samplerate;
  length     = ifo[ifonr]->samplesize;
  //printf("length = %d\n",length);
  
  // double hplusLAL[length+2];
  // double hcrossLAL[length+2];
  //double wave[length+2];
  double *wave = (double*)calloc(length+2,sizeof(double));  //MvdS: should this make a difference? Vivien: no it shouldn't.
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
  LALHpHc12(&status, &waveform, &injParams, &ppnParams, &lengthLAL, par, ifo[ifonr]);  //MvdS: ifonr never used. This routine computes and returns hplusLAL, hcrossLAL, which are never used... However, this information should also be contained in thewaveform
																					//Vivien: ifonr is only used in a commented printf to know which interferometer called the routine. Just for debugging purposes
            																		//Vivien: hplusLAL and hcrossLAL are indeed unecessary (was before I used the structure thewaveform)


  
  // Compute the detector response
  //double delay = LALFpFc(&thewaveform, wave, &lengthLAL, length, par, ifonr);
  double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); //MvdS: lengthLAL never used or set. Uses waveforms in thewaveform to compute the detector response in wave (?)
																	//Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x)
		
  // printf("LALdelay = %10.10f\n", delay);
  
  LALfreedom(&waveform);
  
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
  

}







void LALHpHc12(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parset *par, struct interferometer *ifo) {
  // Compute h_+ and h_x form the parameters in par and interferometer information in ifo. l is a pointer to get the lenght of the waveform computed. this length is also available in waveform->phi->data->length
  
//  static LALStatus    mystatus;
  
//  SimInspiralTable    injParams;
//  PPNParamStruc       ppnParams;
  
  INT4        i;
  int			lengthLAL;
  //REAL8       a1, a2, phi, shift;
  pi=M_PI;
  
  
  ////////////////////////////////////////////////////////////initialisation of memory/////////////////////////////////
  
//  memset( &mystatus, 0, sizeof(LALStatus) );
//  memset( waveform, 0, sizeof(CoherentGW));
//  memset( &injParams, 0, sizeof(SimInspiralTable) );
//  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  
  ////////////////////////////////////////////////////////////conversion between the parameter set of the Apostolatos waveform (parameter stored in par) to the parameter set used in LAL//////////////
  
  
  double pmc       =par->par[0];
  double peta      =par->par[1];
  double ptc       =par->par[2];
  double plogdl    =par->par[3];
  double pspin     =par->par[4];
  double pkappa    =par->par[5];
  double plongi    = fmod(longitude(par->par[6],GMST(ptc))+mtpi,tpi);  //par[6] contains RA
  double psinlati  =par->par[7];
  double pphase    =par->par[8];
  double psinthJ0  =par->par[9];
  double pphiJ0    =par->par[10];
  double palpha    =par->par[11];
  
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
  
  
  //double D_L = exp(plogdl)*Mpcs;                                                                                    //Source luminosity distance, in seconds
  double coslati = sqrt(1.0-psinlati*psinlati);
                                        
  double n_N[3] = {cos(plongi)*coslati,sin(plongi)*coslati,psinlati};										//n_N: Position unit vector = N^
  
  double sthJ0   = psinthJ0;                                                                                        //n_J0: 'total' AM unit vector, J0^  (almost equal to the real J, see Eq.15)
  double cthJ0   = sqrt(1. - sthJ0*sthJ0);
  double n_J0[3] = { cos(pphiJ0)*cthJ0 , sin(pphiJ0)*cthJ0 , sthJ0 };
  
  par->NdJ = dotproduct(n_N,n_J0);																						//Inclination of J_0; only for printing purposes, should be removed from this routine
  
  //Get individual masses from Mch and eta  CHECK: use mceta2masses()
  double root = sqrt(0.25-peta);
  double fraction = (0.5-root) / (0.5+root);
  double inversefraction = 1.0/fraction;
  double Mc = pmc*M0;                                                                                               //Chirp mass in seconds
  x = exp(0.6*log(fraction));
  m1 = Mc * (pow(1.0+fraction,0.2) / x);
  m2 = Mc * (pow(1.0+inversefraction,0.2) * x);
  M = m1+m2;                                                                                                            // Eq.16a
  mu = m1*m2/M;                                                                                                         // Eq.16b
  spin = pspin*m1*m1;
  
  
  double cst5 = spin*sqrt(1.0-pkappa*pkappa);
     
  facvec(n_J0,-sthJ0,tvec1);      
  addvec(n_z,tvec1,tvec2);                                                                                                      //cvec1 = (n_z - J0^*cos(theta_J0))/sin(theta_J0)
  facvec(tvec2,1.0/cthJ0,cvec1);
  
  //Constant vector 2 for the construction of Eq.41e
  crossproduct(n_J0,n_z,tvec1);                                                                                               //cvec2 = (J0^ x z^) / sin(theta_J0)
  facvec(tvec1,1.0/cthJ0,cvec2);
  
  
  //Constant vector 3 for the construction of Eq.12
  facvec(n_N,-dotproduct(normalvec,n_N),tvec1);                                                                          //tvec1 = -N^(z^'.N^)
  addvec(normalvec,tvec1,cvec3);                                                                                         //cvec3 = z^' - N^(z^'.N^)
  
  
  double alpha=0.0;
  double omega_orb=0.0,l_L=0.0,Y=0.0,Gsq=0.0,G=0.0,slamL=0.0,clamL=0.0,LdotN=0.0;
  double cst4=0.0,x1=0.0,x2=0.0,x3=0.0;
  
  omega_orb=pi*f_lower;
  
  l_L = m1*m2*exp(-c3rd*log(omega_orb*M));
  
    
  Y = spin/l_L;                                                                                                    //Y = |S|/|L|, Eq.43
  Gsq = 1.0 + 2.0*pkappa*Y + Y*Y;                                                                                      //G^2, Eq.46
  G   = sqrt(Gsq);
  
  cst4 = l_L+pkappa*spin;
  x = mu*M;
  x1 = x*x*x;
  x = G*l_L;
  x2 = x*x*x;
  x3 = spin*spin*spin;
  alpha = palpha - 5.0/(96.0*x1) * (1.0+0.75*m2/m1) * 
    (2.0*x2 - 3.0*pkappa*spin*cst4*G*l_L - 3.0*pkappa*x3*(1.0-pkappa*pkappa) * asinh(cst4/cst5));                                            //Eq.47
  
  slamL = cst5/(l_L*G);                                                                                            //sin(lambda_L), Eq.48a
  clamL = cst4/(l_L*G);                                                                                            //cos(lambda_L), Eq.48b
  
  //Construct Eq.41e
  facvec(n_J0,clamL,tvec1);                                                                                        //tvec1 = J0^*cos(lambda_L)
  facvec(cvec1,slamL*cos(alpha),tvec4);                                                                            //tvec4 = (n_z - J0^*cos(theta_J0))*sin(lambda_L)*cos(alpha)/sin(theta_J0)
  facvec(cvec2,slamL*sin(alpha),tvec6);                                                                            //tvec6 = (J0^ x z^) * sin(lambds_L)*sin(alpha)/sin(theta_J0)
  addvec(tvec1,tvec4,tvec7);                                                                                       //Construct Eq.41e
  addvec(tvec7,tvec6,n_L);
  //Eq.41e: n_L=L^
  
  LdotN = dotproduct(n_L,n_N);
  
  double r = pow(M/(omega_orb*omega_orb),1.0/3.0);
  double e = (16.0/5.0)*sqrt((M/r)*(M/r)*(M/r)) / ((1.0+(3.0/4.0)*m2/m1)*(1.0+2.0*pkappa*Y+Y*Y));
  double e2 = e*e;
  double e3 = e2*e;
  
  
  double n_S[3];
  
  n_S[0] = -(-n_J0[0] - n_J0[1]*n_L[2]*e + n_L[1]*n_J0[2]*e - n_L[0]*n_L[0]*n_J0[0]*e2 - n_L[0]*n_L[1]*n_J0[1]*e2 - n_L[0]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[0];
  n_S[1] = -(-n_J0[1] + n_J0[0]*n_L[2]*e - n_L[0]*n_J0[2]*e - n_L[0]*n_J0[0]*n_L[1]*e2 - n_L[1]*n_L[1]*n_J0[1]*e2 - n_L[1]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[1];
  n_S[2] = -(-n_J0[2] - n_J0[0]*n_L[1]*e + n_L[0]*n_J0[1]*e - n_L[0]*n_J0[0]*n_L[2]*e2 - n_L[1]*n_J0[1]*n_L[2]*e2 - n_L[2]*n_L[2]*n_J0[2]*e2) / (1.0 + n_L[0]*e3 + n_L[1]*n_L[1]*e2 + n_L[2]*n_L[2]*e2);  //-n_L[2];
  
  n_S[0] = n_S[0]*G*l_L - l_L*n_L[0];
  n_S[1] = n_S[1]*G*l_L - l_L*n_L[1];
  n_S[2] = n_S[2]*G*l_L - l_L*n_L[2];
  
  double xloc[3],yloc[3],zloc[3];                        // coordinates in the global frame (in which N is defined) of the local vectors (e.g. z=N)                                                                          
  for(i=0;i<3;i++) zloc[i] = n_N[i];                                                             
  for(i=0;i<3;i++) yloc[i] = 0.0;
  for(i=0;i<3;i++) xloc[i] = 0.0;
  
  crossproduct(zloc,n_L,yloc);
  normalise(yloc);
  crossproduct(yloc,zloc,xloc);
  
  
  double n_Sloc[3];
  normalise(n_S);
  
  n_Sloc[0] = dotproduct(n_S,xloc);
  n_Sloc[1] = dotproduct(n_S,yloc);
  n_Sloc[2] = dotproduct(n_S,zloc);
  
  for(i=0;i<3;i++) n_S[i] = n_Sloc[i];
  
  
  n_S[0] = pspin*n_S[0];
  n_S[1] = pspin*n_S[1];
  n_S[2] = pspin*n_S[2];
  
  ////////////////////////////////////////////////////////////now we fill the injParam structure with the converted parameters//////////////
  
  injParams->mass1 = (float)(m1/M0);
  injParams->mass2 = (float)(m2/M0);
  
  injParams->f_final = (float)ifo->highCut;  //It seems injParams->f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams->f_lower = (float)f_lower;
  
  snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"SpinTaylorthreePointFivePN");//"SpinTaylortwoPN");
  
  /* this is given in Mpc */    
  injParams->distance = (float)exp(plogdl);//d_L;
  
  injParams->inclination = (float)acos(LdotN);
  
  injParams->spin1x = (float)n_S[0];
  injParams->spin1y = (float)n_S[1];
  injParams->spin1z = (float)n_S[2];
  
  injParams->spin2x = 0.0;
  injParams->spin2y = 0.0;
  injParams->spin2z = 0.0;
  
  // 4 parameters used after the computation of h+,x ********************//
  injParams->coa_phase = (float)pphase;
  injParams->longitude = (float)plongi;
  injParams->latitude = (float)asin(psinlati);
  injParams->polarization = (float)palpha;    
  
  ppnParams->deltaT = inversesamplerate;//1.0 / 4096.0;
  
  
  //Print output:
  //  printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //	 injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //	 injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  /* --- now we can call the injection function --- */
  
  LALGenerateInspiral( status, waveform, injParams, ppnParams );
  if ( status->statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform\n" );
      exit( 1 );
    }
  
  lengthLAL  = waveform->phi->data->length;
  *l = lengthLAL;
  
  //Print output:
  //printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //	 injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //	 injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  ///////////////////////////////////////////////////////at this point the structure waveform is still allocated in memory and will have to be freed. See LALfreedom below//////////
  
}


void templateLAL15(struct parset *par, struct interferometer *ifo[], int ifonr)
//Use the LAL 3.5/2.5 PN spinning waveform, with 2 spinning objects (15 parameters)
{
  int i=0;
  double localtc=0.0,samplerate=0.0,inversesamplerate=0.0;
  int length=0;
  
  samplerate = (double)ifo[ifonr]->samplerate;
  inversesamplerate = 1.0/samplerate;
  length     = ifo[ifonr]->samplesize;
  //printf("length = %d\n",length);
  
  // double hplusLAL[length+2];
  // double hcrossLAL[length+2];
  //double wave[length+2];
  double *wave = (double*)calloc(length+2,sizeof(double));  //MvdS: should this make a difference? Vivien: no it shouldn't.
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
  LALHpHc15(&status, &waveform, &injParams, &ppnParams, &lengthLAL, par, ifo[ifonr]);  //MvdS: ifonr never used. This routine computes and returns hplusLAL, hcrossLAL, which are never used... However, this information should also be contained in thewaveform
																					//Vivien: ifonr is only used in a commented printf to know which interferometer called the routine. Just for debugging purposes
            																		//Vivien: hplusLAL and hcrossLAL are indeed unecessary (was before I used the structure thewaveform)


  
  // Compute the detector response
  //double delay = LALFpFc(&thewaveform, wave, &lengthLAL, length, par, ifonr);
  double delay = LALFpFc(&status, &waveform, &injParams, &ppnParams, wave, length, par, ifo[ifonr], ifonr); //MvdS: lengthLAL never used or set. Uses waveforms in thewaveform to compute the detector response in wave (?)
																	//Vivien: lentghLAL is set in LALinteface.c But is is also availble in the structure thewaveform (which holds h+,x) and the structure wave (which holds F+,x)
		
  // printf("LALdelay = %10.10f\n", delay);
  
  LALfreedom(&waveform);
  
  // printf("localtc = %f\n", localtc);
 
 // localtc = ((par->par[2] - ifo[ifonr]->FTstart) - delay);
  
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
  
}



void LALHpHc15(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, int *l, struct parset *par, struct interferometer *ifo) {
  // Compute h_+ and h_x form the parameters in par and interferometer information in ifo. l is a pointer to get the lenght of the waveform computed. this length is also available in waveform->phi->data->length
  
 // static LALStatus    mystatus;
  
 // SimInspiralTable    injParams;
 // PPNParamStruc       ppnParams;
  
 // INT4        i;
  int			lengthLAL;
  //REAL8       a1, a2, phi, shift;
  pi=M_PI;
  
  
  ////////////////////////////////////////////////////////////initialisation of memory/////////////////////////////////
  
  
 // memset( waveform, 0, sizeof(CoherentGW));
 // memset( &injParams, 0, sizeof(SimInspiralTable) );
 // memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  
  ////////////////////////////////////////////////////////////conversion/////////////////////////////////


  double f_lower=ifo->lowCut;
  
  double samplerate = (double)ifo->samplerate;
  double inversesamplerate = 1.0/samplerate;
  


  //Get masses from Mch and eta  CHECK: use mceta2masses()
  
  double m1,m2;
  mceta2masses(par->par[0],par->par[1],&m1,&m2);
  //printf("%f\t%f\n",m1,m2);
  /*
  double root = sqrt(0.25-par->par[1]);
  double fraction = (0.5-root) / (0.5+root);
  double inversefraction = 1.0/fraction;                                                                                             //Chirp mass in seconds
  double x = exp(0.6*log(fraction));
  double m1 = par->par[0] * (pow(1.0+fraction,0.2) / x);
  double m2 = par->par[0] * (pow(1.0+inversefraction,0.2) * x);*/

  ////////////////////////////////////////////////////////////now we fill the injParam structure with the parameters//////////////
  
  injParams->mass1 = (float)m1;
  injParams->mass2 = (float)m2;
  
  injParams->f_final = (float)ifo->highCut;  //It seems injParams->f_final gets overwritten by LALGenerateInspiral; it's an output parameter rather than input. This will also somewhat affect SNR comparisons with the Apostolatos waveform.
  injParams->f_lower = (float)f_lower;
  
  snprintf(injParams->waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"SpinTayloronePointFivePN");//"SpinTaylortwoPN");
  
  /* this is given in Mpc */    
  injParams->distance = (float)exp(par->par[3]);//d_L;
  
  injParams->inclination = (float)acos(par->par[6]);
  
  double sin1=sqrt(1.0 - par->par[10]*par->par[10]);
  
  injParams->spin1x = (float)(par->par[9]*sin1*cos(par->par[11]));
  injParams->spin1y = (float)(par->par[9]*sin1*sin(par->par[11]));
  injParams->spin1z = (float)(par->par[9]*par->par[10]);
  
  double sin2=sqrt(1.0 - par->par[13]*par->par[13]);
  
  injParams->spin2x = (float)(par->par[12]*sin2*cos(par->par[14]));
  injParams->spin2y = (float)(par->par[12]*sin2*sin(par->par[14]));
  injParams->spin2z = (float)(par->par[12]*par->par[13]);
    
  // 4 parameters used after the computation of h+,x ********************//
  injParams->coa_phase = (float)par->par[7];
  injParams->longitude = (float)fmod(longitude(par->par[4],GMST(par->par[2]))+mtpi,tpi);  //par[4] contains RA 
  injParams->latitude = (float)asin(par->par[5]);
  injParams->polarization = (float)par->par[8];
  
  REAL8 geocent_end_time = par->par[2];
  
  XLALGPSSetREAL8( &(injParams->geocent_end_time), geocent_end_time );
  
   ppnParams->deltaT = inversesamplerate;
    
  
  //Print output:
  //  printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //	 injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //	 injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  /* --- now we can call the injection function --- */
  
  LALGenerateInspiral( status, waveform, injParams, ppnParams );
  if ( status->statusCode )
    {
      fprintf( stderr, "LALSTPNWaveformTest: error generating waveform\n" );
      exit( 1 );
    }
 // printf("ppnParams->tc = %f\n",ppnParams->tc);
 // LALInfo( status, ppnParams.termDescription );
    
  lengthLAL  = waveform->phi->data->length;
  *l = lengthLAL;
  
  //Print output:
  //printf("  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",
  //	 injParams->mass1,injParams->mass2,injParams->f_final,injParams->f_lower,injParams->distance,injParams->inclination,injParams->spin1x,injParams->spin1y,injParams->spin1z,
  //	 injParams->spin2x,injParams->spin2y,injParams->spin2z,injParams->coa_phase,injParams->longitude,injParams->latitude,injParams->polarization);
  
  
  ///////////////////////////////////////////////////////at this point the structure waveform is still allocated in memory and will have to be freed. See LALfreedom below//////////
  
}








// Compute the detector response for a given detector (ifonr) and h_+,h_x. the structure waveform must already hold the computed values of h+,x (or just a1, a2, phi and shift as a function of time)

double LALFpFc(LALStatus *status, CoherentGW *waveform, SimInspiralTable *injParams, PPNParamStruc *ppnParams, double *wave, int length, struct parset *par, struct interferometer *ifo, int ifonr) {
  
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
  detector.ephemerides = NULL;

  /* set up units for the transfer function */
/*	RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LALUnitRaise( &stat, &unit, pair.unitTwo, &negOne );
    pair.unitTwo = &unit;
    LALUnitMultiply( &stat, &(detector.transfer->sampleUnits), &pair );*/

 
  //detector.transfer = NULL;
  //detector.ephemerides = NULL;
  
  
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

  
   
  
  REAL4TimeSeries signal;        // GW signal 
  
  memset( &signal, 0, sizeof(REAL4TimeSeries) );
  
  //REAL4TimeSeries chan;        // channel
  
  //memset( &chan, 0, sizeof(REAL4TimeSeries) );

  
  //signal.epoch.gpsSeconds = (INT4)par->par[2];
  //signal.epoch.gpsNanoSeconds = (INT4)(100000000.0*(par->par[2] - (double)signal.epoch.gpsSeconds));

  //waveform->f->epoch = waveform->phi->epoch = waveform->a->epoch = signal.epoch; 
  INT8 waveformStartTime;
  
 //  waveformStartTime = par->par[2];

/* set the start times for injection */
      //LALINT8toGPS( &stat, &(waveform->a->epoch), &waveformStartTime );
	//  LALFloatToGPS( &stat, &(waveform->a->epoch), &waveformStartTime);

      waveformStartTime = XLALGPSToINT8NS( &(injParams->geocent_end_time) );
	  
	  

      waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams->tc );
	  

	  
	 XLALINT8NSToGPS(&(waveform->a->epoch), waveformStartTime );
	  
      memcpy( &(waveform->f->epoch), &(waveform->a->epoch), sizeof(LIGOTimeGPS) );
      memcpy( &(waveform->phi->epoch), &(waveform->a->epoch), sizeof(LIGOTimeGPS) );


      /* set the start time of the signal vector to the start time of ifo[ifonr]->FTstart */
	  

	  
	  

      /* set the parameters for the signal time series */
      signal.deltaT = waveform->phi->deltaT;
      
      signal.sampleUnits = lalADCCountUnit;
  signal.f0 = 0.0;
  signal.data = NULL;
      /* simulate the detectors response to the inspiral */
      LALSCreateVector( status, &(signal.data), (UINT4)length );
	  XLALGPSSetREAL8( &(signal.epoch), ifo->FTstart);
	  
	        /* set the parameters for the signal time series */
   //   chan.deltaT = waveform->phi->deltaT;
      
   //   chan.sampleUnits = lalADCCountUnit;
 // chan.f0 = 0.0;
 // chan.data = NULL;
      /* simulate the detectors response to the inspiral */
    //  LALSCreateVector( status, &(chan.data), (UINT4)length );
	//  LALFloatToGPS( status, &(chan.epoch), &(ifo->FTstart));
  
waveform->position.system=COORDINATESYSTEM_GEOGRAPHIC;

	LALSimulateCoherentGW( status, &signal, waveform, &detector );//////////////////this is were F+,x are being computed.

//LALFloatToGPS( status, &(chan.epoch), &(ifo->FTstart));
	  





  
 // signal.deltaT = waveform->phi->deltaT;
 // signal.f0 = 0.0;
 // signal.data = NULL;
      
//	  LALSSInjectTimeSeries(status, &chan, &signal );
	  
  for ( i = 0; i < signal.data->length && i < length; i++ ){
  
   //printf("%d\t%10.10e\n", i, chan.data->data[i]);
  
    wave[i] = signal.data->data[i]; // wave is my array of doubles to send back the waveform to the rest of SPINspiral.
  }
  
  /*********TIME DELAY***********/
  
  
  REAL8            delay=0.0;
//  DetTimeAndASource     det1_and_source;
//  LALPlaceAndGPS        det1_and_gps;
  
//  det1_and_gps.p_detector = detector.site;
 // det1_and_gps.p_gps      = &(signal.epoch);
  
 // det1_and_source.p_det_and_time = &det1_and_gps;
//  det1_and_source.p_source       = &(waveform->position);
  
 // LALTimeDelayFromEarthCenter(status, &delay, &det1_and_source);
  
  LALSDestroyVector( status, &( signal.data ) );
 // LALSDestroyVector( status, &( chan.data ) );
  
  // if(waveform->position.system==COORDINATESYSTEM_EQUATORIAL) printf("youpi\n");
  // printf("position=%f,%f\n",waveform->position.longitude,waveform->position.latitude);
  
    if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  return delay;
  
}






void LALfreedom(CoherentGW *waveform) {
  // Free LAL stuff  
  static LALStatus stat;     /* status structure */
  
  memset( &stat, 0, sizeof(LALStatus) );
  
  LALSDestroyVectorSequence(&stat, &( waveform->a->data ));
  LALSDestroyVector(&stat, &( waveform->f->data ));
  LALDDestroyVector(&stat, &( waveform->phi->data ));
  LALSDestroyVector(&stat, &( waveform->shift->data ));
  
  LALFree( waveform->a );
  LALFree( waveform->f ); 
  LALFree( waveform->phi) ;
  LALFree( waveform->shift );
  
}

