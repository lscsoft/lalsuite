/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_templates.c:    waveform templates and related routines
   
   
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



#include <SPINspiral.h>


/**
 * \file SPINspiral_templates.c
 * \brief Contains routines that handle waveform-template generation
 */



// ****************************************************************************************************************************************************  
/**
 * \brief Compute an inspiral waveform
 * 
 * Compute an inspiral waveform template.
 * waveformVersion determines the template to use.
 * injectionWF indicates whether this is an injection waveform (1) or not (0).
 */
// ****************************************************************************************************************************************************  
void waveformTemplate(struct parSet *par, struct interferometer *ifo[], int ifonr, int waveformVersion, int injectionWF, struct runPar run)
{
  
  //CHECK: test - remove this
  /*
    int i=0;
    printf("\n\n\n*** waveformTemplate(): %i %i  ", waveformVersion, injectionWF);
    for(i=0;i<par->nPar;i++) {
    printf(" %i:%6.3lf", i, par->par[i]);
    }
    printf("\n");
  */
  if(waveformVersion==1) {
    templateApostolatos(par, ifo, ifonr, injectionWF, run);  // Apostolatos 12-parameter template
  } else if(waveformVersion==2) {
    templateLAL12(par, ifo, ifonr, injectionWF, run);  // LAL 12-parameter template
  } else if(waveformVersion==3) {
    templateLAL15(par, ifo, ifonr, injectionWF, run);  // LAL 15-parameter template
  } else if(waveformVersion==4) {
    templateLALnonSpinning(par, ifo, ifonr, injectionWF, run);  // LAL non-spinning template
  } else {
    fprintf(stderr,"\n\n   ERROR:  waveformTemplate(): waveformVersion %i not defined!\n\n",waveformVersion);
    exit(1);
  }
  
} // End of waveformTemplate()
// ****************************************************************************************************************************************************  







// ****************************************************************************************************************************************************  
/**
 * \brief Compute an Apostolatos, 12-parameter, single spin, simple-precession waveform
 * 
 * Compute a spinning, 'simple-precession' template in restricted 1.5PN order with 1 spin (Apostolatos et al., 1994, PhRvD..49.6274A).
 * The output vector ifo[ifonr]->FTin is of length ifo[ifonr]->samplesize,  starting at 'tstart'(?) and with resolution ifo[ifonr]->samplerate.
 */
// ****************************************************************************************************************************************************  
void templateApostolatos(struct parSet *par, struct interferometer *ifo[], int ifonr, int injectionWF, struct runPar run)
{
  
  double pMc=0.0,pEta=0.0,pTc=0.0,pSpin1=0.0,pSpCosTh1=0.0,pRA=0.0,pSinDec=0.0,pPhase=0.0,pSinThJ0=0.0,pPhiJ0=0.0,pSpPhi1=0.0;
  double pLongi=0.0,pDl=0.0; //,pLogDl=0.0;
  
  if(injectionWF==1) {                                               // Then this is an injection waveform template:
    pMc       = par->par[run.injRevID[61]];                                           // 61: Mc
    pEta      = par->par[run.injRevID[62]];                                           // 62: eta
    pTc       = par->par[run.injRevID[11]];                                           // 11: t_c
    if(run.injParUse[21]) pDl = exp(log(par->par[run.injRevID[21]])/3.0);             // 21: (d_L)^3 -> d_L
    if(run.injParUse[22]) pDl = exp(par->par[run.injRevID[22]]);                      // 22: log(d_L) -> d_L
    //pLogDl    = par->par[run.injRevID[22]];                                           // 22: log(d_L) 
    pSpin1    = par->par[run.injRevID[71]];                                           // 71: a_spin1
    pSpCosTh1 = par->par[run.injRevID[72]];                                           // 72: cos(theta_spin1)
    pRA       = par->par[run.injRevID[31]];                                           // 31: RA
    pSinDec   = par->par[run.injRevID[32]];                                           // 32: sin(Dec)
    pPhase    = par->par[run.injRevID[41]];                                           // 41: phi_c
    pSinThJ0  = par->par[run.injRevID[53]];                                           // 53: sin(theta_J0)
    pPhiJ0    = par->par[run.injRevID[54]];                                           // 54: phi_J0
    pSpPhi1   = par->par[run.injRevID[73]];                                           // 73: phi_spin1    
    
  } else {                                                           // Then this is an MCMC waveform template:
    pMc       = par->par[run.parRevID[61]];                                           // 61: Mc
    pEta      = par->par[run.parRevID[62]];                                           // 62: eta
    pTc       = par->par[run.parRevID[11]];                                           // 11: t_c
    if(run.mcmcParUse[21]) pDl = exp(log(par->par[run.parRevID[21]]/3.0));             // 21: (d_L)^3 -> d_L
    if(run.mcmcParUse[22]) pDl = exp(par->par[run.parRevID[22]]);                      // 22: log(d_L) -> d_L
    //pLogDl    = par->par[run.parRevID[22]];                                           // 22: log(d_L) 
    pSpin1    = par->par[run.parRevID[71]];                                           // 71: a_spin1            
    pSpCosTh1 = par->par[run.parRevID[72]];                                           // 72: cos(theta_spin1)
    pRA       = par->par[run.parRevID[31]];                                           // 31: RA
    pSinDec   = par->par[run.parRevID[32]];                                           // 32: sin(Dec)     
    pPhase    = par->par[run.parRevID[41]];                                           // 41: phi_c           
    pSinThJ0  = par->par[run.parRevID[53]];                                           // 53: sin(theta_J0)
    pPhiJ0    = par->par[run.parRevID[54]];                                           // 54: phi_J0          
    pSpPhi1   = par->par[run.parRevID[73]];                                           // 73: phi_spin1    
  }
  
  pLongi    = fmod(longitude(pRA, GMST(pTc)) + mtpi, tpi);   // RA -> 'lon'
  
  
  //printf(" Apo WF pars:  injWF: %i, Mc: %f, eta: %f, tc: %f, logD: %f, RA: %f, dec: %f, phi: %f, sin(th_J0): %f, phi_J0: %f, spin1: %f, spcos(th): %f, sp_phi1: %f  \n",
  //injectionWF,pMc, pEta, pTc, pLogDl, pRA, pSinDec, pPhase, pSinThJ0, pPhiJ0, pSpin1, pSpCosTh1, pSpPhi1);
  //printf("  %f\n",log(pDl));
  
  double x=0.0;
  double Mc=0.0,m1=0.0,m2=0.0,Mtot=0.0,mu=0.0;
  double cvec1[3],cvec2[3],cvec3[3];
  double tvec1[3],tvec2[3],tvec4[3],tvec6[3],tvec7[3];
  int i=0,terminate=0;
  double n_L[3];
  double localtc=0.0,altitude=0.0,azimuth=0.0,spin=0.0,samplerate=0.0,inversesamplerate=0.0;
  int length=0;
  for(i=0;i<3;i++) {
    cvec1[i] = 0.0;
    cvec2[i] = 0.0;
    cvec3[i] = 0.0;
    tvec1[i] = 0.0;
    tvec2[i] = 0.0;
    tvec4[i] = 0.0;
    tvec6[i] = 0.0;
    tvec7[i] = 0.0;
    n_L[i] = 0.0;
  }
  
  //for(i=0;i<100;i++) {
  // printf("  %5i  %5i  %5i  %5i  %s\n",i,run.parDef[i],run.mcmcParUse[i],run.injParUse[i],run.parAbrev[i]);
  //}
  
  localtc    = par->loctc[ifonr];
  altitude   = par->localti[ifonr];
  azimuth    = par->locazi[ifonr];
  samplerate = (double)ifo[ifonr]->samplerate;
  inversesamplerate = 1.0/samplerate;
  length     = ifo[ifonr]->samplesize;
  
  
  double n_z[3] = {0.0,0.0,1.0};                                                                                         // North in global coordinates
  double normalvec[3];                                                                                                  
  for(i=0;i<3;i++) normalvec[i] = ifo[ifonr]->normalvec[i];                                                              // Detector position normal vector = local zenith vector z'
  double D_L = pDl*Mpcs;                                                                                                 // Source luminosity distance, in seconds
  //double D_L = exp(pLogDl)*Mpcs;                                                                                                 // Source luminosity distance, in seconds
  double coslati = sqrt(1.0-pSinDec*pSinDec);
  double n_N[3] = { cos(pLongi)*coslati , sin(pLongi)*coslati , pSinDec };                                               // n_N: Position unit vector = N^
  
  double sthJ0   = pSinThJ0;                                                                                             // n_J0: 'total' AM unit vector, J0^  (almost equal to the real J, see Eq.15)
  double cthJ0   = sqrt(1.0 - sthJ0*sthJ0);
  double n_J0[3] = { cos(pPhiJ0)*cthJ0 , sin(pPhiJ0)*cthJ0 , sthJ0 };                                                    // Here, theta_Jo is a latitude-like angle like Dec (-pi/2-pi/2).
  
  par->NdJ = dotProduct(n_N,n_J0);                                                                                       // Inclination of J_0; only for printing purposes, should be removed from this routine
  //printf("  N.J: %lf,  acos(N.J):  %lf\n",par->NdJ,acos(par->NdJ));
  
  //Get individual masses from Mch and eta:
  Mc = pMc*M0;                                                                                                           // Chirp mass in seconds
  McEta2masses(Mc, pEta, &m1, &m2);                                                                                      //Mc,eta->M1,M2; accepts 0.25<eta<0.50
  Mtot = m1+m2;
  if(pEta>0.25) pEta = 0.5 - pEta;
  mu = m1*m2/Mtot;                                                                                                       // Eq.16b
  spin = pSpin1*m1*m1;
  
  
  //if(beVerbose>=1) {
  //printf("Ms: eta: %g  Mc: %g  m1: %g  m2: %g  Mtot: %g  mu: %g  Mo: %g\n",pEta,Mc/M0,m1/M0,m2/M0,Mtot/M0,mu/M0,M0);
  //}
  //printf("  Apo, local parameters:  %d  %lf  %lf  %lf  %lf  %d\n",ifonr,localtc,altitude,azimuth,samplerate,length);
  
  double beta = 1.0/12.0*(113.0*(m1*m1)/(Mtot*Mtot) + 75.0*pEta)*pSpCosTh1*spin/(m1*m1);                                 // Eq.20, for S2=0 or m1=m2,S1=S2:  kappa*spin/(m1*m1) = L^.S/m1^2, see Blanchet et al., PRL 74, 3515, 1995
  
  double cst1 = 743.0/336.0 + 11.0/4.0*pEta;
  double cst2 = (4.0*pi-beta);
  double cst5 = spin*sqrt(1.0-pSpCosTh1*pSpCosTh1);
  
  //Constant vector 1 for the construction of Eq.41e
  facVec(n_J0,-sthJ0,tvec1);                                                                                             //tvec1 = -J0^*cos(theta_J0)   MvdS: theta_J0 is a latitude, not a co-latitude
  addVec(n_z,tvec1,tvec2);                                                                                               //tvec2 = n_z - J0^*cos(theta_J0)
  facVec(tvec2,1.0/cthJ0,cvec1);                                                                                         //cvec1 = (n_z - J0^*cos(theta_J0))/sin(theta_J0)
  
  //Constant vector 2 for the construction of Eq.41e
  crossProduct(n_J0,n_z,tvec1);                                                                                          //tvec1 = J0^ x z^
  facVec(tvec1,1.0/cthJ0,cvec2);                                                                                         //cvec2 = (J0^ x z^) / sin(theta_J0)
  
  //Constant vector 3 for the construction of Eq.12 (local polarisation for F+,x)
  facVec(n_N,-dotProduct(normalvec,n_N),tvec1);                                                                          //tvec1 = -N^(z^'.N^)
  addVec(normalvec,tvec1,cvec3);                                                                                         //cvec3 = z^' - N^(z^'.N^)
  
  //Construct Eq.8ab, needed for F+,Fx
  double cosalti   = cos(altitude);
  double sin2azi   = sin(2.0*azimuth);
  double cos2azi   = cos(2.0*azimuth);
  double cst6  = 0.5*(1.0+cosalti*cosalti)*cos2azi;
  double cst7  = cosalti*sin2azi;
  
  
  double omega_low  = pi*ifo[ifonr]->lowCut;   //30 or 40 Hz, translated from f_gw to omega_orb
  double omega_high = pi*ifo[ifonr]->highCut;  //1600 Hz, translated from f_gw to omega_orb
  //double omega_high = min(pi*ifo[ifonr]->highCut, exp(-1.5*log(cutoff_a) - log(Mtot)) );  //1600 Hz, translated from f_gw to omega_orb, or a/Mtot = cutoff_a, whichever is smaller
  
  /*
    if(beVerbose>=1) {
    double Tcoal = 5.0*pow(8.0*omega_low,-8.0*c3rd)*pow(Mc,-5.0*c3rd) * (1.0 + 4.0*c3rd*cst1*pow(omega_low*Mtot,2.0*c3rd) - 1.6*cst2*(omega_low*Mtot));   //Time between f_low and coalescence
    double t0 = localtc - Tcoal;
    double deltat = (double)length*inversesamplerate;
    printf("Times:  Tcoal: %g,  t0: %g,  localtc: %g,  length: %d,  dt: %g,  dt-ltc: %g\n",Tcoal,t0,localtc,length,deltat,deltat-localtc);
    }
  */
  
  double oldomega = -1.e30;
  double phi_gw=0.0,alpha=0.0;
  
  //To display the number of wave and precession cycles:
  double phi1 = 0.0;
  //double phi2 = 0.0;
  double alpha1 = 0.0;
  //double alpha2 = 0.0;
  int i1=0,i2=0;
  
  double t=0.0,tau=0.0,tau18=0.0,tau28=0.0,tau38=0.0,tau58=0.0,tau_18=0.0,tau_28=0.0,tau_38=0.0,tau_58=0.0,tau_68=0.0;
  double omega_orb=0.0,l_L=0.0,Y=0.0,Gsq=0.0,G=0.0,slamL=0.0,clamL=0.0,LdotN=0.0;
  double hplus=0.0,hcross=0.0,locpolar=0.0,sin2polar=0.0,cos2polar=0.0,Fplus=0.0,Fcross=0.0;
  double cst4=0.0,x1=0.0,x2=0.0,x3=0.0;
  double taperx[length],omegas[length];
  for(i=0;i<length;i++) {
    taperx[i] = 0.0;
    omegas[i] = 0.0;
  }
  
  // Fill ifo[ifonr]->FTin with time-domain template:
  for(i=0; i<length; ++i){
    // Determine time left until coalescence, "(t_c-t)" in (4.17)/(11):
    t = localtc - ((double)i)*inversesamplerate;  // (time to t_c) = "(t_c-t)" in (4.17)
    if(t<0.0) { 
      if(terminate==0) terminate = 1;  //Set to 1 only if it was 0
    } else {
      tau    = pEta/(5.0*Mtot)*t;   //t = localtc-t already
      tau18  = exp(0.125*log(tau));  //tau^(1/8)
      tau28  = tau18*tau18;
      tau38  = tau28*tau18;
      tau58  = tau28*tau38;
      tau_18 = 1.0/tau18;            //tau^(-1/8)
      tau_28 = tau_18*tau_18;
      tau_38 = tau_18*tau_28;
      tau_58 = tau_28*tau_38;
      tau_68 = tau_38*tau_38;
      
      omega_orb = 1.0/(8.0*Mtot) * (tau_38 + 0.125*cst1*tau_58 - 0.075*cst2*tau_68);   // Orbital frequency
      omegas[i] = omega_orb;
      //printf("  t: %lf  tau: %lf  omega_orb: %lf\n",t,tau,omega_orb);
    }
    
    if((omega_orb>=omega_low) && (terminate==0)) {  // After source comes into window, before tc and before frequency reaches its maximum  Careful, if t<0, omega = nan!!!
      
      if(omega_orb < oldomega || omega_orb >= omega_high){  // Frequency starts decreasing, or frequency higher than highCut --> terminate signal
	//if(omega_orb < oldomega || omega_orb >= omega_high || taperx[i]>0.09){  // Frequency starts decreasing, or frequency higher than highCut, or v_orb>0.3c --> terminate signal
        ifo[ifonr]->FTin[i] = 0.0; 
        if(omega_orb < oldomega) terminate = 2;
        if(omega_orb >= omega_high) terminate = 3;
        //if(taperx[i]>0.09) terminate = 4;
	
      } else {             // Frequency still increasing --> keep on computing...
	
        if(i1==0) i1=i;  //Save initial i for tapering the beginning of the signal
        i2 = i;          //Save final i for tapering the end of the signal
        oldomega = omega_orb;
        taperx[i] = exp(2.0*c3rd*log(Mtot*omega_orb));                                                                      // x := (Mtot*w)^(2/3)  =  v_orb^2
	
        //Compute orbital A.M.
        l_L = m1*m2*exp(-c3rd*log(omega_orb*Mtot));
	
        //GW and orbital phase
        phi_gw = pPhase - 2.0/pEta * (tau58 + 0.625*c3rd*cst1*tau38 - 0.1875*cst2*tau28);                                   // GW phase at coalescence
        if(fabs(phi1)<1.e-30) phi1 = phi_gw;   //Save initial phi
        //phi2 = phi_gw;                       //Save final phi
	
        Y = spin/l_L;                                                                                                    //Y = |S|/|L|, Eq.43
        Gsq = 1.0 + 2.0*pSpCosTh1*Y + Y*Y;                                                                                      //G^2, Eq.46
        G   = sqrt(Gsq);
	
        cst4 = l_L + pSpCosTh1*spin;
        x = mu*Mtot;
        x1 = x*x*x;
        x = G*l_L;
        x2 = x*x*x;
        x3 = spin*spin*spin;
        alpha = pSpPhi1 - 5.0/(96.0*x1) * (1.0+0.75*m2/m1) * 
          (2.0*x2 - 3.0*pSpCosTh1*spin*cst4*G*l_L - 3.0*pSpCosTh1*x3*(1.0-pSpCosTh1*pSpCosTh1) * asinh(cst4/cst5));                                            //Eq.47
        if(fabs(alpha1)<1.e-30) alpha1 = alpha;  //Save initial alpha
        //alpha2 = alpha;                         //Save final alpha
	
        slamL = cst5/(l_L*G);                                                                                            //sin(lambda_L), Eq.48a
        clamL = cst4/(l_L*G);                                                                                            //cos(lambda_L), Eq.48b
	
	
        //Construct Eq.41e
        facVec(n_J0,clamL,tvec1);                                                                                        //tvec1 = J0^*cos(lambda_L)
        facVec(cvec1,slamL*cos(alpha),tvec4);                                                                            //tvec4 = (n_z - J0^*cos(theta_J0))*sin(lambda_L)*cos(alpha)/sin(theta_J0)
        facVec(cvec2,slamL*sin(alpha),tvec6);                                                                            //tvec6 = (J0^ x z^) * sin(lambds_L)*sin(alpha)/sin(theta_J0)
        addVec(tvec1,tvec4,tvec7);                                                                                       //Construct Eq.59
        addVec(tvec7,tvec6,n_L);                                                                                         //Eq.59: n_L=L^
	
	
        LdotN  = dotProduct(n_L,n_N);                                                                                    //L^.N^
        x1     = 2.0*exp(5.0*c3rd*log(Mc))/D_L;
        x3     = exp(2.0*c3rd*log(omega_orb));
        hplus  =      x1 * (1.0 + LdotN*LdotN) * x3 * cos(phi_gw);
        hcross = -2.0*x1 * LdotN               * x3 * sin(phi_gw);
	
	
        //Local polarisation vector, F+,Fx:
        crossProduct(n_L,normalvec,tvec1);                                                                              //tvec1 = n_L x z^'
        locpolar  = atan(dotProduct(n_L,cvec3)/dotProduct(n_N,tvec1));                                                  //Eq.12 of Vecchio, result between -pi/2 and pi/2
        sin2polar = sin(2.0*locpolar);
        cos2polar = sqrt(1.0-sin2polar*sin2polar);                                                                      //Since 2*locpolar should be between -pi and pi (?)
        Fplus     =  cst6*cos2polar + cst7*sin2polar;                                                                   //Eq.8a
        Fcross    = -cst6*sin2polar + cst7*cos2polar;                                                                   //Eq.8b
	
        //Detector signal:
        ifo[ifonr]->FTin[i] = Fplus*hplus + Fcross*hcross;                                                              //  (3.10)
	
	
	
        //Print some stuff for diagnostics:
        //printf("i: %8d   t: %10g   f: %10g   h: %10g\n",i,t,omega_orb/pi,ifo[ifonr]->FTin[i]);
        //if((omega_orb/pi<40.002 || fabs(t)<0.2) && beVerbose>=1) {
        //if(beVerbose>=1) {
        //printf("i: %8d   t: %10g   f: %10g   x: %10g\n",i,t,omega_orb/pi,taperx[i]);
        //printf("omg_orb: %10g  phi_orb: %10g  l_L: %10g  S: %10g  k: %10g  Y: %10g  G: %10g  alpha: %10g  slamL: %10g  clamL: %10g \n",  omega_orb,phi_orb,l_L,spin,pSpCosTh1,Y,G,alpha,slamL,clamL);
        //printf("i: %8d   t: %10g  alpha_c: %10g   alpha0: %10g   alpha: %10g  alpha/2pi: %10g\n",  i,t,pSpPhi1,alpha0,alpha,alpha/tpi);
        //}
	
	
      }
    }  //end if((omega_orb>=omega_low) && (terminate==0)) {  // After source comes into window, before tc and before frequency reaches its maximum  Careful, if t<0, omega = nan!!!
    else {
      ifo[ifonr]->FTin[i]   = 0.0;  //  (If before omega_low, after t_c or after termination)
    }
  }  //i
  
  
  
  //Print some stuff for diagnostics
  //if(i1<=1 && pMc>0.02)  fprintf(stderr, "\n ***  Warning: length is too small to fit waveform template, increase before_tc ***\n\n"); //Don't print when doing null-likelihood
  //if(i2>=length-1 && pMc>0.02)  fprintf(stderr, "\n ***  Warning: length is too small to fit waveform template, increase after_tc ***\n\n"); //Don't print when doing null-likelihood
  //if(beVerbose>=1) 
  //printf("%10.2f  %10.2f  %10.2f  %10.1f  %10.2f  %10.2f",pSpin1,acos(pSpCosTh1)*r2d,(double)(i2-i1)*inversesamplerate,(phi2-phi1)/tpi,(alpha2-alpha1)/tpi,pow(Mtot*oldomega,-2.0*c3rd));
  //if(beVerbose>=1) 
  //printf("  term: %d  i1: %d  i2: %d  i2-i1: %d  length: %d  f_gw,old:%9.3lf  f_gw:%9.3lf  f_gw,low:%9.3lf  f_gw,high:%9.3lf  f_gw1:%9.3lf  f_gw2:%9.3lf\n",terminate,i1,i2,i2-i1,length,oldomega/pi,omega_orb/pi,omega_low/pi,omega_high/pi,omegas[i1]/pi,omegas[i2]/pi);
  //Terminate: 1: t>tc, 2: df/dt<0, 3: f>f_high
  
  
  
  //Apply tapering
  //int i1a = i1 + ceil(2.0*samplerate*pi/omegas[i1]);  //pi/omega_orb = 1/f_gw = P_gw.  *samplerate: number of points in first GW cycle
  int i2a = i2 - (int)ceil(2.0*samplerate*pi/omegas[i2]);  //pi/omega_orb = 1/f_gw = P_gw.  *samplerate: number of points in last GW cycle
  for(i=i1;i<=i2;i++) {
    //ifo[ifonr]->FTin[i] *= 0.5*(1.0 - tanh(15000.0*(taperx[i1a]-taperx[i])));  //Taper beginning of template
    ifo[ifonr]->FTin[i] *= 0.5*(1.0 - tanh(100.0*(taperx[i]-taperx[i2a])));  //Taper end of template
  }
  
} // End of templateApostolatos()
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Calculate the local parameters from the global parameters
 * 
 * Calculate the local (i.e. in the detector frame) parameters from the global parameters.
 *    par   :  pointer to parameter set (struct)
 *    ifo   :  pointer to interferometer data (struct)
 */
// ****************************************************************************************************************************************************  
void localPar(struct parSet *par, struct interferometer *ifo[], int networkSize, int injectionWF, struct runPar run)
{
  int ifonr=0,j=0;
  double lineofsight[3], dummyvec[3], scalprod1=0.0, delay=0.0;
  
  double pTc=0.0,pLongi=0.0,pSinDec=0.0;
  if(injectionWF==1) {                                                               // Then this is for an injection waveform template
    pTc       = par->par[run.injRevID[11]];                                            // 11: Tc
    pLongi    = fmod(longitude(par->par[run.injRevID[31]], GMST(pTc)) + mtpi, tpi);    // 31: RA;  RA -> 'longitude'
    pSinDec   = par->par[run.injRevID[32]];                                            // 32: sin(Dec)
  } else {                                                                           // Then this is for an MCMC waveform template
    pTc       = par->par[run.parRevID[11]];                                            // 11: Tc
    pLongi    = fmod(longitude(par->par[run.parRevID[31]], GMST(pTc)) + mtpi, tpi);    // 31: RA;  RA -> 'longitude'
    pSinDec   = par->par[run.parRevID[32]];                                            // 32: sin(Dec)
  }
  
  
  // Determine local coalescence times:
  coord2vec(pSinDec, pLongi, lineofsight);
  for(ifonr=0; ifonr<networkSize; ifonr++) {
    scalprod1 =  ifo[ifonr]->positionvec[0]*lineofsight[0]  +  ifo[ifonr]->positionvec[1]*lineofsight[1]  +  ifo[ifonr]->positionvec[2]*lineofsight[2];  // Project line of sight onto positionvec, scalprod1 is in units of metres
    delay = scalprod1 / c;                                         // Time delay (wrt geocentre) in seconds
    par->loctc[ifonr] = ((pTc - ifo[ifonr]->FTstart) - delay);
  }
  
  
  // Determine local sky position:
  for(ifonr=0; ifonr<networkSize; ifonr++) {
    
    // 'Altitude' in the IFO' frame:
    par->localti[ifonr] = angle(ifo[ifonr]->normalvec, lineofsight);                     // Actually, this is the colatitude in the IFO' frame (i.e. 0deg=zenith, 90deg=horizon)
    
    // 'Azimuth' in the IFO' frame:
    for(j=0; j<3; ++j) dummyvec[j] = lineofsight[j];                                     // Temp vector with line of sight
    orthoProject(dummyvec, ifo[ifonr]->rightvec, ifo[ifonr]->orthoArm);                  // Project line of sight into IFO' arm plane
    par->locazi[ifonr] = angle(dummyvec, ifo[ifonr]->rightvec);                          // The 'true' azimuth (N=0,E=90deg) of the source at the location of the detector is:  pi - (par->locazi[ifonr] + ifo[ifonr]->rightArm) 
    if(!rightHanded(ifo[ifonr]->rightvec, dummyvec, ifo[ifonr]->normalvec)) par->locazi[ifonr] = 2.0*pi - par->locazi[ifonr];
    
    //printf("  localPar:  %d  %lf  %lf  %s\n",ifonr,ifo[ifonr]->lati/pi*180.0,ifo[ifonr]->longi/pi*180.0,ifo[ifonr]->name);
  }
  
} // End of localPar()
// ****************************************************************************************************************************************************  


