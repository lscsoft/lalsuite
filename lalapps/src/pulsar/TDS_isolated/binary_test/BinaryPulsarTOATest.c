/* Matt Pitkin 2/05/05
   Code to test my binary code against actual pulsar observations */

/* $Id$ */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALConstants.h>

#include <lal/BinaryPulsarTiming.h>

int main(int argc, char *argv[]){
  static LALStatus status;

  char *pulsarAndPath=NULL;
  char *pulsarTOAFile=NULL;
  
  FILE *fpin=NULL;
  FILE *fpout1=NULL, *fpout2=NULL, *fpout3=NULL;
  
  char *telescope=NULL;
  char *psr=NULL;
  double radioFreq=0.0, rf[10000];
  double TOA[10000];
  double num1, num3;
  int num2;
  char *phaseStr=NULL;
  int i=0, j=0, k=0, length;
  
  double PPTimeEll1[10000], PPTimeBT[10000], PPTimeDD[10000]; /* Pulsar proper
time - corrected for solar system and
 binary orbit delay times */
  const double D = 4.148808e3; /* dispersion constant MHz^2 pc^-1 cm^3 s */
  
  /* Binary pulsar variables */
  BinaryPulsarParams params;
  BinaryPulsarInput input;
  BinaryPulsarOutput output;

  /* LAL barycentring variables */
  BarycenterInput baryinput;
  EarthState earth;
  EmissionTime  emit;
  EphemerisData *edat=NULL;
  char earthFile[256], sunFile[256];
  
  pulsarAndPath = argv[1];
  pulsarTOAFile = argv[2];
  
  double MJD_tcorr[10000];
  double tcorr[10000];
  char date[15];
  
  double f0, f1, f2, T;

  fpin = fopen(pulsarTOAFile, "r");
  
  /* read in TOA and phase info */
  while(!feof(fpin)){
    fscanf(fpin, "%s%s%lf%lf%f%d%s%f", &telescope, &psr, &radioFreq, &TOA[i],
&num1, &num2, &phaseStr, &num3);
    rf[i] = radioFreq;
    
    i++;
  }
  
  fclose(fpin);
    
  /* read in telescope time corrections from file */
  fpin = fopen("time_bonn.txt", "r");
  
  while(!feof(fpin)){
    fscanf(fpin, "%lf%d%lf%s%s", &MJD_tcorr[j], &num2, &tcorr[j], &telescope,
&date);
    j++;
  }
  length = j;
  
  fclose(fpin);
  
  /* read in binary params from par file */
  LALReadTEMPOParFile(&status, &params, pulsarAndPath);
  
  //params.Tasc -= 51.184;
  
  if(params.Tasc != 0.0){
    params.T0 = params.Tasc;
    params.w0 = 0.0;
    params.e = sqrt(params.eps1*params.eps1 + params.eps2*params.eps2);
  }
  
  /* show pulsar params */
  fprintf(stderr, "f0 = %lf Hz, f1 = %e Hz/s, f2 = %e Hz/s^2.\n", params.f0,
params.f1, params.f2);
  fprintf(stderr, "RA = %lf rads, DEC = %lf rads.\n", params.ra, params.dec);
  fprintf(stderr, "Binary model is %s.\n", params.model);
  fprintf(stderr, "asini = %lf light seconds, period = %lf days.\n", params.x,
params.Pb);
  fprintf(stderr, "eps1 = %.8lf, eps2 = %.8lf.\n", params.eps1, params.eps2);
  fprintf(stderr, "Tasc = %lf.\n", params.Tasc); 
  
  /* set telescope location - for this case it the Effelsberg telescope and the
x,y,z components are got from tempo */
  baryinput.site.location[0] = 4033949.5/LAL_C_SI;
  baryinput.site.location[1] =  486989.4/LAL_C_SI;
  baryinput.site.location[2] = 4900430.8/LAL_C_SI;
  
  /* initialise the solar system ephemerides */
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));
  sprintf(earthFile,
"/home/matthew/lscsoft/lal/packages/pulsar/test/earth00-04.dat");
  sprintf(sunFile,
"/home/matthew/lscsoft/lal/packages/pulsar/test/sun00-04.dat");
  (*edat).ephiles.earthEphemeris = earthFile;
  (*edat).ephiles.sunEphemeris = sunFile;
  LALInitBarycenter(&status, edat);
  
  /* convert TOAs from MJD into GPS */
  
  fpout1 = fopen("deltaT.txt", "w");
  
  for(j=0;j<i-1;j++){
    double t, DM = 9.023345; /* DM for current pulsar - make more general */
    double deltaD_f2;
    REAL8 fe, uasc, Dt; /* see TEMPO tasc2t0.f */

    while(MJD_tcorr[k] < TOA[j])
      k++;

    /* there's a 2 nanosecond UTC - UTC(NIST) correction (for MJD 51599), but I
can't see that making any difference */
    t = (TOA[j]-44244.0)*86400.0 + 13. + tcorr[k]*1.e-6;// (double)input.leapSecs + tcorr[k]*1.e-6;
    
    deltaD_f2 = D * DM/(rf[j]*rf[j]);
    t -= deltaD_f2; /* dedisperse times */
    
    /* set pulsar position */
    baryinput.delta = params.dec + params.pmdec*(t-params.pepoch);
    baryinput.alpha = params.ra + params.pmra*(t-params.pepoch)/cos(baryinput.delta);
    baryinput.dInv = 0.0;  /* no parallax */

    baryinput.tgps.gpsSeconds = (INT4)floor(t);
    baryinput.tgps.gpsNanoSeconds = (INT4)(fmod(t,1.0)*1.e9);
    /* fprintf(stderr, "%.9f = %d.%09d\n", t, baryinput.tgps.gpsSeconds,
baryinput.tgps.gpsNanoSeconds); */
    
    /* calculate solar system barycentre time delay */
    LALBarycenterEarth(&status, &earth, &baryinput.tgps, edat);
    LALBarycenter(&status, &emit, &baryinput, &earth);
    
    /* calculate binary barycentre time delay */
    input.tb = t + (double)emit.deltaT;

    params.model = "ELL1";
    LALBinaryPulsarDeltaT(&status, &output, &input, &params);

    PPTimeEll1[j] = t + ((double)emit.deltaT + output.deltaT);
    
    params.model = "BT";
    LALBinaryPulsarDeltaT(&status, &output, &input, &params);
    
    PPTimeBT[j] = t + ((double)emit.deltaT + output.deltaT);
    
    params.model = "DD";
    LALBinaryPulsarDeltaT(&status, &output, &input, &params);
    
    PPTimeDD[j] = t + ((double)emit.deltaT + output.deltaT);

    fprintf(fpout1, "%.9f\t%lf\t%lf\n", PPTimeEll1[j], output.deltaT,
emit.deltaT);
  }
  
  fclose(fpout1);
  
  fpout1 = fopen("binaryPhaseEll1.txt", "w");
  fpout2 = fopen("binaryPhaseBT.txt", "w");
  fpout3 = fopen("binaryPhaseDD.txt", "w");
  
  /* f0 and f1 at start of observation */
  T = PPTimeEll1[0] - params.pepoch;
  f0 = params.f0 + params.f1*T + 0.5*params.f2*T*T;
  f1 = params.f1 + params.f2*T;
  f2 = params.f2;
  
  /* work out phase of each consecutive time */
  for(j=1;j<i-1;j++){
    double phase1, phase2, phase3;
    double tt0Ell1, tt0BT, tt0DD;

    tt0Ell1 = PPTimeEll1[j] - PPTimeEll1[0];
    tt0BT = PPTimeBT[j] - PPTimeBT[0];
    tt0DD = PPTimeDD[j] - PPTimeDD[0];
  
    /* phase = params.f0*tt0 + 0.5*params.f1*tt0*tt0 +
params.f2*tt0*tt0*tt0/6.0; */
    phase1 = f0*tt0Ell1 + 0.5*f1*tt0Ell1*tt0Ell1 +
f2*tt0Ell1*tt0Ell1*tt0Ell1/6.0;
    phase1 = fmod(phase1+0.5, 1.0);
    
    phase2 = f0*tt0BT + 0.5*f1*tt0BT*tt0BT + f2*tt0BT*tt0BT*tt0BT/6.0;
    phase2 = fmod(phase2+0.5, 1.0);
    
    phase3 = f0*tt0DD + 0.5*f1*tt0DD*tt0DD + f2*tt0DD*tt0DD*tt0DD/6.0;
    phase3 = fmod(phase3+0.5, 1.0);
    
    fprintf(fpout1, "%.9f\t%f\n", tt0Ell1, phase1);
    fprintf(fpout2, "%.9f\t%f\n", tt0BT, phase2);
    fprintf(fpout3, "%.9f\t%f\n", tt0DD, phase3);
  }
  
  fclose(fpout1);
  fclose(fpout2);
  fclose(fpout3);

  return 0;
}
