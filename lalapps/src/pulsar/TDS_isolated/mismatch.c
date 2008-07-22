#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/Random.h>
#include <lal/LALString.h>
#include <lal/SFTutils.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/MatrixUtils.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/LALRCSID.h>

/* #define TIMESPAN 31557600 */
#define TIMESPAN 3155760
#define SIXTH 0.166666666666666666666666666666666666666667L
#define TWENTYFOURTH 0.04166666666666666666666666666666666666666667L
#define EARTH "/Users/matthew/lscsoft/lal/packages/pulsar/test/earth05-09.dat"
#define SUN "/Users/matthew/lscsoft/lal/packages/pulsar/test/sun05-09.dat"

RCSID("$Id$");

REAL8Vector *get_phi( double start, double deltaT, int npoints,
 BinaryPulsarParams params, BarycenterInput bary, EphemerisData *edat );

int main(int argc, char *argv[]){
  static LALStatus status;

  FILE *fp=NULL;

  BinaryPulsarParams params;

  BinaryPulsarParams deltaparams;

  char *parfile=NULL;

  int doParam[20];
  double errVals[20];
  int i=0, count=0, countBin=0, j=1, k=0, N=0;

  int deltaT=10*60;            /* seconds between points */
  int npoints=TIMESPAN/deltaT; /* number of 10 minutes stretches */

  EphemerisData *edat=NULL;
  LALDetector det;
  BarycenterInput baryinput;

  REAL8Vector *phiMean=NULL;
  REAL8 intre=0., intim=0., integral1=0.;
  REAL8 phaseSum1 = 0., phaseSum2 = 0.;

  REAL8 maxmismatch = 1e-100;

  parfile = argv[1];

  XLALReadTEMPOParFile( &params, parfile );
  XLALReadTEMPOParFile( &deltaparams, parfile );

  /* set up ephemerises */
  det = *XLALGetSiteInfo( "H1" ); /* just set site as LHO */
  baryinput.site = det;
  baryinput.site.location[0] /= LAL_C_SI;
  baryinput.site.location[1] /= LAL_C_SI;
  baryinput.site.location[2] /= LAL_C_SI;

  edat = XLALMalloc(sizeof(*edat));

  (*edat).ephiles.earthEphemeris = EARTH;
  (*edat).ephiles.sunEphemeris = SUN;
  LALInitBarycenter(&status, edat);

 /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  /* set the position and frequency epochs if not already set */
  if(deltaparams.pepoch == 0. && deltaparams.posepoch != 0.)
    deltaparams.pepoch = deltaparams.posepoch;
  else if(deltaparams.posepoch == 0. && deltaparams.pepoch != 0.)
    deltaparams.posepoch = deltaparams.pepoch;

  if(deltaparams.pepoch == 0. && deltaparams.posepoch != 0.)
    deltaparams.pepoch = deltaparams.posepoch;
  else if(deltaparams.posepoch == 0. && deltaparams.pepoch != 0.)
    deltaparams.posepoch = deltaparams.pepoch;

  /* calculate the phase every minute for the mean values */
  phiMean = get_phi( 800000000.0, (double)deltaT, npoints, params, baryinput,
    edat );

  /* calculate phase integral over time */
  for( k=0; k<npoints-1; k++ ){

    intre += 2.;
    intim += 0.;

    //phaseSum1 += (phiMean->data[k] - phiMean->data[0]);
  }
  //phaseSum1 += (phiMean->data[k+1] - phiMean->data[0]);

  integral1 = intre*intre + intim*intim;

  /* see which parameters have an associated error */
  if( params.f0Err != 0. ) count++;
  if( params.f1Err != 0. ) count++;
  if( params.f2Err != 0. ) count++;
  if( params.raErr != 0. ) count++;
  if( params.decErr != 0. ) count++;
  if( params.pmraErr != 0. ) count++;
  if( params.pmdecErr != 0. ) count++;

  /* see if binary model */
  if( params.model != NULL ){
    if( params.xErr != 0. ) countBin++;
    if( params.w0Err != 0. ) countBin++;
    if( params.eErr != 0. ) countBin++;
    if( params.PbErr != 0. ) countBin++;
    if( params.T0Err != 0. ) countBin++;
    if( params.TascErr != 0. ) countBin++;
    if( params.eps1Err != 0. ) countBin++;
    if( params.eps2Err != 0. ) countBin++;
    if( params.gammaErr != 0. ) countBin++;
    if( params.wdotErr != 0. ) countBin++;
    if( params.PbdotErr != 0. ) countBin++;
    if( params.xdotErr != 0. ) countBin++;
    if( params.edotErr != 0. ) countBin++;
  }

  fprintf(stderr, "%d frequency and position errors, %d binary system errors\n",
    count, countBin);

  N = (int)pow(2, count+countBin);
  fprintf(stderr, "grid size = %d\n", N);
  /* calculate the max mismatch between the mean value and the mean plus the
     standard deviation */
  for( i=0; i < N; i++ ){
    REAL8Vector *phiOffset=NULL;
    REAL8 integral2=0.;

    j=1;

    //phaseSum2 = 0.;

    if(count > 0){
      if( params.f0Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over f0\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.f0 = params.f0 + params.f0Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.f0 = params.f0 - params.f0Err;
        j++;
      }
      if( params.f1Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over f1\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.f1 = params.f1 + params.f1Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.f1 = params.f1 - params.f1Err;
        j++;
      }
      if( params.f2Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over f2\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.f2 = params.f2 + params.f2Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.f2 = params.f2 - params.f2Err;
        j++;
      }
      if( params.raErr != 0. ){
        if ( i==0 )fprintf(stderr, "Search over ra\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.ra = params.ra + params.raErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.ra = params.ra - params.raErr;
        j++;
      }
      if( params.decErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over dec\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.dec = params.dec + params.decErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.dec = params.dec - params.decErr;
        j++;
      }
      if( params.pmraErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over pmra\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.pmra = deltaparams.pmraErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.pmra = -deltaparams.pmraErr;
        j++;
      }
      if( params.pmdecErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over pmdec\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.pmdec = deltaparams.pmdecErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.pmdec = params.pmdec - params.pmdecErr;
        j++;
      }
    }

    if( countBin > 0 ){
      if( params.xErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over x\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.x = params.x + params.xErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.x = params.x - params.xErr;
        j++;
      }
      if( params.w0Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over w0\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.w0 = params.w0 + params.w0Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.w0 = params.w0 - params.w0Err;
        j++;
      }
      if( params.eErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over e\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.e = params.e + params.eErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.e = params.e - params.eErr;
        j++;
      }
      if( params.PbErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over Pb\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.Pb = params.Pb + params.PbErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.Pb = params.Pb - params.PbErr;
        j++;
      }
      if( params.T0 != 0. ){
        if( i==0 ) fprintf(stderr, "Search over T0\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.T0 = params.T0 + params.T0Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.T0 = params.T0 - params.T0Err;
        j++;
      }
      if( params.Tasc != 0. ){
        if( i==0 ) fprintf(stderr, "Search over Tasc\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.Tasc = params.Tasc + params.TascErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.Tasc = params.Tasc - params.TascErr;
        j++;
      }
      if( params.eps1Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over eps1\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.eps1 = params.eps1 + params.eps1Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.eps1 = params.eps1 - params.eps1Err;
        j++;
      }
      if( params.eps2Err != 0. ){
        if( i==0 ) fprintf(stderr, "Search over eps2\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.eps2 = params.eps2 + params.eps2Err;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.eps2 = params.eps2 - params.eps2Err;
        j++;
      }
      if( params.gammaErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over gamma\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.gamma = params.gamma + params.gammaErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.gamma = params.gamma - params.gammaErr;
        j++;
      }
      if( params.wdotErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over wdot\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.wdot = params.wdot + params.wdotErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.wdot = params.wdot - params.wdotErr;
        j++;
      }
      if( params.PbdotErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over Pbdot\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.Pbdot = params.Pbdot + params.PbdotErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.Pbdot = params.Pbdot - params.PbdotErr;
        j++;
      }
      if( params.xdotErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over xdot\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.xdot = params.xdot + params.xdotErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.xdot = params.xdot - params.xdotErr;
        j++;
      }
      if( params.edotErr != 0. ){
        if( i==0 ) fprintf(stderr, "Search over edot\n");
        if( i%(int)(N/pow(2,j-1)) < N/pow(2,j) )
          deltaparams.edot = params.edot + params.edotErr;
        else if( i%(int)(N/pow(2,j-1)) >= N/pow(2,j) )
          deltaparams.edot = params.edot - params.edotErr;
        j++;
      }
    }

    /* get new phase */
    phiOffset = get_phi( 800000000.0, (double)deltaT, npoints, deltaparams,
      baryinput, edat );

    /* calculate the mismatch 1 - (P(params + delta))/P(params) */
    intre = 0.;
    intim = 0.;

    if(i == 0)
      fp = fopen("model.txt", "w");

    for( k=0; k<npoints-1; k++ ){
      REAL8 phi1=0., phi2=0.;
      phi1 = LAL_TWOPI*fmod(phiOffset->data[k]-phiMean->data[k], 1.);
      phi2 = LAL_TWOPI*fmod(phiOffset->data[k+1]-phiMean->data[k+1], 1.);
      intre += (cos(phi2) + cos(phi1));
      intim += (sin(phi2) + sin(phi1));

      if( i == 0 )
        fprintf(fp, "%lf\n", phi1);
    }

    if( i == 0 )
      fclose(fp);

    integral2 = intre*intre + intim*intim;  /* square value for power */

    //fprintf(stderr, "phase error = %lf\n", phaseSum2 - phaseSum1 );
    fprintf(stderr, "mismatch = %le\n", 1. - integral2/integral1);

    /* work out mismatch */
    if( fabs(1. - integral2/integral1) > fabs(maxmismatch) )
      maxmismatch = 1. - integral2/integral1;

    XLALDestroyREAL8Vector( phiOffset );
  }

  fprintf(stderr, "Maximum mismatch = %lf\n", maxmismatch);

  XLALDestroyREAL8Vector( phiMean );

  return 0;
}

/* function to return a vector of the pulsar phase for each data point */
REAL8Vector *get_phi( double start, double deltaT, int npoints,
 BinaryPulsarParams params, BarycenterInput bary, EphemerisData *edat ){
  static LALStatus status;

  INT4 i=0;

  REAL8 T0=0., DT=0., DTplus=0., deltat=0., deltat2=0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */
  REAL8 time=0.;

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt=0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  REAL8Vector *phis=NULL;

  FILE *fp=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( edat == NULL)
    return NULL;

  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( npoints );

  for( i=0; i<npoints; i++){
    T0 = params.pepoch;
    time = start + deltaT*(double)i;

    DT = time - T0;

    if(time <= 820108813)
      (*edat).leap = 13;
    else
      (*edat).leap = 14;

    /* only call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary.tgps.gpsSeconds = (UINT8)floor(time);
      bary.tgps.gpsNanoSeconds = (UINT8)floor((fmod(time,1.)*1e9));

      bary.delta = params.dec + DT*params.pmdec;
      bary.alpha = params.ra + DT*params.pmra/cos(bary.delta);

      /* call barycentring routines */
      LALBarycenterEarth(&status, &earth, &bary.tgps, edat);
      LALBarycenter(&status, &emit, &bary, &earth);

      /* add 30 minutes (1800secs) to the time */
      DTplus = DT + interptime;
      bary.tgps.gpsSeconds = (UINT8)floor(time+interptime);
      bary.tgps.gpsNanoSeconds = (UINT8)floor((fmod(time+interptime,1.)*1e9));

      /* No point in updating the positions as difference will be tiny */
      LALBarycenterEarth(&status, &earth2, &bary.tgps, edat);
      LALBarycenter(&status, &emit2, &bary, &earth2);
    }

    /* linearly interpolate to get emitdt */
    emitdt = emit.deltaT + (DT - (DTplus - interptime)) *
      (emit2.deltaT - emit.deltaT)/interptime;

    /* check if need to perform binary barycentring */
    if( params.model != NULL ){
      binput.tb = time + emitdt;

      XLALBinaryPulsarDeltaT( &boutput, &binput, &params );

      deltat = DT + emitdt + boutput.deltaT;
    }
    else
      deltat = DT + emitdt;

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = 2.*(params.f0*deltat + 0.5*params.f1*deltat2
      + SIXTH*params.f2*deltat*deltat2 +
        TWENTYFOURTH*params.f3*deltat2*deltat2);
  }

  return phis;
}
