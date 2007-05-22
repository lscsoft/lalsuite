/** \file InspiralInjectionParams.c
 *  \ingroup InspiralInjectionParams
 *  \author D. Brown, J. Creighton, S. Fairhurst, G. Jones, E. Messaritaki
 * 
 *  \brief Functions for generating random distributions of inspiral parameters
 *  for injection purposes
 *
 * $Id$ 
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/InspiralInjectionParams.h>

NRCSID (INSPIRALINJECTIONPARAMSC,"$Id$");

/** Generates the geocent_end_time for an inspiral injection, based on the
 * given startTime and timeWindow */
SimInspiralTable* XLALRandomInspiralTime(
    SimInspiralTable *inj,   /**< injection for which time will be set */
    RandomParams *randParams,/**< random parameter details*/
    LIGOTimeGPS startTime,   /**< the first time that the injection could be */
    REAL4 timeWindow         /**< the time window within which inj must occur*/
    )
{       
  REAL8 timeOffset;
  timeOffset = (REAL8) timeWindow * XLALUniformDeviate( randParams );
  inj->geocent_end_time = *XLALGPSAdd( &startTime, timeOffset );
  inj->end_time_gmst = XLALGreenwichMeanSiderealTime( 
      &(inj->geocent_end_time) );

  return ( inj );
}

/** Generates the distance for an inspiral injection, based on the requested
 * distribution and max/min distances */
SimInspiralTable* XLALRandomInspiralDistance( 
    SimInspiralTable *inj,     /**< injection for which distance will be set */
    RandomParams *randParams,  /**< random parameter details*/
    DistanceDistribution dDist,/**< requested distance distribution */
    REAL4  distMin,            /**< minimum distance (Mpc) */
    REAL4  distMax             /**< maximum distance (Mpc) */
    )    
{       
  if ( dDist == uniformDistance )
  {
    /* uniform distribution in distance */
    inj->distance = distMin + 
      (distMax - distMin) * XLALUniformDeviate( randParams );
  }
  else if ( dDist == uniformLogDistance )
  {
    /* uniform distribution in log(distance) */
    REAL4 lmin = log10(distMin);
    REAL4 lmax = log10(distMax);
    REAL4 deltaL = lmax - lmin;
    REAL4 exponent;
    exponent = lmin + deltaL * XLALUniformDeviate( randParams );
    inj->distance = pow(10.0,(REAL4) exponent);
  }
  else if (dDist == uniformVolume )
  {
    /* uniform volume distribution */
    REAL4 d2min = distMin * distMin ;
    REAL4 d2max = distMax * distMax ;
    REAL4 deltad2 = d2max - d2min ;
    REAL4 d2;
    d2 = d2min + deltad2 * XLALUniformDeviate( randParams );
    inj->distance = pow(d2, 1.0/3.0);     
  }    
  return ( inj );
} 


/** Generates a random sky location (theta, phi) for an inspiral injection */
SimInspiralTable* XLALRandomInspiralSkyLocation( 
    SimInspiralTable *inj,  /**< injection for which sky location will be set*/
    RandomParams *randParams/**< random parameter details*/
    )
{
  inj->latitude = asin( 2.0 * XLALUniformDeviate( randParams ) - 1.0 ) ;
  inj->longitude = LAL_TWOPI * XLALUniformDeviate( randParams ) ;
  return ( inj );
}


/** Generates a random orientation (polarization, inclination, coa_phase)
 * for an inspiral injection.  If inclinationPeak is non-zero, then peak the 
 * inclination around o with width inclinationPeak */
SimInspiralTable* XLALRandomInspiralOrientation( 
    SimInspiralTable *inj,   /**< injection for which orientation will be set*/
    RandomParams *randParams,/**< random parameter details*/
    REAL4   inclinationPeak /**< width of the peak of the inclination */
    )
{       
  if ( !inclinationPeak )
  {
    inj->inclination = acos( 2.0 * XLALUniformDeviate( randParams ) 
        - 1.0 );
  }
  else
  {
    REAL4Vector *normalDev;

    normalDev = XLALCreateREAL4Vector( 1 );
    XLALNormalDeviates(normalDev, randParams);
    inj->inclination = inclinationPeak * normalDev->data[0];
    XLALDestroyVector( normalDev );

  }

  inj->polarization = LAL_TWOPI * XLALUniformDeviate( randParams ) ;
  inj->coa_phase = LAL_TWOPI * XLALUniformDeviate( randParams ) ;

  return ( inj );
}

/** Generates random masses for an inspiral injection. */
SimInspiralTable* XLALRandomInspiralMasses( 
    SimInspiralTable *inj,   /**< injection for which masses will be set*/
    RandomParams *randParams,/**< random parameter details*/
    MassDistribution mDist,  /**< the mass distribution to use */
    REAL4  mass1Min,         /**< minimum mass for first component */
    REAL4  mass1Max,         /**< maximum mass for first component */
    REAL4  mass2Min,         /**< minimum mass for second component */
    REAL4  mass2Max,         /**< maximum mass for second component */
    REAL4  maxTotalMass      /**< maximum total mass of binary */
    )
{       
  REAL4 mTotal = maxTotalMass +1;

  while ( mTotal > maxTotalMass )
  {

    if ( mDist == uniformComponentMass )
    {
      /* uniformly distributed mass1 and uniformly distributed mass2 */
      inj->mass1 = mass1Min + XLALUniformDeviate( randParams ) * 
        (mass1Max - mass1Min);
      inj->mass2 = mass2Min + XLALUniformDeviate( randParams ) * 
        (mass2Max - mass2Min);
      mTotal = inj->mass1 + inj->mass2 ;
    }
    else if ( mDist == logComponentMass )
    {
      /* distributed logarithmically in mass1 and mass2 */
      inj->mass1 = exp( log(mass1Min) + XLALUniformDeviate( randParams ) *
          (log(mass1Max) - log(mass1Min) ) );
      inj->mass2 = exp( log(mass2Min) + XLALUniformDeviate( randParams ) *
          (log(mass2Max) - log(mass2Min) ) );
      mTotal = inj->mass1 + inj->mass2;
    }
    else if ( mDist == uniformTotalMass )
    {
      /*uniformly distributed total mass */
      REAL4 minMass= mass1Min + mass2Min;
      REAL4 maxMass= mass1Max + mass2Max;
      mTotal = minMass + XLALUniformDeviate( randParams ) * (maxMass - minMass);
      inj->mass2 = -1.0;
      while( inj->mass2 > mass2Max || inj->mass2 <= mass2Min )
      {
        inj->mass1 = mass1Min + 
          XLALUniformDeviate( randParams ) * (mass1Max - mass1Min);
        inj->mass2 = mTotal - inj->mass1;
      }
    }
  }

  inj->eta = inj->mass1 * inj->mass2 / ( mTotal * mTotal );
  inj->mchirp = mTotal * pow(inj->eta, 0.6);

  return ( inj );
}  

/** Generates masses for an inspiral injection. Masses are Gaussian distributed
 * with the requested mean and standard deviation. */
SimInspiralTable* XLALGaussianInspiralMasses( 
    SimInspiralTable *inj,   /**< injection for which masses will be set*/
    RandomParams *randParams,/**< random parameter details*/
    REAL4  mass1Mean,        /**< mean value for mass1 */
    REAL4  mass1Std,         /**< standard deviation of mass1 */
    REAL4  mass2Mean,        /**< mean value of mass2 */
    REAL4  mass2Std          /**< standard deviation of mass2 */
    )
{       
  REAL4Vector *normalDev;
  REAL4 mtotal;

  normalDev = XLALCreateREAL4Vector( 2 );
  XLALNormalDeviates(normalDev, randParams);
  inj->mass1 = mass1Mean + mass1Std * normalDev->data[0];
  inj->mass2 = mass1Mean + mass1Std * normalDev->data[1];
  mtotal = inj->mass1 + inj->mass2;
  XLALDestroyVector( normalDev );

  inj->eta = inj->mass1 * inj->mass2 / ( mtotal * mtotal );
  inj->mchirp = mtotal * pow(inj->eta, 0.6);

  return ( inj );
}  

/** Generates spins for an inspiral injection.  Spin magnitudes lie between the
 * specified max and min values.  Orientations are random */
SimInspiralTable* XLALRandomInspiralSpins( 
    SimInspiralTable *inj,   /**< injection for which spins will be set*/
    RandomParams *randParams,/**< random parameter details*/
    REAL4  spin1Min,         /**< minimum magnitude of spin1 */
    REAL4  spin1Max,         /**< maximum magnitude of spin1 */
    REAL4  spin2Min,         /**< minimum magnitude of spin2 */
    REAL4  spin2Max          /**< maximum magnitude of spin2 */ 
    )
{       
  REAL4 spin1Mag;
  REAL4 spin2Mag;
  REAL4 r1;
  REAL4 r2;
  REAL4 phi1;
  REAL4 phi2;

  /* spin1Mag */
  spin1Mag =  spin1Min + XLALUniformDeviate( randParams ) * 
    (spin1Max - spin1Min);

  /* spin1z */
  inj->spin1z = (XLALUniformDeviate( randParams ) - 0.5) * 2 * (spin1Mag);
  r1 = pow( ((spin1Mag * spin1Mag) - (inj->spin1z * inj->spin1z)) , 
      0.5); 

  /* phi1 */
  phi1 = XLALUniformDeviate( randParams ) * LAL_TWOPI;

  /* spin1x and spin1y */  
  inj->spin1x = r1 * cos(phi1);
  inj->spin1y = r1 * sin(phi1);

  /* spin2Mag */
  spin2Mag =  spin2Min + XLALUniformDeviate( randParams ) * 
    (spin2Max - spin2Min);

  /* spin2z */
  inj->spin2z = (XLALUniformDeviate( randParams ) - 0.5) * 2 * (spin2Mag);
  r2 = pow( ((spin2Mag * spin2Mag) - (inj->spin2z * inj->spin2z)) , 
      0.5); 

  /* phi2 */
  phi2 = XLALUniformDeviate( randParams ) * LAL_TWOPI;

  /* spin2x and spin2y */  
  inj->spin2x = r2 * cos(phi2);
  inj->spin2y = r2 * sin(phi2);


  return ( inj );
}

/** Set end time and effective distance of an injection for a detector */
SimInspiralTable *XLALInspiralSiteTimeAndDist( 
    SimInspiralTable  *inj, /**< the injection details */
    LALDetector       *detector, /**< the detector of interest */
    LIGOTimeGPS       *endTime,  /**< the end time to populate */
    REAL4             *effDist   /**< the effective distance to populate */
    )
{

  REAL8                 tDelay, fplus, fcross;
  REAL4                 splus, scross, cosiota;

  /* calculate the detector end time */
  tDelay = XLALTimeDelayFromEarthCenter( detector->location, inj->longitude,
      inj->latitude, &(inj->geocent_end_time) );
  *endTime = *XLALGPSAdd( &(inj->geocent_end_time), tDelay );

  /* initialize distance with real distance and compute splus and scross */
  *effDist = 2.0 * inj->distance;
  cosiota = cos( inj->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;


  /* calculate the detector response */
  XLALComputeDetAMResponse(&fplus, &fcross, detector->response, inj->longitude, 
      inj->latitude, inj->polarization, inj->end_time_gmst);

  /* compute the effective distance */
  *effDist /= sqrt( 
      splus*splus*fplus*fplus + scross*scross*fcross*fcross );

  return ( inj );
}

/** Set the end time and effective distance for all detectors for this 
 * injection */
SimInspiralTable *XLALPopulateSimInspiralSiteInfo(
    SimInspiralTable           *inj /**< the injection */ 
    )
{
  LALDetector           detector; 
  REAL4                *eff_dist;
  LIGOTimeGPS          *end_time;

  /* LIGO Hanford observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  end_time = &(inj->h_end_time);
  eff_dist = &(inj->eff_dist_h);
  inj = XLALInspiralSiteTimeAndDist(inj, &detector, end_time, eff_dist);

  /* LIGO Livingston observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  end_time = &(inj->l_end_time);
  eff_dist = &(inj->eff_dist_l);
  inj = XLALInspiralSiteTimeAndDist(inj, &detector, end_time, eff_dist);

  /* GEO observatory*/
  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  end_time = &(inj->g_end_time);
  eff_dist = &(inj->eff_dist_g);
  inj = XLALInspiralSiteTimeAndDist(inj, &detector, end_time, eff_dist);

  /* TAMA observatory*/
  detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  end_time = &(inj->t_end_time);
  eff_dist = &(inj->eff_dist_t);
  inj = XLALInspiralSiteTimeAndDist(inj, &detector, end_time, eff_dist);

  /* Virgo observatory*/
  detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  end_time = &(inj->v_end_time);
  eff_dist = &(inj->eff_dist_v);
  inj = XLALInspiralSiteTimeAndDist(inj, &detector, end_time, eff_dist);

  return ( inj );

}

/** Populate a frequency series with the actuation response.  Here, we just use
 * the pendulum part of the actuation function */
COMPLEX8FrequencySeries *generateActuation( 
    COMPLEX8FrequencySeries *resp, /* the frequency series to be populated */
    REAL4                    ETMcal,/* the ETM calibration */
    REAL4                    pendF, /* the pendulum frequency */
    REAL4                    pendQ  /* the pendulum Q */ 
    )
{
  INT4   k;
  REAL4  fNorm = 0;
  COMPLEX8Vector *num = NULL;  
  COMPLEX8Vector *denom = NULL;  

  num = XLALCreateCOMPLEX8Vector( resp->data->length );
  denom = XLALCreateCOMPLEX8Vector( resp->data->length );

  /* the response function is
   *
   *                  ETMcal * pendF^2
   * R(f) =  ---------------------------------------
   *         pendF^2 - freq^2 - i freq (pendF/pendQ)
   */

  for ( k = 0; k < resp->data->length; k++ )
  {
    fNorm = k * resp->deltaF / pendF;
    denom->data[k].re = ( 1 - fNorm * fNorm );
    denom->data[k].im = - fNorm / pendQ;
    num->data[k].re = 1.0 * ETMcal;
    num->data[k].im = 0.0;
  }

  XLALCCVectorDivide( resp->data, num, denom);
  XLALDestroyCOMPLEX8Vector( num );
  XLALDestroyCOMPLEX8Vector( denom );
  return( resp );

}

