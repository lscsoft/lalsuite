/*
*  Copyright (C) 2007 Drew Keppel, Gareth Jones, Patrick Brady, Thomas Cokelaer
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
 * Author : 
 * 
 * Purpose : generate xml file for binary injections (spinning case)
 *
 * */
#include <stdio.h>
#include <stdlib.h>
#include <config.h>

#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/LALErrno.h>
#include <LALAppsVCSInfo.h>


#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "spininj"

#define USAGE \
"lalapps_spininj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                print mass and galactocentic cartesian coordinates\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --waveform WVF           set the injection waveform to WVF\n"\
"                           (EOB, TaylorT1, TaylorT3,PadeT1; followed by the\n"\
"                           order: newtonian, onePN, onePointFivePN, twoPN,\n"\
"                           twoPointFivePN, threePN) (default: EOBtwoPN)\n"\
"  --fLower                 set the lower cutoff frequency of the isnpiral waveform (40)\n\n"\
"  --dist-distr DDISTR      set method of simulated source distance distr to DDISTR \n"\
"                           SPININJ_distance    : uniform distr of sources in distance \n"\
"                           SPININJ_logDistance : uniform distr of sources in log(distance) \n"\
"                           SPININJ_volume      : uniform distr of sources in volume \n"\
"                           Default is SPININJ_logDistance \n"\
"  --distance-min           set minimal value of simulated sources distance in kpc \n"\
"  --distance-max           set maximal value of simulated sources distance in kpc \n\n"\
"  --coa-phase-min          set minimal value of the initial orbital angle coa-phase  (0)\n"\
"  --coa-phase-max          set maximal value of the initial orbital angle coa-phase  (2 pi)\n"\
"  --coa-phase-range        set range of the initial orbital angle coa-phase  (0 -  2 pi)\n\n"\
"  --inclination-min        set minimal value of the initial inclination  (>=0.01 \n"\
"  --inclination-max        set maximal value of the initial inclination  (<=pi)\n"\
"  --inclination-range      set range of the initial inclination (0.01 -  pi)\n\n"\
"  --spin1-min              set minimal value of the initial spin  (>=0)\n"\
"  --spin1-max              set maximal value of the initial spin  (<=1)\n"\
"  --spin1-range            set range of the initial spin (0 -  1)\n\n"\
"  --spin2-min              set minimal value of the initial spin (>=0)\n"\
"  --spin2-max              set maximal value of the initial spin (<=1)\n"\
"  --spin2-range            set range of the initial spin (0 - 1)\n\n"\
"  --mass-distr MDISTR      set method of simulated source mass distr to MDISTR \n"\
"                           SPININJ_m1Andm2     : uniform distr of sources component mass \n"\
"                           SPININJ_totalMass   : uniform distr of sources total mass \n"\
"                           Default is SPININJ_totalMass \n"\
"  --mass1-range            set range of the first mass (3 - 20 )\n"\
"  --mass2-range            set range of the second mass (3 - 20)\n\n"\
"\n"

const INT4            S2StartTime = 729273613; /* Feb 14 2003 16:00:00 UTC */
const INT4            S2StopTime  = 734367613; /* Apr 14 2003 15:00:00 UTC */


typedef enum{ 
  SPININJ_m1Andm2, 
  SPININJ_totalMass 
} massEnum;

typedef enum{
  SPININJ_distance,
  SPININJ_logDistance, 
  SPININJ_volume
} distributionEnum;

typedef struct{
  double  min;
  double  max;
} SPININJrange;


typedef struct {
  REAL4           fLower;
  LIGOTimeGPS     gpsStartTime;
  LIGOTimeGPS     gpsEndTime;
  UINT4           randSeed;
  SPININJrange         distance;
  SPININJrange         m1, m2;
  REAL8           timeInterval;
  REAL8           meanTimeStep;
  CHAR            waveform[LIGOMETA_WAVEFORM_MAX];  
  CHAR            *userTag ;  
  massEnum         mdistr;
  distributionEnum ddistr;
  REAL4           theta0, phi0; 
  SPININJrange         spin1,  spin2;
  SPININJrange         coa_phase;
  SPININJrange         inclination;
} InspiralInjectionParameters;



void LALSimInspiralTablePopulate(LALStatus        *status, 
                                 MetadataTable     injections, 
                                 SimInspiralTable *this_inj);



void LALSetIndividualMasses(LALStatus             *status,
                            InspiralInjectionParameters params,
                            SimInspiralTable *this_inj);



void LALSetDistance(LALStatus *status, 
                    InspiralInjectionParameters,
                    SimInspiralTable *this_inj);

void LALSetSpin(LALStatus                   *status,
                InspiralInjectionParameters params,
                SimInspiralTable            *this_inj  );

void LALSetSpatialDistribution(LALStatus *status,
                               InspiralInjectionParameters params,
                               SimInspiralTable *this_inj);

void LALSetSiteParameters(LALStatus *status,
                          SimInspiralTable *this_inj);





void LALSetGeoCentricEndTime(LALStatus *status,        
                             InspiralInjectionParameters params,
                             SimInspiralTable *this_inj);

void LALParserInspiralInjection(LALStatus *status,
                                int,
                                char**,
                                InspiralInjectionParameters *);


void XLALCheckInspiralInjectionParameters(InspiralInjectionParameters params);



static ProcessParamsTable *next_process_param( const char *name, const char *type,
    const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );
  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  snprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}


extern int vrbflg;
LIGOLwXMLStream       xmlfp;
CHAR                  fname[256];
RandomParams *randParams = NULL;


int main( int argc, char *argv[] )
{

  static  LALStatus             status ;

  /* command line options */

  InspiralInjectionParameters         paramsIn;     
  MetadataTable                       injections;
  SimInspiralTable                    *this_inj = NULL;

    /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  
  LAL_CALL( LALParserInspiralInjection(&status, argc, argv, &paramsIn), 
           &status);
  
  LAL_CALL( LALCreateRandomParams( &status, &randParams, paramsIn.randSeed ), 
            &status );  
  
  /* null out the head of the linked list */
  injections.simInspiralTable = NULL;
  
  /*
   *
   * loop over duration of desired output times
   *
   */

  while ( XLALGPSCmp(&(paramsIn.gpsStartTime), &(paramsIn.gpsEndTime)) < 0 )
  {

    /* rho, z and lGal are the galactocentric galactic axial coordinates */
    /* r and phi are the geocentric galactic spherical coordinates       */

    /* create the sim_inspiral table */
    if ( injections.simInspiralTable )
    {
      this_inj = this_inj->next = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }
    else
    {
      injections.simInspiralTable = this_inj = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }

    LAL_CALL(LALSetGeoCentricEndTime(&status, paramsIn, this_inj), 
             &status);
    LAL_CALL( LALSetIndividualMasses(&status, paramsIn, this_inj), 
              &status);   
    LAL_CALL( LALSetSpin(&status, paramsIn, this_inj), 
              &status); 
    LAL_CALL( LALSetDistance(&status, paramsIn, this_inj), 
              &status);
    LAL_CALL( LALSetSpatialDistribution(&status, paramsIn, this_inj),
              &status);  
    LAL_CALL( LALSetSiteParameters(&status, this_inj), 
              &status);

    /* increment the injection time */
    XLALGPSAdd( &(paramsIn.gpsStartTime), paramsIn.meanTimeStep );

    this_inj->f_lower  = paramsIn.fLower;
    memcpy (this_inj->waveform, paramsIn.waveform, sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX);  
  } /* end loop over injection times */

  LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );

  /*
   *
   * write output to LIGO_LW XML file
   *
   */

  /* write the sim_inspiral table */
  if ( injections.simInspiralTable ){
    LAL_CALL(LALSimInspiralTablePopulate(&status, 
                                         injections, 
                                         this_inj), &status);
  }
  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}



/*################################################*/
void LALSimInspiralTablePopulate(LALStatus        *status, 
                                 MetadataTable     injections, 
                                 SimInspiralTable *this_inj)
{
  
  LAL_CALL( LALBeginLIGOLwXMLTable( status, &xmlfp, sim_inspiral_table ), 
            status );
  LAL_CALL( LALWriteLIGOLwXMLTable( status, &xmlfp, injections, 
                                    sim_inspiral_table ), status );
  LAL_CALL( LALEndLIGOLwXMLTable ( status, &xmlfp ), status );
  
  while ( injections.simInspiralTable )
    {
      this_inj = injections.simInspiralTable;
      injections.simInspiralTable = injections.simInspiralTable->next;
      LALFree( this_inj );
    }
}




/* TODO check uniformity of totalmass and (2) set mchirp*/

void LALSetIndividualMasses(LALStatus                   *status,
                            InspiralInjectionParameters params,
                            SimInspiralTable            *this_inj  )
{ 
      
  REAL4 deltaM1 = params.m1.max - params.m1.min;
  REAL4 deltaM2 = params.m2.max - params.m2.min;

  
  REAL4 u, mtotal; 

  switch(params.mdistr){
  case  SPININJ_m1Andm2:
    LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
    this_inj->mass1 = params.m1.min + u * deltaM1;
    LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
    this_inj->mass2 = params.m2.min + u * deltaM2;
    mtotal = this_inj->mass1 + this_inj->mass2 ;
    this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
    break;
  case SPININJ_totalMass:
    mtotal = params.m1.min + params.m2.min;
    LAL_CALL( LALUniformDeviate( status, &u, randParams ), status);
    mtotal += u* (deltaM1+deltaM2);

    LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
    this_inj->mass1 =  params.m1.min + u * deltaM1;
    this_inj->mass2 = mtotal - this_inj->mass1;
    while (this_inj->mass2 >= params.m2.max 
           || this_inj->mass2 <=  params.m2.min )
      {
        LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
        this_inj->mass1 =  params.m1.min + u * deltaM1;
        this_inj->mass2 = mtotal - this_inj->mass1;
      }
    this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
    break;
  }

  this_inj->mchirp = pow( this_inj->eta, 0.6) * 
      (this_inj->mass1 + this_inj->mass2);

}


void LALSetSpin(LALStatus                   *status,
                InspiralInjectionParameters params,
                SimInspiralTable            *this_inj  )
{       
  REAL4 u; 
  REAL4 spin1Mag;
  REAL4 spin2Mag;
  REAL4 r1;
  REAL4 r2;
  REAL4 phi1;
  REAL4 phi2;

 /* spin1Mag */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 spin1Mag =  params.spin1.min +  u * (params.spin1.max - params.spin1.min);

 /* spin1z */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 this_inj->spin1z = (u - 0.5) * 2 * (spin1Mag);

 r1 = pow( ((spin1Mag * spin1Mag) - (this_inj->spin1z * this_inj->spin1z)) , 0.5); 

 /* phi1 */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status ); 
 phi1 = u * LAL_TWOPI;
 
 /* spin1x and spin1y */  
 this_inj->spin1x = r1 * cos(phi1);
 this_inj->spin1y = r1 * sin(phi1);

  /* spin2Mag */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 spin2Mag = params.spin2.min + u * (params.spin2.max - params.spin2.min);

 /* spin2z */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 this_inj->spin2z = (u - 0.5) * 2 * (spin2Mag);

 r2 = pow( ((spin2Mag * spin2Mag) - (this_inj->spin2z * this_inj->spin2z)) , 0.5);

 /* phi2 */
 LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 phi2 = u * LAL_TWOPI;

 /* spin2x and spin2y */
 this_inj->spin2x = r2 * cos(phi2);
 this_inj->spin2y = r2 * sin(phi2);

}

/* TODO to check */
void LALSetDistance(LALStatus *status, 
                    InspiralInjectionParameters params,
                    SimInspiralTable *this_inj)
{
  REAL4 u,  deltaL, lmin,lmax, d3min,d3max,deltad3,d3, exponent;

  switch (params.ddistr){
  case SPININJ_distance:
    LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
    this_inj->distance = params.distance.min + u * (params.distance.max - params.distance.min);

    break;
    
  case SPININJ_logDistance:
    lmin = log10(params.distance.min);
    lmax = log10(params.distance.max);
    deltaL = lmax - lmin;
    LAL_CALL(  LALUniformDeviate(status,&u,randParams),status );
    exponent = lmin + deltaL * u;
    this_inj->distance = pow(10.0,(REAL4) exponent);
    break;
    
  case SPININJ_volume:
    d3min = params.distance.min * params.distance.min * params.distance.min;
    d3max = params.distance.max * params.distance.max * params.distance.max;
    deltad3 = d3max - d3min ;
    LAL_CALL(  LALUniformDeviate(status,&u,randParams),status );
    d3 = d3min + u * deltad3 ;
    this_inj->distance = pow(d3, 1.0/3.0);
    break;
  }
    
  this_inj->distance = this_inj->distance / 1000.0; /*convert to Mpc */



}


void LALSetSpatialDistribution(LALStatus *status,
                               InspiralInjectionParameters params,
                               SimInspiralTable *this_inj
                               )
{
  REAL4 u; 

  /*TO check */  
  LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
 
  /* we are setting min value of inclination to 0.01
   * inclination = 0 causes initial orb ang mmtm to coincide with 
   * a coord singularity in the injection code. 
   */ 

  this_inj->inclination = 
        acos( cos(params.inclination.min) 
        + u *(cos(params.inclination.max) - cos(params.inclination.min)));
  
  /* provide an input argument ? */
  LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
  this_inj->polarization = LAL_TWOPI * u ;
  
  LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
  this_inj->coa_phase =  params.coa_phase.min + u * (params.coa_phase.max - params.coa_phase.min);
  
  /* compute random longitude and latitude ; TODO input arguments ? */
  LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
  this_inj->longitude = LAL_TWOPI * u;
  LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
  this_inj->latitude = asin(2.0 * u - 1);  
  
}    


void LALSetSiteParameters(LALStatus *status,
                          SimInspiralTable *this_inj)
                       
{
  SkyPosition           skyPos;
  LALSource             source;
  LALDetAndSource       detAndSource;
  REAL4 cosiota, splus, scross;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAMResponse      resp;
  REAL8                 time_diff;
  
  /* set up params for the site end times and detector response */
  memset( &skyPos, 0, sizeof(SkyPosition) );
  memset( &source, 0, sizeof(LALSource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );
  
  skyPos.longitude = this_inj->longitude;
  skyPos.latitude  = this_inj->latitude;
  skyPos.system    = COORDINATESYSTEM_EQUATORIAL;
  
  source.equatorialCoords = skyPos;
  source.orientation      = this_inj->polarization;

  detAndSource.pSource = &source;

  /*
   * compute site end times
   * (copied from SnglInspiralUtils.c)
   */
  
  /* initialize end times with geocentric value */
  this_inj->h_end_time = this_inj->l_end_time = this_inj->geocent_end_time;
  
  /* lho */
  time_diff = XLALTimeDelayFromEarthCenter( lho.location, this_inj->longitude, this_inj->latitude, &(this_inj->geocent_end_time) );
  XLALGPSAdd( &(this_inj->h_end_time), time_diff );

  /* llo */
  time_diff = XLALTimeDelayFromEarthCenter( llo.location, this_inj->longitude, this_inj->latitude, &(this_inj->geocent_end_time) );
  XLALGPSAdd( &(this_inj->l_end_time), time_diff );

  /* temporarily, populate the fields for the */
  /* GEO, TAMA and VIRGO times                */
  
  memset( &(this_inj->g_end_time), 0,sizeof((this_inj->l_end_time)) );
  memset( &(this_inj->t_end_time), 0,sizeof((this_inj->l_end_time)) );
  memset( &(this_inj->v_end_time), 0,sizeof((this_inj->l_end_time)) );
  
  /*
   * compute the effective distance of the inspiral
   * (copied from SnglInspiralUtils.c)
   */

  
  /* initialize distances with real distance and compute splus and scross*/
  this_inj->eff_dist_h = this_inj->eff_dist_l = 2.0 * this_inj->distance;
  cosiota = cos( this_inj->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;


  /* compute the response of the LHO detectors */   
  
  detAndSource.pDetector = &lho;
  LAL_CALL( LALComputeDetAMResponse( status, &resp, &detAndSource,
                                     &this_inj->geocent_end_time ), status );
  

  /* compute the effective distance for LHO */
  this_inj->eff_dist_h /= sqrt( splus*splus*resp.plus*resp.plus +
                                scross*scross*resp.cross*resp.cross );

  
  /* compute the response of the LLO detector */
  detAndSource.pDetector = &llo;
  LAL_CALL( LALComputeDetAMResponse( status, &resp, &detAndSource,
                                     &this_inj->geocent_end_time ), status);

  
  /* compute the effective distance for LLO */
  this_inj->eff_dist_l /= sqrt( splus*splus*resp.plus*resp.plus 
                                + scross*scross*resp.cross*resp.cross );
  
  /* temporarily, populate the fields for the */
  /* GEO, TAMA and VIRGO effective distances  */

  
  memset( &(this_inj->eff_dist_g), 0, sizeof((this_inj->eff_dist_g)) );
  memset( &(this_inj->eff_dist_t), 0, sizeof((this_inj->eff_dist_g)) );
  memset( &(this_inj->eff_dist_v), 0, sizeof((this_inj->eff_dist_g)) );
  
  


}



void LALSetGeoCentricEndTime(LALStatus *status,
                             InspiralInjectionParameters params,
                             SimInspiralTable *this_inj)
{
  REAL4 u;
  
  /* set the geocentric end time of the injection */
  /* XXX CHECK XXX */
  this_inj->geocent_end_time = params.gpsStartTime;
  if ( params.timeInterval )
    {
      LAL_CALL( LALUniformDeviate( status, &u, randParams ), status );
      XLALGPSAdd( &(this_inj->geocent_end_time), u * params.timeInterval );
    }
  
  /* set gmst */
  this_inj->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
      &this_inj->geocent_end_time), LAL_TWOPI) * 24.0 / LAL_TWOPI; /* hours */
  ASSERT( !XLAL_IS_REAL8_FAIL_NAN(this_inj->end_time_gmst), status, LAL_FAIL_ERR, LAL_FAIL_MSG );
}



void LALParserInspiralInjection(LALStatus *status, 
                                int argc, 
                                char **argv, 
                                InspiralInjectionParameters *params)
{
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  int i;
  


  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );

  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  
  /* clear the waveform field */
  memset( params->waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  
  /* Default values */
  params->theta0                      = 0;    /* variable not used for spinning inj */
  params->inclination.min             = 0.01; /* can not be set to zero */
  params->inclination.max             = LAL_PI;
  params->phi0                        = 0.0; /* set =0 with no less of generality */
  params->spin1.min                   = 0.;
  params->spin2.min                   = 0.;
  params->spin1.max                   = 1.;
  params->spin2.max                   = 1.;
  params->coa_phase.min               = 0.;
  params->coa_phase.max               = LAL_TWOPI;
  params->fLower                      = 40; 
  params->gpsStartTime.gpsSeconds     = S2StartTime;
  params->gpsStartTime.gpsNanoSeconds = 0;
  params->gpsEndTime.gpsSeconds       = S2StopTime;  
  params->gpsEndTime.gpsNanoSeconds   = 0;
  params->randSeed                    = 1;
  params->m1.min                      = 3.;
  params->m1.max                      = 20.;
  params->m2.min                      = 3.;
  params->m2.max                      = 20.;
  params->timeInterval                = 1000.;
  params->meanTimeStep                = 2630 / LAL_PI;  
  params->distance.min                = 1.0;
  params->distance.max                = 20000.0;
  params->mdistr                      = SPININJ_totalMass;
  params->ddistr                      = SPININJ_logDistance;
  params->userTag                     = NULL;
  snprintf( params->waveform, LIGOMETA_WAVEFORM_MAX, "EOBtwoPN");  
  snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
               params->randSeed, params->gpsStartTime.gpsSeconds, 
               params->gpsEndTime.gpsSeconds - params->gpsStartTime.gpsSeconds );
  
  /* Actual parsing is here; we do not check validity of the input arguments 
     for the time being*/
  i = 1;

  while (i < argc)
    {
      if ( (strcmp(argv[i], "--help") == 0) ||  (strcmp(argv[i], "-h") == 0) ){
        fprintf( stderr, USAGE );
        exit(1);
      }
      else if ( strcmp(argv[i], "--user-tag") == 0){
        /* create storage for the usertag */
        params->userTag = (CHAR *) calloc( strlen( argv[i+1] ) + 1, sizeof(CHAR) );
        memcpy( params->userTag, argv[i+1], strlen( argv[i+1] ) + 1 );
        this_proc_param = this_proc_param->next = 
          next_process_param( "user-tag",  "string", "%s",  argv[i+1] );
        i++;
      }
      else if  ( strcmp(argv[i], "--gps-start-time") == 0){               
        params->gpsStartTime.gpsSeconds = atol(argv[++i]);
        this_proc_param = this_proc_param->next = 
          next_process_param( "gps-start-time", "int",  "%ld", params->gpsStartTime.gpsSeconds);
      }
      else if  ( strcmp(argv[i], "--gps-end-time") == 0){
        params->gpsEndTime.gpsSeconds = atol(argv[++i]);
        this_proc_param = this_proc_param->next = 
          next_process_param( "gps-end-time", "int",  "%ld", params->gpsEndTime.gpsSeconds );
      }
      else if ( strcmp(argv[i], "--seed") == 0){
        params->randSeed = atoi( argv[++i]);
        this_proc_param = this_proc_param->next = 
          next_process_param( "seed", "int", "%d", params->randSeed );
      }      
      else if ( strcmp(argv[i] , "--time-interval") == 0 ){
        params->timeInterval = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "time-interval", "float", "%le", params->timeInterval );
      }
      else if ( strcmp(argv[i] , "--time-step") == 0 ){
        params->meanTimeStep = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "time-step", "float",  "%le", params->meanTimeStep );
      }
      else if ( strcmp(argv[i] , "--coa-phase-min") == 0 ){
        params->coa_phase.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "coa-phase-min", "float", "%le", params->coa_phase.min );
      }  
      else if ( strcmp(argv[i] , "--coa-phase-max") == 0 ){
        params->coa_phase.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "coa-phase-max", "float", "%le", params->coa_phase.max );
      }
      else if ( strcmp(argv[i] , "--coa-phase-range") == 0 ){
        params->coa_phase.min = atof( argv[++i] );
        params->coa_phase.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "coa-phase-min", "float", "%le", params->coa_phase.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "coa-phase-max", "float", "%le", params->coa_phase.max );
      } 
      else if ( strcmp(argv[i] , "--mass1-min") == 0 ){
        params->m1.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass1-min", "float", "%le", params->m1.min );
      }  
      else if ( strcmp(argv[i] , "--mass1-max") == 0 ){
        params->m1.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass1-max", "float", "%le", params->m1.max );
      }
      else if ( strcmp(argv[i] , "--mass1-range") == 0 ){
        params->m1.min = atof( argv[++i] );
        params->m1.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass1-min", "float", "%le", params->m1.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass1-max", "float", "%le", params->m1.max );
      } 
      else if ( strcmp(argv[i] , "--fLower") == 0 ){
        params->fLower = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "fLower", "float", "%le", params->fLower );
      } 
      else if ( strcmp(argv[i] , "--mass2-min") == 0 ){
        params->m2.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass2-min", "float", "%le", params->m2.min );
      }  
      else if ( strcmp(argv[i] , "--mass2-max") == 0 ){
        params->m2.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass2-max", "float", "%le", params->m2.max );
      }
      else if ( strcmp(argv[i] , "--mass2-range") == 0 ){
        params->m2.min = atof( argv[++i] );
        params->m2.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass2-min", "float", "%le", params->m2.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "mass2-max", "float", "%le", params->m2.max );
      }
      else if ( strcmp(argv[i] , "--distance-min") == 0 ){
        params->distance.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "distance-min", "float", "%le", params->distance.min );
      } 
      else if ( strcmp(argv[i] , "--distance-max") == 0 ){
        params->distance.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "distance-max", "float", "%le", params->distance.max );
      }  
      else if ( strcmp(argv[i] , "--distance-range") == 0 ){
        params->distance.min = atof( argv[++i] );
        params->distance.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "distance-min", "float", "%le", params->distance.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "distance-max", "float", "%le", params->distance.max );
      } 
      else if ( strcmp(argv[i] , "--spin1-min") == 0 ){
        params->spin1.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin1-min", "float", "%le", params->spin1.min );
      } 
      else if ( strcmp(argv[i] , "--spin2-min") == 0 ){
        params->spin2.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin2-min", "float", "%le", params->spin2.min );
      }
      else if ( strcmp(argv[i] , "--spin1-max") == 0 ){
        params->spin1.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin1-max", "float", "%le", params->spin1.max );
      } 
      else if ( strcmp(argv[i] , "--spin2-max") == 0 ){
        params->spin2.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin2-max", "float", "%le", params->spin2.max );
      }
      else if ( strcmp(argv[i] , "--spin1-range") == 0 ){
        params->spin1.min = atof( argv[++i] );
        params->spin1.max = atof( argv[++i] );
        
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin1-min", "float", "%le", params->spin1.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin1-max", "float", "%le", params->spin1.max );
      }  
      else if ( strcmp(argv[i] , "--spin2-range") == 0 ){
        params->spin2.min = atof( argv[++i] );
        params->spin2.max = atof( argv[++i] );

        this_proc_param = this_proc_param->next = 
          next_process_param( "spin2-min", "float", "%le", params->spin2.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "spin2-max", "float", "%le", params->spin2.max );
      } 
      else if ( strcmp(argv[i] , "--inclination-min") == 0 ){
        params->inclination.min = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "inclination-min", "float", "%le", params->inclination.min );
      } 
      else if ( strcmp(argv[i] , "--inclination-max") == 0 ){
        params->inclination.max = atof( argv[++i] );
        this_proc_param = this_proc_param->next = 
          next_process_param( "inclination-max", "float", "%le", params->inclination.max );
      } 
      else if ( strcmp(argv[i] , "--inclination-range") == 0 ){
        params->inclination.min = atof( argv[++i] );
        params->inclination.max = atof( argv[++i] );

        this_proc_param = this_proc_param->next = 
          next_process_param( "inclination-min", "float", "%le", params->inclination.min );
        this_proc_param = this_proc_param->next = 
          next_process_param( "inclination-max", "float", "%le", params->inclination.max );
      } 
      else if ( strcmp(argv[i] , "--waveform") == 0 ){
        snprintf( params->waveform, LIGOMETA_WAVEFORM_MAX, "%s", argv[++i]);
        this_proc_param = this_proc_param->next =
          next_process_param( "waveform", "string",
                              "%s",argv[i] );        
      }
      else if ( strcmp(argv[i] , "--mass-distr") == 0 ){
         if ( ! strcmp( "SPININJ_m1Andm2", argv[++i] ) )
          {
            params->mdistr = SPININJ_m1Andm2;
          }
          else if ( ! strcmp( "SPININJ_totalMass", argv[i] ) )
          {
            params->mdistr = SPININJ_totalMass;
          }
          else
          {
            fprintf( stderr, "invalid arg to --mass-distr \n");
            exit(1); 
          }
        this_proc_param = this_proc_param->next =
          next_process_param( "mass-distr", "string",
                              "%s",argv[i] );        
      }
       else if ( strcmp(argv[i] , "--dist-distr") == 0 ){
         if ( ! strcmp( "SPININJ_distance", argv[++i] ) )
          {
            params->ddistr = SPININJ_distance;
          }
          else if ( ! strcmp( "SPININJ_logDistance", argv[i] ) )
          {
            params->ddistr = SPININJ_logDistance;
          }
          else if ( ! strcmp( "SPININJ_volume", argv[i] ) )
          {
            params->ddistr = SPININJ_volume;
          }
          else
          {
            fprintf( stderr, "invalid arg to --dist-distr \n");
            exit(1);
          }
        this_proc_param = this_proc_param->next =
          next_process_param( "dist-distr", "string",
                              "%s",argv[i] );
      }
      else {
        fprintf(stderr, "option %s unknown !!! abort\n", argv[i]);
        exit(1);
      }
          
    
      i++;
    }
  
  if (params->userTag){
    snprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
              params->randSeed, params->userTag, params->gpsStartTime.gpsSeconds, 
              params->gpsEndTime.gpsSeconds - params->gpsStartTime.gpsSeconds );
  }

  /* Let us check now the validity of the arguments */
  XLALCheckInspiralInjectionParameters(*params);


  
  /* and finally stored the arguments in the process table */

  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( status, &xmlfp, fname), status );
  
  /* write the process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "H1H2L1" );
  XLALGPSTimeNow(&(proctable.processTable->end_time));

  
  LAL_CALL(  LALBeginLIGOLwXMLTable( status, &xmlfp, process_table ),status), 

  LAL_CALL(  LALWriteLIGOLwXMLTable( status, &xmlfp, proctable, 
                                     process_table ), status );  

  
  LAL_CALL(  LALEndLIGOLwXMLTable ( status, &xmlfp ), status ); 


  free( proctable.processTable );
  
 this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

     /* write the process params table */
  if ( procparams.processParamsTable )
  {
      LAL_CALL( LALBeginLIGOLwXMLTable( status, &xmlfp, process_params_table ), 
        status );
      LAL_CALL( LALWriteLIGOLwXMLTable( status, &xmlfp, procparams, 
          process_params_table ), status );
      LAL_CALL( LALEndLIGOLwXMLTable ( status, &xmlfp ), status );
    while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
  }
  /* We dont close right now. it will be done within the main function 
     after having written the injection data*/

  /*  DETATCHSTATUSPTR(status);
      RETURN (status);  */
}


void XLALCheckInspiralInjectionParameters(InspiralInjectionParameters params)
{ 
  if ( params.fLower <= 5 || params.fLower >=1000)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "fLower should not be less than 5Hz or greater than 1000Hz " 
               "(%f specified)\n",
               "fLower", params.fLower );
      exit( 1 );
    } 
  if ( params.gpsStartTime.gpsSeconds < 441417609 )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "GPS start time is prior to " 
               "Jan 01, 1994  00:00:00 UTC:\n"
               "(%d specified)\n",
               "gps-start-time", params.gpsStartTime.gpsSeconds );
      exit( 1 );
    } 
  /* check that the start time is before the end time */
  if ( params.meanTimeStep <= 0 )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "time step must be > 0: (%e seconds specified)\n",
               "mean-time", params.meanTimeStep );
      exit( 1 );
    }
  if ( params.timeInterval < 0 )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "time interval must be >= 0: (%e seconds specified)\n",
               "time-interval", params.timeInterval );
      exit( 1 );
    }
  if ( params.m1.min < 0 || params.m2.max < 0 )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "miniumum component mass must be > 0: \n",
               "mass-range"  );

      exit( 1 );
    } 
  if ( params.m1.min > params.m1.max  )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "minimum component (%f) mass must be < maximum component (%f) : \n",
               "mass-range",  params.m1.min, params.m1.max);

      exit( 1 );
    }
  if ( params.m2.min > params.m2.max  )
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "minimum component (%f) mass must be < maximum component (%f) : \n",
               "mass-range",  params.m2.min, params.m2.max);
      
      exit( 1 );
    }
  if ( params.distance.min <= 0)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "minimum distance must be > 0: "
               "(%f kpc specified)\n", "distance-range",params.distance.min);
      
      exit( 1 );
    }
  if (  params.distance.max < params.distance.min)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "maximum distance max (%f) must be > distance min (%f)",
               "distance-range",params.distance.max, params.distance.min);
      
      exit( 1 );
    }
  if (  params.coa_phase.max < params.coa_phase.min)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "maximum coa_phase max (%f) must be > coa_phase min (%f)",
               "coa-phase-range",params.coa_phase.max, params.coa_phase.min);
      
      exit( 1 );
    }


  if (  params.spin1.min < 0 || params.spin1.max < 0)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "spin1-min (%f) and spin1-max (%f) must be > 0 \n",
               "spin1-min or spin1-max",params.spin1.min, params.spin1.max);

      exit( 1 );
    }

  if (  params.spin1.max < params.spin1.min)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "spin1-min (%f) must be < spin1-max (%f) \n",
               "spin1-min or spin1-max",params.spin1.min, params.spin1.max);

      exit( 1 );
    }

  if (  params.spin2.min < 0 || params.spin2.max < 0)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "spin2-min (%f) and spin2-max (%f) must be > 0 \n",
               "spin2-min or spin2-max",params.spin2.min, params.spin2.max);

      exit( 1 );
    }

  if (  params.spin2.max < params.spin2.min)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "spin2-min (%f) must be < spin2-max (%f) \n",
               "spin2-min or spin2-max",params.spin2.min, params.spin2.max);

      exit( 1 );
    }

  if (  params.inclination.max > LAL_PI)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "inclination-max (%f) must be < pi \n",
               "inclination-max", params.inclination.max);

      exit( 1 );
    }

  if (  params.inclination.min < 0.01)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "inclination-min (%f) must be > 0.01 \n",
               "inclination-min", params.inclination.min);

      exit( 1 );
    }

  if (  params.inclination.max < params.inclination.min)
    {
      fprintf( stderr, "invalid argument to --%s:\n"
               "inclination-min (%f) must be < inclination-max (%f) \n",
               "inclination-min or inclination-max",params.inclination.min, params.inclination.max);

      exit( 1 );
    }


}

