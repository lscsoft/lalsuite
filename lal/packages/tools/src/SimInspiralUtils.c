/*----------------------------------------------------------------------- 
 * 
 * File Name: SimInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SimInspiralUtilsCV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 
#endif

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

NRCSID( SIMINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SimInspiralUtils.c}}

Provides a set of utilities for manipulating \texttt{simInspiralTable}s.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SimInspiralUtilsCP}
\idx{LALGalacticInspiralParamsToSimInspiralTable()}
\idx{LALInspiralSiteTimeAndDist()}
\idx{LALPopulateSimInspiralSiteInfo()}
  
\subsubsection*{Description}

The function \texttt{LALInspiralSiteTimeAndDist()} calculates detector end
time (\texttt{endTime}) and effective distance (\texttt{effDist}) for an
inspiral signal from a specific location in the sky (\texttt{skyPos}) assumed
to be given in equatorial coordinates.  The detector end time is obtained by
using \texttt{LALTimeDelayFromEarthCenter()}, while the effective distance
requires calculation of the detector response, calculated using
\texttt{LALComputeDetAMResponse()}. 

The function \texttt{LALPopulateSimInspiralSiteInfo()} populates the end time
and effective distance for each of the interferometer sites.  The sky location
(in equatorial coordinates) is assumed to be already contained in the input
\texttt{SimInspiralTable}.  The end time and effective distance for each site
is calculated by calling \texttt{LALInspiralSiteTimeAndDist()} once for each
of the detectors, and setting the \texttt{detector} appropriately.
    
    

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent LALGetInspiralParams, LALGPStoGMST1, LALTimeDelayFromEarthCenter, 
  LALAddFloatToGPS, LALComputeDetAMResponse.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SimInspiralUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
LALGalacticInspiralParamsToSimInspiralTable(
    LALStatus                  *status,
    SimInspiralTable           *output,
    GalacticInspiralParamStruc *input,
    RandomParams               *params
    )
/* </lalVerbatim> */
{
  PPNParamStruc         ppnParams;
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };
  LALGPSandAcc          gpsAndAcc;
  SkyPosition           skyPos;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8     time_diff_ns;
  REAL4                 splus, scross, cosiota;

  INITSTATUS( status, "LALGalacticParamsToSimInspiral", SIMINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( params, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  /*
   *
   * compute sky position and inspiral params
   *
   */


  /* generate the ppn inspiral params */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
  LALGetInspiralParams( status->statusPtr, &ppnParams, input, params );
  CHECKSTATUSPTR( status );

  if ( ppnParams.position.system != COORDINATESYSTEM_EQUATORIAL )
  {
    ABORT( status, LIGOMETADATAUTILSH_ECOOR, LIGOMETADATAUTILSH_MSGECOOR );
  }

  /* copy the inspiral data into sim_inspiral table */
  output->mass1        = input->m1;
  output->mass2        = input->m2;
  output->eta          = ppnParams.eta;
  output->distance     = ppnParams.d / (1.0e6 * LAL_PC_SI); /* Mpc */
  output->longitude    = ppnParams.position.longitude;
  output->latitude     = ppnParams.position.latitude;
  output->inclination  = ppnParams.inc;
  output->coa_phase    = ppnParams.phi;
  output->polarization = ppnParams.psi;

  /* populate geocentric end time */
  output->geocent_end_time = input->geocentEndTime;

  /* populate gmst field */
  LALGPStoGMST1( status->statusPtr, &(output->end_time_gmst), 
      &(output->geocent_end_time), &gmstUnits );
  CHECKSTATUSPTR( status );

  /* set up params for the site end times and detector response */
  memset( &skyPos, 0, sizeof(SkyPosition) );
  memset( &source, 0, sizeof(LALSource) );
  memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
  memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );

  skyPos.longitude = output->longitude;
  skyPos.latitude  = output->latitude;
  skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

  source.equatorialCoords = skyPos;
  source.orientation      = output->polarization;

  placeAndGPS.p_gps = &(output->geocent_end_time);

  detTimeAndSource.p_det_and_time = &placeAndGPS;
  detTimeAndSource.p_source = &skyPos;

  detAndSource.pSource = &source;

  gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
  gpsAndAcc.gps = output->geocent_end_time;


  /*
   *
   * compute site end times
   *
   */


  /* initialize end times with geocentric value */
  output->h_end_time = output->l_end_time = input->geocentEndTime;

  /* ligo hanford observatory */
  placeAndGPS.p_detector = &lho;
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, &(output->h_end_time),
      &(output->h_end_time), time_diff_ns );
  CHECKSTATUSPTR( status );

  /* ligo livingston observatory */
  placeAndGPS.p_detector = &llo;
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, &(output->l_end_time), 
      &(output->l_end_time), time_diff_ns );
  CHECKSTATUSPTR( status );


  /*
   *
   * compute the effective distance of the inspiral
   *
   */


  /* initialize distances with real distance and compute splus and scross */
  output->eff_dist_h = output->eff_dist_l = 2.0 * output->distance;
  cosiota = cos( output->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;

  /* compute the response of the LHO detectors */
  detAndSource.pDetector = &lho;
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance for LHO */
  output->eff_dist_h /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );

  /* compute the response of the LLO detector */
  detAndSource.pDetector = &llo;
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance for LLO */
  output->eff_dist_l /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );


  /* 
   *
   * normal exit
   *
   */


  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
LALInspiralSiteTimeAndDist( 
    LALStatus         *status,
    SimInspiralTable  *output,
    LALDetector       *detector,
    LIGOTimeGPS       *endTime,
    REAL4             *effDist,
    SkyPosition       *skyPos
    )
/* </lalVerbatim> */
{
  LALGPSandAcc          gpsAndAcc;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8                 time_diff_ns;
  REAL4                 splus, scross, cosiota;

  INITSTATUS( status, "LALInspiralSiteTimeAndDist", SIMINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that the arguments are not null */
  ASSERT( output, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( detector, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( endTime, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( effDist, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( skyPos, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  memset( &source, 0, sizeof(LALSource) );
  memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
  memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
  memset( &detAndSource, 0, sizeof(LALDetAndSource) );


  source.equatorialCoords = *skyPos;
  source.orientation      = output->polarization;

  placeAndGPS.p_gps = &(output->geocent_end_time);

  detTimeAndSource.p_det_and_time = &placeAndGPS;
  detTimeAndSource.p_source = skyPos;
  detTimeAndSource.p_det_and_time->p_detector = detector;

  detAndSource.pSource = &source;
  detAndSource.pDetector = detector;

  gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
  gpsAndAcc.gps = output->geocent_end_time;

  /* initialize end time with geocentric value */
  *endTime = output->geocent_end_time;

  /* calculate the detector end time */
  LALTimeDelayFromEarthCenter( status->statusPtr, &time_diff_ns, 
      &detTimeAndSource );
  CHECKSTATUSPTR( status );
  LALAddFloatToGPS( status->statusPtr, endTime,
      endTime, time_diff_ns );
  CHECKSTATUSPTR( status );

  /* initialize distance with real distance and compute splus and scross */
  *effDist = 2.0 * output->distance;
  cosiota = cos( output->inclination );
  splus = -( 1.0 + cosiota * cosiota );
  scross = -2.0 * cosiota;

  /* compute the response of the detector */
  LALComputeDetAMResponse( status->statusPtr, &resp, &detAndSource,
      &gpsAndAcc );
  CHECKSTATUSPTR( status );

  /* compute the effective distance */
  *effDist /= sqrt( 
      splus*splus*resp.plus*resp.plus + scross*scross*resp.cross*resp.cross );

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
LALPopulateSimInspiralSiteInfo(
    LALStatus                  *status,
    SimInspiralTable           *output
    )
/* </lalVerbatim> */
{
  SkyPosition           skyPos;
  LALDetector           detector; 
  REAL4                *eff_dist;
  LIGOTimeGPS          *end_time;


  INITSTATUS( status, "LALPopulateSimInspiralSiteInfo", SIMINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* set up params for the geocent end time and source location */
  memset( &skyPos, 0, sizeof(SkyPosition) );

  skyPos.longitude = output->longitude;
  skyPos.latitude  = output->latitude;
  skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

  /* LIGO Hanford observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  end_time = &(output->h_end_time);
  eff_dist = &(output->eff_dist_h);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time, 
      eff_dist, &skyPos);

  /* LIGO Livingston observatory*/
  detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  end_time = &(output->l_end_time);
  eff_dist = &(output->eff_dist_l);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time, 
      eff_dist, &skyPos);

  /* GEO observatory*/
  detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
  end_time = &(output->g_end_time);
  eff_dist = &(output->eff_dist_g);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time, 
      eff_dist, &skyPos);

  /* TAMA observatory*/
  detector = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
  end_time = &(output->t_end_time);
  eff_dist = &(output->eff_dist_t);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time, 
      eff_dist, &skyPos);

  /* Virgo observatory*/
  detector = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
  end_time = &(output->v_end_time);
  eff_dist = &(output->eff_dist_v);
  LALInspiralSiteTimeAndDist(status->statusPtr, output, &detector, end_time, 
      eff_dist, &skyPos);

  /* 
   *
   * normal exit
   *
   */


  DETATCHSTATUSPTR (status);
  RETURN (status);
}
