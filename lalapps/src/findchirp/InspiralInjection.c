 /* *****************************************
  * Author: Thomas Cokelaer.
  * [PURPOSE] Inject inspiral waveform test
  * [USAGE]   Just type the executable name
  * [INPUT]   An xml file with injection data
  * [OUTPUT]  Write a file "injection.dat"
  * 
  * Documentation in lalapps.pdf in ./doc
  * **************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <FrameL.h>

#include <lalapps.h>
#include <series.h>
#include <processtable.h>
#include <lalappsfrutils.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>

#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>

#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/FindChirp.h>
#include <lal/PrintFTSeries.h>


RCSID(  "$Id$");


/* Temporary version ; haven't taken care of units, time and so on 
 * for the time being.  */
int
main()
{
  
 
  /* input xml data */
  CHAR  	*injectionFile	= NULL;         

  /* some variables to create the data */
  UINT4 	startTime;
  UINT4		endTime;
  UINT4 	k;
  UINT4 	numPoints 	= 524288 * 2 / 16;

  REAL8 	sampling	= 2048.;

  LALStatus 	status 		= blank_status;

  /* injection structure */
  SimInspiralTable *injections 	= NULL;

  /* the data */
  REAL4TimeSeries               ts;						/* A time series to store the injection */
  COMPLEX8FrequencySeries       fs;						/* A freq series to store the psd 	*/

  /* output data to check injection */
  FILE		*output;



  /* --- MAIN --- */
  lalDebugLevel = 1;

  injectionFile = "injection.xml"; 						/* an xml file with the injection	*/
  output 	= fopen("injection.dat","w");					/* an ascii file for results		*/
  
  memset(&ts, 0, sizeof(REAL4TimeSeries));					/* allocate memory 			*/
  LAL_CALL( LALSCreateVector( &status, &(ts.data), numPoints), &status);	/* and ts null				*/
  memset( &fs, 0, sizeof(COMPLEX8FrequencySeries));				/* idem for fs				*/
  LAL_CALL( LALCCreateVector( &status, &(fs.data), numPoints / 2 + 1 ), 
    &status );

  ts.epoch.gpsSeconds 	= 729273610;						/* gps time of the time series		*/
  startTime 		= 729273610;						/* gps start time and end time of ..	*/	
  endTime   		= startTime + 100;					/* ..injection; should be in agreement..*/
  										/* ..with the xml file			*/ 
  ts.sampleUnits 	= lalADCCountUnit;					/*  UNITY ?? 				*/
  ts.deltaT 		= 1./sampling;						/* sampling				*/
  fs.deltaF 		= sampling / numPoints;					/* idem for fs				*/
  fs.sampleUnits 	= lalADCCountUnit;

  /* --- the psd is flat for simplicity --- */		
  for( k = 0 ; k< numPoints/2; k++){
      fs.data->data[k].re = 1.;
      fs.data->data[k].im = 0.;
    }
  /* --- and time series is zero --- */
  for( k = 0 ; k< numPoints; k++){
      ts.data->data[k] = 0.;
    }

  /* --- injection is here --- */
  SimInspiralTableFromLIGOLw( &injections, 
			      injectionFile,
			      startTime,
			      endTime);

  LAL_CALL( LALFindChirpInjectSignals( &status, &ts, injections, &fs ), 
	&status );

  /* --- let's print the results in the file ---*/
  for (k=0; k<numPoints; k++){
      fprintf(output,"%i %e\n", k, ts.data->data[k]);
    }
  fclose(output);

  return 0;
}


/* An example of xml file to be used with this file */




#if 0
<?xml version='1.0' encoding='utf-8' ?>
<!DOCTYPE LIGO_LW [
<!ELEMENT LIGO_LW ((LIGO_LW|Comment|Param|Table|Array|Stream|IGWDFrame|AdcData|AdcInterval|Time|Detector)*)>
<!ATTLIST LIGO_LW
         Name CDATA #IMPLIED
         Type CDATA #IMPLIED>

<!ELEMENT Comment (#PCDATA)>

<!ELEMENT Param (#PCDATA|Comment)*>
<!ATTLIST Param 
         Name CDATA #IMPLIED
         Type CDATA #IMPLIED
         Start CDATA #IMPLIED
         Scale CDATA #IMPLIED
         Unit CDATA #IMPLIED
         DataUnit CDATA #IMPLIED>

<!ELEMENT Table (Comment?,Column*,Stream?)>
<!ATTLIST Table 
         Name CDATA #IMPLIED
         Type CDATA #IMPLIED>

<!ELEMENT Column EMPTY>
<!ATTLIST Column
         Name CDATA #IMPLIED
         Type CDATA #IMPLIED
         Unit CDATA #IMPLIED>

<!ELEMENT Array (Dim*,Stream?)>
<!ATTLIST Array 
         Name CDATA #IMPLIED
         Type CDATA #IMPLIED
         Unit CDATA #IMPLIED>

<!ELEMENT Dim (#PCDATA)>
<!ATTLIST Dim 
         Name  CDATA #IMPLIED
         Unit CDATA #IMPLIED
         Start CDATA #IMPLIED
         Scale CDATA #IMPLIED>

<!ELEMENT Stream (#PCDATA)>
<!ATTLIST Stream 
         Name      CDATA #IMPLIED
         Type      (Remote|Local) "Local"
         Delimiter CDATA ","
         Encoding  CDATA #IMPLIED
         Content   CDATA #IMPLIED>

<!ELEMENT IGWDFrame ((Comment|Param|Time|Detector|AdcData|LIGO_LW|Stream?|Array|IGWDFrame)*)>
<!ATTLIST IGWDFrame 
         Name CDATA #IMPLIED>

<!ELEMENT Detector ((Comment|Param|LIGO_LW)*)>
<!ATTLIST Detector 
         Name CDATA #IMPLIED>

<!ELEMENT AdcData ((AdcData|Comment|Param|Time|LIGO_LW|Array)*)>
<!ATTLIST AdcData 
         Name CDATA #IMPLIED>

<!ELEMENT AdcInterval ((AdcData|Comment|Time)*)>
<!ATTLIST AdcInterval 
         Name CDATA #IMPLIED
         StartTime CDATA #IMPLIED
         DeltaT CDATA #IMPLIED>

<!ELEMENT Time (#PCDATA)>
<!ATTLIST Time 
         Name CDATA #IMPLIED
         Type (GPS|Unix|ISO-8601) "ISO-8601">
]>

<LIGO_LW>
   <Comment>metadata</Comment>
   <Table Name="processgroup:process:table">
      <Column Name="processgroup:process:program" Type="lstring"/>
      <Column Name="processgroup:process:version" Type="lstring"/>
      <Column Name="processgroup:process:cvs_repository" Type="lstring"/>
      <Column Name="processgroup:process:cvs_entry_time" Type="int_4s"/>
      <Column Name="processgroup:process:comment" Type="lstring"/>
      <Column Name="processgroup:process:is_online" Type="int_4s"/>
      <Column Name="processgroup:process:node" Type="lstring"/>
      <Column Name="processgroup:process:username" Type="lstring"/>
      <Column Name="processgroup:process:unix_procid" Type="int_4s"/>
      <Column Name="processgroup:process:start_time" Type="int_4s"/>
      <Column Name="processgroup:process:end_time" Type="int_4s"/>
      <Column Name="processgroup:process:jobid" Type="int_4s"/>
      <Column Name="processgroup:process:domain" Type="lstring"/>
      <Column Name="processgroup:process:ifos" Type="lstring"/>
      <Column Name="processgroup:process:process_id" Type="ilwd:char"/>
      <Stream Name="processgroup:process:table" Type="Local" Delimiter=",">
         "inspinj","1.19","/usr/local/cvs/lal/lalapps/src/inspiral/inspinj.c\,v",759519007," ",0,"cokelaerlt.astro.cf.ac.uk","cokelaer",14350,763339866,763339867,0,"lalapps","","process:process_id:0"
      </Stream>
   </Table>
   <Table Name="process_paramsgroup:process_params:table">
      <Column Name="process_paramsgroup:process_params:program" Type="lstring"/>
      <Column Name="process_paramsgroup:process_params:process_id" Type="ilwd:char"/>
      <Column Name="process_paramsgroup:process_params:param" Type="lstring"/>
      <Column Name="process_paramsgroup:process_params:type" Type="lstring"/>
      <Column Name="process_paramsgroup:process_params:value" Type="lstring"/>
      <Stream Name="process_paramsgroup:process_params:table" Type="Local" Delimiter=",">
         "inspinj","process:process_id:0","--mass-file","string","BBHMasses.dat",
         "inspinj","process:process_id:0","--source-file","string","inspsrcs.dat"
      </Stream>
   </Table>
   <Table Name="sim_inspiralgroup:sim_inspiral:table">
      <Column Name="sim_inspiralgroup:sim_inspiral:process_id" Type="ilwd:char"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:waveform" Type="lstring"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:geocent_end_time" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:geocent_end_time_ns" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:h_end_time" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:h_end_time_ns" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:l_end_time" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:l_end_time_ns" Type="int_4s"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:end_time_gmst" Type="real_8"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:source" Type="lstring"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:mass1" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:mass2" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:eta" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:distance" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:longitude" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:latitude" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:inclination" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:coa_phase" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:polarization" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:eff_dist_h" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:eff_dist_l" Type="real_4"/>
      <Column Name="sim_inspiralgroup:sim_inspiral:simulation_id" Type="ilwd:char"/>
      <Stream Name="sim_inspiralgroup:sim_inspiral:table" Type="Local" Delimiter=",">
         "process:process_id:0","TaylorT2twoPN",729273613,0,729273612,998227760,729273612,994513400,1.6175260586695899e+00,"M33",1.000000e+01,1.000000e+01,2.500000e-01,8.400000e-01,4.097160e-01,5.349434e-01,1.791432e-01,2.040999e+00,8.858610e-01,1.717248e+00,1.666386e+00,"sim_inspiral:simulation_id:0",
      </Stream>
   </Table>
</LIGO_LW>
#endif
