/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes  
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
 * \file
 * \ingroup pulsarApps
 * \author Badri Krishnan, Berit Behnke
 * \brief Little Helper code for getting velocities from v2 SFTs
*/


/* lalapps includes */
#include <lalapps.h>

/* globals, constants and defaults */



/* boolean global variables for controlling output */

#define EARTHEPHEMERIS "/home/badkri/lscsoft/share/lal/earth05-09.dat"
#define SUNEPHEMERIS "/home/badkri/lscsoft/share/lal/sun05-09.dat"


#define TRUE (1==1)
#define FALSE (1==0)


#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <errno.h> 

#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/DetectorStates.h>
#include <lal/LALDatatypes.h>
#include <lal/Velocity.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>

#include <lalapps.h>




/******************************************/

int main(int argc, char *argv[]){

  /* LALStatus pointer */
  static LALStatus  status;  
  
  /* time and velocity  */
  static LIGOTimeGPSVector    timeV;
  REAL8 timeBase;


  LALDetector *det;  
  DetectorStateSeries *detStates=NULL;

  /* ephemeris */
  EphemerisData    *edat=NULL;

  /* user input variables */
  BOOLEAN  uvar_help; 

  CHAR     *uvar_earthEphemeris=NULL;
  CHAR     *uvar_sunEphemeris=NULL;
  CHAR     *uvar_sftDir=NULL;

  /* Set up the default parameters */

  /* LAL error-handler */
  lal_errhandler = LAL_ERR_EXIT;
  
  uvar_help = FALSE;

  uvar_earthEphemeris = (CHAR *)LALCalloc( 512 , sizeof(CHAR));
  strcpy(uvar_earthEphemeris,EARTHEPHEMERIS);

  uvar_sunEphemeris = (CHAR *)LALCalloc( 512 , sizeof(CHAR));
  strcpy(uvar_sunEphemeris,SUNEPHEMERIS);

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar( &status, "help",             'h',  UVAR_HELP,     "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "earthEphemeris",  'E', UVAR_REQUIRED, "Earth Ephemeris file",  &uvar_earthEphemeris),  &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sunEphemeris",    'S', UVAR_REQUIRED, "Sun Ephemeris file", &uvar_sunEphemeris), &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "sftDir",          'D', UVAR_REQUIRED, "SFT filename pattern", &uvar_sftDir), &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /*  get ephemeris  */
  edat = (EphemerisData *)LALCalloc(1, sizeof(EphemerisData));
  (*edat).ephiles.earthEphemeris = uvar_earthEphemeris;
  (*edat).ephiles.sunEphemeris = uvar_sunEphemeris;

  LAL_CALL( LALInitBarycenter( &status, edat), &status);
  
  /* read sft Files and set up weights and nstar vector */
  {
    /* new SFT I/O data types */
    SFTCatalog *catalog = NULL;
    static SFTConstraints constraints;
    UINT4 k;

    /* set detector constraint */
    constraints.detector = NULL;
    
    /* get sft catalog */
    LAL_CALL( LALSFTdataFind( &status, &catalog, uvar_sftDir, &constraints), &status);
    if ( (catalog == NULL) || (catalog->length == 0) ) {
      fprintf (stderr,"Unable to match any SFTs with pattern '%s'\n", uvar_sftDir );
      exit(1);
    }

    for (k = 0; k < catalog->length; k++) {

      det = XLALGetSiteInfo ( catalog->data[k].header.name);
      timeBase =  1./catalog->data[k].header.deltaF;

      timeV.length = 1;
      timeV.deltaT = timeBase;
      timeV.data = &(catalog->data[k].header.epoch);

      LAL_CALL( LALGetDetectorStates ( &status, &detStates, &timeV, det, edat, 0.5*timeBase), &status);


      fprintf(stdout, "%g  %g   %g\n", detStates->data[0].vDetector[0], detStates->data[0].vDetector[1],
	      detStates->data[0].vDetector[2]);
      
      LALDestroyDetectorStateSeries (&status, &detStates );
      LALFree (det);

    }   /* end loop over sfts */

    LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);  	

  } /* end of sft reading block */

  LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);

  LAL_CALL (LALDestroyUserVars(&status), &status);
  
  LALCheckMemoryLeaks();

  if ( lalDebugLevel )
    REPORTSTATUS ( &status);
  
  return status.statusCode;

}



  







