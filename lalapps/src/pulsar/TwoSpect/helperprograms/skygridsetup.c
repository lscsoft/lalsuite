/*
 *  Copyright (C) 2011 Evan Goetz
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


#include <stdio.h>
#include <math.h>

#include <lal/DopplerScan.h>

#include "cmdline_skygridsetup.h"
#include "../antenna.h"

//Main program
int main(int argc, char *argv[])
{
   
   const CHAR *fn = __func__;
   
   CHAR *earth_ephemeris = NULL, *sun_ephemeris = NULL;
   FILE *OUTPUT;
   char s[20000];
   LALDetector det;
   LALStatus status;          //LALStatus structure
   status.statusPtr = NULL;   //Set statuspointer to NULL
   
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();
   configparams->initialize = 0;
   if ( cmdline_parser(argc, argv, &args_info) ) {
      fprintf(stderr, "%s: cmdline_parser() failed.\n", fn);
      return -1;
   }
   if ( args_info.config_given ) {
      if ( cmdline_parser_config_file(args_info.config_arg, &args_info, configparams) ) {
         fprintf(stderr, "%s: cmdline_parser_config_file() failed.\n", fn);
         return -1;
      }
   }
   
   snprintf(s, 20000, "%s", args_info.outfilename_arg);
   OUTPUT = fopen(s,"w");
   if (OUTPUT==NULL) {
      fprintf(stderr, "%s: Output file could not be opened.\n", fn);
      return -1;
   }
   
   //Allocate memory for files
   earth_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*earth_ephemeris));
   sun_ephemeris = XLALCalloc(strlen(args_info.ephemDir_arg)+25, sizeof(*sun_ephemeris));
   if (earth_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*earth_ephemeris));
      return -1;
   } else if (sun_ephemeris==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sun_ephemeris));
      return -1;
   }
   sprintf(earth_ephemeris, "%s/earth%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);
   sprintf(sun_ephemeris, "%s/sun%s.dat", args_info.ephemDir_arg, args_info.ephemYear_arg);
   
   //Interferometer
   CHAR *IFO = XLALCalloc(strlen(args_info.IFO_arg)+1, sizeof(*IFO));
   if (IFO==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*IFO));
      XLAL_ERROR(XLAL_ENOMEM);
   }
   sprintf(IFO, "%s", args_info.IFO_arg);
   if (strcmp("L1", IFO)==0) {
      fprintf(stderr,"IFO = %s\n", IFO);
      det = lalCachedDetectors[LALDetectorIndexLLODIFF]; //L1
   } else if (strcmp("H1", IFO)==0) {
      fprintf(stderr,"IFO = %s\n", IFO);
      det = lalCachedDetectors[LALDetectorIndexLHODIFF]; //H1
   } else {
      fprintf(stderr, "%s: Not using valid interferometer! Expected 'H1' or 'L1' not %s.\n", fn, IFO);
      return -1;
   }
   XLALFree((CHAR*)IFO);
   
   //Parameters for the sky-grid
   CHAR *sky = XLALCalloc(strlen(args_info.skyRegion_arg)+1, sizeof(*sky));
   if (sky==NULL) {
      fprintf(stderr, "%s: XLALCalloc(%zu) failed.\n", fn, sizeof(*sky));
      return -1;
   }
   sprintf(sky, "%s", args_info.skyRegion_arg);
   fprintf(stderr, "Sky region = %s\n", sky);
   DopplerSkyScanInit scanInit = empty_DopplerSkyScanInit;
   DopplerSkyScanState scan = empty_DopplerSkyScanState;
   PulsarDopplerParams dopplerpos;
   scanInit.gridType = GRID_ISOTROPIC;     //Default value for an approximate-isotropic grid
   scanInit.skyRegionString = sky;      //"allsky" = Default value for all-sky search
   scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
   scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
   
   //Initialize ephemeris data structure
   EphemerisData *edat = XLALInitBarycenter(earth_ephemeris, sun_ephemeris);
   if (edat==NULL) {
      fprintf(stderr, "%s: XLALInitBarycenter() failed.\n", fn);
      return -1;
   }
   
   //Maximum orbital earth speed in units of c from start of S6 TwoSpect data for 104 weeks total time
   REAL4 detectorVmax = CompDetectorVmax(931081500.0+args_info.SFToverlap_arg, args_info.Tcoh_arg, args_info.SFToverlap_arg, 62899200.0-args_info.SFToverlap_arg, det, edat);
   if (xlalErrno!=0) {
      fprintf(stderr, "%s: CompDetectorVmax() failed.\n", fn);
      return -1;
   }
   
   //Initialize the sky-grid
   scanInit.dAlpha = 0.5/((args_info.fmin_arg+0.5*args_info.fspan_arg) * args_info.Tcoh_arg * detectorVmax);
   scanInit.dDelta = scanInit.dAlpha;
   InitDopplerSkyScan(&status, &scan, &scanInit);
   if (status.statusCode!=0) {
      fprintf(stderr, "%s: InitDopplerSkyScan() failed.\n", fn);
      return -1;
   }
   
   //Start at first location
   if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
      fprintf(stderr, "%s: XLALNextDopplerSkyPos() failed.\n", fn);
      return -1;
   }
   
   //loop through and output to the specified file
   while (scan.state != STATE_FINISHED) {
      fprintf(OUTPUT, "%.6f %.6f\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      //Iterate to next sky location
      if ((XLALNextDopplerSkyPos(&dopplerpos, &scan))!=0) {
         fprintf(stderr,"%s: XLALNextDopplerSkyPos() failed.\n", fn);
         return -1;
      }
   }
   
   
   //Destroy
   cmdline_parser_free(&args_info);
   XLALFree((CHAR*)earth_ephemeris);
   XLALFree((CHAR*)sun_ephemeris);
   XLALFree((CHAR*)sky);
   XLALDestroyEphemerisData(edat);
   fclose(OUTPUT);
   
   return 0;
   
}

