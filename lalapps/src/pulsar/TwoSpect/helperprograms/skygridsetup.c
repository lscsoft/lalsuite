/*
 *  Copyright (C) 2011, 2013 Evan Goetz
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
   
   FILE *OUTPUT;
   char s[20000];
   LALDetector det;
   LALStatus XLAL_INIT_DECL(status);
   
   struct gengetopt_args_info args_info;
   struct cmdline_parser_params *configparams;
   configparams = cmdline_parser_params_create();
   configparams->check_required = 0;  //don't check for required values at the step
   cmdline_parser_ext(argc, argv, &args_info, configparams);
   configparams->initialize = 0;  //don't reinitialize the parameters structure
   if (args_info.config_given) cmdline_parser_config_file(args_info.config_arg, &args_info, configparams);
   cmdline_parser_required(&args_info, argv[0]);

   if (args_info.v2_given && !args_info.Tobs_given) {
      fprintf(stderr, "%s: Tobs must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   if (args_info.v2_given && !args_info.t0_given) {
      fprintf(stderr, "%s: Tobs must be specified.\n", __func__);
      XLAL_ERROR(XLAL_FAILURE);
   }
   
   snprintf(s, 20000, "%s", args_info.outfilename_arg);
   XLAL_CHECK( (OUTPUT = fopen(s,"w")) != NULL, XLAL_EIO, "Output file %s could not be opened\n", s );
   
   //Interferometer
   if (strcmp("L1", args_info.IFO_arg)==0) {
      fprintf(stderr,"IFO = %s\n", args_info.IFO_arg);
      det = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; //L1
   } else if (strcmp("H1", args_info.IFO_arg)==0) {
      fprintf(stderr,"IFO = %s\n", args_info.IFO_arg);
      det = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; //H1
   } else if (strcmp("V1", args_info.IFO_arg)==0) {
      fprintf(stderr,"IFO = %s\n", args_info.IFO_arg);
      det = lalCachedDetectors[LAL_VIRGO_DETECTOR]; //V1
   } else {
      fprintf(stderr, "%s: Not using valid interferometer! Expected 'H1', 'L1', or 'V1' not %s.\n", __func__, args_info.IFO_arg);
      XLAL_ERROR( XLAL_EFUNC );
   }
   
   //Parameters for the sky-grid
   fprintf(stderr, "Sky region = %s\n", args_info.skyRegion_arg);
   DopplerSkyScanInit XLAL_INIT_DECL(scanInit);
   DopplerSkyScanState XLAL_INIT_DECL(scan);
   PulsarDopplerParams dopplerpos;
   scanInit.gridType = GRID_ISOTROPIC;     //Default value for an approximate-isotropic grid
   scanInit.skyRegionString = args_info.skyRegion_arg;      //"allsky" = Default value for all-sky search
   scanInit.numSkyPartitions = 1;   //Default value so sky is not broken into chunks
   scanInit.Freq = args_info.fmin_arg+0.5*args_info.fspan_arg;  //Mid-point of the frequency band
   
   //Initialize ephemeris data structure
   EphemerisData *edat = NULL;
   XLAL_CHECK( (edat = XLALInitBarycenter(args_info.ephemEarth_arg, args_info.ephemSun_arg)) != NULL, XLAL_EFUNC );
   
   //Maximum orbital earth speed in units of c from start of S6 TwoSpect data for 104 weeks total time
   REAL4 detectorVmax = 0.0;
   if (!args_info.v2_given) detectorVmax = CompDetectorVmax(931081500.0+args_info.SFToverlap_arg, args_info.Tcoh_arg, args_info.SFToverlap_arg, 62899200.0-args_info.SFToverlap_arg, det, edat);
   else detectorVmax = CompDetectorVmax(args_info.t0_arg-args_info.Tcoh_arg, args_info.Tcoh_arg, args_info.SFToverlap_arg, args_info.Tobs_arg+2.0*args_info.Tcoh_arg, det, edat);
   XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC, "CompDetectorVmax() failed\n" );
   
   //Initialize the sky-grid
   scanInit.dAlpha = 0.5/((args_info.fmin_arg+0.5*args_info.fspan_arg) * args_info.Tcoh_arg * detectorVmax);
   scanInit.dDelta = scanInit.dAlpha;
   InitDopplerSkyScan(&status, &scan, &scanInit);
   XLAL_CHECK( status.statusCode == 0, XLAL_EFUNC );
   
   //Start at first location
   XLAL_CHECK( XLALNextDopplerSkyPos(&dopplerpos, &scan) == XLAL_SUCCESS, XLAL_EFUNC );
   
   //loop through and output to the specified file
   while (scan.state != STATE_FINISHED) {
      fprintf(OUTPUT, "%.6f %.6f\n", dopplerpos.Alpha, dopplerpos.Delta);
      
      //Iterate to next sky location
      XLAL_CHECK( XLALNextDopplerSkyPos(&dopplerpos, &scan) == XLAL_SUCCESS, XLAL_EFUNC );
   }

   //Destroy
   cmdline_parser_free(&args_info);
   XLALDestroyEphemerisData(edat);
   fclose(OUTPUT);
   
   return 0;
   
}

