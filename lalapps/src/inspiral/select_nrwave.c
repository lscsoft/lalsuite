/*
 * Copyright (C) 2007 Badri Krishnan, Chad Hanna, Lucia Santamaria Lara, Robert Adam Mercer, Stephen Fairhurst
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Revision: $Id$
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include <lalapps.h>
#include <FrameL.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/NRWaveIO.h>
#include <lal/NRWaveInject.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Inject.h>
#include <lal/FileIO.h>
#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/LALDetectors.h>
#include <lal/LALFrameIO.h>
#include <lal/UserInput.h>
#include <lalappsfrutils.h>
#include <lal/FrameStream.h>

RCSID( "$Id$" );

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "nr_wave"

/* function prototypes */


#define TRUE (1==1)
#define FALSE (1==0)

/* verbose flag */
extern int vrbflg;


/* main program entry */
INT4 main( INT4 argc, CHAR *argv[] )
{
  LALStatus status = blank_status;

  FrCache *frGlobCache = NULL;
  FrCache *frInCache = NULL;
  FrCacheSieve  sieve;
  FrStream  *frStream = NULL;

  NRWaveCatalog nrCatalog;

  UINT4 k;

  /* user input variables */
  BOOLEAN  uvar_help = FALSE;
  CHAR *uvar_nrDir=NULL;
  INT4 uvar_minEta=0, uvar_maxEta=0;

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;

  lalDebugLevel = 0;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);


  LAL_CALL( LALRegisterBOOLUserVar( &status, "help", 'h', UVAR_HELP, "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "nrDir", 'D', UVAR_REQUIRED, "Directory with NR data", &uvar_nrDir), &status);
  LAL_CALL( LALRegisterINTUserVar( &status, "minEta", 0, UVAR_OPTIONAL, "Minimum symm. mass ratio", &uvar_minEta),  &status);
  LAL_CALL( LALRegisterINTUserVar( &status, "maxEta", 0, UVAR_OPTIONAL, "Minimum symm. mass ratio", &uvar_maxEta),  &status);


  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  
  /* create a frame cache by globbing *.gwf in specified dir */
  LAL_CALL( LALFrCacheGenerate( &status, &frGlobCache, uvar_nrDir, NULL ), 
	    &status );

  memset( &sieve, 0, sizeof(FrCacheSieve) );
  LAL_CALL( LALFrCacheSieve( &status, &frInCache, frGlobCache, &sieve ), 
	    &status );
  LAL_CALL( LALDestroyFrCache( &status, &frGlobCache ), &status );  
  
  /* check we globbed at least one frame file */
  if ( ! frInCache->numFrameFiles )  {
    fprintf( stderr, "error: no numrel frame files found\n");
    exit(1);
  }

  LAL_CALL( LALFrCacheOpen( &status, &frStream, frInCache ), &status );  

  /* loop over frame files and select the ones with nr-params in the right range */
  for (k = 0; k < frInCache->numFrameFiles; k++) {

    
  }

  LAL_CALL( LALFrClose( &status, &frStream ), &status );
  LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );  


  LAL_CALL (LALDestroyUserVars(&status), &status);
  LALCheckMemoryLeaks();

  exit(0);
}


