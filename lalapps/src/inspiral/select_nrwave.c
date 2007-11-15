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
#include <lal/LIGOLwXML.h>
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

typedef struct {
  REAL8 massRatioMin; 
  REAL8 massRatioMax;
  REAL8 sx1Min; 
  REAL8 sx1Max;
  REAL8 sx2Min; 
  REAL8 sx2Max;
  REAL8 sy1Min; 
  REAL8 sy1Max;
  REAL8 sy2Min; 
  REAL8 sy2Max;
  REAL8 sz1Min; 
  REAL8 sz1Max;
  REAL8 sz2Min; 
  REAL8 sz2Max; 
} NrParRange;


void GetNRMetaDataFromFrameHistory(NRWaveMetaData *data,
				   FrHistory *history);

void ParseThisComment(NRWaveMetaData *data, CHAR *comment);

REAL8 GetNumericValue(CHAR *comment);

REAL8 MetaDataInRange(NRWaveMetaData *data, NrParRange *range);


/* main program entry */
INT4 main( INT4 argc, CHAR *argv[] )
{
  LALStatus status = blank_status;

  /* frame file stuff */
  FrCache *frGlobCache = NULL;
  FrCache *frInCache = NULL;
  FrCacheSieve  sieve;
  FrameH *frame=NULL;
  FrFile *frFile=NULL;

  /* inspiral table stuff */
  SimInspiralTable  *this_inj = NULL;
  LIGOLwXMLStream   xmlfp;
  MetadataTable  injections;

  /* nrwave stuff */
  NRWaveMetaData metaData;
  NrParRange range;

  UINT4 k;

  /* user input variables */
  BOOLEAN  uvar_help = FALSE;
  CHAR *uvar_nrDir=NULL;
  CHAR *uvar_nrGroup=NULL;
  REAL8 uvar_minMassRatio=1, uvar_maxMassRatio=0;
  REAL8 uvar_minSx1=-1, uvar_minSx2=-1, uvar_maxSx1=1, uvar_maxSx2=1;
  REAL8 uvar_minSy1=-1, uvar_minSy2=-1, uvar_maxSy1=1, uvar_maxSy2=1;
  REAL8 uvar_minSz1=-1, uvar_minSz2=-1, uvar_maxSz1=1, uvar_maxSz2=1;

  /* default debug level */
  lal_errhandler = LAL_ERR_EXIT;

  lalDebugLevel = 0;
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  LAL_CALL( LALRegisterBOOLUserVar( &status, "help", 'h', UVAR_HELP, "Print this message", &uvar_help), &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "nrDir", 'D', UVAR_REQUIRED, "Directory with NR data", &uvar_nrDir), &status);

  LAL_CALL( LALRegisterREALUserVar( &status, "minMassRatio", 0, UVAR_OPTIONAL, "Minimum mass ratio", &uvar_minMassRatio),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxMassRatio", 0, UVAR_OPTIONAL, "Maximum mass ratio", &uvar_maxMassRatio),  &status);

  LAL_CALL( LALRegisterREALUserVar( &status, "minSx1", 0, UVAR_OPTIONAL, "Minimum x-spin of first BH", &uvar_minSx1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "minSx2", 0, UVAR_OPTIONAL, "Minimum x-Spin of second BH", &uvar_minSx2),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSx1", 0, UVAR_OPTIONAL, "Maximum x-spin of first BH", &uvar_maxSx1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSx2", 0, UVAR_OPTIONAL, "Maximum x-spin of second BH", &uvar_maxSx2),  &status);

  LAL_CALL( LALRegisterREALUserVar( &status, "minSy1", 0, UVAR_OPTIONAL, "Minimum y-spin of first BH", &uvar_minSy1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "minSy2", 0, UVAR_OPTIONAL, "Minimum y-Spin of second BH", &uvar_minSy2),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSy1", 0, UVAR_OPTIONAL, "Maximum y-spin of first BH", &uvar_maxSy1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSy2", 0, UVAR_OPTIONAL, "Maximum y-spin of second BH", &uvar_maxSy2),  &status);

  LAL_CALL( LALRegisterREALUserVar( &status, "minSz1", 0, UVAR_OPTIONAL, "Minimum z-spin of first BH", &uvar_minSz1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "minSz2", 0, UVAR_OPTIONAL, "Minimum z-Spin of second BH", &uvar_minSz2),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSz1", 0, UVAR_OPTIONAL, "Maximum z-spin of first BH", &uvar_maxSz1),  &status);
  LAL_CALL( LALRegisterREALUserVar( &status, "maxSz2", 0, UVAR_OPTIONAL, "Maximum z-spin of second BH", &uvar_maxSz2),  &status);

  LAL_CALL( LALRegisterSTRINGUserVar( &status, "nrGroup", 0, UVAR_OPTIONAL, "NR group", &uvar_nrGroup), &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 


  range.massRatioMin = uvar_minMassRatio;
  range.massRatioMax = uvar_maxMassRatio;

  range.sx1Min = uvar_minSx1;
  range.sx1Max = uvar_maxSx1;

  range.sx2Min = uvar_minSx2;
  range.sx2Max = uvar_maxSx2;

  range.sy1Min = uvar_minSy1;
  range.sy1Max = uvar_maxSy1;

  range.sy2Min = uvar_minSy2;
  range.sy2Max = uvar_maxSy2;

  range.sz1Min = uvar_minSz1;
  range.sz1Max = uvar_maxSz1;

  range.sz2Min = uvar_minSz2;
  range.sz2Max = uvar_maxSz2;

  
  /* create a frame cache by globbing *.gwf in specified dir */
  LAL_CALL( LALFrCacheGenerate( &status, &frGlobCache, uvar_nrDir, NULL ), 
	    &status );


  memset( &sieve, 0, sizeof(FrCacheSieve) );
  /* sieve doesn't actually do anything yet */
  LAL_CALL( LALFrCacheSieve( &status, &frInCache, frGlobCache, &sieve ), 
	    &status );

  LAL_CALL( LALDestroyFrCache( &status, &frGlobCache ), &status );  
  
  /* check we globbed at least one frame file */
  if ( !frInCache->numFrameFiles )  {
    fprintf( stderr, "error: no numrel frame files found\n");
    exit(1);
  }

  /* initialize head of inspiral table linked list to null */
  injections.simInspiralTable = NULL;

  /* loop over frame files and select the ones with nr-params in the right range */
  for (k = 0; k < frInCache->numFrameFiles; k++) {

    frFile =  XLALFrOpenURL( frInCache->frameFiles[k].url);

    frame = FrameRead (frFile);
    
    memset(&metaData, 0, sizeof(NRWaveMetaData));    
    LAL_CALL( GetNRMetaDataFromFrameHistory( &metaData, frame->history), &status);

    /* if we find parameters in range then write to output */
    if ( MetaDataInRange(&metaData, &range)) {

      REAL8 tmp;

      /* need to print other pars as well */
      fprintf(stdout,"%s\n",frInCache->frameFiles[k].url);

      /* set up the inspiral table linked list*/
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

      tmp = sqrt(metaData.massRatio) + 1.0/sqrt(metaData.massRatio);
      this_inj->eta = 1.0/( tmp * tmp );

      this_inj->spin1x = metaData.spin1[0];
      this_inj->spin1y = metaData.spin1[1];
      this_inj->spin1z = metaData.spin1[2];

      this_inj->spin2x = metaData.spin2[0];
      this_inj->spin2y = metaData.spin2[1];
      this_inj->spin2z = metaData.spin2[2];

    }    


  } /* end loop over framefiles */



  /* open the xml file */
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, "nrout.xml" ), &status );


  if ( injections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
				      sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }

  /* close the linked list */
  while ( injections.simInspiralTable )
  {
    this_inj = injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    LALFree( this_inj );
  }

  /* close cache */
  /*   LAL_CALL( LALFrClose( &status, &frStream ), &status ); */
  LAL_CALL( LALDestroyFrCache( &status, &frInCache ), &status );  

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  LAL_CALL (LALDestroyUserVars(&status), &status);
  LALCheckMemoryLeaks();

  exit(0);
}


/* metadata is stored in the history field comment 
   -- this function parses the comment to fill the metadata struct */
void GetNRMetaDataFromFrameHistory(NRWaveMetaData *data,
				   FrHistory *history)
{

  UINT4 strlen=128;
  CHAR *comment=NULL; /* the comments string */
  FrHistory *localhist;

  comment = LALMalloc(strlen*sizeof(CHAR));

  localhist = history;
  while (localhist) {

    strcpy(comment,localhist->comment);
    ParseThisComment(data, comment);

    localhist = localhist->next;
  }

  LALFree(comment);
  return;
}


void ParseThisComment(NRWaveMetaData *data,
		      CHAR *comment)
{

  if (strstr(comment,"spin1x"))
    data->spin1[0] = GetNumericValue(comment);
  else if (strstr(comment,"spin1y"))
    data->spin1[1] = GetNumericValue(comment);
  else if (strstr(comment,"spin1z"))
    data->spin1[2] = GetNumericValue(comment);
  else if (strstr(comment,"spin2x"))
    data->spin2[0] = GetNumericValue(comment);
  else if (strstr(comment,"spin2y"))
    data->spin2[1] = GetNumericValue(comment);
  else if (strstr(comment,"spin2z"))
    data->spin2[2] = GetNumericValue(comment);
  else if (strstr(comment,"mass-ratio"))
    data->massRatio = GetNumericValue(comment);
  
  return;

}

REAL8 GetNumericValue(CHAR *comment)
{

  CHAR *tmp;
  tmp = strstr(comment, ":");
  tmp++;

  /* todo: might want to skip whitespace before calling atof */

  return(atof(tmp));

}


/* typedef struct  */
/* { */
/*   REAL8 massRatio; /\**< Mass ratio m1/m2 where we assume m1 >= m2*\/ */
/*   REAL8 spin1[3];  /\**< Spin of m1 *\/ */
/*   REAL8 spin2[3];  /\**< Spin of m2 *\/ */
/*   INT4  mode[2];   /\**< l and m values *\/ */
/*   CHAR  filename[LALNameLength]; /\**< filename where data is stored *\/ */
/* }  */

/* typedef struct { */
/*   REAL8 massRatioMin, massRatioMax; */
/*   REAL8 sx1Min, sx1Max; */
/*   REAL8 sx2Min, sx2Max; */
/*   REAL8 sy1Min, sy1Max; */
/*   REAL8 sy2Min, sy2Max; */
/*   REAL8 sz1Min, sz1Max; */
/*   REAL8 sz2Min, sz2Max;  */
/* } NrParRange; */

REAL8 MetaDataInRange(NRWaveMetaData *data, NrParRange *range)
{

  REAL8 ret;
  BOOLEAN flag = FALSE;

  flag = (data->massRatio >= range->massRatioMin) && (data->massRatio <= range->massRatioMax);
  flag = flag && (data->spin1[0] >= range->sx1Min) && (data->spin1[0] <= range->sx1Max);
  flag = flag && (data->spin2[0] >= range->sx2Min) && (data->spin2[0] <= range->sx2Max);  
  flag = flag && (data->spin1[1] >= range->sy1Min) && (data->spin1[1] <= range->sy1Max);
  flag = flag && (data->spin2[1] >= range->sy2Min) && (data->spin2[1] <= range->sy2Max);
  flag = flag && (data->spin1[2] >= range->sz1Min) && (data->spin1[2] <= range->sz1Max);
  flag = flag && (data->spin2[2] >= range->sz2Min) && (data->spin2[2] <= range->sz2Max);

  if (flag) 
    ret = 1;
  else
    ret = 0; 

  return(ret);

}
