 /*
  * Copyright (C) 2004, 2005 Cristina V. Torres
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

#include "tracksearch.h"
#include "tracksearchAverager.h"



#define PROGRAM_NAME "TSaverager"
typedef struct
{
  INT4    argc;
  CHAR**  argv;
}LALInitSearchParams;

/* Code Identifying information */
NRCSID( TRACKSEARCHAVERAGERC, "tracksearchAverager $Id$");
RCSID( "tracksearchAverager $Id$");
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */

/* Usage format string. */
#define USAGE "Still in flux"

/*
 * Begin MAIN
 */
int main(int argc, char *argv[])
{
  LALStatus       status=blank_status;
  TSAparams       params;
  UINT4           i=0;
  TSAcache       *multiSetCache=NULL;
  /*
   *Sleep for Attaching DDD 
   */
  unsigned int doze = 0;
  pid_t myPID;
  myPID = getpid( );
  fprintf( stdout, "pid %d sleeping for %d seconds\n", myPID, doze );
  fflush( stdout );
  sleep( doze );
  fprintf( stdout, "pid %d awake\n", myPID );
  fflush( stdout );

  /* SET LAL DEBUG STUFF */
  set_debug_level("ERROR | WARNING");
  lal_errhandler = LAL_ERR_RTRN;

  /*set_debug_level("ALLDBG");*/

  LALappsTSAInitialize(&status,
		       argc,
		       argv,
		       &params);
  /*
   * Load up and process multiple caches if requested
   */
  if (params.multiCacheFilename)
    {
      LALappsTSALoadCacheFile(&status,
			      params.multiCacheFilename,
			      &multiSetCache);
      /*
       * Loop through each entry to process that map cache
       */
      for (i=0;i<multiSetCache->numMapFilenames;i++)
	{
	  params.cacheFilename=XLALCreateCHARVector(maxFilenameLength);
	  strcpy(params.cacheFilename->data,multiSetCache->filename[i]->data);
	  LALappsTSAProcessSingleCache(&status,
				       &params);
	  if (params.cacheFilename)
	    XLALDestroyCHARVector(params.cacheFilename);
	}

      if (multiSetCache)
	LALappsTSADestroyCache(&status,
			       &multiSetCache);
    }
  else
    {
      /*
       * Else just execute a single cache job
       */
      LALappsTSAProcessSingleCache(&status,
				   &params);
    }
  
  if (params.cacheFilename)
    XLALDestroyCHARVector(params.cacheFilename);
  
  if (params.multiCacheFilename)
    XLALDestroyCHARVector(params.multiCacheFilename);

  LALCheckMemoryLeaks();
  return 0;
}
/*
 * End MAIN
 */


/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/*
 * Averager Routines
 */
void
LALappsTSAInitialize(LALStatus        *status,
		     int               argc,
		     char*             argv[],
		     TSAparams        *params)
{
  int C;

  struct option long_options[] =
    {
      {"map_cache",            required_argument,  0,   'a'},
      {"new_t_bins",           required_argument,  0,   'b'},
      {"new_f_bins",           required_argument,  0,   'c'},
      {"merge_maps",           no_argument,        0,   'd'},
      {"multi_cache",          required_argument,  0,   'e'},
      {0,                      0,                  0,     0}
    };

  /*
   * Default Arguements
   */
  params->cacheFilename=NULL;
  params->multiCacheFilename=NULL;
  params->operation=Undefined;
  params->colParams.averageTBins=-1;
  params->colParams.averageFBins=-1;
  params->colParams.newTDim=-1;
  params->colParams.newFDim=-1;
  params->verbosity=printFiles;

  if (argc < TRACKSEARCHAVERAGERC_NARGS ) /* Not enough arguments to run */
    {
      fprintf(stderr,TRACKSEARCHAVERAGERC_MSGEARGS);
      fprintf(stderr,"\n");
      exit(TRACKSEARCHAVERAGERC_EARGS);
    }

  while (TRUE)
    {
      int option_index=0;
      C = getopt_long_only(argc,
			   argv,
			   "a:b:c:d",
			   long_options,
			   &option_index);
      if ( C == -1)
	{
	  break;
	}
      switch( C )
	{
	case 0 :
	  if ( long_options[option_index].flag != 0 )
	    {
	      break;
	    }
	  else
	    {
	      fprintf(stderr, "error parsing option %s with argument %s\n",
		      long_options[option_index].name, optarg );
	      exit( 1 );
	    }
	  break;

	case 'a':
	  {
	    params->cacheFilename=XLALCreateCHARVector(maxFilenameLength);
	    strcpy(params->cacheFilename->data,optarg);
	  }
	  break;

	case 'b':
	  {
	    params->colParams.newTDim=(UINT4) atoi(optarg);
	  }
	  break;

	case 'c':
	  {
	    params->colParams.newFDim=(UINT4) atoi(optarg);
	  }
	  break;

	case 'd':
	  {
	    params->operation=Merge;
	  }
	  break;

	case 'e':
	  {
	    /*
	     * If we are sending a txt cache of cache's
	     */
	    params->multiCacheFilename=XLALCreateCHARVector(maxFilenameLength);
	    strcpy(params->multiCacheFilename->data,optarg);
	  }
	  break;

	default :
	  {
	    fprintf(stderr,TRACKSEARCHAVERAGERC_MSGEMISC);
	    exit(TRACKSEARCHAVERAGERC_EMISC);
	  }
	  break;
	}
    }
  /*
   * Simple argument consistency checking 
   */
  LALappsTSassert((
		   (params->cacheFilename != NULL)
		   ||
		   (params->multiCacheFilename != NULL)
		   ),
		  TRACKSEARCHAVERAGERC_EVAL,
		  TRACKSEARCHAVERAGERC_MSGEVAL);
    
  return;
}
/*
 * End LALappsTSAInitialize
 */
/*JUST A THOUGHT OR IMPOSE THAT COLLAPSE VALUE MOD LENGTH =0
 * We should add all the maps together first to avoid losing n points
 *	  off the end of every map to be merged together
 *        ie 256 point collapsing every 10 points loses 6 points
 *        resulting map is 25 pixels long.  In the event of adding 10 maps
 *        we loose a total of 60 time bins
 */
void
LALappsTSACollapseMap(LALStatus         *status,
		      TSAMap            **mapDptr,
		      TSACollapseParams   cparams)
{
  /*tCol     -> TimeBins -> Height  -> PGM Freq Axis (X)*/
  /*fRow/2+1 -> FreqBins -> Width   -> PGM Time Axis (Y)*/

  INT4    h=0;
  INT4    i=0;
  INT4    j=0;
  INT4    k=0;
  INT4    pixelCount=0;
  TSAMap  *mapPrime=NULL;/*Map collapsed in T and F directions*/
  TSAMap  *mapPrimeF=NULL;/*Intermediate map collapsed in F only*/
  TSAMap  *map=NULL;
  REAL4    tmpSum=0;
  CreateTimeFreqIn  mapPrimeCreateParams;
  CreateTimeFreqIn  mapPrimeFreqCreateParams;
  TrackSearchMapMarkingParams mapPrimeMarkerParams;
  TrackSearchMapMarkingParams mapPrimeFreqMarkerParams;

  /*
   * Are any inputs NULL
   */
  LALappsTSassert(((*mapDptr)!=NULL),
		  TRACKSEARCHAVERAGERC_ENPTR,
		  TRACKSEARCHAVERAGERC_MSGENPTR);
  /*
   * Use pointer on double pointer for manipulations
   */
  map=*mapDptr;
  /*
   *Check that we can "evenly" collapse the map via the input
   *  collapse bins specified
   */
  cparams.averageFBins=abs(cparams.averageFBins);
  cparams.averageTBins=abs(cparams.averageTBins);
  LALappsTSassert((map->imageRep->fRow/2+1 >= cparams.averageFBins),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);
  LALappsTSassert((map->imageRep->tCol >= cparams.averageTBins),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);
  LALappsTSassert((map->imageRep->fRow/2+1 >= cparams.newFDim),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);
  LALappsTSassert((map->imageRep->tCol >= cparams.newTDim),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);

  /*
   * Set averaging field of cparams from user input 
   */
  /*Rounding via floor(x.y+.5) = x if y <0.5 x+1 if y > 0.5*/
  if (cparams.newTDim > 0)
    cparams.averageTBins=
      (INT4)floor((map->imageRep->tCol/cparams.newTDim)+0.5);
  else
    cparams.averageTBins=1;
  if (cparams.newFDim > 0)
    cparams.averageFBins=
      (INT4)floor((map->imageRep->fRow/2+1/cparams.newFDim)+0.5);
  else
    cparams.averageFBins=1;

  /*
   * Determine size of new collapsed map and create 
   *  map structure to send back include appropriate metadata
   */
  mapPrimeCreateParams=map->imageCreateParams;
  mapPrimeCreateParams.fRow=map->imageCreateParams.fRow/cparams.averageFBins;
  mapPrimeCreateParams.tCol=map->imageCreateParams.tCol/cparams.averageTBins;
  /*
   * We don't touch the imageBorder.deltaT information
   * we need to figure out how though later!
   */
  mapPrimeMarkerParams=map->imageBorders;
  mapPrimeMarkerParams.mapTimeBins=mapPrimeCreateParams.tCol;
  mapPrimeMarkerParams.mapFreqBins=mapPrimeCreateParams.fRow/2+1;
  LALappsTSACreateMap(status,		      &mapPrime,
		      &mapPrimeMarkerParams,
		      &mapPrimeCreateParams);

  mapPrimeFreqCreateParams=map->imageCreateParams;
  mapPrimeFreqCreateParams.fRow=mapPrimeCreateParams.fRow;
  mapPrimeFreqMarkerParams=map->imageBorders;
  mapPrimeFreqMarkerParams.mapFreqBins=mapPrimeFreqCreateParams.fRow/2+1;
  
  LALappsTSACreateMap(status,
		      &mapPrimeF,
		      &mapPrimeFreqMarkerParams,
		      &mapPrimeFreqCreateParams);
  /*
   * Begin the map collapse
   * Collapse fBins first into mapPrimeF
   */
  /* 
   * Note: Logic flaw on map collapse!
   * Consider a pixel of avg power P in each pixel.  If the curve is
   * 40 pixels long the total power is 40*P.  If we collapse that map
   * in the time direction by 10 columns -> 1 column.  The new curve
   * would be 4 pixels long. It should have the same total power!  We
   * can not average the pixels and keep the average from the 40 pixel
   * map.  Instead the 4 pixels must have an average power per 
   * pixel of 10 times the original map pixel average. ie P'=10P
   * We were doing P'=10P/10=P but we know P'>P it has to be!
   * If we average we keep the curvature dynamic range down!
   */ 
  /*
   * Cristina  Tue-Oct-30-2007:200710301503 
   * To address my comment above we will remove the averaging. So that
   * the power in 10 pixels collapsed into 1 is equal.  The new 1
   * pixel has the same power as original sum in 10 pixels
   */
  for (h=0;h<map->imageRep->tCol;h++)
    {
      for (i=0,k=0;i<(map->imageRep->fRow/2+1);i=i+cparams.averageFBins,k++)
	{
	  for (j=0,tmpSum=0,pixelCount=0;
	       ((j<cparams.averageFBins) && ((i+j)<
					     map->imageRep->fRow/2+1));
	       j++,pixelCount++)
	    {
	      tmpSum=tmpSum+map->imageRep->map[h][i+j];
	    }
	  /*mapPrimeF->imageRep->map[h][k]=tmpSum/pixelCount;*/
	  mapPrimeF->imageRep->map[h][k]=tmpSum;
	}
    }
  /*
   * Now collapse mapPrimeF into mapPrime 
   * We will collapse in the tCol direction now
   */
  for (h=0;(h<mapPrimeF->imageRep->fRow/2+1);h++)
    {
      for (i=0,k=0;(i<mapPrimeF->imageRep->tCol);i=i+cparams.averageTBins,k++)
	{
	  for(j=0,tmpSum=0,pixelCount=0;
	      ((j<cparams.averageTBins) &&
	       ((i+j)<mapPrimeF->imageRep->tCol));
	      j++,pixelCount++)
	    {
	      /* [time][freq] */
	      tmpSum=tmpSum+mapPrimeF->imageRep->map[(i+j)][h];
	    }
	  /*mapPrime->imageRep->map[k][h]=tmpSum/pixelCount;*/
	  mapPrime->imageRep->map[k][h]=tmpSum;
	}
    }
  /*
   * NonCritical Tue-Oct-30-2007:200710301618 
   * Add in manipulation of freqBin and timeInstant entries after the 
   * collapse or merge of the map.  
   */
  /*
   * Free original map and swap pointer to mapPrime
   * Also free intermediate mapPrimeF
   */
  LALappsTSADestroyMap(status,&mapPrimeF);
  LALappsTSADestroyMap(status,&map);
  *mapDptr=mapPrime;
  mapPrime=NULL;
  return;
}
/*
 * End LALappsTSACollapseMap
 */

void
LALappsTSAMergeMap(LALStatus  *status,
		   TSAMap    **output,
		   TSAMap      mapA,
		   TSAMap      mapB)
{
  TSAMap   *inputA=NULL;
  TSAMap   *inputB=NULL;
  TSAMap   *inputC=NULL;
  TSAMap   *outputPtr=NULL;
  REAL8     timeA=0;
  REAL8     timeB=0;
  REAL8     timeAstop=0;
  REAL8     timeBstop=0;
  REAL8     timeResA=0;
  REAL8     timeResB=0;
  REAL8     avgTimeRes=0;
  const REAL8     resTolPercent=0.01;/*map time resolution tol*/
  BOOLEAN   linkMaps=0;
  INT4     i=0;
  INT4     j=0;
  INT4     k=0;
  INT4     splitPixelA=0;
  INT4     splitPixelB=0;
  INT4     relativePixelMarkerB=0;
  INT4     pixelOverlap=0;
  INT4     outputTimeBins=0;
  TrackSearchMapMarkingParams    imageBordersOut;
  CreateTimeFreqIn               imageCreateParamsOut;

  /* 
   * User pointers to place input maps in time order
   */
  inputA=&mapA;
  inputB=&mapB;
  timeA = XLALGPSGetREAL8(&(inputA->imageBorders.mapStartGPS));
  timeB = XLALGPSGetREAL8(&(inputB->imageBorders.mapStartGPS));
  LALappsTSassert((*output==NULL),
		  TRACKSEARCHAVERAGERC_ENPTR,
		  TRACKSEARCHAVERAGERC_MSGENPTR);
  /* TimeB comes before TimeA swap pointers so A comes before B*/
  if (timeB < timeA)
    {
      inputC=inputA;
      inputA=NULL;
      inputA=inputB;
      inputB=NULL;
      inputB=inputC;
      inputC=NULL;
    }
  /* 
   * Check inputs A and B are same dimensions in frequency
   */
  LALappsTSassert((inputA->imageRep->fRow==inputB->imageRep->fRow),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);
  /*
   *  LALappsTSassert((inputA->imageRep->tCol==inputB->imageRep->tCol),
   *	  TRACKSEARCHAVERAGERC_EDIMS,
   *	  TRACKSEARCHAVERAGERC_MSGEDIMS);
   */
  timeA = XLALGPSGetREAL8(&(inputA->imageBorders.mapStartGPS));
  timeB = XLALGPSGetREAL8(&(inputB->imageBorders.mapStartGPS));
  timeAstop = XLALGPSGetREAL8(&(inputA->imageBorders.mapStopGPS));
  timeBstop = XLALGPSGetREAL8(&(inputB->imageBorders.mapStopGPS));
  /*
   * InputA should end before inputB ends and after equal B starts
   */
  LALappsTSassert(((timeAstop<=timeBstop) && (timeAstop >= timeB)),
		  TRACKSEARCHAVERAGERC_EMAPO,
		  TRACKSEARCHAVERAGERC_MSGEMAPO);
  /*
   * Check to see that mapB starts before or at end of mapA and starts
   * after mapA starts
   */
  LALappsTSassert(((timeAstop>=timeB) && (timeA<=timeB)),
		  TRACKSEARCHAVERAGERC_EMAPO,
		  TRACKSEARCHAVERAGERC_MSGEMAPO);
  /*
   * Check that both maps don't cover matching time stretches
   */
  LALappsTSassert(((timeA!=timeB)||(timeAstop!=timeBstop)),
		  TRACKSEARCHAVERAGERC_EMAPO,
		  TRACKSEARCHAVERAGERC_MSGEMAPO);
  /* 
   *Determine the number of time bins in the output structure
   */
  timeResA=(REAL8)((timeAstop-timeA)/(REAL8)inputA->imageCreateParams.tCol);
  timeResB=(REAL8)((timeBstop-timeB)/(REAL8)inputB->imageCreateParams.tCol);
  /*
   * Check that the resolutions of the maps match
   */
  avgTimeRes=(timeResA+timeResB)/2;
  LALappsTSassert((fabs((timeResA-timeResB)/avgTimeRes)<resTolPercent),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);
  outputTimeBins=floor((timeBstop-timeA)/avgTimeRes);
  /*
   * Check for reasonable outputTimeBins must be less that sum of 
   * the two input maps timebins and greate that 2 time cols
   */
  LALappsTSassert((outputTimeBins<=(inputA->imageCreateParams.tCol +
				    inputB->imageCreateParams.tCol)),
		  TRACKSEARCHAVERAGERC_ESUB,
		  TRACKSEARCHAVERAGERC_MSGESUB);
  LALappsTSassert((outputTimeBins>=2),
		  TRACKSEARCHAVERAGERC_ESUB,
		  TRACKSEARCHAVERAGERC_MSGESUB);
  /*
   * Allocate an output map equal to the two maps combined in time
   * we account for the overlaping sections such that our new output
   * maps spans the time inputA..gpsStart until inputB..gpsStop
   */
  imageBordersOut=inputA->imageBorders;
  imageBordersOut.mapTimeBins=outputTimeBins;
  imageBordersOut.mapStopGPS=inputB->imageBorders.mapStopGPS;
  imageCreateParamsOut=inputA->imageCreateParams;
  imageCreateParamsOut.tCol=outputTimeBins;
  /*
   * Allocate output map struct
   */
  LALappsTSACreateMap(status,
		      output,
		      &imageBordersOut,
		      &imageCreateParamsOut);
  /*
   * Handle map "blending" if sets are overlapped
   * or handle a straight joining if possible
   */
  outputPtr=*output;
  if (((outputPtr)->imageRep->tCol==
       (inputB->imageRep->tCol + inputA->imageRep->tCol)))
    linkMaps=1;

  if (linkMaps)
    {
      fprintf(stdout,"Doing a direct map joining\n");
      /*
       * Copy over the map matrix
       */
      outputPtr->clippedWith=0;
      /*
       * Start appending in time
       */
      for (i=0;i<inputA->imageRep->fRow/2+1;i++)
	{
	  for (j=0,k=0;j<inputA->imageRep->tCol;j++,k++)
	    {
	      /*[time][freq]*/
	      outputPtr->imageRep->map[k][i]=inputA->imageRep->map[j][i];
	    }
	}
      for (i=0;i<inputB->imageRep->fRow/2+1;i++)
	{
	  for (j=0,k=inputA->imageRep->tCol;j<inputB->imageRep->tCol;j++,k++)
	    {
	      /*[time][freq]*/
	      outputPtr->imageRep->map[k][i]=inputB->imageRep->map[j][i];
	    }
	}
    }
  else
    {
      fprintf(stdout,"Doing a mergee/clip of the two maps\n");
      /* 
       * Figure out how to merge the two maps accordingly
       */
      /*
       * Since we are "clipping inputB" we need to set the clipping field
       */
      outputPtr->clippedWith=1;
      outputPtr->clipperMapStart=inputB->imageBorders.mapStartGPS;
      /* 
       * We assume that outputPtr->imageRep has been correctly allocated
       * outside of this function via a call to LALCreateTimeFreqRep
       */
      /*
       * We must determine where to mark the maps for the truncations
       * and then the appending the marked sections into outputPtr
       */
      relativePixelMarkerB=(timeB-timeA)/avgTimeRes;
      pixelOverlap=inputA->imageRep->tCol-relativePixelMarkerB;
      splitPixelA=
	inputA->imageRep->tCol-(pixelOverlap/2);
      splitPixelB=(pixelOverlap/2+1);
      /*
       *Stop inputA copy at splitPixelA
       *Start inputB copy at splitPixelB
       */
      for (i=0;i<inputA->imageRep->fRow/2+1;i++)
	{
	  for (j=0;j<splitPixelA;j++)
	    {
	      /*[time][freq]*/
	      outputPtr->imageRep->map[j][i]=inputA->imageRep->map[j][i];
	    }
	}
      for (i=0;i<inputB->imageRep->fRow/2+1;i++)
	{
	  for (k=splitPixelA,j=splitPixelB;j<inputB->imageRep->tCol;j++,k++)
	    {
	      /*[time][freq]*/
	      outputPtr->imageRep->map[k][i]=inputB->imageRep->map[j][i];
	    }
	}
    }
  /*
   * NonCritical Tue-Oct-30-2007:200710301618 
   * Add in manipulation of freqBin and timeInstant entries after the 
   * collapse or merge of the map.  
   */
  return;
}
/*
 * End LALappsTSAMergeMap
 */


void
LALappsTSAProcessSingleCache(LALStatus    *status,
			     TSAparams    *params)
{
  TSAcache       *cache=NULL;
  TSAMap        **mapArray=NULL;
  TSAMap         *mergeResultMap=NULL;
  TSAMap         *mergeMapTmp=NULL;
  UINT4           i=0;

  /*
   * Read in the cache file to memory
   */
  if (params->cacheFilename == NULL)
    {
      fprintf(stderr,"No map cache file specified!\n");
      LALappsTSassert((params->cacheFilename != NULL),
		      TRACKSEARCHAVERAGERC_ENULL,
		      TRACKSEARCHAVERAGERC_MSGENULL);
    }
  else
    {
      /*
       * Load up the cache file
       */
      if (params->verbosity > quiet)
	fprintf(stdout,"Loading map cache file\n");

      LALappsTSALoadCacheFile(status,
			      params->cacheFilename,
			      &cache);
      /*
       * Sort cache structure chrono
       logically
      */
      if (params->verbosity > quiet)
	fprintf(stdout,"Sorting cache file entries chronologically\n");
      LALappsTSASortCache(status,
			  cache,
			  0);
      if (params->verbosity > quiet)
	{
	  /*
	   * Write out sorted list of map files to screen
	   */
	  for (i=0;i<cache->numMapFilenames;i++)
	    fprintf(stdout,"Time:%f Map name %s\n",
		    cache->mapStartTime[i],
		    cache->filename[i]->data);
	}
    }
  /*
   * Create and Load array of TSA maps from disk
   */
  mapArray=(TSAMap**)LALMalloc(sizeof(TSAMap*)*cache->numMapFilenames);
  for (i=0;i<cache->numMapFilenames;i++)
    {
      mapArray[i]=NULL;
      LALappsTSAReadMapFile(status,
			    &(mapArray[i]),
			    cache->filename[i]);
      LALappsTSassert((mapArray[i] != NULL),
		      TRACKSEARCHAVERAGERC_EMEM,
		      TRACKSEARCHAVERAGERC_MSGEMEM);
    }
  /*
   * Begin merging the files into one output TSAMap
   */
  if (params->operation == Merge)
    {
      if (params->verbosity > quiet)
	fprintf(stdout,"Merging maps\n");
      LALappsTSassert((cache->numMapFilenames >= 2),
		      TRACKSEARCHAVERAGERC_EVAL,
		      TRACKSEARCHAVERAGERC_MSGEVAL);
      /*
       * Merge first two maps "by hand"
       */
      if (params->verbosity > quiet)
	fprintf(stdout,"Merging first two maps by hand\n");
      LALappsTSAMergeMap(status,
			 &mergeResultMap,
			 *(mapArray[0]),
			 *(mapArray[1]));
      /*
       * Loop through the rest here
       */
      if (cache->numMapFilenames > 2)
	{
	  for (i=2;i<cache->numMapFilenames;i++)
	    {
	      if (params->verbosity > quiet)
		fprintf(stdout,"Looping maps %i \n",i);
	      mergeMapTmp=mergeResultMap;
	      mergeResultMap=NULL;
	      LALappsTSAMergeMap(status,
				 &mergeResultMap,
				 *mergeMapTmp,
				 *(mapArray[i]));
	      if (mergeMapTmp != NULL)
		{
		  if (params->verbosity > quiet)
		    fprintf(stdout,"Deleting Struct Entry %i\n",i);
		  LALappsTSADestroyMap(status,
				       &mergeMapTmp);
		}
	    }
	}
      /*
       * Collapse each map if the user has requested it
       */
      if (
	  (params->colParams.newTDim > -1)
	  ||
	  (params->colParams.newFDim > -1)
	  )
	{
	  if (params->verbosity > quiet)
	    fprintf(stdout,"Collapsing individual maps\n");
	  for (i=0;i<cache->numMapFilenames;i++)
	    LALappsTSACollapseMap(status,
				  &(mergeResultMap),
				  params->colParams);
	}
      /*
       * Output the merged map
       */
      if (params->verbosity > quiet)
	fprintf(stdout,"Writing Merged Map\n");

      LALappsTSAWriteMapFile(status,
			     mergeResultMap,
			     NULL);
      /*
       * Free RAM for mergeResultMap
       */
      if (mergeResultMap)
	LALappsTSADestroyMap(status,
			     &mergeResultMap);
    }
  else
    {
      /*ONLY COLLAPSING*/
      /*
       * Collapse each map if the user has requested it
       */
      if (
	  (params->colParams.newTDim > -1)
	  ||
	  (params->colParams.newFDim > -1)
	  )
	{
	  if (params->verbosity > quiet)
	    fprintf(stdout,"Collapsing individual maps\n");
	  for (i=0;i<cache->numMapFilenames;i++)
	    LALappsTSACollapseMap(status,
				  &(mapArray[i]),
				  params->colParams);
	}
      /* 
       * Dump out the array of collapse maps if merge not requested
       */

      if (params->verbosity > quiet)
	fprintf(stdout,"Maps not merged writing maps to disk");

      for (i=0;i<cache->numMapFilenames;i++)
	{
	  LALappsTSAWriteMapFile(status,
				 mapArray[i],
				 NULL);
	}
    }
  /*
   * Deallocate the array of TSA maps
   */
  if (mapArray)
    {
      for (i=0;i<cache->numMapFilenames;i++)
	LALappsTSADestroyMap(status,
			     &(mapArray[i]));
      LALFree(mapArray);
    }
  /*
   * Deallocate the cache file structure
   */
  LALappsTSADestroyCache(status,
			 &cache);
  return;
}
/*
 * End LALappsTSAProcessSingleCache
 */

void
tsaTest(
	LALStatus      *status
	)
{
  TSAMap        *testMapMerge=NULL;
  TSAMap        *testMapA=NULL;
  TSAMap        *testMapB=NULL;
  TSACollapseParams collapseParams;
  /*  LALStatus      status=blank_status;*/
  CHARVector    *fileNameVector=NULL;
  TrackSearchMapMarkingParams  markerParams;
  CreateTimeFreqIn             createParams;
  UINT4          slipPoint=0;


  if (1) /*Run a small section of test code*/
    {
      /*Do Something*/
      fprintf(stderr,"Hi There\n");
      /*
       * Create blanked test map structure
       */
      createParams.type=Undefined;
      /*
       * Number of F rows from -f -> f
       */
      createParams.fRow=64;
      createParams.tCol=128;
      createParams.wlengthT=32;
      createParams.wlengthF=0;
      markerParams.deltaT=1;
      markerParams.dataDeltaT=1;
      markerParams.mapStartGPS.gpsSeconds=1;
      markerParams.mapStartGPS.gpsNanoSeconds=1;
      markerParams.mapStopGPS.gpsSeconds=128;
      markerParams.mapStopGPS.gpsNanoSeconds=0;
      markerParams.mapTimeBins=createParams.tCol;
      /*
       * Only plotable real F bins are tracted
       */
      markerParams.mapFreqBins=createParams.fRow/2+1;
      fprintf(stderr,"Allocating\n");
      LALappsTSACreateMap(status,&testMapA,&markerParams,&createParams);
      fprintf(stderr,"Writing to disk\n");
      fileNameVector=XLALCreateCHARVector(64);
      strcpy(fileNameVector->data,"TestFile.dat");
      LALappsTSAWriteMapFile(status,testMapA,fileNameVector);
      fprintf(stderr,"Reading from disk\n");
      LALappsTSAReadMapFile(status,&testMapB,fileNameVector);
      /*
       * Artifically change time markers of B
       */
      slipPoint=128;
      testMapB->imageBorders.mapStartGPS.gpsSeconds=slipPoint;
      testMapB->imageBorders.mapStartGPS.gpsNanoSeconds=0;
      testMapB->imageBorders.mapStopGPS.gpsSeconds=slipPoint+128;
      testMapB->imageBorders.mapStopGPS.gpsNanoSeconds=1;
      fprintf(stderr,"Writing the READ map back to disk!\n");
      LALappsTSAWriteMapFile(status,testMapB,NULL);
      fprintf(stderr,"Merging both maps structs\n");
      LALappsTSAMergeMap(status,
			 &testMapMerge,
			 *testMapA,
			 *testMapB);
      fprintf(stderr,"Writing the MERGED map back to disk!\n");
      LALappsTSAWriteMapFile(status,testMapMerge,NULL);
      fprintf(stderr,"Collapsing Output test map\n");
      /*
       * Should collapse every two time bins and no 
       * collapse on freq bins
       * Makeing it just like the A and B inputs
       */
      collapseParams.newTDim=(INT4)floor(createParams.tCol/2);
      collapseParams.newFDim=createParams.fRow;
      LALappsTSACollapseMap(status,
			    &testMapMerge,
			    collapseParams);
      fprintf(stderr,"Writing COLLAPSED map back to disk!\n");
      LALappsTSAWriteMapFile(status,testMapMerge,NULL);
      fprintf(stderr,"Deallocating\n");
      LALappsTSADestroyMap(status,&testMapA);
      LALappsTSADestroyMap(status,&testMapB);
      LALappsTSADestroyMap(status,&testMapMerge);
      XLALDestroyCHARVector(fileNameVector);
      fprintf(stderr,"Done\n");
    }
  return;
}
/*
 * End tsaTest Routine
 */
 
