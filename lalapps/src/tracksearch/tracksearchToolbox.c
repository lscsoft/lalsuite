/*
 * Author: Torres C. (Univ of TX at Brownsville)
 */
#include "tracksearch.h"
#include "tracksearchToolbox.h"
#include "unistd.h"

/* Code Identifying information */
NRCSID( TRACKSEARCHC, "tracksearch $Id$");
RCSID( "tracksearch $Id$");
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


void LALappsTSAReadMapFile( LALStatus         *status,
			    TSAMap           **tfImage,
			    CHARVector        *fileNameVec
			    )
{
  INT4 file;
  INT4 dataRead;
  TrackSearchMapMarkingParams tmpImageBorders;
  CreateTimeFreqIn            tmpImageCreateParams;
  UINT4                        cW=0;
  INT4                        i=0;
  INT4                        j=0;
  LIGOTimeGPS                 cMS;
  /*
   * Read in first element of struture so we
   * can properly allocate/read in second and third elements
   */
  if ((file = open(fileNameVec->data,O_RDONLY,0)) == -1)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  if ((dataRead=read(file,
		     (void*) &tmpImageCreateParams,
		     sizeof(tmpImageCreateParams))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  if ((dataRead=read(file,
		     (void*) &cW,
		     sizeof(cW))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);
    
  if ((dataRead=read(file,
		     (void*) &cMS,
		     sizeof(cMS))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  if ((dataRead=read(file,
		     (void*) &tmpImageBorders,
		     sizeof(tmpImageBorders))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);
  /* 
   * Allocate the struture
   */
  LALappsTSACreateMap(status,
		      tfImage,
		      &tmpImageBorders,
		      &tmpImageCreateParams);
  (*tfImage)->clippedWith=cW;
  (*tfImage)->clipperMapStart=cMS;
  /*
   * We need to read it into our allocated memory sequentially
   */
  if ((dataRead=read(file,
		     (void*) &((*tfImage)->imageRep->type),
		     sizeof((*tfImage)->imageRep->type))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  if ((dataRead=read(file,
		     (void*) &((*tfImage)->imageRep->fRow),
		     sizeof((*tfImage)->imageRep->fRow))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  if ((dataRead=read(file,
		     (void*) &((*tfImage)->imageRep->tCol),
		     sizeof((*tfImage)->imageRep->tCol))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EREAD,
		    TRACKSEARCHAVERAGERC_MSGEREAD);

  for (i=0;i<(*tfImage)->imageRep->fRow/2+1;i++)
      if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->freqBin[i]),
			 sizeof((*tfImage)->imageRep->freqBin[i]))) <= 0)
	LALappsTSassert(0,
			TRACKSEARCHAVERAGERC_EREAD,
			TRACKSEARCHAVERAGERC_MSGEREAD);
  

  for (i=0;i<(*tfImage)->imageRep->tCol;i++)
        if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->timeInstant[i]),
			 sizeof((*tfImage)->imageRep->timeInstant[i]))) <= 0)
	  LALappsTSassert(0,
			  TRACKSEARCHAVERAGERC_EREAD,
			  TRACKSEARCHAVERAGERC_MSGEREAD);
	
  for (i=0;i<(*tfImage)->imageRep->tCol;i++)
    for (j=0;j<(*tfImage)->imageRep->fRow/2+1;j++)
	 if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->map[i][j]),
			 sizeof((*tfImage)->imageRep->map[i][j]))) <= 0)
	   LALappsTSassert(0,
			   TRACKSEARCHAVERAGERC_EREAD,
			   TRACKSEARCHAVERAGERC_MSGEREAD);
  /*
   * Close the file
   */
  close(file);
  return;
}
/*
 * End LALappsTSAReadMapFile
 */

void LALappsTSAWriteMapFile( LALStatus         *status,
			     TSAMap            *tfImage,
			     CHARVector        *fileNameVec
			     )
{
  INT4 file;
  CHARVector  *thisFilename=NULL;
  UINT4        dflag=0;
  INT4 i=0;
  INT4 j=0;
  UINT4 dataWrite=0;

  if (fileNameVec == NULL)
    {
      dflag=1;
      LAL_CALL(LALCHARCreateVector(status,&thisFilename,256),status);
      sprintf(thisFilename->data,
	      "MAP:Start:%i,%i:Stop:%i,%i:TF:%i,%i:.dat",
	      tfImage->imageBorders.mapStartGPS.gpsSeconds,
	      tfImage->imageBorders.mapStartGPS.gpsNanoSeconds,
	      tfImage->imageBorders.mapStopGPS.gpsSeconds,
	      tfImage->imageBorders.mapStopGPS.gpsNanoSeconds,
	      tfImage->imageBorders.mapTimeBins,
	      tfImage->imageBorders.mapFreqBins);

    }
  else
    thisFilename=fileNameVec;

  file=creat(thisFilename->data,TSAPERMS);
  /*
   * Always write the field
   * imageCreateParams
   * first to so that when reading we know the
   * dimensions of the fields that we need to alloc
   */
  /*
   * Error checking
   */
  LALappsTSassert(tfImage->imageRep != NULL,
		  TRACKSEARCHAVERAGERC_EVAL,
		  TRACKSEARCHAVERAGERC_MSGEVAL);
  /*
   * Write first element of struct
   */
  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageCreateParams),
	    sizeof(tfImage->imageCreateParams))) != 
      sizeof(tfImage->imageCreateParams)
      )
    LALappsTSassert(0,
	    TRACKSEARCHAVERAGERC_EWRITE,
	    TRACKSEARCHAVERAGERC_MSGEWRITE);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->clippedWith),
	    sizeof(tfImage->clippedWith))) != sizeof(tfImage->clippedWith)
      )
    LALappsTSassert(0,
	    TRACKSEARCHAVERAGERC_EWRITE,
	    TRACKSEARCHAVERAGERC_MSGEWRITE);

    if ((dataWrite=write(file,
	    (void*) &(tfImage->clipperMapStart),
	    sizeof(tfImage->clipperMapStart))) != 
	sizeof(tfImage->clipperMapStart)
      )
    LALappsTSassert(0,
	    TRACKSEARCHAVERAGERC_EWRITE,
	    TRACKSEARCHAVERAGERC_MSGEWRITE);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageBorders),
	    sizeof(tfImage->imageBorders))) != sizeof(tfImage->imageBorders)
      )
    LALappsTSassert(0,
	    TRACKSEARCHAVERAGERC_EWRITE,
	    TRACKSEARCHAVERAGERC_MSGEWRITE);

  /*
   * Write the (TimeFreqRep) structure out as atomic variables in order
   */
  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->type),
	    sizeof(tfImage->imageRep->type))) != 
      sizeof(tfImage->imageRep->type)
      )
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EWRITE,
		    TRACKSEARCHAVERAGERC_MSGEWRITE);
  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->fRow),
	    sizeof(tfImage->imageRep->fRow))) != 
      sizeof(tfImage->imageRep->fRow)
      )
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EWRITE,
		    TRACKSEARCHAVERAGERC_MSGEWRITE);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->tCol),
	    sizeof(tfImage->imageRep->tCol))) != 
      sizeof(tfImage->imageRep->tCol)
      )
    LALappsTSassert(0,
		    TRACKSEARCHAVERAGERC_EWRITE,
		    TRACKSEARCHAVERAGERC_MSGEWRITE);

  for (i=0;i<tfImage->imageRep->fRow/2+1;i++)
      if ((dataWrite=write(file,
		(void*) &(tfImage->imageRep->freqBin[i]),
		sizeof(tfImage->imageRep->freqBin[i]))) != 
	  sizeof(tfImage->imageRep->freqBin[i])
	  )
	LALappsTSassert(0,
			TRACKSEARCHAVERAGERC_EWRITE,
			TRACKSEARCHAVERAGERC_MSGEWRITE);

  for (i=0;i<tfImage->imageRep->tCol;i++)
      if ((dataWrite=write(file,
		(void*) &(tfImage->imageRep->timeInstant[i]),
		sizeof(tfImage->imageRep->timeInstant[i]))) != 
	  sizeof(tfImage->imageRep->timeInstant[i])
	  )
	LALappsTSassert(0,
			TRACKSEARCHAVERAGERC_EWRITE,
			TRACKSEARCHAVERAGERC_MSGEWRITE);

  for (i=0;i<tfImage->imageRep->tCol;i++)
    for (j=0;j<tfImage->imageRep->fRow/2+1;j++)
	if ((dataWrite=write(file,
		  (void*) &(tfImage->imageRep->map[i][j]),
		  sizeof(tfImage->imageRep->map[i][j]))) != 
	    sizeof(tfImage->imageRep->map[i][j])
	    )
	  LALappsTSassert(0,
			  TRACKSEARCHAVERAGERC_EWRITE,
			  TRACKSEARCHAVERAGERC_MSGEWRITE);
  /*
   * Close the file
   */
  /*  fsync(file);*/
  close(file);
  /*
   * Deallocation of local ram
   */
  if (dflag==1)
    LAL_CALL(LALCHARDestroyVector(status,&thisFilename),status);
  return;
}
/* 
 * End LALappsTSAWriteMapFile
 */


void
LALappsTSACreateMap(LALStatus                    *status,
		    TSAMap                      **map,
		    TrackSearchMapMarkingParams  *imageBorders,
		    CreateTimeFreqIn             *createParams)
{
  INT4 i=0;
  INT4 j=0;
  /*
   * Error checking
   */
  LALappsTSassert((status != NULL),
		  TRACKSEARCHAVERAGERC_EARGS,
		  TRACKSEARCHAVERAGERC_MSGEARGS);
  LALappsTSassert(((*map) == NULL),
		  TRACKSEARCHAVERAGERC_ENULL,
		  TRACKSEARCHAVERAGERC_MSGENULL);
  LALappsTSassert(((createParams) != NULL),
		  TRACKSEARCHAVERAGERC_EARGS,
		  TRACKSEARCHAVERAGERC_MSGEARGS);

  /*
   * Allocate the memory for TSAMap variable
   */
  *map = (TSAMap*)LALMalloc(sizeof(TSAMap));
  LALappsTSassert((*map)!=NULL,
		  TRACKSEARCHAVERAGERC_EMEM,
		  TRACKSEARCHAVERAGERC_MSGEMEM);
      

  /*
   * Set image border field and creation field
   */
  (*map)->imageCreateParams=*createParams;
  (*map)->imageBorders=*imageBorders;
  (*map)->clippedWith=0;
  (*map)->clipperMapStart.gpsSeconds=0;
  (*map)->clipperMapStart.gpsNanoSeconds=0;
  (*map)->imageRep=NULL;
    LAL_CALL( LALCreateTimeFreqRep(status,
				   &((*map)->imageRep),
				   createParams),
	      status);
    /*
     * Initialize the struct elements to zeros
     */
    /*
     * REMEMBER: Since we fft a real tseries the freq components -f
     * the same for f so we only need to specify the freq components
     * F_display=fRow/2+1 SEE: CreateTimeFreqRep.c Line:107
     */
    for (i=0;i<(*map)->imageRep->tCol;i++)
    for (j=0;j<(*map)->imageRep->fRow/2+1;j++)
      (*map)->imageRep->map[i][j]=0;
    return;
}
/*
 * End LALappTSACreateMap
 */

void
LALappsTSADestroyMap(LALStatus *status,
		     TSAMap   **map)
{
  LAL_CALL( LALDestroyTimeFreqRep(status,
				  &((*map)->imageRep)),
	    status);
  LALFree(*map);
  *map=NULL;
  return;
}
/*
 * End LALappTSADestroyMap
 */


void
LALappsTSAWritePGM(LALStatus  *status,
		   TSAMap     *map,
		   CHARVector *overrideMask)
{
  INT4       i;
  INT4       j;
  FILE      *fp;
  REAL4      maxval=0;
  REAL4      maxnum;
  REAL4      minval=0;
  INT4       pgmval;
  REAL4      temppgm;
  REAL4      currentval;
  REAL4      currentabsval;
  REAL4      temppoint;
  INT4       scale=255;
  INT4       width=0;
  INT4       height=0;
  INT4       killNeg=0;
  CHARVector   *pgmFile=NULL;
  CHARVector   *auxFile=NULL;
  CHARVector   *rawFile=NULL;
  CHARVector   *basicMask=NULL;
  /*
   * Build filenames
   */
  LAL_CALL(LALCHARCreateVector(status,&pgmFile,maxFilenameLength),status);
  LAL_CALL(LALCHARCreateVector(status,&auxFile,maxFilenameLength),status);
  LAL_CALL(LALCHARCreateVector(status,&rawFile,maxFilenameLength),status);
  LAL_CALL(LALCHARCreateVector(status,&basicMask,maxFilenameLength),status);
  /*
   * Setup filename mask
   */
  if (overrideMask != NULL)
    strcpy(basicMask->data,overrideMask->data);
  else
  sprintf(basicMask->data,"MAP:Start:%i,%i:Stop:%i,%i:TF:%i,%i:",
	  map->imageBorders.mapStartGPS.gpsSeconds,
	  map->imageBorders.mapStartGPS.gpsNanoSeconds,
	  map->imageBorders.mapStopGPS.gpsSeconds,
	  map->imageBorders.mapStopGPS.gpsNanoSeconds,
	  map->imageBorders.mapTimeBins,
	  map->imageBorders.mapFreqBins);
  /*
   * Set filenames
   */
  sprintf(pgmFile->data,
	  "%s.pgm",
	  basicMask->data);
  sprintf(rawFile->data,
	  "%s.raw",
	  basicMask->data);
  sprintf(auxFile->data,
	  "%s.aux",
	  basicMask->data);
  
  /* Getting max image value to normalize to scale */
  width=map->imageRep->fRow/2+1;
  height=map->imageRep->tCol;
  killNeg=0;
  /* Write Space Delimited table */
  temppoint=0;
  fp = fopen(rawFile->data,"w");
  maxnum = map->imageRep->map[0][0]; /* Max Value no abs taken on value */
  for (i=(width-1);i>-1;i--)
    { 
      for (j=0;j<height;j++)
	{
	  temppoint = map->imageRep->map[j][i];
	  if (killNeg)
	    {
	      if (map->imageRep->map[j][i] < 0)
		{
		  temppoint = 0;
		}
	    };
	  currentval = temppoint;
	  currentabsval = fabs(temppoint);
	  if (maxval < currentabsval)
	    {
	      maxval = currentabsval;
	    }
	  if (maxnum < currentval)
	    {
	      maxnum = currentval;
	    }
	  if (minval > currentval)
	    {
	      minval = currentval;
	    }
	  fprintf(fp,"%6.18f ",currentval);
	}
      fprintf(fp,"\n");
    }
  fclose(fp);
  /* PGM File Creation */
  fp = fopen(pgmFile->data,"w");
  fprintf(fp,"P2\n");
  fprintf(fp,"#Written by ImageDump\n");
  fprintf(fp,"%i %i\n",height,width);
  fprintf(fp,"%i\n",scale);
  for (i=(width-1);i>-1;i--)
    {
      for (j=0;j<height;j++)
	{
	  temppoint = map->imageRep->map[j][i];
	  if (killNeg)
	    {
	      if (map->imageRep->map[j][i] < 0)
		{
		  temppoint = 0;
		}
	    };
	  currentval = temppoint;
	  temppgm = ((currentval-minval)*(scale/(maxnum-minval)));
	  pgmval = floor(temppgm);
	  fprintf(fp,"%i\n",pgmval);
	}
    }
  fclose(fp);
  /* Writing Aux file for information purpose only */
  fp = fopen(auxFile->data,"w");
  fprintf(fp,"Aux Data\n");
  fprintf(fp,"TF Type %i\n",map->imageCreateParams.type);
  fprintf(fp,"fRow %i\n",map->imageCreateParams.fRow);
  fprintf(fp,"tCol %i\n",map->imageCreateParams.tCol);
  fprintf(fp,"wlengthT %i\n",map->imageCreateParams.wlengthT);
  fprintf(fp,"wlengthF %i\n",map->imageCreateParams.wlengthF);
  fprintf(fp,"\n");
  fprintf(fp,"clippedWith %i\n",map->clippedWith);
  fprintf(fp,"%i_%i\n",
	  map->clipperMapStart.gpsSeconds,
	  map->clipperMapStart.gpsNanoSeconds);
  fprintf(fp,"deltaT %f\n",map->imageBorders.deltaT);
  fprintf(fp,"Start GPS\n");
  fprintf(fp,"%i_%i\n",
	  map->imageBorders.mapStartGPS.gpsSeconds,
	  map->imageBorders.mapStartGPS.gpsNanoSeconds);
  fprintf(fp,"Stop GPS\n");
  fprintf(fp,"%i_%i\n",
	  map->imageBorders.mapStopGPS.gpsSeconds,
	  map->imageBorders.mapStopGPS.gpsNanoSeconds);
  fprintf(fp,"TimeBins %i\n",map->imageBorders.mapTimeBins);
  fprintf(fp,"FreqBins %i\n",map->imageBorders.mapFreqBins);
  fclose(fp);
  if (auxFile)
    LAL_CALL(LALCHARDestroyVector(status,&auxFile),
	     status);
  if (rawFile)
    LAL_CALL(LALCHARDestroyVector(status,&rawFile),
	     status);
  if (pgmFile)
    LAL_CALL(LALCHARDestroyVector(status,&pgmFile),
	     status);
  if (basicMask)
    LAL_CALL(LALCHARDestroyVector(status,&basicMask),
	     status);
  return;
}
/* 
 * End utility function dat to pgm/aux
 */

void
LALappsTSASortCache(LALStatus   *status,
		    TSAcache    *inputCache,
		    UINT4        ignoreMissingFiles)
{
  TSAMap      *tempMap=NULL;
  UINT4        i=0;
  UINT4        index=0;
  UINT4        lastIndex=0;
  UINT4        STOP=0;
  REAL8        tmpTime=0;
  CHARVector  *tmpFilename=NULL;

  /*
   * Value of -1 in mapStartTime means file IO problem
   */
  for (i=0;i<inputCache->numMapFilenames;i++)
    {
      /* 
       * Try to open the ith file
       */
          LALappsTSAReadMapFile(status,&tempMap,inputCache->filename[i]);
      /*
       * Check Error Code
       */
      if (status->statusCode != 0)
	inputCache->mapStartTime[i]=-1;
      else
	LAL_CALL(
		 LALGPStoFloat(status,
			       &(inputCache->mapStartTime[i]),
			       &(tempMap->imageBorders.mapStartGPS)),
		 status);
      LALappsTSADestroyMap(status,
			   &tempMap);
    }

  LAL_CALL(LALCHARCreateVector(status,
			       &tmpFilename,
			       inputCache->filename[0]->length),
	   status);
  
  /*
   * Perform the bubble sort to put the maps in cronological order
   * We aren't considering map gaps here the map joiner does that!
   */
  STOP=0;
  while (STOP < inputCache->numMapFilenames-1)
    {
      STOP=0;
      while (index < inputCache->numMapFilenames-1)
	{
	  if (inputCache->mapStartTime[index] >
	      inputCache->mapStartTime[index+1])
	    {
	      /*
	       * Swap the entries
	       */
	      tmpTime=inputCache->mapStartTime[index];
	      strcpy(tmpFilename->data,inputCache->filename[index]->data);
	      inputCache->mapStartTime[index]=inputCache->mapStartTime[index+1];
	      strcpy(inputCache->filename[index]->data,
		     inputCache->filename[index+1]->data);
	      inputCache->mapStartTime[index+1]=tmpTime;
	      strcpy(inputCache->filename[index+1]->data,
		     tmpFilename->data);
	    }
	  else
	    STOP++;
	  index++;
	}
      index=0;
    }
  LAL_CALL(LALCHARDestroyVector(status,
				&tmpFilename),
	   status);
  return;
}
/*
 * End LALappsTSASortCache
 */

void
LALappsTSALoadCacheFile(LALStatus    *status,
			CHARVector   *filename,
			TSAcache    **mapCache)
{
  UINT4       i=0;
  UINT4       lineCount=0;
  INT4       errCode=0;
  CHAR        scanString[512];
  FILE       *inputFilePtr=NULL; 

  /*
   * Assume the formatting is one filename and path per line
   */
  inputFilePtr=fopen(filename->data,"r");
  lineCount=0;
  /*
   * Estimating number of lines in the file
   */
  while (errCode != EOF)
    {
      errCode=fscanf(inputFilePtr,"%s\n",scanString);
      /*
       * Ignore blank lines
       */
      if (strcmp("",scanString) == 0)
	lineCount--;
      lineCount++;
    }
  lineCount=lineCount-1;
  fclose(inputFilePtr);
  /*
   * Allocate RAM for cache file
   */
  (*mapCache)=LALMalloc(sizeof(TSAcache));
  (*mapCache)->numMapFilenames=lineCount;
  (*mapCache)->mapStartTime=LALMalloc(sizeof(REAL8)*
					   (*mapCache)->numMapFilenames);
  (*mapCache)->filename=(CHARVector**) LALMalloc((*mapCache)->numMapFilenames*sizeof(CHARVector*));
  for (i=0;i<(*mapCache)->numMapFilenames;i++)
    {
      (*mapCache)->mapStartTime[i]=-1;
      (*mapCache)->filename[i]=NULL;
      LAL_CALL(LALCHARCreateVector(status,
				   &((*mapCache)->filename[i]),
				   maxFilenameLength),
	       status);
    }
  /*
   * 2nd Pass reopen the file and read in the entries
   */
    inputFilePtr=fopen(filename->data,"r");
    for (i=0;i<(*mapCache)->numMapFilenames;i++)
      {
	fscanf(inputFilePtr,"%s\n",(*mapCache)->filename[i]->data);
	LALappsTSassert((inputFilePtr != EOF),
			TRACKSEARCHAVERAGERC_EREAD,
			TRACKSEARCHAVERAGERC_MSGEREAD);
			
	/*
	 * Ignore blank lines
	 */
	if (strcmp("",scanString) == 0)
	  i--;
      }
    fclose(inputFilePtr);
}
/*
 * End LALappsTSALoadCacheFile
 */

void
LALappsTSADestroyCache(LALStatus       *status,
		       TSAcache       **cache)
{
  LALappsTSassert(((*cache) != NULL),
		  TRACKSEARCHTOOLBOXC_EFAIL,
		  TRACKSEARCHTOOLBOXC_EMSGFAIL);
  UINT4   i=0;

  for (i=0;i<(*cache)->numMapFilenames;i++)
    LAL_CALL(
	     LALCHARDestroyVector(status,
				  &((*cache)->filename[i])),
	     status);
  LALFree((*cache)->filename);
  LALFree((*cache)->mapStartTime);
  (*cache)->numMapFilenames=0;
  LALFree(*cache);
  return;
}
/*
 * End LALappsTSADestroyCache
 */

/*
 * Misc shorthand function
 */

void
LALappsTSassert(UINT4               err,
		INT4               code,
		const CHAR         *msg)

{
  if (err != 1)
    {
      fprintf(stderr,"Error %i: %s\n",code,msg);
      exit(code);
    }
  return;
}
