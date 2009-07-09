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


void
LALappsTSACropMap(
		  LALStatus          *status,
		  TSAMap            **inputMap,
		  UINT4               binsToCrop)
{
  TSAMap                       *tmpMap=NULL;
  CreateTimeFreqIn              tmpCreateParams;
  TrackSearchMapMarkingParams   tmpMarkingParams;
  INT4                         i=0;
  INT4                         j=0;
  REAL8                         startTime=0;
  REAL8                         stopTime=0;
  CHARVector                    *name=NULL;
  LALappsTSassert(((*inputMap)!=NULL),
		  TRACKSEARCHAVERAGERC_ENPTR,
		  TRACKSEARCHAVERAGERC_MSGENPTR);

  LALappsTSassert(((*inputMap)->imageRep->tCol > (INT4)(2*binsToCrop)),
		  TRACKSEARCHAVERAGERC_EDIMS,
		  TRACKSEARCHAVERAGERC_MSGEDIMS);

  name=XLALCreateCHARVector(1024);
/*    strcpy(name->data,"Uncroppedmap");  */
/*    LALappsTSAWritePGM(status,(*inputMap),name);  */

  /*
   * New 'cropped' map to be create and swapped into old maps place.
   */
  memcpy(&tmpCreateParams,&((*inputMap)->imageCreateParams),sizeof(CreateTimeFreqIn));
  tmpCreateParams.tCol=(*inputMap)->imageCreateParams.tCol-(2*binsToCrop);
  /**/
  memcpy(&tmpMarkingParams,&((*inputMap)->imageBorders),sizeof(TrackSearchMapMarkingParams));
  tmpMarkingParams.mapTimeBins=(*inputMap)->imageBorders.mapTimeBins-(2*binsToCrop);
  /**/
  /*
   * Adjust times on mapStopGPS and mapStartGPS to be times
   * after cropping the MAP
   */
  /* FIXME:  loss of precision; consider
  XLALGPSAdd(&tmpMarkingParams.mapStartGPS, tmpMarkingParams.deltaT*binsToCrop);
  XLALGPSAdd(&tmpMarkingParams.mapStopGPS, -tmpMarkingParams.deltaT*binsToCrop);
  */
  startTime = XLALGPSGetREAL8(&tmpMarkingParams.mapStartGPS);

  startTime=startTime+(tmpMarkingParams.deltaT*binsToCrop);

  XLALGPSSetREAL8(&tmpMarkingParams.mapStartGPS, startTime);

  stopTime = XLALGPSGetREAL8(&tmpMarkingParams.mapStopGPS);

  stopTime=stopTime-(tmpMarkingParams.deltaT*binsToCrop);

  XLALGPSSetREAL8(&tmpMarkingParams.mapStopGPS, stopTime);

  LALappsTSACreateMap(status,
		      &tmpMap,
		      &tmpMarkingParams,
		      &tmpCreateParams);
  /* 
   * Copy over the data that is to be 'kept'
   */	/*[time][freq]*/
  for (i=0;(i<tmpMap->imageRep->fRow/2+1);i++)
    for(j=0;(j<tmpMap->imageRep->tCol);j++)
      tmpMap->imageRep->map[j][i]=(*inputMap)->imageRep->map[j+binsToCrop][i];

  for (i=0;(i<tmpMap->imageRep->fRow/2+1);i++)
    tmpMap->imageRep->freqBin[i]=(*inputMap)->imageRep->freqBin[i];
  for (i=0;(i<tmpMap->imageRep->tCol);i++)
    tmpMap->imageRep->timeInstant[i]=(*inputMap)->imageRep->timeInstant[i+binsToCrop];

  /*
   * Copy the newly created information into the structure.
   * This should return a restructed 'structure'
   */
  memcpy(&((*inputMap)->imageCreateParams),&(tmpMap->imageCreateParams),sizeof(CreateTimeFreqIn));
  /**/
  (*inputMap)->clippedWith=tmpMap->clippedWith;
  memcpy(&((*inputMap)->clipperMapStart),&(tmpMap->clipperMapStart),sizeof(LIGOTimeGPS));
  /**/
  memcpy(&((*inputMap)->imageBorders),&(tmpMap->imageBorders),sizeof(TrackSearchMapMarkingParams));
  /**/
  LAL_CALL(
	   LALDestroyTimeFreqRep(status,
				 &((*inputMap)->imageRep)),
	   status);
  LAL_CALL(
	   LALCreateTimeFreqRep(status,
				&((*inputMap)->imageRep),
				&((tmpMap)->imageCreateParams)),
	   status);
  /*Copy the elements now*/
  for (i=0;(i<(*inputMap)->imageRep->fRow/2+1);i++)
    for(j=0;(j<(*inputMap)->imageRep->tCol);j++)
      (*inputMap)->imageRep->map[j][i]=tmpMap->imageRep->map[j][i];

  if (name)
    XLALDestroyCHARVector(name);

  if (tmpMap)
    LALappsTSADestroyMap(status,&tmpMap);
  return;
}


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
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);

  if ((dataRead=read(file,
		     (void*) &tmpImageCreateParams,
		     sizeof(tmpImageCreateParams))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);

  if ((dataRead=read(file,
		     (void*) &cW,
		     sizeof(cW))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);
    
  if ((dataRead=read(file,
		     (void*) &cMS,
		     sizeof(cMS))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);

  if ((dataRead=read(file,
		     (void*) &tmpImageBorders,
		     sizeof(tmpImageBorders))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);
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
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);
  if ((dataRead=read(file,
		     (void*) &((*tfImage)->imageRep->fRow),
		     sizeof((*tfImage)->imageRep->fRow))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);

  if ((dataRead=read(file,
		     (void*) &((*tfImage)->imageRep->tCol),
		     sizeof((*tfImage)->imageRep->tCol))) <= 0)
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_ERHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGRHEAD);

  for (i=0;i<(*tfImage)->imageRep->fRow/2+1;i++)
      if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->freqBin[i]),
			 sizeof((*tfImage)->imageRep->freqBin[i]))) <= 0)
	LALappsTSassert(0,
			TRACKSEARCHTOOLBOXC_ERDATA,
			TRACKSEARCHTOOLBOXC_EMSGRDATA);

  for (i=0;i<(*tfImage)->imageRep->tCol;i++)
        if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->timeInstant[i]),
			 sizeof((*tfImage)->imageRep->timeInstant[i]))) <= 0)
	  LALappsTSassert(0,
			  TRACKSEARCHTOOLBOXC_ERDATA,
			  TRACKSEARCHTOOLBOXC_EMSGRDATA);
	
  for (i=0;i<(*tfImage)->imageRep->tCol;i++)
    for (j=0;j<(*tfImage)->imageRep->fRow/2+1;j++)
	 if ((dataRead=read(file,
			 (void*) &((*tfImage)->imageRep->map[i][j]),
			 sizeof((*tfImage)->imageRep->map[i][j]))) <= 0)
	  LALappsTSassert(0,
			  TRACKSEARCHTOOLBOXC_ERDATA,
			  TRACKSEARCHTOOLBOXC_EMSGRDATA);
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
      LALappsDetermineFilename(status,
			       tfImage->imageBorders,
			       &thisFilename,
			       ".dat");
    }
  else
    thisFilename=fileNameVec;

  if ((file=creat(thisFilename->data,TSAPERMS)) == -1)
    {
      /* A quick pause to see if we can try again to create the file.*/
      sleep(1);
      if ((file=creat(thisFilename->data,TSAPERMS)) == -1)
	{
	  fprintf(stderr,"Error creating:%s\n",thisFilename->data);
	  LALappsTSassert(0,
			  TRACKSEARCHTOOLBOXC_EWRITE,
			  TRACKSEARCHTOOLBOXC_EMSGWRITE);
	}
    }
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
		    TRACKSEARCHTOOLBOXC_EHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGHEAD);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->clippedWith),
	    sizeof(tfImage->clippedWith))) != sizeof(tfImage->clippedWith)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGHEAD);

    if ((dataWrite=write(file,
	    (void*) &(tfImage->clipperMapStart),
	    sizeof(tfImage->clipperMapStart))) != 
	sizeof(tfImage->clipperMapStart)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGHEAD);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageBorders),
	    sizeof(tfImage->imageBorders))) != sizeof(tfImage->imageBorders)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EHEAD,
		    TRACKSEARCHTOOLBOXC_EMSGHEAD);
  /*
   * Write the (TimeFreqRep) structure out as atomic variables in order
   */
  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->type),
	    sizeof(tfImage->imageRep->type))) != 
      sizeof(tfImage->imageRep->type)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EDATA,
		    TRACKSEARCHTOOLBOXC_EMSGDATA);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->fRow),
	    sizeof(tfImage->imageRep->fRow))) != 
      sizeof(tfImage->imageRep->fRow)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EDATA,
		    TRACKSEARCHTOOLBOXC_EMSGDATA);

  if ((dataWrite=write(file,
	    (void*) &(tfImage->imageRep->tCol),
	    sizeof(tfImage->imageRep->tCol))) != 
      sizeof(tfImage->imageRep->tCol)
      )
    LALappsTSassert(0,
		    TRACKSEARCHTOOLBOXC_EDATA,
		    TRACKSEARCHTOOLBOXC_EMSGDATA);

  for (i=0;i<tfImage->imageRep->fRow/2+1;i++)
      if ((dataWrite=write(file,
		(void*) &(tfImage->imageRep->freqBin[i]),
		sizeof(tfImage->imageRep->freqBin[i]))) != 
	  sizeof(tfImage->imageRep->freqBin[i])
	  )
	LALappsTSassert(0,
			TRACKSEARCHTOOLBOXC_EDATA,
			TRACKSEARCHTOOLBOXC_EMSGDATA);

  for (i=0;i<tfImage->imageRep->tCol;i++)
      if ((dataWrite=write(file,
		(void*) &(tfImage->imageRep->timeInstant[i]),
		sizeof(tfImage->imageRep->timeInstant[i]))) != 
	  sizeof(tfImage->imageRep->timeInstant[i])
	  )
	LALappsTSassert(0,
			TRACKSEARCHTOOLBOXC_EDATA,
			TRACKSEARCHTOOLBOXC_EMSGDATA);

  for (i=0;i<tfImage->imageRep->tCol;i++)
    for (j=0;j<tfImage->imageRep->fRow/2+1;j++)
	if ((dataWrite=write(file,
		  (void*) &(tfImage->imageRep->map[i][j]),
		  sizeof(tfImage->imageRep->map[i][j]))) != 
	    sizeof(tfImage->imageRep->map[i][j])
	    )
	  LALappsTSassert(0,
			  TRACKSEARCHTOOLBOXC_EDATA,
			  TRACKSEARCHTOOLBOXC_EMSGDATA);
  /*
   * Close the file
   */
  /*  fsync(file);*/
  close(file);
  /*
   * Deallocation of local ram
   */
  if (dflag==1)
    XLALDestroyCHARVector(thisFilename);

  return;
}
/* 
 * End LALappsTSAWriteMapFile
 */

void
LALappsDetermineFilename(LALStatus                   *status,
			 TrackSearchMapMarkingParams  imageBorders,
			 CHARVector                 **thisFilename,
			 const CHAR*                  myExt)
{
  LALappsTSassert((thisFilename != NULL),
		  TRACKSEARCHTOOLBOXC_EFAIL,
		  TRACKSEARCHTOOLBOXC_EMSGFAIL);
  LALappsTSassert((*thisFilename == NULL),
		  TRACKSEARCHTOOLBOXC_EFAIL,
		  TRACKSEARCHTOOLBOXC_EMSGFAIL);
  *thisFilename=XLALCreateCHARVector(maxFilenameLength);

  sprintf((*thisFilename)->data,
	  "MAP:Start:%i,%i:Stop:%i,%i:TF:%i,%i:%s",
	  imageBorders.mapStartGPS.gpsSeconds,
	  imageBorders.mapStartGPS.gpsNanoSeconds,
	  imageBorders.mapStopGPS.gpsSeconds,
	  imageBorders.mapStopGPS.gpsNanoSeconds,
	  imageBorders.mapTimeBins,
	  imageBorders.mapFreqBins,
	  myExt);
}
/*
 * End LALappsDetermineFilename
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
		  TRACKSEARCHTOOLBOXC_EALLOC,
		  TRACKSEARCHTOOLBOXC_EMSGALLOC);

  LALappsTSassert(((*map) == NULL),
		  TRACKSEARCHTOOLBOXC_EALLOC,
		  TRACKSEARCHTOOLBOXC_EMSGALLOC);

  LALappsTSassert(((createParams) != NULL),
		  TRACKSEARCHTOOLBOXC_EALLOC,
		  TRACKSEARCHTOOLBOXC_EMSGALLOC);
  /*
   * Allocate the memory for TSAMap variable
   */
  *map = (TSAMap*) XLALMalloc(sizeof(TSAMap));
  LALappsTSassert((*map)!=NULL,
		  TRACKSEARCHTOOLBOXC_EALLOC,
		  TRACKSEARCHTOOLBOXC_EMSGALLOC);
  

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
  pgmFile=XLALCreateCHARVector(maxFilenameLength);
  auxFile=XLALCreateCHARVector(maxFilenameLength);
  rawFile=XLALCreateCHARVector(maxFilenameLength);
  /*
   * Setup filename mask
   */
  if (overrideMask != NULL)
    {
      basicMask=XLALCreateCHARVector(maxFilenameLength);
      strcpy(basicMask->data,overrideMask->data);
    }
  else
    LALappsDetermineFilename(status,
			     map->imageBorders,
			     &basicMask,
			     "");
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
    XLALDestroyCHARVector(auxFile);

  if (rawFile)
    XLALDestroyCHARVector(rawFile);

  if (pgmFile)
    XLALDestroyCHARVector(pgmFile);

  if (basicMask)
    XLALDestroyCHARVector(basicMask);

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
        inputCache->mapStartTime[i] = XLALGPSGetREAL8(&(tempMap->imageBorders.mapStartGPS));
      LALappsTSADestroyMap(status,
			   &tempMap);
    }

  tmpFilename=XLALCreateCHARVector(inputCache->filename[0]->length);
  
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
  XLALDestroyCHARVector(tmpFilename);

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
  INT4        errCode=0;
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
  (*mapCache)=XLALMalloc(sizeof(TSAcache));
  (*mapCache)->numMapFilenames=lineCount;
  (*mapCache)->mapStartTime=LALMalloc(sizeof(REAL8)*
					   (*mapCache)->numMapFilenames);
  (*mapCache)->filename=(CHARVector**) XLALMalloc((*mapCache)->numMapFilenames*sizeof(CHARVector*));
  for (i=0;i<(*mapCache)->numMapFilenames;i++)
    {
      (*mapCache)->mapStartTime[i]=-1;
      (*mapCache)->filename[i]=NULL;
      (*mapCache)->filename[i]=XLALCreateCHARVector(maxFilenameLength);
    }
  /*
   * 2nd Pass reopen the file and read in the entries
   */
    inputFilePtr=fopen(filename->data,"r");
    for (i=0;i<(*mapCache)->numMapFilenames;i++)
      {
	fscanf(inputFilePtr,"%s\n",(*mapCache)->filename[i]->data);
	LALappsTSassert((inputFilePtr != EOF),
			TRACKSEARCHTOOLBOXC_EALLOC,
			TRACKSEARCHTOOLBOXC_EMSGALLOC);
		
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
  UINT4   i=0;

  LALappsTSassert(((*cache) != NULL),
		  TRACKSEARCHTOOLBOXC_EFAIL,
		  TRACKSEARCHTOOLBOXC_EMSGFAIL);

  for (i=0;i<(*cache)->numMapFilenames;i++)
    XLALDestroyCHARVector(((*cache)->filename[i]));

  XLALFree((*cache)->filename);
  XLALFree((*cache)->mapStartTime);
  (*cache)->numMapFilenames=0;
  XLALFree(*cache);
  return;
}
/*
 * End LALappsTSADestroyCache
 */

void LALappsCreateR4FromR8TimeSeries(LALStatus             *status,
				     REAL4TimeSeries      **R4TS,
				     REAL8TimeSeries       *R8TS)

{
  UINT4 i;
  CHAR txtNumber[32];
  *R4TS = XLALCreateREAL4TimeSeries(
    R8TS->name,
    &(R8TS->epoch),
    R8TS->f0,
    R8TS->deltaT,
    &R8TS->sampleUnits,
    R8TS->data->length
  );
  if(!*R4TS)
    {
    /* ERROR */
    }
  for (i=0;i<(*R4TS)->data->length;i++)
    /*    (*R4TS)->data->data[i]=((REAL4) R8TS->data->data[i]);*/
    {
      sprintf(txtNumber,"%3.6e",R8TS->data->data[i]);
      (*R4TS)->data->data[i]=atof(txtNumber);
    }
}
/* End LALappsCreateR4FromR8TimeSeries */

void LALappsPSD_Check(REAL8TimeSeries       *dataIn)
{
  UINT4 size;
  REAL4TimeSeries       *R4TimeSeries=NULL;
  REAL4FrequencySeries  *avgPSDR4=NULL;
  REAL4FFTPlan          *fftplanR4=NULL;
  REAL8FrequencySeries  *avgPSDR8=NULL;
  REAL8FFTPlan          *fftplanR8=NULL;
  LALStatus             status;
  const LIGOTimeGPS        gps_zero = LIGOTIMEGPSZERO;

  memset(&status, 0, sizeof(status));
  print_real8tseries(dataIn,"PreCast_REAL8_Tseries.diag");
  LALappsCreateR4FromR8TimeSeries(&status,&R4TimeSeries,dataIn);
  print_real4tseries(R4TimeSeries,"PostCast_REAL4_Tseries.diag");

  size=dataIn->data->length/2 +1;
  avgPSDR8 = XLALCreateREAL8FrequencySeries(
    "psd",
    &(gps_zero),
    0,
    1/(dataIn->deltaT*dataIn->data->length),
    &lalADCCountUnit,
    size
  );
  if(!avgPSDR8)
    {
    /* ERROR */
    }

  fftplanR8=XLALCreateForwardREAL8FFTPlan(dataIn->data->length,0);
  XLALREAL8PowerSpectrum(avgPSDR8->data,(dataIn->data),fftplanR8);

  print_real8fseries(avgPSDR8,"CastREAL4_REAL8_PSD.diag");

  avgPSDR4 = XLALCreateREAL4FrequencySeries(
					    "psd",
					    &(gps_zero),
					    0,
					    1/(dataIn->deltaT*dataIn->data->length),
					    &lalADCCountUnit,
					    size
					    );
  if(!avgPSDR4)
    {
    /* ERROR */
    }
  fftplanR4=XLALCreateForwardREAL4FFTPlan(R4TimeSeries->data->length,0);
  XLALREAL4PowerSpectrum(avgPSDR4->data,
			 R4TimeSeries->data,
			 fftplanR4);

 print_real4fseries(avgPSDR4,"CastREAL4_from_REAL8_PSD.diag");
 /****** FROM HERE ******* Edit above to REAL4 */
 
  XLALDestroyREAL4TimeSeries(R4TimeSeries);
  XLALDestroyREAL8FrequencySeries(avgPSDR8);
  XLALDestroyREAL4FrequencySeries(avgPSDR4);
 if (fftplanR8 != NULL)
   XLALDestroyREAL8FFTPlan(fftplanR8);

 if (fftplanR4 != NULL)
   XLALDestroyREAL4FFTPlan(fftplanR4);

}
/* LALappsPSD_CHECK */


/* TEMPORARY */
/* Non Compliant code taken from EPSearch.c */
void print_real4tseries(const REAL4TimeSeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALSPrintTimeSeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  LALStatus   status=blank_status;
  REAL8   timeT;
  size_t i;
  timeT = XLALGPSGetREAL8(&(fseries->epoch));
  if(fp) 
    {
      for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%3.6e\n", (i * fseries->deltaT)+timeT, fseries->data->data[i]);
      fclose(fp);
    }
#endif
}

void print_real8tseries(const REAL8TimeSeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALDPrintTimeSeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  LALStatus   status=blank_status;
  REAL8   timeT;
  size_t i;
  timeT = XLALGPSGetREAL8(&(fseries->epoch));
  if(fp) 
    {
      for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%3.15e\n", (i * fseries->deltaT)+timeT, fseries->data->data[i]);
      fclose(fp);
    }
#endif
}

void print_real4fseries(const REAL4FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%3.6e\n", (i * fseries->deltaF), fseries->data->data[i]);
    fclose(fp);
  }
#endif
}

void print_real8fseries(const REAL8FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
	fprintf(fp, "%f\t%3.15e\n", (i * fseries->deltaF), fseries->data->data[i]);
    fclose(fp);
  }
#endif
}

void print_complex8fseries(const COMPLEX8FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
      fprintf(fp, "%f\t%g\n", i * fseries->deltaF, sqrt(fseries->data->data[i].re * fseries->data->data[i].re + fseries->data->data[i].im * fseries->data->data[i].im));
    fclose(fp);
  }
#endif
}

void print_complex8_RandC_fseries(const COMPLEX8FrequencySeries *fseries, const char *file)
{
#if 0
  /* FIXME: why can't the linker find this function? */
  LALCPrintFrequencySeries(fseries, file);
#else
  FILE *fp = fopen(file, "w");
  size_t i;

  if(fp) {
    for(i = 0; i < fseries->data->length; i++)
      fprintf(fp, "%f\t%g\t%g\n", i * fseries->deltaF, 
	      fseries->data->data[i].re,
	      fseries->data->data[i].im);
    fclose(fp);
  }

#endif
}


void print_lalUnit(LALUnit unit,const char *file)
{
  FILE *fp = fopen(file,"w");
  CHARVector *unitString=NULL;
  LALStatus  status=blank_status;

  unitString=XLALCreateCHARVector(maxFilenameLength);
  XLALUnitAsString(unitString->data,unitString->length,&unit);
  if (fp)
    fprintf(fp,"%s\n",unitString->data);
  fclose(fp);
  XLALDestroyCHARVector(unitString);
}
/*
 * End diagnostic code
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
 
