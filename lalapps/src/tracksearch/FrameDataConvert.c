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
  *
  */

 /*
  * Author: Cristina Valeria Torres (LIGO Livingston)
  */

#include "tracksearch.h"
#include "tracksearchToolbox.h"
#include "unistd.h"
NRCSID( TRACKSEARCHC, "tracksearch $Id$");
RCSID( "tracksearch $Id$");
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


















int XLALFrGetREAL8FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, FrStream *stream)
{
REAL8TimeSeries *tmpData=NULL;
REAL8TimeSeries *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;

tmpData=XLALCreateREAL8TimeSeries (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
		     inputSeries->data->length);
errcode=XLALFrGetREAL8TimeSeriesMetadata (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"XLALFrGetREAL8TimeSeriesMetadata : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=XLALCreateREAL8TimeSeries (tmpData->name,
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
                      loadPoints);
if (tmpData)
   XLALDestroyREAL8TimeSeries (tmpData);
   
errcode=XLALFrGetREAL8TimeSeries (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrGetREAL8TimeSeries : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,20));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-20;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   XLALDestroyREAL8TimeSeries (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %f Hz to %f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}



















int XLALFrGetREAL4FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, FrStream *stream)
{
REAL4TimeSeries *tmpData=NULL;
REAL4TimeSeries *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;

tmpData=XLALCreateREAL4TimeSeries (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
		     inputSeries->data->length);
errcode=XLALFrGetREAL4TimeSeriesMetadata (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"XLALFrGetREAL4TimeSeriesMetadata : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=XLALCreateREAL4TimeSeries (tmpData->name,
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
                      loadPoints);
if (tmpData)
   XLALDestroyREAL4TimeSeries (tmpData);
   
errcode=XLALFrGetREAL4TimeSeries (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrGetREAL4TimeSeries : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,0));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-0;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   XLALDestroyREAL4TimeSeries (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %f Hz to %f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}



















int XLALFrGetINT2FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, FrStream *stream)
{
INT2TimeSeries *tmpData=NULL;
INT2TimeSeries *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;

tmpData=XLALCreateINT2TimeSeries (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
		     inputSeries->data->length);
errcode=XLALFrGetINT2TimeSeriesMetadata (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"XLALFrGetINT2TimeSeriesMetadata : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=XLALCreateINT2TimeSeries (tmpData->name,
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
                      loadPoints);
if (tmpData)
   XLALDestroyINT2TimeSeries (tmpData);
   
errcode=XLALFrGetINT2TimeSeries (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrGetINT2TimeSeries : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,0));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-0;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   XLALDestroyINT2TimeSeries (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %f Hz to %f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}



















int XLALFrGetINT4FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, FrStream *stream)
{
INT4TimeSeries *tmpData=NULL;
INT4TimeSeries *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;

tmpData=XLALCreateINT4TimeSeries (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
		     inputSeries->data->length);
errcode=XLALFrGetINT4TimeSeriesMetadata (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"XLALFrGetINT4TimeSeriesMetadata : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=XLALCreateINT4TimeSeries (tmpData->name,
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
                      loadPoints);
if (tmpData)
   XLALDestroyINT4TimeSeries (tmpData);
   
errcode=XLALFrGetINT4TimeSeries (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrGetINT4TimeSeries : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,0));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-0;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   XLALDestroyINT4TimeSeries (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %f Hz to %f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}



















int XLALFrGetINT8FrameConvertToREAL4TimeSeries (REAL4TimeSeries *inputSeries, FrStream *stream)
{
INT8TimeSeries *tmpData=NULL;
INT8TimeSeries *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;

tmpData=XLALCreateINT8TimeSeries (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
		     inputSeries->data->length);
errcode=XLALFrGetINT8TimeSeriesMetadata (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"XLALFrGetINT8TimeSeriesMetadata : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=XLALCreateINT8TimeSeries (tmpData->name,
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
                      loadPoints);
if (tmpData)
   XLALDestroyINT8TimeSeries (tmpData);
   
errcode=XLALFrGetINT8TimeSeries (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrGetINT8TimeSeries : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,0));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-0;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   XLALDestroyINT8TimeSeries (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %f Hz to %f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}

