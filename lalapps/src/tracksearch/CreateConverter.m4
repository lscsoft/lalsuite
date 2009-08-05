define(`SCALE',0)
ifelse(TYPECODE,`D', `define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S', `define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`',  `define(`TYPE',`REAL4')')
ifelse(TYPECODE,`D', `define(`SCALE',20)')

define(`XFUNC',`format(`XLALFrGet%sFrameConvertToREAL4TimeSeries',TYPE)')
define(`CREATESERIES',`format(`XLALCreate%sTimeSeries',TYPE)')
define(`GETMETA',`format(`XLALFrGet%sTimeSeriesMetaData',TYPE)')
define('GETDATA',`format(`XLALFrGet%sTimeSeriesData',TYPE)')
define(`DESTROYSERIES',`format(`XLALDestroy%sTimeSeries',TYPE)')
define(`VARTYPE',`format(`%sTimeSeries',TYPE)')

int XFUNC (REAL4TimeSeries *inputSeries FrStream *stream)
{
VARTYPE *tmpData=NULL;
VARTYPE *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

INT4  i=0;
INT4  errcode=0;
UINT4  loadPoints=0;
REAL8 originalSampleRate=0;

tmpData=CREATESERIES (inputSeries->name,
		     inputSeries->epoch,
		     inputSeries->deltaT,
		     inputSeries->sampleUnits,
		     inputSeries->data->length);
errcode=GETMETA (tmpData,stream);
if (errcode!=0)
   {
    fprintf(stderr,"GETMETA : Metadata on %s unavailable in stream.\n",tmpData->name);
    fflush(stderr);
    return errcode;
   }
loadPoints=(1/tmpData->deltaT)*(inputSeries->deltaT*inputSeries->data->length);
tmpData2=CREATESERIES (tmpData->name,
		      tmpData->epoch,
                      tmpData->deltaT,
		      tmpData->sampleUnits,
                      loadPoints);
if (tmpData)
   DESTROYSERIES (tmpData);
   
errcode=GETDATA (tmpData2,stream);
if (errcode!=0)
   {
   fprintf(stderr,"GETDATA : Error getting data from frame stream!\n");
   fflush(stderr);
   return errcode;
   }
errcode=XLALFrSeek(stream,tmpData2->epoch);
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
for (i=0;i<tmpData2->data->length;i++)
tmpData2->data->data[i]=tmpData2->data->data[i]*(1/pow(10,SCALE));

tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-SCALE;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   tmpData2->epoch,
				   tmpData2->deltaT,
				   tmpData2->sampleUnits,
				   tmpData2->data->length);
for (i=0;i<tmpData3->data->length;i++)
    tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]);

if (tmpData2)
   DESTROYSERIES (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampleing from %f Hz to %f Hz\n.",(1/tmpData2->deltaT,(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData2->data->data[i];

return 0;
}
