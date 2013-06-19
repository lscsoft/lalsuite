#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define XFUNC CONCAT3(XLALFrGet,TYPE,FrameConvertToREAL4TimeSeries)
#define CREATESERIES CONCAT3(XLALCreate,TYPE,TimeSeries)
#define GETMETA CONCAT3(XLALFrStreamGet,TYPE,TimeSeriesMetadata)
#define GETDATA CONCAT3(XLALFrStreamGet,TYPE,TimeSeries)
#define DESTROYSERIES CONCAT3(XLALDestroy,TYPE,TimeSeries)
#define VARTYPE CONCAT2(TYPE,TimeSeries)

int XFUNC (REAL4TimeSeries *inputSeries, LALFrStream *stream)
{
VARTYPE *tmpData=NULL;
VARTYPE *tmpData2=NULL;
REAL4TimeSeries *tmpData3=NULL;

UINT4  i=0;
INT4   errcode=0;
UINT4  loadPoints=0;
REAL8  factor=1;

tmpData=CREATESERIES (inputSeries->name,
		     &(inputSeries->epoch),
		     0,
		     inputSeries->deltaT,
		     &(inputSeries->sampleUnits),
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
		      &(tmpData->epoch),
		      0,
                      tmpData->deltaT,
		      &(tmpData->sampleUnits),
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
errcode=XLALFrStreamSeek(stream,&(tmpData2->epoch));
if (errcode!=0)
   {
   fprintf(stderr,"XLALFrStreamSeek : Can not seek frame to start epoch?\n");
   fflush(stderr);
   return errcode;
   }
if (SCALE > 0)
  tmpData2->sampleUnits.powerOfTen=tmpData2->sampleUnits.powerOfTen-SCALE;

tmpData3=XLALCreateREAL4TimeSeries (tmpData2->name,
				   &(tmpData2->epoch),
				   0,
				   tmpData2->deltaT,
				   &(tmpData2->sampleUnits),
				   tmpData2->data->length);
factor=(1/pow(10,SCALE));
for (i=0;i<tmpData3->data->length;i++)
  tmpData3->data->data[i]=((REAL4) tmpData2->data->data[i]*factor);

if (tmpData2)
   DESTROYSERIES (tmpData2);

errcode=XLALResampleREAL4TimeSeries (tmpData3,inputSeries->deltaT);
if (errcode!=0)
   {
   fprintf(stderr,"XLALResampleREAL4TimeSeries : Error resampling from %6.2f Hz to %6.2f Hz\n.",(1/tmpData2->deltaT),(1/inputSeries->deltaT));
   fflush(stderr);
   return errcode;
   }

for (i=0;i<inputSeries->data->length;i++)
    inputSeries->data->data[i]=tmpData3->data->data[i];

return 0;
}
