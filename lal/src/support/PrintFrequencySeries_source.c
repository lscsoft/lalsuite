#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,FrequencySeries)
#define VTYPE CONCAT2(TYPE,Sequence)
#define FUNC CONCAT3(LAL,TYPECODE,PrintFrequencySeries)


void FUNC ( STYPE *series, const CHAR *filename ) 
 
{
  REAL8 f;
  TYPE *data;
  FILE *fp;
  UINT4 i;
  CHAR unitString[LALUnitTextSize];

  if (series==NULL) return;

  /* *(series->data) is a VTYPE */
  /* series->data->data is a pointer to TYPE */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->epoch.gpsSeconds && series->epoch.gpsNanoSeconds) {
    fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
            series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  }
  else {
    fprintf(fp,"# \n");
  }
  XLALUnitAsString(unitString, LALUnitTextSize, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString);
  fprintf(fp,HEADER);
  for ( i = 0; i < series->data->length; ++i )
  {
    f = series->f0 + i * series->deltaF;
    data = &(series->data->data[i]);
    fprintf(fp,FMT,f,ARG);
  }

  LALFclose(fp);

  return;
}
