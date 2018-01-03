#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,TimeSeries)
#define VTYPE CONCAT2(TYPE,Sequence)
#define FUNC CONCAT3(LAL,TYPECODE,PrintTimeSeries)


void FUNC ( STYPE *series, const CHAR *filename ) 
{ 
  REAL8 t;
  TYPE *data;
  FILE *fp;
  UINT4 i;
  CHAR unitString[LALUnitTextSize];

  if (series==NULL) return;

  /* *(series->data) is a VTYPE */
  /* series->data->data is a pointer to TYPE */

  /* Make a TYPE pointer which points to the first memory address not
   * belonging to the sequence
   */

  /* open output file */
  fp=LALFopen(filename,"w");
  fprintf(fp,"# %s\n",series->name);
  if (series->f0) {
     fprintf(fp,"# Heterodyned at %g Hz\n",series->f0);
  }
  else {
    fprintf(fp,"# \n");
  }
  fprintf(fp,"# Epoch is %d seconds, %d nanoseconds\n",
          series->epoch.gpsSeconds,series->epoch.gpsNanoSeconds);
  XLALUnitAsString(unitString, LALUnitTextSize, &(series->sampleUnits));
  fprintf(fp,"# Units are (%s)\n",unitString);
  fprintf(fp,HEADER);
  for ( i = 0; i < series->data->length; ++i )
  {
    t = i * series->deltaT;
    data = &(series->data->data[i]);
    fprintf(fp,FMT,t,ARG);
  }	

  LALFclose(fp);

  return;
}
