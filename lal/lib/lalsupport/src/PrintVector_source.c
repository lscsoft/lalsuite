#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRINGx(a) #a
#define STRING(a) STRINGx(a)

#define VTYPE CONCAT2(TYPE,Vector)
#define FUNC CONCAT3(LAL,TYPECODE,PrintVector)


void FUNC ( VTYPE *vector ) 
{ 
  int i;
  static int filenum=0;
  FILE *fp;
  char fname[FILENAME_MAX];


  if (vector==NULL) return;

  /* open output file */
  snprintf(fname, FILENAME_MAX, "%sPrintVector.%03d", STRING(TYPECODE), filenum++);
  fp=LALFopen(fname,"w");

  for (i=0;i<(int)vector->length;i++)
    fprintf(fp,FMT,i,ARG);

  LALFclose(fp);

  return;
}
