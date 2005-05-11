/************************************************************************************/
/* This code reads in results data from a search and reads in results from          */
/* injecting many signals into the data.  By comparing the loudest event in the     */
/* real results with the distribution of events using injections we can calculate   */
/* an upperlimit confidence as a function of injected signal amplitude.             */
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2005                            */
/************************************************************************************/

#include "ComputeUL_v1.h"
#include "GenerateBinaryMesh_v1.h"
#include "FindCoincidence_v1.h" 


/* clargs */
char outfile[256],freqmeshfile[256],h0file[256],maxfile[256],coresultsdir[256],injectiondir[256];
REAL8 conf;

/* needs to be globally visible */
REAL8 *tempULdataconf;

int ReadCommandLine(int argc,char *argv[]);
int ReadFreqMeshFile(char *,FreqMeshes **);
int Readh0File(char *,INT4 *,REAL8 **);
int FindLoudest(char *,FreqMeshes *,char *,Loudest **);
int ComputeConfidence(char *,FreqMesh,REAL8,REAL8,REAL8 *,REAL8 *);
int ExtrapolateUL(ULData,REAL8,REAL8 **);
int OutputConfidence(char *,ULData,REAL8,REAL8 *);
int compare(const void *, const void *);

int main(int argc, char **argv){
  
  INT4 i,j;
  FreqMeshes *freqmeshdata;
  INT4 Nh0;
  REAL8 *h0data;
  Loudest *loudest;
  REAL8 *h0conf=NULL;
  ULData ULdata;
  REAL8 temp_conf,temp_conf_err;

   /* read the command line arguments */
  if (ReadCommandLine(argc,argv)) return 1;

  /* read in freq-mesh file */
  if (ReadFreqMeshFile(freqmeshfile,&freqmeshdata)) return 2;

  /* Read in the h0 values used for the injections */
  if (Readh0File(h0file,&Nh0,&h0data)) return 3;

  /* find the loudest events on a band by band basis */
  if (FindLoudest(coresultsdir,freqmeshdata,maxfile,&loudest)) return 4;
      
  /* allocate memory to the UL data structure */
  ULdata.UL=(SingleULData *)LALMalloc(freqmeshdata->Nheaders*sizeof(SingleULData));
  ULdata.N_h0=Nh0;
  ULdata.N_UL=freqmeshdata->Nheaders;
  for (i=0;i<freqmeshdata->Nheaders;i++) {
    ULdata.UL[i].h0=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    ULdata.UL[i].conf=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    ULdata.UL[i].conf_err=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
  }
  if (Nh0>3) h0conf=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));

  /* loop over each upper-limit band */
  for (i=0;i<freqmeshdata->Nheaders;i++) {
    
    /* fill in the UL data */
    ULdata.UL[i].f_min=freqmeshdata->freqmesh[i].f_min;
    ULdata.UL[i].f_max=freqmeshdata->freqmesh[i].f_max;
    ULdata.UL[i].p_loudestsig=loudest[i].p_sig;
    ULdata.UL[i].s_loudestsig=loudest[i].s_sig;
    ULdata.UL[i].co_loudestsig=loudest[i].co_sig;

    /* loop over the h0 values */
    for (j=0;j<Nh0;j++) {
      
      /* fill in the UL data */
      ULdata.UL[i].h0[j]=h0data[j];

      /* find confidence */
      if (ComputeConfidence(injectiondir,freqmeshdata->freqmesh[i],ULdata.UL[i].co_loudestsig,h0data[j],&temp_conf,&temp_conf_err)) return 5;
      ULdata.UL[i].conf[j]=temp_conf;
      ULdata.UL[i].conf_err[j]=temp_conf_err;

    }

  }

  /* extrapolate the requested confidence from the descrete values */
  if (h0conf!=NULL) {

    if (ExtrapolateUL(ULdata,conf,&h0conf)) return 6;

  }


    /* output UL results to file */
  if (OutputConfidence(outfile,ULdata,conf,h0conf)) return 7;

  /* free all memory */
  

  return 0;
  
}

/******************************************************************************/

int ExtrapolateUL(ULData ULdata,REAL8 confidence, REAL8 **h0conf)
{

  /* this function finds the three closest points in confidence to the requested */
  /* confidence and fits a quadratic to the points to estimate the location of */
  /* the requested confidence */

  INT4 i,j,q;
  REAL8 *diff;
  REAL8 c1,c2,c3,h1,h2,h3,r,s,t;
  INT4 *indexes;

  /* allocate space for diff */
  tempULdataconf=(REAL8 *)LALMalloc(ULdata.N_h0*sizeof(REAL8));

  /* create an array of indexes */
  indexes=(INT4 *)LALMalloc(sizeof(INT4)*ULdata.N_h0);
  

  /* loop over the frequency bands */
  for (i=0;i<ULdata.N_UL;i++) {
    
    /* first find the three closest points to the requested confidence */
    for (j=0;j<ULdata.N_h0;j++) {
      
      tempULdataconf[j]=fabs(confidence-ULdata.UL[i].conf[j]);

    }

    /* populate it */
    for (q=0;q<ULdata.N_h0;q++) {
      indexes[q]=q;
    }
    
    /* sort array of indexes */
    qsort((void *)indexes, (size_t)ULdata.N_h0, sizeof(INT4), compare);

    /* need to pick three points that are close and either side of the target in confidence */
    c1=ULdata.UL[i].conf[indexes[0]];
    c2=ULdata.UL[i].conf[indexes[1]];
    c3=-1.0;
    h1=log10(ULdata.UL[i].h0[indexes[0]]);
    h2=log10(ULdata.UL[i].h0[indexes[1]]);
    h3=-1.0;

    if ((c1<conf)&&(c2<conf)) {

      q=0;
      while (q<ULdata.N_h0) {
	if (ULdata.UL[i].conf[indexes[q]]>conf) {
	  c3=ULdata.UL[i].conf[indexes[q]];
	  h3=log10(ULdata.UL[i].h0[indexes[q]]);
	  q=ULdata.N_h0+1;
	}
	q++;
      }

    }
    else if ((c1>conf)&&(c2>conf)) {

      q=0;
      while (q<ULdata.N_h0) {
	if (ULdata.UL[i].conf[indexes[q]]<conf) {
	  c3=ULdata.UL[i].conf[indexes[q]];
	  h3=log10(ULdata.UL[i].h0[indexes[q]]);
	  q=ULdata.N_h0+1;
	}
	q++;
      }

    }
    else if (((c1>conf)&&(c2<conf))||((c1<conf)&&(c2>conf))){

      c3=ULdata.UL[i].conf[indexes[3]];
      h3=log10(ULdata.UL[i].h0[indexes[3]]);

    }
   
    if (c3>0.0) {

      
      /*printf("three closest that span conf are %f %e\n",c1,h1);
      printf("three closest that span conf are %f %e\n",c2,h2);
      printf("three closest that span conf are %f %e\n",c3,h3);*/
      
            
      s=( ((c2-c3)*(c2+c3)*(h1-h2)) - ((h2-h3)*(c1-c2)*(c1+c2)) ) / ( ((c1-c2)*(c2-c3)*(c2+c3)) - ((c2-c3)*(c1-c2)*(c1+c2)) );
      r=( (h2-h3) - (s*(c2-c3)) ) / ((c2-c3)*(c2+c3));
      t=h1-(r*c1*c1)-s*c1;
      
            
      (*h0conf)[i]=pow(10,r*conf*conf+s*conf+t);
      
      /*printf("done the conf it is %e\n",(*h0conf)[i]);*/
   
    }
    else {
      
      (*h0conf)[i]=0.0;
      
    }
     
  }
 
  LALFree(tempULdataconf);
  LALFree(indexes);

  

  
  return 0;

}
/*******************************************************************************/

/* Sorting function to sort into DECREASING order */
int compare(const void *ip, const void *jp)
{
  
  REAL8 di,dj;

  di=tempULdataconf[*(const int *)ip];
  dj=tempULdataconf[*(const int *)jp];

  
  if (di>dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/******************************************************************************/

int ComputeConfidence(char *inputdir,FreqMesh freqmesh,REAL8 loudsig,REAL8 h0,REAL8 *confidence,REAL8 *conf_err)
{

  /* this function takes a frequency band and value of h0 and computes the */
  /* upperlimit confidence */

  char filename[512];
  char line[1024];
  FILE *fp;
  INT4 tot,count,i;
  INT4 fileno=0;
  REAL8 sig;
  char **filelist;
  glob_t globbuf;
  INT4 nfiles;
  REAL8 a;
  INT4 A;
  REAL8 bum,arse;

  /* make the template filename for searching */
  sprintf(filename,"%s/injection_%.3f-%.3f_%.3e_*.data",inputdir,freqmesh.f_min,freqmesh.f_max,h0);
      
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(filename, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any files of this name in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo files in directory %s ... Exiting.\n", inputdir);
      *confidence=-1.0;
      *conf_err=-1.0;
      return 0;
    }
  
  /* allocate memory for the pathnames */
  filelist=(char **)LALMalloc(globbuf.gl_pathc*sizeof(char *));
  for (i=0;i<(INT4)globbuf.gl_pathc;i++) filelist[i]=(char *)LALMalloc(256*sizeof(char));
  
  /* read all file names into memory */
  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"\nToo many files in directory! Exiting... \n");
	  exit(1);
	}
    }
  globfree(&globbuf);
  nfiles=fileno;
  
 
  tot=0;
  count=0;
  /* loop over these files */
  for (i=0;i<nfiles;i++) {
    
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : cannot open the file %s\n",filelist[i]);
      exit(1);
    }
    
    while (fgets(line,1023,fp)!=NULL) {
    

      sscanf(line,"%lf %lf%lf%lf%lf%lf %d%d %lf%lf %d %lf%lf%lf%lf %lf%lf%lf%lf%lf%d%d%lf%lf%d%lf%lf%lf%lf%lf",&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&sig);  
      /* printf("sig is %f\n",sig); */
      if (sig<=loudsig) {
	/* printf("sig read as %f loudest is %f -> counting as louder !!! [count = %d]\n",sig,loudsig,count); */
	count++;
      }
      /* else {
	printf("sig read as %f loudest is %f -> counting as quieter !!!\n",sig,loudsig);
	}*/
      tot++;
    }
   
    fclose(fp);
    
  }
  
  *confidence=(REAL8)((REAL8)count/(REAL8)tot);
  *conf_err=(REAL8)(1.0/sqrt((REAL8)tot));
  
  
  return 0;

}

/******************************************************************************/

int ReadFreqMeshFile(char *filename, FreqMeshes **freqmesh)
{

  /* this function reads in the freq-mesh info file that tells the code */
  /* which template bank was used for which frequency band */

  FILE *fp,*fp_pmesh,*fp_smesh;
  INT4 i;
  REAL8 min_f,max_f,band;
  char p_mesh[256],s_mesh[256];
  BinaryMeshFileHeader *test;

  fp=fopen(filename,"r");
  if (fp==NULL) {
    printf("ERROR : cannot open file %s\n",filename);
    exit(1);
  }

  i=0;
  while (fscanf(fp,"%lf %lf %lf %s %s\n",&min_f,&max_f,&band,p_mesh,s_mesh)!=EOF) {
    i++;
  }
  fclose(fp);
  

  fp=fopen(filename,"r");
  if (fp==NULL) {
    printf("ERROR : cannot open file %s\n",filename);
    exit(1);
  }

  /* allocate some memory */
  (*freqmesh)=(FreqMeshes *)LALMalloc(sizeof(FreqMeshes));
  (*freqmesh)->freqmesh=(FreqMesh *)LALMalloc(i*sizeof(FreqMesh));

  (*freqmesh)->Nheaders=i;


  i=0;

  while (fscanf(fp,"%lf %lf %lf %s %s",&min_f,&max_f,&band,p_mesh,s_mesh)!=EOF) {
    

    /* open up mesh files */
    fp_pmesh=fopen(p_mesh,"r");
    if (fp_pmesh==NULL) {
      printf("ERROR : could not open mesh file %s\n",p_mesh);
      exit(1);
    }
    fp_smesh=fopen(s_mesh,"r");
    if (fp_smesh==NULL) {
      printf("ERROR : could not open mesh file %s\n",s_mesh);
      exit(1);
    }


    (*freqmesh)->freqmesh[i].f_min=min_f;
    (*freqmesh)->freqmesh[i].f_max=max_f;
    (*freqmesh)->freqmesh[i].f_band=band;


    test = &((*freqmesh)->freqmesh[i].p_header);

    if (ReadMeshFileHeader(fp_pmesh,test)) return 1;

    if (ReadMeshFileHeader(fp_smesh,&((*freqmesh)->freqmesh[i]).s_header)) return 1;
    


    fclose(fp_pmesh);
    fclose(fp_smesh);


    i++;
  }

  fclose(fp);


  return 0;

}

/****************************************************************************/

int Readh0File(char *filename, INT4 *N, REAL8 **h0)
{

  /* this function reads in the values of h0 from a file */
 
  FILE *fp;
  INT4 i;
  REAL8 dummy;

  fp=fopen(filename,"r");
  if (fp==NULL) {
    printf("ERROR : cannot open file %s\n",filename);
    exit(1);
  }

  i=0;
  while (fscanf(fp,"%lf\n",&dummy)!=EOF) i++;
  fclose(fp);
  (*N)=i;

  /* allocate some memory */
  (*h0)=(REAL8 *)LALMalloc((*N)*sizeof(REAL8));

  i=0;
  fp=fopen(filename,"r");
  while (fscanf(fp,"%lf\n",&dummy)!=EOF) {
    (*h0)[i]=dummy;
    i++;
  }

  fclose(fp);


  return 0;

}



/****************************************************************************/

int FindLoudest(char *resultsdir,FreqMeshes *freqmeshes,char *maxoutfile,Loudest **loudest)
{

  INT4 fileno=0;
  FILE *fp;
  char **filelist;
  char command[512];
  char line[512];
  char **loudestline;
  glob_t globbuf;
  INT4 i,j;
  INT4 nfiles;
  REAL8 a,b,c,d,e,f,g,h,l,m,n,o,p,q,r,s,t,u;
  INT4 A,B,C,D;
  REAL8 p_sig,s_sig,co_sig;
  FILE *fpmax;

  /* set up the datadir name */
  strcpy(command, resultsdir);
  strcat(command,"/coincidence_*");
    
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(command, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any SFT's in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo files in directory %s ... Exiting.\n", resultsdir);
      exit(1);
    }
  
  /* allocate memory for the pathnames */
  filelist=(char **)LALMalloc(globbuf.gl_pathc*sizeof(char *));
  for (i=0;i<(INT4)globbuf.gl_pathc;i++) filelist[i]=(char *)LALMalloc(256*sizeof(char));
  
  /* read all file names into memory */
  while ((UINT4)fileno < globbuf.gl_pathc) 
    {
      strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);
      fileno++;
      if (fileno > MAXFILES)
	{
	  fprintf(stderr,"\nToo many files in directory! Exiting... \n");
	  exit(1);
	}
    }
  globfree(&globbuf);
  nfiles=fileno;

  /* allocate memory for loudest events */
  (*loudest)=(Loudest *)LALMalloc(freqmeshes->Nheaders*sizeof(Loudest));
  loudestline=(char **)LALMalloc(freqmeshes->Nheaders*sizeof(char *));
  for (i=0;i<freqmeshes->Nheaders;i++) loudestline[i]=(char *)LALMalloc(1024*sizeof(char));

  for (j=0;j<freqmeshes->Nheaders;j++) (*loudest)[j].co_sig=999999999;
  
  /* so lets loop over every file in the results directory */
  for (i=0;i<fileno;i++) {

    /* open this file */
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : Strange can't open file %s\n",filelist[i]);
      exit(1);
    }
    
    /* loop over each result */
    while (fgets(line,511,fp)!=NULL) {
    
      sscanf(line,"%lf%lf%lf%lf%lf %d%d %lf%lf %d %lf%lf%lf%lf %lf%lf%lf%lf%lf %d%d %lf%lf %d %lf%lf%lf%lf%lf\n",
	     &f,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&p_sig,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&s_sig,&co_sig);

      /* now we need to loop over each upperlimit injection band */
      for (j=0;j<freqmeshes->Nheaders;j++) {


	/* if the result is in the right band */
	if ((f>=freqmeshes->freqmesh[j].f_min)&&(f<freqmeshes->freqmesh[j].f_max)) {
	  

	  /* if the result is louder than the previous loudest */
	  if (co_sig<(*loudest)[j].co_sig) {
	    (*loudest)[j].co_sig=co_sig;
	    (*loudest)[j].p_sig=p_sig;
	    (*loudest)[j].s_sig=s_sig;
	    strcpy(loudestline[j],line);
	    
	  }

	}

      }

    }
      
    fclose(fp);

  }
  
  /* open the max out file */
  fpmax=fopen(maxoutfile,"w");
  if (fpmax==NULL) {
    printf("ERROR : Strange can't open file %s\n",maxoutfile);
    exit(1);
  }
  
  /* output loudest events info */
  for (j=0;j<freqmeshes->Nheaders;j++) {
    fprintf(fpmax,"%s\n",loudestline[j]);
  }
  
  fclose(fpmax);
  

  return 0;  
}

/*********************************************************************************************/

 int OutputConfidence(char *outputfile,ULData ULdata,REAL8 confidence,REAL8 *h0conf) 
{

  /* this function outputs the final upperlimit results to a single file */

  FILE *fp;
  INT4 i,j;
  char outputfilename[256];

  /* output a series of files for plotting */
  for (i=0;i<ULdata.N_UL;i++) {
    
    sprintf(outputfilename,"%s_%.3f-%.3f.data",outputfile,ULdata.UL[i].f_min,ULdata.UL[i].f_max);
    fp=fopen(outputfilename,"w");
    if (fp==NULL) {
      printf("ERROR : cannot open file %s\n",outputfilename);
      exit(1);
    }

    fprintf(fp,"*** Upperlimit results file for freq = [%.3f - %.3f] ***\n***\n",ULdata.UL[i].f_min,ULdata.UL[i].f_max);
    for (j=0;j<ULdata.N_h0;j++) fprintf(fp,"%e\t%f\t%f\n",ULdata.UL[i].h0[j],ULdata.UL[i].conf[j],ULdata.UL[i].conf_err[j]);

    fclose(fp);

}

  /* output main full upperlimit results file */
  sprintf(outputfilename,"%s_FULL.data",outputfile);
  fp=fopen(outputfilename,"w");
  if (fp==NULL) {
    printf("ERROR : cannot open file %s\n",outputfilename);
    exit(1);
  }

  fprintf(fp,"*** Upperlimit results file ***\n***\n");
  fprintf(fp,"injection results = %s\n",injectiondir);
  fprintf(fp,"coincidence results = %s\n",coresultsdir);
  fprintf(fp,"h0 injection values file = %s\n***\n",h0file);
  fprintf(fp,"f_min\t\tf_max\t\tp_loudest_log10(1-sig)\ts_loudest_log10(1-sig)\tco_loudest_log10(1-sig)\t");
  for (i=0;i<ULdata.N_h0;i++) fprintf(fp,"%e\t",ULdata.UL[0].h0[i]);
  fprintf(fp,"UL_(%f)\n--------------------------------------------------------------------------------------------------------------------\n",confidence);
  
  /* loop over each band */
  for (i=0;i<ULdata.N_UL;i++) {
    fprintf(fp,"%f\t%f\t%f\t\t%f\t\t%f\t\t",ULdata.UL[i].f_min,ULdata.UL[i].f_max,ULdata.UL[i].p_loudestsig,ULdata.UL[i].s_loudestsig,ULdata.UL[i].co_loudestsig);
    for (j=0;j<ULdata.N_h0;j++) fprintf(fp,"%f\t",ULdata.UL[i].conf[j]);
    fprintf(fp,"%e\n",h0conf[i]);
  }

  fclose(fp);

  /* output to file ready for plotting freq vs UL */
  /* output main full upperlimit results file */
  sprintf(outputfilename,"%s_freq_vs_UL.data",outputfile);
  fp=fopen(outputfilename,"w");
  if (fp==NULL) {
    printf("ERROR : cannot open file %s\n",outputfilename);
    exit(1);
  }

  fprintf(fp,"*** Upperlimit results file for freq vs UL ***\n***\n");
  fprintf(fp,"%.3f %e\n",ULdata.UL[0].f_max,h0conf[0]);
  for (i=1;i<ULdata.N_UL;i++) {
    fprintf(fp,"%.3f %e\n",ULdata.UL[i-1].f_min,h0conf[i-1]);
    fprintf(fp,"%.3f %e\n",ULdata.UL[i].f_max,h0conf[i]);
  }
  fprintf(fp,"%.3f %e\n",ULdata.UL[ULdata.N_UL-1].f_min,h0conf[ULdata.N_UL-1]);
  
  fclose(fp);

  return 0;

}

/*********************************************************************************************/

 int ReadCommandLine(int argc,char *argv[]) 
{
  INT4 c, errflg = 0;
  CHAR *temp;
  optarg = NULL;
  
  /* Initialize default values */
  sprintf(injectiondir," ");
  sprintf(coresultsdir," ");
  sprintf(freqmeshfile," ");
  sprintf(outfile," ");
  sprintf(maxfile," ");
  conf=0.95;

   
  {
    int option_index = 0;
    static struct option long_options[] = {
      {"injectiondir", required_argument, 0, 'i'},
      {"coresultsdir", required_argument, 0, 'c'},
      {"freqmeshfile", required_argument, 0, 'f'},
      {"gwamplitudefile", required_argument, 0, 'H'},
      {"uloutfile", required_argument, 0, 'U'},
      {"maxoutfile", required_argument, 0, 'x'},
      {"confidence", required_argument, 0, 'C'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hi:c:f:H:U:x:C:",long_options, &option_index)))!=-1)
      switch (c) {
      case 'i':
	temp=optarg;
	sprintf(injectiondir,temp);
	break;
      case 'c':
	temp=optarg;
	sprintf(coresultsdir,temp);
	break;
      case 'f':
	temp=optarg;
	sprintf(freqmeshfile,temp);
	break;
      case 'H':
	temp=optarg;
	sprintf(h0file,temp);
	break;
      case 'U':
	temp=optarg;
	sprintf(outfile,temp);
	break;
      case 'x':
	temp=optarg;
	sprintf(maxfile,temp);
	break;
      case 'C':
	conf=atof(optarg);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--injectiondir       STRING\t Name of the directory where injection data is stored [DEFAULT= ]\n");
	fprintf(stdout,"\t--coresultsdir       STRING\t Name of the directory where coincident results are stored [DEFAULT= ]\n");
	fprintf(stdout,"\t--freqmeshfile       STRING\t Name of the file containing freq-mesh information [DEFAULT= ]\n");
	fprintf(stdout,"\t--gwamplpitudefile   STRING\t Name of the file containing h0 injection values [DEFAULT = ]\n");
	fprintf(stdout,"\t--uloutfile          STRING\t Name of the upperlimit output file [DEFAULT = ]\n");
	fprintf(stdout,"\t--maxoutfile         STRING\t Name of the loudest event output file [DEFAULT = 0.0]\n");
	fprintf(stdout,"\t--confidence         REAL8\t Confidence for which h0_UL(conf) will be calculated [DEFAULT = 0.0]\n");
	exit(0);
	break;
      default:
	/* unrecognized option */
	errflg++;
	fprintf(stderr,"Unrecognized option argument %c\n",c);
	exit(1);
	break;
      }
  }
  /* Need to add some CLA error checking here */
  
  /* update global variable and return */
  return errflg;
}

/************************************************************************/

