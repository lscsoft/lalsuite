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

static LALStatus status;
/* clargs */
char outfile[256],freqmeshfile[256],h0file[256],maxfile[256],coresultsdir[256],injectiondir[256];
REAL8 conf;
INT4 nsplit;

/* needs to be globally visible */
REAL8 *tempULdataconf_lower;
REAL8 *tempULdataconf_higher;

int ReadCommandLine(int argc,char *argv[]);
int ReadFreqMeshFile(char *,FreqMeshes **);
int Readh0File(char *,INT4 *,REAL8 **);
int FindLoudest(char *,FreqMeshes *,char *,Loudest **);
int ComputeConfidence(REAL8Vector *,REAL8,REAL8 *);
int ExtrapolateUL(ULData,REAL8,REAL8 **);
int compare_lower(const void *, const void *);
int compare_higher(const void *, const void *);
int ReadInjections(char *,FreqMesh,REAL8,REAL8Vector **);
int Interpolateh0(SingleSplitULData *,INT4,REAL8 *,REAL8);
int CalculateStats(SingleULData *,INT4);
int OutputConfidence(char *,ULData,REAL8); 

int main(int argc, char **argv){
  
  INT4 i,j,k=0,q;
  FreqMeshes *freqmeshdata;
  INT4 Nh0;
  INT4 start;
  REAL8 *h0data;
  Loudest *loudest;
  ULData ULdata;
  REAL8Vector *results=NULL;
  REAL8Vector *tempresults=NULL;
  REAL8 confsum;

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
  ULdata.N_split=nsplit;
  ULdata.N_UL=freqmeshdata->Nheaders;
  for (i=0;i<freqmeshdata->Nheaders;i++) {
    /* allocate memory for each subgroup of injections */
    ULdata.UL[i].SSULdata=(SingleSplitULData *)LALMalloc(ULdata.N_split*sizeof(SingleSplitULData));
    /* allocate memory for each final confidence */
    ULdata.UL[i].final_conf=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    ULdata.UL[i].final_conf_err=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    ULdata.UL[i].h0=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    for (j=0;j<ULdata.N_split;j++) {
      /* allocate memory for each value of conf */
      ULdata.UL[i].SSULdata[j].conf=(REAL8 *)LALMalloc(Nh0*sizeof(REAL8));
    }
  }
  printf("*** allocated memory for the new UL structure\n");
 
  /* loop over each upper-limit band */
  for (i=0;i<freqmeshdata->Nheaders;i++) {
    
    /* fill in the UL data */
    ULdata.UL[i].f_min=freqmeshdata->freqmesh[i].f_min;
    ULdata.UL[i].f_max=freqmeshdata->freqmesh[i].f_max;
    ULdata.UL[i].p_loudestsig=loudest[i].p_sig;
    ULdata.UL[i].s_loudestsig=loudest[i].s_sig;
    ULdata.UL[i].co_loudestsig=loudest[i].co_sig;
    confsum=0.0;

    /* loop over the chunks that we have split the dataset into */
    for (j=0;j<Nh0;j++) {
      
      /* read in injections results into memory */
      if (ReadInjections(injectiondir,freqmeshdata->freqmesh[i],h0data[j],&results)) return 4;
      /*printf("*** read injections\n");*/

      /* calculate the value of N_per_split */
      ULdata.UL[i].N_per_split=floor(results->length/(REAL8)ULdata.N_split);
      /* printf("N_per_slit is %d\n",ULdata.UL[i].N_per_split);*/

      /* fill in the UL data */
      ULdata.UL[i].h0[j]=h0data[j];

      /* loop over the h0 values */
      for (k=0;k<ULdata.N_split;k++) {
   
	/* extract the required chunk from the results */
	start=k*ULdata.UL[i].N_per_split;
	LALDCreateVector(&status,&tempresults,ULdata.UL[i].N_per_split);
	for (q=0;q<ULdata.UL[i].N_per_split;q++) {
	  tempresults->data[q]=results->data[q+start];
	}
	
	/* find confidence */
	if (ComputeConfidence(tempresults,ULdata.UL[i].co_loudestsig,&ULdata.UL[i].SSULdata[k].conf[j])) return 5;
	/* free vector */
	LALDDestroyVector(&status,&tempresults);
	
      } /* end of loop over N_split */

    } /* end loop over h0 */

  } /* end loop over frequency bands */

  /* loop over frequency again */
  for (i=0;i<ULdata.N_UL;i++) {
    
    /* initialise the confidence results */
    for (k=0;k<Nh0;k++) {
      ULdata.UL[i].final_conf[k]=0.0;
      ULdata.UL[i].final_conf_err[k]=0.0;
    }

    /* loop over N_split */
    for (j=0;j<ULdata.N_split;j++) {

      /* compute interpolated best h0 value for confidence */
      if (Interpolateh0(&ULdata.UL[i].SSULdata[j],Nh0,ULdata.UL[i].h0,conf)) return 8;
      /* printf("interpolated result = %e\n",ULdata.UL[i].SSULdata[j].est_h0); */

      /* loop over h0 values */
      for (k=0;k<Nh0;k++) {
	
	/* compute total confidences */
	ULdata.UL[i].final_conf[k]+=ULdata.UL[i].SSULdata[j].conf[k]/ULdata.N_split;
	ULdata.UL[i].final_conf_err[k]=1.0/sqrt(ULdata.N_split*ULdata.UL[i].N_per_split);
	/*printf("final conf = %f +/- %f\n",ULdata.UL[i].SSULdata[j].final_conf,ULdata.UL[i].SSULdata[j].final_conf_err);*/

      }

    } /* end loop over N_split */

    /* calculate the average interpolated result + errors */
    if (CalculateStats(&ULdata.UL[i],ULdata.N_split)) return 9;

    
  } /* end loop over frequency bands */

  /* output the results */
  if (OutputConfidence(outfile,ULdata,conf)) return 10;

  /* free the memory */
  for (i=0;i<ULdata.N_UL;i++) {
    for(j=0;j<ULdata.N_split;j++) LALFree(ULdata.UL[i].SSULdata[j].conf);
    LALFree(ULdata.UL[i].SSULdata);
    LALFree(ULdata.UL[i].final_conf);
    LALFree(ULdata.UL[i].final_conf_err);
    LALFree(ULdata.UL[i].h0);
  }
  LALFree(ULdata.UL);

  return 0;
  
}

/*********************************************************************************************/
int OutputConfidence(char *outputfile,ULData ULdata,REAL8 confidence) 
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
    printf("*** outputting to %s\n",outputfilename);

    fprintf(fp,"*** Upperlimit results file for freq = [%.3f - %.3f] ***\n***\n",ULdata.UL[i].f_min,ULdata.UL[i].f_max);
    printf("*** Upperlimit results file for freq = [%.3f - %.3f] ***\n***\n",ULdata.UL[i].f_min,ULdata.UL[i].f_max);
    for (j=0;j<ULdata.N_h0;j++) {
      fprintf(fp,"%e\t%f\t%f\n",ULdata.UL[i].h0[j],ULdata.UL[i].final_conf[j],ULdata.UL[i].final_conf_err[j]);
      printf("%e\t%f\t%f\n",ULdata.UL[i].h0[j],ULdata.UL[i].final_conf[j],ULdata.UL[i].final_conf_err[j]);
    }

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
  for (i=0;i<ULdata.N_h0;i++) fprintf(fp,"%e\t\t",ULdata.UL[0].h0[i]);
  fprintf(fp,"UL_(%f)\t\tUL_(%f)_err\n--------------------------------------------------------------------------------------------------------------------\n",confidence,confidence);
  
  /* loop over each band */
  for (i=0;i<ULdata.N_UL;i++) {
    fprintf(fp,"%f\t%f\t%f\t\t%f\t\t%f\t\t",ULdata.UL[i].f_min,ULdata.UL[i].f_max,ULdata.UL[i].p_loudestsig,ULdata.UL[i].s_loudestsig,ULdata.UL[i].co_loudestsig);
    for (j=0;j<ULdata.N_h0;j++) fprintf(fp,"%f\t\t",ULdata.UL[i].final_conf[j]);
    fprintf(fp,"%e\t\t%e\n",ULdata.UL[i].final_h0,ULdata.UL[i].final_h0_err);
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
  fprintf(fp,"%.3f %e\n",ULdata.UL[0].f_max,ULdata.UL[0].final_h0);
  for (i=1;i<ULdata.N_UL;i++) {
    fprintf(fp,"%.3f %e\n",ULdata.UL[i-1].f_min,ULdata.UL[i-1].final_h0);
    fprintf(fp,"%.3f %e\n",ULdata.UL[i].f_max,ULdata.UL[i].final_h0);
  }
  fprintf(fp,"%.3f %e\n",ULdata.UL[ULdata.N_UL-1].f_min,ULdata.UL[ULdata.N_UL-1].final_h0);
  
  fclose(fp);

  printf("*** finished outputting results\n");

  return 0;

}

/*******************************************************************************/
int CalculateStats(SingleULData *UL,INT4 N_split)
{

  REAL8 av=0.0;
  REAL8 var=0.0;
  REAL8 std;
  INT4 i;
  INT4 count;

  /* loop over the number of injection subgroups to get average */
  count=0;
  for (i=0;i<N_split;i++) {
    if (UL->SSULdata[i].est_h0>0.0) {
      av+=UL->SSULdata[i].est_h0;
      count++;
    }
  }

  /* if we have any results above confidence */
  if (count>0) {
    
    av/=(REAL8)count;
    
    /* calculate variance */
    for (i=0;i<N_split;i++) {
      var+=(av-UL->SSULdata[i].est_h0)*(av-UL->SSULdata[i].est_h0);
    }
    
    /* calculate std */
    var/=(REAL8)(count-1);
    std=sqrt(var);
    
    UL->final_h0=av;
    UL->final_h0_err=std/sqrt((REAL8)count-1.0);
    
  }
  else {
     UL->final_h0=0.0;
     UL->final_h0_err=0.0;
  }    

  
  printf("final h0 = %e +/- %e\n",UL->final_h0,UL->final_h0_err);
    
  return 0;
  
}
  
/*******************************************************************************/
int Interpolateh0(SingleSplitULData *UL,INT4 N_h0,REAL8 *h0,REAL8 confidence)
{

  INT4 j,q,r;
  REAL8 c1,c2,h1,h2;
  INT4 *indexes;
  REAL8 grad;
  
  
  /* allocate space for diff */
  tempULdataconf_lower=(REAL8 *)LALMalloc(N_h0*sizeof(REAL8));
  tempULdataconf_higher=(REAL8 *)LALMalloc(N_h0*sizeof(REAL8));

  /* create an array of indexes */
  indexes=(INT4 *)LALMalloc(sizeof(INT4)*N_h0);
  
  /* fill in the temporary data structures */
  q=0;
  r=0;
  for (j=0;j<N_h0;j++) {
    if (conf>=(*UL).conf[j]) {
      tempULdataconf_lower[j]=fabs(confidence-UL->conf[j]);
      tempULdataconf_higher[j]=1.1;
      q++;
    }
    else if (conf<(*UL).conf[j]) {
      tempULdataconf_higher[j]=fabs(confidence-UL->conf[j]);
      tempULdataconf_lower[j]=1.1;
      r++;
    }
  }
  
  /* if we have a result either side of the confidence */
  if ((q>0)&&(r>0)) {

    /* populate indexes */
    for (q=0;q<N_h0;q++) {
      indexes[q]=q;
    }
    
    /* sort array of indexes */
    qsort((void *)indexes, (size_t)N_h0, sizeof(INT4), compare_lower);
    
    /* need to pick point closest to the confidence but lower than it */
    c1=UL->conf[indexes[0]];
    h1=h0[indexes[0]];
    
    /* populate indexes */
    for (q=0;q<N_h0;q++) {
      indexes[q]=q;
    }
    
    /* sort array of indexes */
    qsort((void *)indexes, (size_t)N_h0, sizeof(INT4), compare_higher);
    
    /* need to pick point closest to the confidence but higher than it */
    c2=UL->conf[indexes[0]];
    h2=h0[indexes[0]];
    
    /* approximate a straight line through these points */
    grad=(c2-c1)/(h2-h1);
    (*UL).est_h0=h1+(confidence-c1)/grad;
    

  }
  else (*UL).est_h0=0.0;

  /* free some memory */
  LALFree(tempULdataconf_lower);
  LALFree(tempULdataconf_higher);
  LALFree(indexes);

  return 0;

}
/*******************************************************************************/

/* Sorting function to sort into DECREASING order */
int compare_lower(const void *ip, const void *jp)
{
  
  REAL8 di,dj;

  di=tempULdataconf_lower[*(const int *)ip];
  dj=tempULdataconf_lower[*(const int *)jp];

  
  if (di>dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/*******************************************************************************/

/* Sorting function to sort into DECREASING order */
int compare_higher(const void *ip, const void *jp)
{
  
  REAL8 di,dj;

  di=tempULdataconf_higher[*(const int *)ip];
  dj=tempULdataconf_higher[*(const int *)jp];

  
  if (di>dj)
    return 1;
  
  if (di==dj)
    return 0;

  return -1;
}

/******************************************************************************/
int ReadInjections(char *inputdir,FreqMesh freqmesh,REAL8 h0,REAL8Vector **results)
{

  /* this function takes a frequency band and value of h0 and computes the */
  /* upperlimit confidence */

  char filename[512];
  char line[1024];
  FILE *fp;
  INT4 count,i;
  INT4 fileno=0;
  REAL8 sig;
  char **filelist;
  glob_t globbuf;
  INT4 nfiles;
  REAL8 a;
  INT4 A;
 
  /* make the template filename for searching */
  sprintf(filename,"%s/injection_%.3f-%.3f_%.3e_*.data",inputdir,freqmesh.f_min,freqmesh.f_max,h0);
      
  /* set up some glob stuff */
  globbuf.gl_offs = 1;
  glob(filename, GLOB_ERR, NULL, &globbuf);
  
  /* check if there are any files of this name in the directory */
  if(globbuf.gl_pathc==0)
    {
      fprintf (stderr,"\nNo files in directory %s ... Exiting.\n", inputdir);
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

  count=0;
  /* loop over these files to ocount results first*/
  for (i=0;i<nfiles;i++) {
    
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : cannot open the file %s\n",filelist[i]);
      exit(1);
    }
    
    while (fgets(line,1023,fp)!=NULL) {
    
      sscanf(line,"%lf %lf%lf%lf%lf%lf %d%d %lf%lf %d %lf%lf%lf%lf %lf%lf%lf%lf%lf%d%d%lf%lf%d%lf%lf%lf%lf%lf",&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&sig);  
      count++;
    }
   
    fclose(fp);
    
  }

  /* allocate memory */
  LALDCreateVector(&status,results,count);

  count=0;
  /* loop over these files to record results */
  for (i=0;i<nfiles;i++) {
    
    fp=fopen(filelist[i],"r");
    if (fp==NULL) {
      printf("ERROR : cannot open the file %s\n",filelist[i]);
      exit(1);
    }
    
    while (fgets(line,1023,fp)!=NULL) {
    
      sscanf(line,"%lf %lf%lf%lf%lf%lf %d%d %lf%lf %d %lf%lf%lf%lf %lf%lf%lf%lf%lf%d%d%lf%lf%d%lf%lf%lf%lf%lf",&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&a,&a,&a,&a,&a,&A,&A,&a,&a,&A,&a,&a,&a,&a,&sig);  
      (*results)->data[count]=sig;
      count++;
    }
   
    fclose(fp);
    
  }
  
 
  for (i=0;i<nfiles;i++) LALFree(filelist[i]);
  LALFree(filelist);

  return 0;

}
/******************************************************************************/
int ComputeConfidence(REAL8Vector *results,REAL8 loudsig,REAL8 *confidence)
{

  /* this function takes a frequency band and value of h0 and computes the */
  /* upperlimit confidence */

  INT4 tot=0,count=0;
  UINT4 i; 
 
  for (i=0;i<results->length;i++) {  
    if (results->data[i]<=loudsig) count++;
    tot++;
  }
  
  *confidence=(REAL8)((REAL8)count/(REAL8)tot);
  
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
  REAL8 a,f;
  INT4 A;
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
  nsplit=10;

   
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
      {"nsplit", required_argument, 0, 'n'},
      {"help", no_argument, 0, 'h'}
    };
    /* Scan through list of command line arguments */
    while (!errflg && ((c = getopt_long (argc, argv,"hi:c:f:H:U:x:C:n:",long_options, &option_index)))!=-1)
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
      case 'n':
	nsplit=atoi(optarg);
	break;
      case 'h':
	/* print usage/help message */
	fprintf(stdout,"Arguments are:\n");
	fprintf(stdout,"\t--injectiondir       STRING\t Name of the directory where injection data is stored [DEFAULT= ]\n");
	fprintf(stdout,"\t--coresultsdir       STRING\t Name of the directory where coincident results are stored [DEFAULT= ]\n");
	fprintf(stdout,"\t--freqmeshfile       STRING\t Name of the file containing freq-mesh information [DEFAULT= ]\n");
	fprintf(stdout,"\t--gwamplitudefile    STRING\t Name of the file containing h0 injection values [DEFAULT = ]\n");
	fprintf(stdout,"\t--uloutfile          STRING\t Name of the upperlimit output file [DEFAULT = ]\n");
	fprintf(stdout,"\t--maxoutfile         STRING\t Name of the loudest event output file [DEFAULT = 0.0]\n");
	fprintf(stdout,"\t--confidence         REAL8\t Confidence for which h0_UL(conf) will be calculated [DEFAULT = 0.0]\n");
	fprintf(stdout,"\t--nsplit             INT4\t Number of subsets to split each injection set into for error estimation [DEFAULT = 10]\n");
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

