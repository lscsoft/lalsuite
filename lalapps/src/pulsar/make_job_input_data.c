#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* max number of lines in FILE1 (list of locked stretches) */
#define N 100000

/* max number of lines in FILE2 (sorted list of data) */
#define FN 200000

/* number of SFTs for each job file to produce */
#define SFTPERJOB 50

/* number of seconds in a frame */
#define FRAMELEN 16

/* Maximum length of the path to a frame file */
#define FILENAMEMAX 66

/* Number of nodes on cluster (approximate number of jobs desired) */
#define NODES 296

/* If nonzero, you want to have a fixed number of segments per job */
#define SEGSPERJOB 3

char filenames[FN][FILENAMEMAX+1],printed[FN];
int starttimes[FN];
int nseg[N],tseg[N],tbase[N];

/* First argument is file containing # segs and start times */
/* Second argument is ordered file names */
int main(int argc, char* argv[]) {
  FILE *fp,*fp1=NULL,*fp2=NULL;
  int *nstart,*tstart,*bstart,code,totalsegs=0;
  int locksegs=0,i,j=0,k,firstframe,lastframe,nframes;
  int jobno=0,segno=0,filepointer1=0,filepointer2=0,fileno=0,thistime;
  char bigbuff[1024];

  if (argc!=4) {
	fprintf(stderr,"Syntax:\n%s FILE1 FILE2 DIRNAME\n",argv[0]);
	fprintf(stderr,"where FILE1 contains segment numbers, start times, and SFT length\n");
	fprintf(stderr,"and FILE2  contains ordered frame file names\n");
	fprintf(stderr,"and DIRNAME is a directory number to place output job files.\n");
	return 1;
  }
 
  /* read file containing list of segment numbers and start times */
  if (!(fp=fopen(argv[1],"r"))){
    perror("Could not open file");
    fprintf(stderr,"Unable to open file %s for reading\n",argv[1]);
    exit(1);
  }
  nstart=nseg;
  tstart=tseg;
  bstart=tbase;
  while (3==(code=fscanf(fp,"%d%d%d",nstart++,tstart++,bstart++))){
    totalsegs+=*(nstart-1);
  }

  if (code!=EOF){
    fprintf(stderr,"File %s terminated unexpectedly at line %d\n", argv[1], nstart-nseg);
    exit(1); 
  }
  fclose(fp);
  locksegs=nstart-nseg-1;

  printf("Total number of segments is %d, in %d locked stretches\n", totalsegs, locksegs);

  /* open read file containing file names "in order".  */
  if (!(fp=fopen(argv[2],"r"))){
    perror("Could not open file");
    fprintf(stderr,"Unable to open file %s for reading\n",argv[2]);
    exit(1);
  }
  
  /* parse file containing files names in order Assumes end of
     filename is of the form TGPS-XX.gwf where TGPS is 9 digits.  Note
     that this string is 9+3+4=16 digits long */
  while (1==(code=fscanf(fp,"%s",bigbuff))){
    /* check that file name is not too long for array */
    int namelen;
    bigbuff[1023]='\0';
    namelen=strlen(bigbuff);
    if (namelen>FILENAMEMAX){
      fprintf(stderr,"File name %s length %d too large.  Recompile %s with bigger FILENAMEMAX\n",
	      bigbuff, namelen, argv[0]);
      exit(1);
    }
    
    /* save file name in array */
    strcpy(filenames[fileno], bigbuff);
    
    /* the magic 16 in the next line follows from the file-naming convention above! */
    starttimes[fileno]=atoi(filenames[fileno]+strlen(filenames[fileno])-16);
    fileno++;
  };

  /* Check that file of frame file names was in expected format */
  if (code!=EOF){
    fprintf(stderr,"File %s terminated unexpectedly at line %d\n", argv[2], fileno);
    exit(1); 
  }
  fclose(fp);
  
  /* clear array to keep track of what is printed */
  for (i=0;i<fileno;i++)
    printed[i]=0;
  
  /* check that files names/times are properly ordered */
  for (i=0;i<fileno-1;i++){
    if (starttimes[i]>=starttimes[i+1]){
      fprintf(stderr,"Problem with file time stamps at line %d %s %s\n",i,filenames[i],filenames[i+1]);
      return 1;
    }
  }

  /* start writing job description files loop over all locked segments
     spread the pain equally over the different nodes. */
  for (i=0;i<locksegs;i++){
    for (k=0;k<nseg[i];k++){
      int startat=(totalsegs*jobno)/NODES;
      
      if (
#if (SEGSPERJOB)
	  (segno%SEGSPERJOB)==0
#else
	  segno==startat
#endif
	  ){
	/* clear any lists of printed files first... */
	int j0,j1;
	char fname[256];
	
	j0=j-2;
	if (j0<0)
	  j0=0;
	j1=j+2;
	if (j1>fileno)
	  j1=fileno;
	for (j=j0;j<j1;j++)
	  printed[j]=0;

	/* Need to close old files and open new files */
	if (jobno) {
	  fclose(fp);
	  fclose(fp1);
	  if (fp2) {
	    fclose(fp2);
	    fp2=NULL;
	  }
	}
	/* open new files */
	sprintf(fname,"%s/jobtimes.%05d",argv[3],jobno);
	if (!(fp=fopen(fname,"w"))){
	  perror("Could not write file");
	  fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	  return 1;
	};
	sprintf(fname,"%s/jobdata.%05d.ffl",argv[3],jobno++);
	if (!(fp1=fopen(fname,"w"))){
	  perror("Could not write file");
	  fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	  return 1;
	};
      }
      
      /* time at which next output segment starts */
      thistime=tseg[i]+k*tbase[i];
      /* fprintf(fp,"%d\n",thistime); */
      
      /* How many frame files do we need? */
      firstframe=thistime-thistime%FRAMELEN;
      lastframe =(thistime+tbase[i]-1)-(thistime+tbase[i]-1)%FRAMELEN;
      nframes=1+(lastframe-firstframe)/FRAMELEN;
      
      /* find correct data files */
      while (filepointer1<fileno && starttimes[filepointer1]<firstframe)
	filepointer1++;
      while (filepointer2<fileno && starttimes[filepointer2]<lastframe)
	filepointer2++;
      
      /* see if data files are there */
      if (starttimes[filepointer1]==firstframe &&
	  starttimes[filepointer2]==lastframe &&
	  filepointer2-filepointer1==nframes-1){
      
	/* We have another segment! */
	segno++;
	
	fprintf(fp,"%d %d\n",thistime, tbase[i]);
	for (j=filepointer1;j<=filepointer2;j++){
	  /* don't print a file name twice, please.  Note though */
	  /* that at the end of one input file, and the begining of */
	  /* the next, we haven't printed anything out yet, so we'd */
	  /* better clear printed at that point... */
	  if (!printed[j]){
	    /* The 16 that appears is how far into the filename the detector is identified */
	    /* fprintf(fp1,"%c R %d 16 file:%s\n",filenames[j][15],starttimes[j],filenames[j]); */
	    fprintf(fp1,"%s\n",filenames[j]);
	    printed[j]=1;
	  }
	}
      }
      else { 
	/* missing data */
	if (fp2==NULL){
	  char fname[256];
	  sprintf(fname,"%s/joberror.%05d",argv[3],jobno);
	  if (!(fp2=fopen(fname,"w"))){
	    perror("Could not write file");
	    fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	    return 1;
	  }; 
	}
        fprintf(fp2,"Data missing for SFT starting at %d\n",thistime);
	fprintf(fp2,"Missing data GPS frames in range %d to %d\n",firstframe,lastframe);
      }
    }
  }
  
  fclose(fp);
  fclose(fp1);
  if (fp2)
    fclose(fp2);
 
  return 0;
}
