#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* max number of lines in FILE1 */
#define N 100000
/* max number of lines in FILE2 */
#define FN 200000

#define SFTPERJOB 50 /* number of SFTs for each job file to produce */
#define FRAMELEN 16  /* number of seconds in a frame */

char filenames[FN][75],printed[FN];
int starttimes[FN];
int nseg[N],tseg[N],tbase[N];

/* First argument is file containing # segs and start times */
/* Second argument is ordered file names */
int main(int argc, char* argv[]) {
  FILE *fp,*fp1=NULL,*fp2=NULL;
  int *nstart,*tstart,*bstart;
  int locksegs=0,i,j=0,k,firstframe,lastframe,nframes;
  int jobno=0,segno=0,filepointer1=0,filepointer2=0,fileno=0,thistime;

  if (argc!=4) {
	fprintf(stderr,"Syntax:\n%s FILE1 FILE2 DIRNAME\n",argv[0]);
	fprintf(stderr,"where FILE1 contains segment numbers, start times, and SFT length\n");
	fprintf(stderr,"and FILE2  contains ordered frame file names\n");
	fprintf(stderr,"and DIRNAME is a directory number to place output job files.\n");
	return 1;
  }
 
  /* read file containing list of segment numbers and start times */
  fp=fopen(argv[1],"r");
  nstart=nseg;
  tstart=tseg;
  bstart=tbase;
  while (3==fscanf(fp,"%d%d%d",nstart++,tstart++,bstart++));
  fclose(fp);
  locksegs=nstart-nseg-1;

  /* read file containing file names "in order".  Assumes end of */
  /* filename is of the form TGPS-XX.gwf where TGPS is 9 digits */
  fp=fopen(argv[2],"r");
  while (1==fscanf(fp,"%s",filenames[fileno])){
    /* the magic 16 in the next line follows from the file-naming convention above! */
    starttimes[fileno]=atoi(filenames[fileno]+strlen(filenames[fileno])-16);
    fileno++;
  };
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

  /* start writing job description files */
  /* loop over all locked segments */
  for (i=0;i<locksegs;i++){
    for (k=0;k<nseg[i];k++){
      if (segno%SFTPERJOB==0){
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
	  fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	  return 1;
	};
	sprintf(fname,"%s/jobdata.%05d.ffl",argv[3],jobno++);
	if (!(fp1=fopen(fname,"w"))){
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
