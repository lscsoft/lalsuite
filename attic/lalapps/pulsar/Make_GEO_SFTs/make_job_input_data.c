/*
*  Copyright (C) 2007 Bruce Allen
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
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Bruce Allen
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* max number of lines in FILE1 (list of locked stretches) */
#define N 100000

/* max number of lines in FILE2 (sorted list of data) */
#define FN 400000

/* number of SFTs for each job file to produce */
#define SFTPERJOB 50

/* Maximum length of the path to a frame file */
#define FILENAMEMAX 128

/* Number of nodes on cluster (approximate number of jobs desired) */
#define NODES 10

/* If nonzero, you want to have a fixed number of segments per job */
#define SEGSPERJOB 0

char filenames[FN][FILENAMEMAX+1],printed[FN];
int starttimes[FN];
int nseg[N],tseg[N],tbase[N];

/* First argument is file containing # segs and start times */
/* Second argument is ordered file names */
int main(int argc, char* argv[]) {
  FILE *fp=NULL,*fp1=NULL,*fp2=NULL;
  int *nstart,*tstart,*bstart,code,totalsegs=0,framesec=-1;
  int locksegs=0,j=0,k,firstframe,lastframe,nframes;
  int jobno=0,segno=0,filepointer1=0,filepointer2=0,filenum=0,thistime;
  char bigbuff[1024];

  if (argc!=4) {
	fprintf(stderr,"Syntax:\n%s FILE1 FILE2 DIRNAME\n",argv[0]);
	fprintf(stderr,"where FILE1 contains segment numbers, start times, and SFT length\n");
	fprintf(stderr,"and FILE2  contains ordered frame file names\n");
	fprintf(stderr,"and DIRNAME is a directory number to place output job files.\n");
	exit(1);
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
    fprintf(stderr,"File %s terminated unexpectedly at line %d\n", argv[1], (int)(nstart-nseg));
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
    int namelen, i, deltat;

    bigbuff[1023]='\0';
    namelen=strlen(bigbuff);
    if (namelen>FILENAMEMAX){
      fprintf(stderr,"File name %s length %d too large.  Recompile %s with bigger FILENAMEMAX\n",
	      bigbuff, namelen, argv[0]);
      exit(1);
    }

    /* check that there is still space left */
    if (filenum>=FN){
      fprintf(stderr, "Too many files: %d.  Please recompile with bigger FN=%d\n", filenum, FN);
      exit(1);
    }

    /* save file name in array */
    strcpy(filenames[filenum], bigbuff);

    /* locate the start time and length of the data */
    for (i=strlen(filenames[filenum])-1; i>=0; i--)
      if (filenames[filenum][i]=='-')
	break;

    /* does dash appear as separator in the right place */
    if (i<9){
      fprintf(stderr,"Filename:\n%s\nat line: %d of file: %s\ndoesn't have '-' separator in right place\n",
	      filenames[filenum], filenum+1, argv[2]);
      exit(1);
    }

    /* get start time and duration of file from name */
    starttimes[filenum]=atoi(filenames[filenum]+i-9);
    deltat=atoi(filenames[filenum]+i+1);
    if (deltat<0)
      deltat=0;

    /* record framelength the first time through */
    if (framesec<0)
      framesec=deltat;

    /* check that duration is what's expected */
    if (deltat!=framesec){
            fprintf(stderr,"Filename:\n%s\nat line: %d of file: %s\n has length %d (framesec=%d)\n",
	      filenames[filenum], filenum+1, argv[2], deltat, framesec);
      exit(1);
    }

    /* see if start time same as the previous file (NDAS data can do this!) */
    if (filenum==0 || starttimes[filenum]>starttimes[filenum-1])
      filenum++;
    else if (starttimes[filenum]==starttimes[filenum-1]){
      /* NOTE THAT WE DON'T EXIT -- THIS IS A WARNING ONLY */
#if (0)
      fprintf(stderr,"Identical timestamps:\n%d: %s\n%d: %s\n",
	      starttimes[filenum], filenames[filenum], starttimes[filenum-1], filenames[filenum-1]);
#endif
    }
    else {
      fprintf(stderr,"Problem with file time stamps at line %d of file: %s:\n%d: %s\n%d: %s\n",
	      filenum+1, argv[2], starttimes[filenum], filenames[filenum],
	      starttimes[filenum-1], filenames[filenum-1]);
      exit(1);
    }
  };

  /* Check that file of frame file names was in expected format */
  if (code!=EOF){
    fprintf(stderr,"File %s terminated unexpectedly at line %d\n", argv[2], filenum);
    exit(1);
  }
  fclose(fp);

  /* clear array to keep track of what is printed */
  int i;
  for (i=0;i<filenum;i++)
    printed[i]=0;

  /* check that files names/times are properly ordered */
  for (i=0;i<filenum-1;i++){
    if (starttimes[i]>=starttimes[i+1]){
      fprintf(stderr,"Problem with file time stamps at line %d of file: %s:\n%d: %s\n%d: %s\n",
	      i+1, argv[2], starttimes[i], filenames[i], starttimes[i+1], filenames[i+1]);
      return 1;
    }
  }

  /* start writing job description files loop over all locked segments
     spread the pain equally over the different nodes. */
  fp=fp1=fp2=NULL;
  for (i=0;i<locksegs;i++){
    for (k=0;k<nseg[i];k++){
      int startat;

      /* Segment that the job should start with */
      startat=(totalsegs*jobno)/NODES;

      if (
#if (SEGSPERJOB)
	  /* predetermined fixed number of segments per file */
	  (segno%SEGSPERJOB)==0
#else
	  /* spread segments evenly between nodes */
	  segno==startat
#endif
	  ){
	int j0,j1;
	char fname[256];

	/* clear any lists of printed files first... */
	j0=j-2;
	if (j0<0) j0=0;
	j1=j+2;
	if (j1>filenum) j1=filenum;
	for (j=j0;j<j1;j++)
	  printed[j]=0;

	/* Need to close old files and open new files */
	if (fp)  fclose(fp);
	if (fp1) fclose(fp1);
	if (fp2) fclose(fp2);
	fp=fp1=fp2=NULL;

	/* open new files */
	sprintf(fname,"%s/jobtimes.%05d",argv[3],jobno);
	if (!(fp=fopen(fname,"w"))){
	  perror("Could not write file");
	  fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	  return 1;
	};
	sprintf(fname,"%s/jobdata.%05d.ffl",argv[3],jobno);
	if (!(fp1=fopen(fname,"w"))){
	  perror("Could not write file");
	  fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	  return 1;
	};
	jobno++;
      }

      /* time at which next output segment starts */
      thistime=tseg[i]+k*tbase[i];
      /* fprintf(fp,"%d\n",thistime); */

      /* How many frame files do we need? */
      firstframe=thistime-thistime%framesec;
      lastframe =(thistime+tbase[i]-1)-(thistime+tbase[i]-1)%framesec;
      nframes=1+(lastframe-firstframe)/framesec;

      /* find correct data files */
      while (filepointer1<filenum && starttimes[filepointer1]<firstframe)
	filepointer1++;
      while (filepointer2<filenum && starttimes[filepointer2]<lastframe)
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
	int shouldbethere,my1=filepointer1;

	/* missing data */
	if (fp2==NULL){
	  char fname[256];
	  sprintf(fname,"%s/joberror.%05d",argv[3],jobno-1);
	  if (!(fp2=fopen(fname,"w"))){
	    perror("Could not write file");
	    fprintf(stderr,"Unable to open file %s for writing.\n",fname);
	    return 1;
	  };
	}
        fprintf(fp2,"Data missing for %d-second SFT starting at %d\n",tbase[i],thistime);
	fprintf(fp2,"First frame GPS start: %d\n", firstframe);
	fprintf(fp2,"Last  frame GPS start: %d\n", lastframe);

	/* to list the missing frames, first step through all the
	   frames that SHOULD be there */
	for (shouldbethere=firstframe; shouldbethere<=lastframe; shouldbethere+=framesec){
	  /* for each of these, step through the list are frames that ARE there */
	  int arethere, foundit=0;
	  for (arethere=my1; arethere<=filepointer2; arethere++)
	    /* if we find the frame in the list, then break out of this inner loop */
	    if (starttimes[arethere]==shouldbethere){
	      my1=arethere;
	      foundit=1;
	      break;
	    }
	  /* see if we found it */
	  if (!foundit)
	    fprintf(fp2,"  Missing frame has start time: %d\n", shouldbethere);
	}
      }
    }
  }

  fclose(fp);
  fclose(fp1);
  if (fp2)
    fclose(fp2);

  return 0;}
