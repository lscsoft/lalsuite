#include <BurstProcess.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


int main(int argc, char *argv[]) {

  char **files;
  int nfiles;

  double twin = 0.1;
  double toff = 0.5;

  char *outfile=NULL;
  int dotwin=0, dooutfile=0;

  BurstSegInfo *ptr, *optr;

  unsigned int nback, nfore;
  BurstSegInfo backData, foreData;
  BackFunNumberOfUnclusteredEventsParams backParams;
  ForeFunIsDetectedParams foreParams;

  int ttype = 0;

  if(argc < 2) {
    fprintf(stderr,"GetEffnBac [-t twin toff ttype] [-o outfile] files\n");
    return 1;
  }

  if(!strcmp(argv[1],"-t")) {
    twin = atof(argv[2]);
    toff = atof(argv[3]);
    ttype = atoi(argv[4]);
    dotwin = 1;

    if(!strcmp(argv[5],"-o")) {
      outfile = argv[6];
      dooutfile=1;
    }
  }  

  if(!strcmp(argv[1],"-o")) {
    outfile = argv[2];
    dooutfile=1;
  
    if(!strcmp(argv[3],"-t")) {
      twin = atof(argv[4]);
      toff = atof(argv[5]);
      ttype = atoi(argv[6]);
      dotwin = 1;
    }
  }


  if(dotwin) {
    argv += 4;
    argc -= 4;
  }

  if(dooutfile) {
    argv += 2;
    argc -= 2;
  }

  files = argv+1;
  nfiles = argc-1;

  nback = nfore = 0;
  bzero(&backData, sizeof(BurstSegInfo));
  bzero(&foreData, sizeof(BurstSegInfo));
  bzero(&backParams, sizeof(BackFunNumberOfUnclusteredEventsParams));
  bzero(&foreParams, sizeof(ForeFunIsDetectedParams));

  foreParams.twin = twin;
  foreParams.toff = toff;
  foreParams.dType = ttype;

  if(BurstProcess(&nback, &backData, BackFunNumberOfUnclusteredEvents, &backParams, sizeof(backParams),
		  &nfore, &foreData, ForeFunIsDetected, &foreParams, sizeof(foreParams),
		  files, nfiles)) {
    fprintf(stderr,"Error in BurstProcess\n");
    return 1;
  }

  if(dooutfile) {
    FILE *out;
    char fn[32000];

    sprintf(fn,"%s.eff",outfile);
    if(!(out=fopen(fn,"w"))) {
      fprintf(stderr,"Can't open file %s\n",fn);
      return 1;
    }

    ptr = foreData.next;
    while(ptr) {
      char *obuf = (char *)calloc(1+strlen(ptr->params),sizeof(char));
      strcpy(obuf,ptr->params);
      if(obuf[strlen(obuf)-1] == '\n') {
	obuf[strlen(obuf)-1] = 0;
      }

      fprintf(out,"%s,%s\n",obuf,ptr->info);
      ptr = ptr->next;

      free(obuf);
    }

    ptr = foreData.next;
    while(ptr) {
      optr = ptr;
      ptr = ptr->next;
      free(optr->params);
      free(optr->info);
      free(optr);
    }

    fclose(out);

    sprintf(fn,"%s.bac",outfile);
    if(!(out=fopen(fn,"w"))) {
      fprintf(stderr,"Can't open file %s\n",fn);
      return 1;
    }

    ptr = backData.next;
    while(ptr) {
      char *obuf = (char *)calloc(1+strlen(ptr->params),sizeof(char));
      strcpy(obuf,ptr->params);
      if(obuf[strlen(obuf)-1] == '\n') {
	obuf[strlen(obuf)-1] = 0;
      }

      fprintf(out,"%s,%s\n",obuf,ptr->info);
      ptr = ptr->next;

      free(obuf);
    }

    ptr = backData.next;
    while(ptr) {
      optr = ptr;
      ptr = ptr->next;
      free(optr->params);
      free(optr->info);
      free(optr);
    }

    fclose(out);    
  } else {
    printf("Efficiencies:\n");
    ptr = foreData.next;
    while(ptr) {
      printf("%s,%s\n",ptr->params,ptr->info);
      ptr = ptr->next;
    }

    ptr = foreData.next;
    while(ptr) {
      optr = ptr;
      ptr = ptr->next;
      free(optr->params);
      free(optr->info);
      free(optr);
    }

    printf("Backgrounds:\n");
    ptr = backData.next;
    while(ptr) {
      printf("%s,%s\n",ptr->params,ptr->info);
      ptr = ptr->next;
    }

    ptr = backData.next;
    while(ptr) {
      optr = ptr;
      ptr = ptr->next;
      free(optr->params);
      free(optr->info);
      free(optr);
    }
  }

  return 0;
}
