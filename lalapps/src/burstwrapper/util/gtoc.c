#include <stdio.h>
#include <stdlib.h>
#include <metaio.h>
#include <string.h>
#include <unistd.h>

#define LIGOMETA_IFO_MAX 3
#define LIGOMETA_SEARCH_MAX 25
#define LIGOMETA_CHANNEL_MAX 65


#define BUFSIZE 65536
#define MAXFILE 100

int gettable(char *, char *, char *);

#ifdef linux
static void endian_swap(char * pdata, int dsize, int nelements);
#endif

char *getparams(char *bbuf);

int main(int argc, char *argv[]) {

  FILE *fout;
  char *file;
  int mei, md, mid, status, re, mt0s, mt0ns;
  struct MetaioParseEnvironment parseEnv1;
  const MetaioParseEnv env1 = &parseEnv1;
  char **pgrp = NULL;

  unsigned int i,ii,ngrp = 0;
  unsigned int **igrp = NULL;
  unsigned int *Ngrp = NULL;
  unsigned int *NN = NULL;
  int t0s, t0ns;
  int first = 1;

  char sp[256], sm[256];
  pid_t pid = getpid();
  char buf[BUFSIZE];

  char *toctarget = NULL;

  int Zip = 0;

  size_t size1 = (LIGOMETA_IFO_MAX + LIGOMETA_SEARCH_MAX + LIGOMETA_CHANNEL_MAX) * sizeof(char) + 2*sizeof(int) + 6 * sizeof(float);

  if(argc < 2) {
    fprintf(stderr,"sblit [toc] file\n");
    return 1;
  }

  if(argc == 2) {
    file = argv[1];
  } else {
    file = argv[2];
    toctarget = argv[1];
  }

  if(!strcmp(file+strlen(file)-3,".gz")) {
    char buf[65000];
    
    Zip = 1;

    sprintf(buf,"gzip -df %s",file);
    system(buf);

    *(file+strlen(file)-3) = 0;
  }

  if(!strcmp(file+strlen(file)-4,".bin")) {
    
    FILE *fid, *out;
    unsigned int ui;
    size_t bsiz=0;
    char *bbuf=NULL;
    char fname[1024];

    if(!(fid=fopen(file,"r"))) {
      fprintf(stderr,"Can't open %s\n",file);
      perror("System error");
      return 1;
    }

    i = 0;

    while(fread(buf,sizeof(unsigned int),3,fid) == 3) {

#ifdef linux
      endian_swap(buf, sizeof(unsigned int), 3);
#endif

      memcpy(&t0s,buf,sizeof(unsigned int));
      memcpy(&t0ns,buf+sizeof(unsigned int),sizeof(unsigned int));
      memcpy(&ui,buf+2*sizeof(unsigned int),sizeof(unsigned int));

      ui -= 3*sizeof(unsigned int);

      if(ui > bsiz) {
	bsiz = ui + 1024;
	if(!(bbuf=(char *)realloc(bbuf,bsiz))) {
	  fprintf(stderr,"memory error:%u\n",bsiz);
	  return 1;
	}
      }

      if(fread(bbuf,1,ui,fid) != ui) {
	fprintf(stderr,"File error!\n");
	return 1;
      }

      if(!toctarget ||
	 !strcmp(bbuf,toctarget)) {

	if(i==0) {
	  printf("%i %i 0\n",t0s,t0ns);
	}

	printf("%s %u\n",bbuf,(ui - strlen(bbuf) - 1)/size1);

      }

    }

    free(bbuf);
    fclose(fid);

  } else {

  /* split into tables */
  sprintf(sp,"/tmp/sb.%lu.tmp",pid);
  sprintf(sm,"/tmp/sm.%lu.tmp",pid);
  
  re = gettable(file,sp,"sngl_burst");
  if(re) {

    if(re == -1) {

      /* binary case */
      sprintf(buf,"rm -f %s",sp);
      system(buf);

      if(gettable(file,sm,"sngl_transdata")) {
	return 1;
      }

      status = MetaioOpen( env1, sm );
      if ( status != 0 ) {
	printf( "3-Error opening file\n" );
	MetaioClose( env1 );
	return 1;
      }

      md = MetaioFindColumn(env1, "transdata");
      mei = MetaioFindColumn(env1, "transdata_length");
      mid = MetaioFindColumn(env1, "event_id");
      mt0s = MetaioFindColumn(env1, "x_bins");
      mt0ns = MetaioFindColumn(env1, "y_bins");

      i = 0;
      while((status = MetaioGetRow(env1))==1) {
	unsigned char *dat = env1->ligo_lw.table.elt[md].data.ilwd_char_u.data;
	int dlen = env1->ligo_lw.table.elt[mei].data.int_4s;



	int dim;

	if(first) {
	  t0s = env1->ligo_lw.table.elt[mt0s].data.int_4s;
	  t0ns = env1->ligo_lw.table.elt[mt0ns].data.int_4s;

	  printf("%i %i 0\n",t0s,t0ns);
	  first = 0;
	}

#ifdef linux
	endian_swap((char *)(dat), sizeof(char), strlen(dat));
#endif

	memcpy(&dim, dat + strlen(dat) + 1, sizeof(unsigned int));
#ifdef linux
	endian_swap((char *)(&dim), sizeof(int), 1);
#endif

	if(dim != (dlen-sizeof(int) - strlen(dat) - 1)/size1) {
	  printf("Big datasets not yet supported:%u:%u:%u\n",dim,dlen,(dlen-sizeof(int) - strlen(dat) - 1)/size1);
	  return 1;
	}

	sprintf(buf,"%s.%u",file,i);
	if(!(fout = fopen(buf,"w"))) {
	  fprintf(stderr,"can't write to %s\n",buf);
	  return 1;
	}
	fwrite(dat + strlen(dat) + 1, sizeof(char), (dlen - strlen(dat) - 1)*sizeof(char),fout);
	fclose(fout);

	printf("%s %s %u\n",dat,buf,(dlen-sizeof(int))/size1);

	i++;

      }

      MetaioClose(env1);

      sprintf(buf,"rm -f %s",sm);
      system(buf);
 
      return 0;
    }

    return 1;
  }

  if(gettable(file,sm,"sngl_mime")) {
    return 1;
  }

  /* read mime data */
  status = MetaioOpen( env1, sm );
  if ( status != 0 ) {
    printf( "Error opening file\n" );
    MetaioClose( env1 );
    return 1;
  }

  mei = MetaioFindColumn(env1, "event_id");
  md = MetaioFindColumn(env1, "mimedata");
  /*
  mdl = MetaioFindColumn(env1, "mimedata_length");
  */

  while((status = MetaioGetRow(env1))==1) {
    char *dat = (char *)env1->ligo_lw.table.elt[md].data.ilwd_char_u.data;
    char *eid = env1->ligo_lw.table.elt[mei].data.ilwd_char.data;

    /*
    printf("%s,%u,%s\n",
	   env1->ligo_lw.table.elt[mei].data.ilwd_char.data,
	   env1->ligo_lw.table.elt[mdl].data.int_4s,
	   env1->ligo_lw.table.elt[md].data.ilwd_char_u.data
	   );
    */

    for(i=0;i<ngrp;i++) {
      if(!strcmp(pgrp[i],dat)) {
	break;
      }
    }

    if(i==ngrp) {
      ngrp++;
      pgrp = (char **)realloc(pgrp, ngrp * sizeof(char *));
      pgrp[ngrp-1] = (char *)calloc(1+strlen(dat), sizeof(char));
      strcpy(pgrp[ngrp-1], dat);

      igrp = (unsigned int **)realloc(igrp, ngrp*sizeof(unsigned int *));
      Ngrp = (unsigned int *)realloc(Ngrp, ngrp*sizeof(unsigned int));
      Ngrp[ngrp-1] = 1;
      igrp[i] = (unsigned int *)malloc(sizeof(unsigned int));

    } else {
      (Ngrp[i])++;
      igrp[i] = (unsigned int *)realloc(igrp[i], Ngrp[i]*sizeof(unsigned int));
    }

    sscanf(eid,"sngl_burstgroup:sngl_burst:%u",igrp[i]+Ngrp[i]-1);

  }

  MetaioClose(env1);


  NN = (unsigned int *)calloc(ngrp, sizeof(unsigned int));

  /* write output */
  for(ii=0;ii<ngrp;ii+=MAXFILE)
  {
    FILE *in;
    FILE **out;
    char nam[256] = "sngl_burstgroup:sngl_burst:";
    unsigned int m;

    out = (FILE **)calloc(ngrp, sizeof(FILE *));

    for(i=ii;i<ngrp && i-ii<MAXFILE;i++) {
      sprintf(buf,"%s.%u",file,i);
      if(!(out[i] = fopen(buf,"w"))) {
	fprintf(stderr,"2-Unable to open %s\n",buf);
	return 1;
      }
    }


    if(!(in = fopen(sp,"r"))) {
      fprintf(stderr,"2-Unable to open %s\n",sp);
      return 1;
    }


    while(fgets(buf,BUFSIZE,in)) {

      for(i=ii;i<ngrp && i-ii<MAXFILE;i++) {
	fputs(buf,out[i]);
      }

      if(strstr(buf,"<Stream Delimiter=\",\" Name=\"sngl_burstgroup:sngl_burst:table\" Type=\"Local\">")) {

	while(fgets(buf,BUFSIZE,in)) {
	  char *ptr, *p2;
	  char op2;
	  unsigned int j,n;

	  if(strstr(buf,"</Stream>")) {
	    for(i=ii;i<ngrp && i-ii<MAXFILE;i++) {
	      fseek(out[i],-2,SEEK_CUR);
	      fputs(" \n",out[i]);
	      fputs(buf,out[i]);
	    }
	    break;
	  } 

	  if(!(ptr = strstr(buf,nam))) {
	    fprintf(stderr,"wrong format!!\n");
	    return 1;
	  }

	  ptr += strlen(nam);

	  if(!(p2 = strstr(ptr,"\""))) {
	    fprintf(stderr,"wrong format!!\n");
	    return 1;
	  }

	  op2 = *p2;
	  *p2 = 0;

	  j = (unsigned int)atoi(ptr);

	  *p2 = op2;

	  for(m=0;m<ngrp;m++) {
	    for(n=0;n<Ngrp[m];n++) {
	      if(igrp[m][n] == j) {
		break;
	      }
	    }
	    if(n<Ngrp[m]) {
	      break;
	    }
	  }

	  if(m>=ii && m-ii<MAXFILE) {
	    (NN[m])++;

	    if(buf[strlen(buf)-2] != ',') {
	      char *p = buf + strlen(buf) - 1;
	      *p = ',';
	      *(p+1) = '\n';
	      *(p+2) = 0;
	    }
	  
	    fputs(buf,out[m]);
	  }
	}

      }
    }    

    for(i=ii;i<ngrp && i-ii<MAXFILE;i++) {
      fclose(out[i]);
    }

    fclose(in);

  }

  sprintf(buf,"rm -f %s",sp);
  system(buf);

  sprintf(buf,"rm -f %s",sm);
  system(buf);

  /* print TOC */
  printf("%i %i 0\n",t0s,t0ns);

  for(i=0;i<ngrp;i++) {
    char bu[BUFSIZE];
    sprintf(bu,"%s.%i",file,i);
    printf("%s %s %u\n",pgrp[i],bu,NN[i]);
  }
  }


  if(Zip) {
    char buf[65000];
    sprintf(buf,"gzip -f %s",file);
    system(buf);
  }

  return 0;
}



int gettable(char *file, char *ofile, char *table) {

  FILE *fid, *out;
  int head = 1, notable = 1;
  char buf[BUFSIZE];

  if(!(fid = fopen(file,"r"))) {
    fprintf(stderr,"Unable to open %s\n",file);
    return 1;
  }

  if(!(out = fopen(ofile,"w"))) {
    fprintf(stderr,"2-Unable to open %s\n",ofile);
    return 1;
  }

  while(fgets(buf,BUFSIZE,fid)) {
    if(strstr(buf,"<Table Name")) {

      head = 0;

      if(!strstr(buf,table)) {
	while(!strstr(buf,"</Table>")) {
	  if(!fgets(buf,BUFSIZE,fid)) {
	    fprintf(stderr,"missing /table\n");
	    return 1;
	  }
	}
      } else {
	notable = 0;
	while(!strstr(buf,"</Table>")) {
	  fputs(buf,out);
	  if(!fgets(buf,BUFSIZE,fid)) {
	    fprintf(stderr,"2-Wrong file format\n");
	    return 1;
	  }
	}
	fputs(buf,out);
      }

    } else {
      if(head) {
	fputs(buf,out);
      }
    }
  }

  fprintf(out,"</LIGO_LW>\n");

  fclose(fid);
  fclose(out);

  if(notable) {
    return -1;
  }

  return 0;

}


#ifdef linux
static void endian_swap(char * pdata, int dsize, int nelements)

{

        int i,j,indx;
        char tempbyte;

        if (dsize <= 1) return;

        for (i=0; i<nelements; i++)
        {
                indx = dsize;
                for (j=0; j<dsize/2; j++)
                {
                        tempbyte = pdata[j];
                        indx = indx - 1;
                        pdata[j] = pdata[indx];
                        pdata[indx] = tempbyte;
                }

                pdata = pdata + dsize;
        }

        return;

}
#endif
