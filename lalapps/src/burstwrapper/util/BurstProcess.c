#include <stdio.h>
#include <stdlib.h>
#include <metaio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <BurstProcess.h>
/*
#include <libgen.h>
*/
#include <sys/types.h>
#include <regex.h>

#define LIGOMETA_IFO_MAX 3
#define LIGOMETA_SEARCH_MAX 25
#define LIGOMETA_CHANNEL_MAX 65

#define BUFSIZE 65536

#ifdef linux
static void endian_swap(char * pdata, int dsize, int nelements);
#endif

/******************************************************************/
int BackFunNumberOfUnclusteredEvents(
				     char **info,
				     SnglBurstTableC *input,
				     void *parameters,
				     BurstSegParams *bparams
				     ) {
  BackFunNumberOfUnclusteredEventsParams *params = (BackFunNumberOfUnclusteredEventsParams *)parameters;

  SnglBurstTableC *ptr;

  int i;

  (params->nSegments)++;

  if(params->Nbands == 0) {

    ptr = input->next;
    while(ptr) {
      (params->nBursts)++;
      ptr = ptr->next;
    }

    if(!(*info)) {
      *info = (char *)calloc(256,sizeof(char));
    }

    sprintf(*info,"%i,%i",params->nBursts,params->nSegments);
  } else {

    ptr = input->next;
    while(ptr) {

      for(i=0;i<params->Nbands;i++) {
	if(ptr->central_freq - ptr->bandwidth/2.0 <= params->Fbands[i+1] &&
	   ptr->central_freq + ptr->bandwidth/2.0 >= params->Fbands[i]) {   
	  (params->nBurstsBands[i])++;
	}
      }

      ptr = ptr->next;

    }

    if(!(*info)) {
      *info = (char *)calloc(1024,sizeof(char));
    }

    sprintf(*info,"%i",params->nBurstsBands[0]);

    for(i=1;i<params->Nbands;i++) {
      sprintf(*info,"%s,%i",*info,params->nBurstsBands[i]);
    }
    sprintf(*info,"%s,%i",*info,params->nSegments);

  }

  return 0;
}

/******************************************************************/
#define TSIZ 65000

int ForeFunIsDetected(
		      char **info,
		      SnglBurstTableC *input,
		      void *parameters,
		      BurstSegParams *bparams
		      ) {

  static char timstr[TSIZ];

  ForeFunIsDetectedParams *params = (ForeFunIsDetectedParams *)parameters;

  int i,po=0, nc=0, iind, nout, nin;
  char *str, *str2;

  char *bbuf = bparams->params;
  double twin = params->twin;
  double toff = params->toff;
  int t0s = bparams->t0s;
  int t0ns = bparams->t0ns;

  int dType = params->dType;

  SnglBurstTableC *ptr;
  
  size_t len = strlen(bbuf);

  for(i=0;i<len && nc<5;i++) {
    if(bbuf[i] == '(') {
      po = 1;
    } else if(bbuf[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf[i] == ',') {
      nc++;
    }    
  }

  str = bbuf + i + 1;
  str2 = strchr(str,')');

  strncpy(timstr,str,str2-str);

  str = strtok(timstr,",");

  nout = 0;
  nin = 0;

  while(str) {
    
    nin++;

    iind = atoi(str);

    ptr = input->next;

    while(ptr) {

      if(params->doFCut ==0 ||
	 (ptr->central_freq - 0.5*ptr->bandwidth <= params->f1 &&
	  ptr->central_freq + 0.5*ptr->bandwidth >= params->f0)) {

	int t,tn;
	double dur;

	t = ptr->start_time.gpsSeconds;
	tn = ptr->start_time.gpsNanoSeconds;
	dur = ptr->duration;

	if(dType == 0) {
	  if(fabs((double)(t-t0s)+1E-9*(double)(tn-t0ns)-(double)iind/16384.0-toff) < twin) {
	    nout++;
	    break;
	  }
	} else {
	  if((double)(t-t0s)+1E-9*(double)(tn-t0ns)-toff-twin <= (double)iind/16384.0 &&
	     (double)(t-t0s)+1E-9*(double)(tn-t0ns)-toff+twin+dur >= (double)iind/16384.0) {
	    nout++;
	    break;
	  }
	}
      }
	
      ptr = ptr->next;

    }

    str = strtok(NULL,",");
  }

  (params->nDetected) += nout;
  (params->nInjected) += nin;

  if(!(*info)) {
    *info = (char *)calloc(256,sizeof(char));
  }

  sprintf(*info,"%i,%i",params->nDetected,params->nInjected);

  return 0;
}

/******************************************************************/
#define TSIZ 65000

int EstimationErrors(
		     char **info,
		     SnglBurstTableC *input,
		     void *parameters,
		     BurstSegParams *bparams
		     ) {

  static char timstr[TSIZ];

  EstimationErrorsParams *params = (EstimationErrorsParams *)parameters;

  int i,po=0, nc=0, iind, nin, f1=1;
  char *str, *str2;

  FILE *fid = NULL;
  char waveName[256];
  char abuf[64]; 

  char *bbuf = bparams->params;
  double toff = params->toff;
  int t0s = bparams->t0s;
  int t0ns = bparams->t0ns;

  SnglBurstTableC *ptr;
  SnglBurstTableC *ptre = NULL;

  size_t len = strlen(bbuf);

  double f0, dur0, bw0;
  double h0 = 0.0; 

  for(i=0;i<len && nc<5;i++) {
    if(bbuf[i] == '(') {
      po = 1;
    } else if(bbuf[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf[i] == ',') {
      if(f1) {
	char *ptr;

	strncpy(waveName, bbuf, i);
	waveName[i] = 0;
	f1 = 0;


	ptr = strchr(bbuf+i+1,',');
	if(bbuf[i+1] != '(') {
	  strncpy(abuf, bbuf+i+1, ptr-bbuf-i-1);
	  abuf[ptr-bbuf-i-1] = 0;
	} else {
	  strncpy(abuf, bbuf+i+2, ptr-bbuf-i-2);
	  abuf[ptr-bbuf-i-2] = 0;
	}

	h0 = atof(abuf);
	sprintf(waveName,"%s.%g",waveName,h0);

      } 

      nc++;
    }    
  }

  if(params->outfile) {
    fid = fopen(waveName,"a");
  }

  if(params->f0 <= 0.0 || params->bw0 <= 0.0 || params->dur0 <= 0.0)
    {
      double wavefreq;
      char *fptr = strpbrk(waveName,"0123456789");
      
      if(fptr) {
	wavefreq = atof(fptr);
	
	if(strstr(waveName,"SG")) {
	  f0 = wavefreq;
	  bw0 = 0.1499 * f0; /* 50% power, Q=9 */
	  dur0 = 1.43613 / f0; /* for Q=9, 50% power */
	} else {
	  /*
	  fprintf(stderr,"Unknown waveform!\n");
	  return 1;
	  */
	  f0 = bw0 = dur0 = 0.0;
	}
      } else {
	fprintf(stderr,"Can't get frequency out of wavefrom name!\n");
	return 1;
      }
    } else {
      f0 = params->f0;
      bw0 = params->bw0;
      dur0 = params->dur0;
    }

  nin = 0;

  str = bbuf + i + 1;
  str2 = strchr(str,')');

  strncpy(timstr,str,str2-str);

  str = strtok(timstr,",");

  while(str) { /* loop over times in toc */
    
    double mte = 1E30;
    double mtes = mte;

    nin++;

    iind = atoi(str);

    ptr = input->next;
    while(ptr) { /* loop over all events */

      int t,tn;
      double te, tes;

      t = ptr->start_time.gpsSeconds;
      tn = ptr->start_time.gpsNanoSeconds;

      /* NOTE: iind is time index of injection with respect to earth barycentre, i.e. is the same for all IFOs, even if delays are added due to sky position. */
      tes = (double)(t-t0s)+1E-9*(double)(tn-t0ns)-(double)iind/16384.0-toff;
      te = fabs(tes);

      if(te < mte) {
	mte = te;
	mtes = tes;
	ptre = ptr;
      }

      ptr = ptr->next;

    }

    if(mte <= params->twin) { /* got a detection */

      if(params->outfile) {
	fprintf(fid,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",mtes,ptre->central_freq,ptre->bandwidth,ptre->duration,ptre->amplitude,ptre->confidence,ptre->snr);
      }

      (params->nDetected)++;
      (params->dt) += mte;
      (params->dt2) += mte * mte;
      (params->df) += (ptre->central_freq - f0);
      (params->df2) += pow(ptre->central_freq - f0, 2.0);
      (params->dur) += (ptre->duration - dur0);
      (params->dur2) += pow(ptre->duration - dur0, 2.0);
      (params->bw) += (ptre->bandwidth - bw0);
      (params->bw2) += pow(ptre->bandwidth - bw0, 2.0);
      (params->h) += (ptre->amplitude);
      (params->h2) += pow(ptre->amplitude, 2.0);
      (params->conf) += (ptre->confidence);
      (params->conf2) += pow(ptre->confidence, 2.0);
      (params->snr) += (ptre->snr);
      (params->snr2) += pow(ptre->snr, 2.0);

    } else if(params->reportAll) {
      fprintf(fid,"-1\t-1\t-1\t-1\t-1\t-1\t-1\n");
    }

    str = strtok(NULL,",");
  }

  (params->nInjected) += nin;

  if(!(*info)) {
    *info = (char *)calloc(1024,sizeof(char));
  }

  sprintf(*info,"%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%i,%i",params->dt,params->dt2,params->df,params->df2,params->dur,params->dur2,params->bw,params->bw2,params->h,params->h2,params->conf,params->conf2,params->snr,params->snr2,params->nDetected,params->nInjected);

  if(params->outfile) {
    fclose(fid);
  }

  return 0;
}

/******************************************************************/
int BackFunNumberOfUnclusteredEventsC2(
				       char **info, 
				       SnglBurstTableC *input1, 
				       SnglBurstTableC *input2, 
				       void *parameters, 
				       BurstSegParams *bparams1, 
				       BurstSegParams *bparams2
				       ) {

  BackFunNumberOfUnclusteredEventsC2Params *params = (BackFunNumberOfUnclusteredEventsC2Params *)parameters;

  SnglBurstTableC coin;

  SnglBurstTableC *ptr, *optr;

  double dt;

  int nt,nb;

  if(bparams1->t0s != bparams2->t0s ||
     bparams1->t0ns != bparams2->t0ns) {
    fprintf(stderr,"Error: trying to compare two segments of different time origine\n");
    fprintf(stderr,"%i:%i\t\t%i:%i\n",bparams1->t0s,bparams1->t0ns,bparams2->t0s,bparams2->t0ns);
    return 1;
  }

  for(nt=nb=0, dt = params->dtmin; dt <= params->dtmax; dt += params->dt) {
    
    nt++;
    
    ptr = input2->next;
    while(ptr) {
      long long ti = (int)floor((double)(ptr->start_time.gpsNanoSeconds) + 1e9*dt);
      if(ti >= 1000000000) {
	ptr->start_time.gpsSeconds += ti / 1000000000;
	ptr->start_time.gpsNanoSeconds = ti - 1000000000*(ti / 1000000000);
      }
      if(ptr->start_time.gpsNanoSeconds < 0) {
	ti = llabs(ti);
	ptr->start_time.gpsSeconds -= 1 + ti / 1000000000;
	ti -= 1000000000*(ti / 1000000000);
	ptr->start_time.gpsNanoSeconds = 1000000000 - ti;
      }

      ptr = ptr->next;
    }

    bzero(&coin, sizeof(SnglBurstTableC));

    if(Get2Coincidences(&coin,input1,input2,params->cparams)) {
      fprintf(stderr,"Can't get coincidences!\n");
      return 1;
    }

    ptr = coin.next;
    while(ptr) {
      nb++;
      optr = ptr;
      ptr = ptr->next;
      free(optr);
    }
  }

  (params->nSegments)++;
  (params->nBursts)+= (double)nb / (double)nt;

  if(!(*info)) {
    *info = (char *)calloc(256,sizeof(char));
  }

  sprintf(*info,"%g,%i",params->nBursts,params->nSegments);


  return 0;
}

/******************************************************************/
int ForeFunIsDetectedC2(
			char **info,
			SnglBurstTableC *input1,
			SnglBurstTableC *input2,
			void *parameters,
			BurstSegParams *bparams1,
			BurstSegParams *bparams2
			) {

  static char timstr1[TSIZ];
  /*  static char timstr2[TSIZ]; */

  ForeFunIsDetectedC2Params *params = (ForeFunIsDetectedC2Params *)parameters;

  SnglBurstTableC coin;

  int i,po=0, nc=0, iind, nout, nin;
  char *str, *str2;

  double twin = params->twin;
  double toff = params->toff;

  char *bbuf1 = bparams1->params;
  /*  char *bbuf2 = bparams2->params; */

  int t0s = bparams1->t0s;
  int t0ns = bparams1->t0ns;

  SnglBurstTableC *ptr, *optr;

  size_t len1; /* , len2;*/

  if(bparams1->t0s != bparams2->t0s ||
     bparams1->t0ns != bparams2->t0ns) {
    fprintf(stderr,"Error: trying to compare two segments of different time origine\n");
    return 1;
  }

  len1 = strlen(bbuf1);

  po=nc=0;
  for(i=0;i<len1 && nc<5;i++) {
    if(bbuf1[i] == '(') {
      po = 1;
    } else if(bbuf1[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf1[i] == ',') {
      nc++;
    }    
  }
  str = bbuf1 + i + 1;
  str2 = strchr(str,')');
  strncpy(timstr1,str,str2-str);

  /*
  len2=strlen(bbuf2);
  po=nc=0;
  for(i=0;i<len2 && nc<5;i++) {
    if(bbuf2[i] == '(') {
      po = 1;
    } else if(bbuf2[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf2[i] == ',') {
      nc++;
    }    
  }
  str = bbuf2 + i + 1;
  str2 = strchr(str,')');
  strncpy(timstr2,str,str2-str);

  if(strcmp(timstr1,timstr2)) {
    fprintf(stderr,"Error: trying to compare segments with different injection times\n");
    return 1;
  }
  */

  bzero(&coin, sizeof(SnglBurstTableC));

  if(Get2Coincidences(&coin,input1,input2,params->cparams)) {
    fprintf(stderr,"Can't get coincidences!\n");
    return 1;
  }

  str = strtok(timstr1,",");
  nout = 0;
  nin = 0;

  while(str) {
    
    nin++;

    iind = atoi(str);

    ptr = coin.next;

    while(ptr) {

      int t,tn;

      t = ptr->start_time.gpsSeconds;
      tn = ptr->start_time.gpsNanoSeconds;

      if(fabs((double)(t-t0s)+1E-9*(double)(tn-t0ns)-(double)iind/16384.0-toff) < twin) {
	nout++;
	break;
      }

      ptr = ptr->next;

    }

    str = strtok(NULL,",");
  }

  (params->nDetected) += nout;
  (params->nInjected) += nin;

  if(!(*info)) {
    *info = (char *)calloc(256,sizeof(char));
  }

  sprintf(*info,"%i,%i",params->nDetected,params->nInjected);

  ptr = coin.next;
  while(ptr) {
    optr = ptr;
    ptr = ptr->next;
    free(optr);
  }

  return 0;
}

/******************************************************************/
#define SnglBurstTableCCopy(a,b) a=b
/*
 memcpy(a.ifo, b.ifo, LIGOMETA_IFO_MAX * sizeof(char)); \
 memcpy(a.search, b.search, LIGOMETA_SEARCH_MAX * sizeof(char)); \
 memcpy(a.channel, b.channel, LIGOMETA_CHANNEL_MAX * sizeof(char))
*/

void hpsort2(unsigned long n, int ra[], SnglBurstTableC rb[]) 
{ 
  unsigned long i,ir,j,l; 
  int rra; 
  SnglBurstTableC rrb;

  ra--;
  rb--;

  if (n < 2) return; 

  l=(n >> 1)+1; 
  ir=n; 

  for (;;) { 
    if (l > 1) {
      rra=ra[--l]; 
      SnglBurstTableCCopy(rrb,rb[l]);
    } else { 
      rra=ra[ir];
      SnglBurstTableCCopy(rrb,rb[ir]);

      ra[ir]=ra[1]; 
      SnglBurstTableCCopy(rb[ir],rb[1]);

      if (--ir == 1) { 
	ra[1]=rra; 
	SnglBurstTableCCopy(rb[1],rrb);
	break; } 
    } 
    
    i=l; 
    j=l+l; 

    while (j <= ir) { 
      if (j < ir && ra[j] < ra[j+1]) 
	j++; 
      if (rra < ra[j]) { 
	ra[i]=ra[j];
	SnglBurstTableCCopy(rb[i],rb[j]);

	i=j; 
	j <<= 1; 
      } else break; 
    } 
    ra[i]=rra; 
    SnglBurstTableCCopy(rb[i],rrb); 
  }
}

/******************************************************************/
int AppendEvents(SnglBurstTableC *events, char *bbuf) {

  char *bbuf1;
  unsigned int i, N;
  SnglBurstTableC *ev, *eve;
  int *tim;

  size_t size1 = (LIGOMETA_IFO_MAX + LIGOMETA_SEARCH_MAX + LIGOMETA_CHANNEL_MAX) * sizeof(char) + 2*sizeof(int) + 6 * sizeof(float);

  bbuf1 = bbuf + strlen(bbuf) + 1;

  memcpy(&N,bbuf1,sizeof(unsigned int));

  bbuf1 += sizeof(unsigned int);

#ifdef linux
  endian_swap((char *)(&N), sizeof(unsigned int), 1);
#endif

  if(N>0) {
    tim = (int *)malloc(N*sizeof(int));
    eve = (SnglBurstTableC *)calloc(N, sizeof(SnglBurstTableC));
  
    for(i=0;i<N;i++) {
      char *p1 = bbuf1 + i*size1, *pst, *pr4;

      ev = eve + i;

      /*
      memcpy(ev->ifo, p1, LIGOMETA_IFO_MAX * sizeof(char));
      */

      p1 += LIGOMETA_IFO_MAX * sizeof(char);
      /*
      memcpy(ev->search, p1,LIGOMETA_SEARCH_MAX * sizeof(char));
      */

      p1 += LIGOMETA_SEARCH_MAX * sizeof(char);
      /*
      memcpy(ev->channel, p1, LIGOMETA_CHANNEL_MAX * sizeof(char));
      */

      pst = p1 + LIGOMETA_CHANNEL_MAX * sizeof(char);

      memcpy(&(ev->start_time.gpsSeconds), pst, sizeof(int));
      memcpy(&(ev->start_time.gpsNanoSeconds), pst+sizeof(int), sizeof(int));

#ifdef linux
      endian_swap((char *)(&(ev->start_time.gpsSeconds)), sizeof(int), 1);
      endian_swap((char *)(&(ev->start_time.gpsNanoSeconds)), sizeof(int), 1);
#endif

      tim[i] = ev->start_time.gpsSeconds;

      pr4 = pst + 2*sizeof(int);

      memcpy(&(ev->duration), pr4, sizeof(float));

#ifdef linux
      endian_swap((char *)(&(ev->duration)), sizeof(float), 1);
#endif


      p1 = pr4 + sizeof(float);
      memcpy(&(ev->central_freq), p1, sizeof(float));

#ifdef linux
      endian_swap((char *)(&(ev->central_freq)), sizeof(float), 1);
#endif

      p1 += sizeof(float);
      memcpy(&(ev->bandwidth), p1, sizeof(float));
      
#ifdef linux
      endian_swap((char *)(&(ev->bandwidth)), sizeof(float), 1);
#endif

      p1 += sizeof(float);
      memcpy(&(ev->amplitude), p1, sizeof(float));
    
#ifdef linux
      endian_swap((char *)(&(ev->amplitude)), sizeof(float), 1);
#endif

      p1 += sizeof(float);
      memcpy(&(ev->snr), p1, sizeof(float));

#ifdef linux
      endian_swap((char *)(&(ev->snr)), sizeof(float), 1);
#endif

      p1 += sizeof(float);
      memcpy(&(ev->confidence), p1, sizeof(float));
    
#ifdef linux
      endian_swap((char *)(&(ev->confidence)), sizeof(float), 1);
#endif

   }

    /* now sort events in time (seconds only) */
  
    hpsort2(N,tim,eve);

    ev = events;
    for(i=0;i<N;i++) {
      ev->next = (SnglBurstTableC *)calloc(1,sizeof(SnglBurstTableC));
      ev = ev->next;

      SnglBurstTableCCopy((*ev),(eve[i]));
    }
    
    free(eve);
    free(tim);
  }

  return 0;
}

/******************************************************************/

char *getparams(char *bbuf) {

  int i,po=0, nc=0;
  size_t len = strlen(bbuf);

  for(i=0;i<len && nc<6;i++) {
    if(bbuf[i] == '(') {
      po = 1;
    } else if(bbuf[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf[i] == ',') {
      nc++;
    }    
  }

  return bbuf+i;

}

/******************************************************************/

char *getallparams(char *bbuf) {

  char *str;
  int i,j,po=0, nc=0, peo=0;
  size_t len = strlen(bbuf);
  
  str = calloc(2+len,1);

  for(i=j=0;i<len && nc<5;i++) {
    if(bbuf[i] == '(') {
      peo = po = 1;
    } else if(bbuf[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf[i] == ',') {
      nc++;
    }    
    
    if(!peo) {
      str[j] = bbuf[i];
      j++;
    } 

    if(peo && !po) {
      peo = 0;
    }

  }

  strcat(str,getparams(bbuf));

  return str;

}

/******************************************************************/

int paramsdiff(char *bbuf, char *toc, int compETG) {

  int i, i1, oi, oi1, po=0, nc=0, nc1=0, peo, peo1;

  /* first compare ETG parameters */
  if(compETG) {
    if(strcmp(getparams(bbuf),getparams(toc))) {
      return 1;
    }
  }

  /* now compare injection parameters. Thing in parentheses are not compared */
  oi = oi1 = 0;

  while(nc<5) {

    size_t len = strlen(bbuf);

    for(peo=po=0,i=oi;i<len;i++) {
      if(bbuf[i] == '(') {
	peo = po = 1;
      } else if(bbuf[i] == ')') {
	po = 0;
      } else if(po==0 && bbuf[i] == ',') {
	nc++;
	break;
      }    
    }

    len = strlen(toc);
    for(peo1=po=0,i1=oi1;i1<len;i1++) {
      if(toc[i1] == '(') {
	peo1 = po = 1;
      } else if(toc[i1] == ')') {
	po = 0;
      } else if(po==0 && toc[i1] == ',') {
	nc1++;
	break;
      }    
    }

    if(!peo && !peo1) {
      if(strncmp(bbuf+oi,toc+oi1,i-oi)) {
	return 1;
      }
    }

    oi = i+1;
    oi1 = i1+1;

  }

  return 0;
}

/******************************************************************/
int paramsdiffcp(char *bbuf, char *toc, int ETGcomp) {

  int i, po=0, nc=0;
  size_t len;

  if(ETGcomp) {
    if(strcmp(getparams(bbuf),getparams(toc))) {
      return 1;
    }
  }

  len = strlen(bbuf);

  for(i=0;i<len && nc<5;i++) {
    char bbufi = bbuf[i];

    if(bbufi != toc[i]) {
      return 1;
    }

    if(bbufi == '(') {
      po = 1;
    } else if(bbufi == ')') {
      po = 0;
    } else if(po==0 && bbufi == ',') {
      nc++;
    }    
  }

  /*
  if(strncmp(bbuf,toc,i)) {
    return 1;
  }
  */

  return 0;
}

/******************************************************************/
int paramsdiffIFO(char *s1, char *s2, char *IFO1, char *IFO2) {

  size_t off1, off2;
  size_t l1,l2;

  if(!strcmp(s1,s2)) {
    return 0;
  }

  off1 = off2 = 0;

  l1 = strlen(s1);
  l2 = strlen(s2);

  while(off1 < l1 &&
	off2 < l2) {

    while(s1[off1] == s2[off2]) {
      off1++;
      off2++;

      if(off1 >= l1 || off2>=l2) {
	return 0;
      }
    }

    if(l1<off1+2 || l2<off2+2) {
      return 1;
    }

    if(strncmp(s1+off1,IFO1,2) && strncmp(s1+off1,IFO2,2)) {
      return 1;
    }

    if(strncmp(s2+off2,IFO1,2) && strncmp(s2+off2,IFO2,2)) {
      return 1;
    }

    off1 += 2;
    off2 += 2;

  }

  return 0;
}

/******************************************************************/
int isbackground(char *bbuf) {

  int i,po=0, nc=0;
  size_t len = strlen(bbuf);

   for(i=0;i<len && nc<5;i++) {
    if(bbuf[i] == '(') {
      po = 1;
    } else if(bbuf[i] == ')') {
      po = 0;
    } else if(po==0 && bbuf[i] == ',') {
      nc++;
    }    
   }

   if(i==5) {
     return 1;
   } 
     
   return 0;
}


/******************************************************************/
int RemoveRepetitions(char *bbuf, size_t ui) {

  char *ptr;

  ptr = strchr(bbuf,',');

  if(ptr && ptr[1]=='(') {
    char *pend, *buf, *ph, *pin, *pout;
    double h0;

    pend = strchr(ptr+2,')');
    buf = (char *)malloc(1+strlen(bbuf));
    bzero(buf,1+strlen(bbuf));
    strncpy(buf,ptr+2,pend-ptr-1);
    
    ph = strtok(buf,",)");

    h0 = atof(ph);

    ph = strtok(NULL,",)");

    while(ph) {
      double h;
      h = atof(ph);

      if(fabs(h-h0)/fabs(h+h0) > 1e-6) {
	free(buf);
	return 0;
      }
      ph = strtok(NULL,",)");
    }
  
    free(buf);
    
    buf = (char *)malloc(1+ui);
    memcpy(buf,bbuf,ui);
    bzero(bbuf,ui);

    pin = buf;
    pout = bbuf;

    pend = strchr(pin,',');
    memcpy(pout,pin,pend-pin+1);
    pout += pend-pin+1;
    pin = pend+2; /* put us on char after ( */

    pend = strchr(pin,',');
    memcpy(pout,pin,pend-pin+1);
    pout += pend-pin+1;

    pin = strchr(pin,'(') + 1;
    pend = strchr(pin,',');
    memcpy(pout,pin,pend-pin+1);
    pout += pend-pin+1; /* alpha */

    pin = strchr(pin,'(') + 1;
    pend = strchr(pin,',');
    memcpy(pout,pin,pend-pin+1);
    pout += pend-pin+1; /* delta */

    pin = strchr(pin,'(') + 1;
    pend = strchr(pin,',');
    memcpy(pout,pin,pend-pin+1); 
    pout += pend-pin+1; /* psi */

    pin = strchr(pin,'(');
    memcpy(pout,pin,1+ui-(pin-buf));

    free(buf);
  }

  return 0;
}

/******************************************************************/
int ProcessFilesSub(
		 unsigned int *nback, /* number of backgnd parameters */
		 unsigned int *nfore, /* number of foregnd parameters */
		 char **files,        /* list of input files */
		 int nfiles,          /* number of files */
		 int *nSeg,           /* number of jobs in all files */
		 int **isBack,        /* whether each job is backgnd */
		 int **jInd,          /* index in tocs of all jobs */
		 int *ntocs,          /* number of tocs entries */
		 char ***tocs,        /* indep. parameters list */
		 BurstSegParams **bparams, /* t0 & params for all jobs */
		 SnglBurstTableC ***events, /* events[tocs Id][job order] */
		 int **jeId,          /* second index in events of all jobs */
		 int **Nevents,       /* number of events per tocs entry */
		 int **backInd,       /* index in tocs of back params */
		 int **foreInd,        /* index in tocs of fore params */
		 regex_t *rexp
		 ) {

  FILE *fid;
  unsigned int j, fi, t0s, t0ns, ui, new1;
  int Zip = 0;

  clock_t tic, toc;

  char buf[BUFSIZE];

  size_t bsiz=0;
  char *bbuf=NULL;
  char *str2;

  char *file;

  char tmpfil[2048];

  *nback = 0;
  *nfore = 0;
  *ntocs = 0;
  *nSeg = 0;

  *tocs = NULL;
  *events = NULL;
  *jeId = NULL;
  *Nevents = NULL;
  *backInd = NULL;
  *foreInd = NULL;
  *isBack = NULL;
  *jInd = NULL;
  *bparams = NULL;

  for(fi=0;fi<nfiles;fi++) {

    file = files[fi];

    tic = clock();

    printf("Processing %s...",file);
    fflush(stdout);

    if(!strcmp(file+strlen(file)-3,".gz")) {
      char buf[65000];
      
      Zip = 1;
      
      sprintf(tmpfil,"/var/tmp/BurstProcessTmpFileXXXXXX");

      mkstemp(tmpfil);

      sprintf(buf,"cp %s %s.gz",file,tmpfil);
      system(buf);

      sprintf(buf,"gzip -df %s.gz",tmpfil);
      system(buf);

      /*
      *(file+strlen(file)-3) = 0;
      */
    } else {
      strcpy(tmpfil,file);
    }
    
    if(!(fid=fopen(tmpfil,"r"))) {
      fprintf(stderr,"Can't open %s\n",tmpfil);
      perror("System error");
      return 1;
    }

    while(fread(buf,sizeof(unsigned int),3,fid) == 3) {
#ifdef linux
      endian_swap(buf, sizeof(unsigned int), 3);
#endif


      memcpy(&t0s,buf,sizeof(unsigned int));
      memcpy(&t0ns,buf + sizeof(unsigned int),sizeof(unsigned int));
      memcpy(&ui,buf + 2*sizeof(unsigned int),sizeof(unsigned int));
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

      if(!rexp ||
	 regex(rexp,bbuf)) {

      (*nSeg)++;

      (*isBack) = (int *)realloc((*isBack), (*nSeg) * sizeof(int));
      (*jInd) = (int *)realloc((*jInd), (*nSeg) * sizeof(int));
      (*jeId) = (int *)realloc(*jeId, *nSeg * sizeof(int));

      (*bparams) = (BurstSegParams *)realloc((*bparams), (*nSeg) * sizeof(BurstSegParams));
      (*bparams)[(*nSeg)-1].t0s = t0s;
      (*bparams)[(*nSeg)-1].t0ns = t0ns;

      RemoveRepetitions(bbuf,ui);

      (*bparams)[(*nSeg)-1].params = (char *)calloc(1+strlen(bbuf),sizeof(char));
      strcpy((*bparams)[(*nSeg)-1].params, bbuf);

      str2 = getallparams(bbuf);

      for(j=0;j<(*ntocs);j++) {
	if(!paramsdiff(str2,(*tocs)[j],1)) break;
      }

      (*jInd)[(*nSeg)-1] = j;
      
      if(j==(*ntocs)) {
	new1 = 1;
	(*ntocs)++;
	(*tocs) = (char **)realloc((*tocs), (*ntocs) * sizeof(char *));
	(*tocs)[j] = str2;
	
	(*events) = (SnglBurstTableC **)realloc((*events), (*ntocs) * sizeof(SnglBurstTableC *));
	(*events)[j] = NULL;

	(*Nevents) = (int *)realloc((*Nevents), (*ntocs)*sizeof(int));
	(*Nevents)[j] = 0;

      } else {
	new1 = 0;
	free(str2);
      }
      
      ((*Nevents)[j])++;
      (*events)[j] = (SnglBurstTableC *)realloc((*events)[j], ((*Nevents)[j]) * sizeof(SnglBurstTableC));
      bzero((*events)[j] + (*Nevents)[j] - 1, sizeof(SnglBurstTableC));
      AppendEvents((*events)[j] + (*Nevents)[j] - 1,bbuf);
      (*jeId)[*nSeg-1] = (*Nevents)[j] - 1;

      if(isbackground(bbuf)) {

	(*isBack)[(*nSeg)-1] = 1;

	if(new1) {
	  (*nback)++;	
	  (*backInd) = (int *)realloc((*backInd), (*nback)*sizeof(int));
	  (*backInd)[(*nback)-1] = j;
	}

      } else {

	(*isBack)[(*nSeg)-1] = 0;

	if(new1) {
	  (*nfore)++;
	  (*foreInd) = (int *)realloc((*foreInd), (*nfore)*sizeof(int));
	  (*foreInd)[(*nfore)-1] = j;
	}
      }
      }      
    }
    
    fclose(fid);

    if(Zip) {
      char buf[65000];
      sprintf(buf,"rm %s",tmpfil);
      system(buf);
    }

    toc = clock();
    
    printf("done (%g s)\n",(double)(toc-tic)/(double)CLOCKS_PER_SEC);
    
  }

  free(bbuf);

  return 0;
}

int ProcessFiles(
		 unsigned int *nback, /* number of backgnd parameters */
		 unsigned int *nfore, /* number of foregnd parameters */
		 char **files,        /* list of input files */
		 int nfiles,          /* number of files */
		 int *nSeg,           /* number of jobs in all files */
		 int **isBack,        /* whether each job is backgnd */
		 int **jInd,          /* index in tocs of all jobs */
		 int *ntocs,          /* number of tocs entries */
		 char ***tocs,        /* indep. parameters list */
		 BurstSegParams **bparams, /* t0 & params for all jobs */
		 SnglBurstTableC ***events, /* events[tocs Id][job order] */
		 int **jeId,          /* second index in events of all jobs */
		 int **Nevents,       /* number of events per tocs entry */
		 int **backInd,       /* index in tocs of back params */
		 int **foreInd        /* index in tocs of fore params */
		 ) {
  return ProcessFilesSub(nback,nfore,files,nfiles,nSeg,isBack,jInd,ntocs,tocs,bparams,events,jeId,Nevents,backInd,foreInd,NULL);
}


/******************************************************************/
int BurstProcess(
		 unsigned int *nback, /* number of background segments */
		 BurstSegInfo *backData, /* params,info for each segment */ 
		 int (*backFun)(char **, SnglBurstTableC *, void *, BurstSegParams *), /* function to process each segment */
		 void *backParams0,
		 size_t backSize,
		 unsigned int *nfore,
		 BurstSegInfo *foreData,
		 int (*foreFun)(char **, SnglBurstTableC *, void *, BurstSegParams *),
		 void *foreParams0,
		 size_t foreSize,
		 char **files,
		 int nfiles
		 ) {
  return BurstProcessSub(nback,backData,backFun,backParams0,backSize,nfore,foreData,foreFun,foreParams0,foreSize,files,nfiles,NULL);
}

int BurstProcessSub(
		 unsigned int *nback, /* number of background segments */
		 BurstSegInfo *backData, /* params,info for each segment */ 
		 int (*backFun)(char **, SnglBurstTableC *, void *, BurstSegParams *), /* function to process each segment */
		 void *backParams0,
		 size_t backSize,
		 unsigned int *nfore,
		 BurstSegInfo *foreData,
		 int (*foreFun)(char **, SnglBurstTableC *, void *, BurstSegParams *),
		 void *foreParams0,
		 size_t foreSize,
		 char **files,
		 int nfiles,
		 char *matchexp
		 ) {

  regex_t *rexp = NULL;

  int i,j;

  char **tocs = NULL;
  char **infos = NULL;
  int *Nevents = NULL;
  SnglBurstTableC **events = NULL;
  int ntocs = 0;
  
  int *backInd = NULL;
  int *foreInd = NULL;

  int *isBack = NULL;
  int *jInd = NULL;
  BurstSegParams *bparams = NULL;
  int nSeg = 0;

  char **backParams = NULL, 
    **foreParams = NULL;

  int *jeId;

  int *eId;

  *nback = 0;
  *nfore = 0;

  if(matchexp) {
    regcomp(rexp,matchexp,REG_NOSUB);
  }

  if(ProcessFilesSub(
		  nback,
		  nfore,
		  files,
		  nfiles,
		  &nSeg,
		  &isBack,
		  &jInd,
		  &ntocs,
		  &tocs,
		  &bparams,
		  &events,
		  &jeId,
		  &Nevents,
		  &backInd,
		  &foreInd,
		  rexp
		  )) {
    fprintf(stderr,"Error in ProcessFiles\n");
    return 1;
  }

  if(rexp) {
    regfree(rexp);
  }

  infos = (char **)calloc(ntocs, sizeof(char *));

  eId = (int *)calloc(ntocs, sizeof(int));

  backParams = (char **)malloc(ntocs * sizeof(char *));
  for(i=0;i<ntocs;i++) {
    backParams[i] = (char *)malloc(backSize);
    memcpy(backParams[i], backParams0, backSize);
  }

  foreParams = (char **)malloc(ntocs * sizeof(char *));
  for(i=0;i<ntocs;i++) {
    foreParams[i] = (char *)malloc(foreSize);
    memcpy(foreParams[i], foreParams0, foreSize);
  }

   
  for(i=0;i<nSeg;i++) {

    j = jInd[i];

    if(isBack[i]) {

      if((*backFun)(infos + j, events[j] + eId[j], backParams[j], bparams + i)) {
	  fprintf(stderr,"Error in backFun\n");
	  return 1;
	}

      } else {

	if((*foreFun)(infos + j, events[j] + eId[j], foreParams[j], bparams + i)) {
	  fprintf(stderr,"Error in foreFun\n");
	  return 1;
	}

      }
   
    (eId[j])++;
   
  }
    
  {
    BurstSegInfo *pbsi = backData;
      
    for(j=0;j<*nback;j++) {
      pbsi->next = (BurstSegInfo *)calloc(1,sizeof(BurstSegInfo));
      pbsi = pbsi->next;
      
      pbsi->params = tocs[backInd[j]];
      pbsi->info = infos[backInd[j]];
    }
  }

  {
    BurstSegInfo *pbsi = foreData;
    
    for(j=0;j<*nfore;j++) {
      pbsi->next = (BurstSegInfo *)calloc(1,sizeof(BurstSegInfo));
      pbsi = pbsi->next;
      
      pbsi->params = tocs[foreInd[j]];
      pbsi->info = infos[foreInd[j]];
    }
  }

  for(j=0;j<ntocs;j++) {
    for(i=0;i<Nevents[j];i++) {
      SnglBurstTableC *eptr = (events[j]+i)->next, *oeptr;
      while(eptr) {
	oeptr = eptr->next;
	free(eptr);
	eptr = oeptr;
      }
    }
    free(events[j]);

    free(backParams[j]);
    free(foreParams[j]);
  }

  free(jeId);
  free(events);
  free(Nevents);
  free(backParams);
  free(foreParams);
  
  free(tocs);
  free(infos);
  free(foreInd);
  free(backInd);

  free(jInd);
  free(isBack);

  for(i=0;i<nSeg;i++) {
    free(bparams[i].params);
  }
  free(bparams);

  free(eId);

  return 0;
}

/******************************************************************/

int BurstProcess2(
		  unsigned int *nback, /* number of background segments */
		  BurstSegInfo *backData, /* params,info for each segment */ 
		  int (*backFun)(char **, SnglBurstTableC *, SnglBurstTableC *, void *, BurstSegParams *, BurstSegParams *), /* function to process each segment */
		  void *backParams0,
		  size_t backSize,
		  unsigned int *nfore,
		  BurstSegInfo *foreData,
		  int (*foreFun)(char **, SnglBurstTableC *, SnglBurstTableC *, void *, BurstSegParams *, BurstSegParams *),
		  void *foreParams0,
		  size_t foreSize,
		  char **files1,
		  int nfiles1,
		  char **files2,
		  int nfiles2,
		  int mixETGParams,
		  char *IFO1,
		  char *IFO2
		  ) {

  int i,i1,i2,j,j1,j2,kk,k1,k2;

  int **jInd1Inv, **jInd2Inv;
  int *NjInd1Inv, *NjInd2Inv;

  char **infosFore = NULL,
    **infosBack = NULL;
  int *j1Fore = NULL,
    *j2Fore = NULL,
    *j1Back = NULL,
    *j2Back = NULL;

  char **backParams = NULL, 
    **foreParams = NULL;

  int *foreNew1 = NULL,
    *backNew1 = NULL;

  char **tocs1 = NULL;
  int *Nevents1 = NULL;
  SnglBurstTableC **events1 = NULL;
  int ntocs1 = 0;
  int *backInd1 = NULL;
  int *foreInd1 = NULL;
  int *isBack1 = NULL;
  int *jInd1 = NULL;
  BurstSegParams *bparams1 = NULL;
  int nSeg1 = 0;

  char **tocs2 = NULL;
  int *Nevents2 = NULL;
  SnglBurstTableC **events2 = NULL;
  int ntocs2 = 0;
  int *backInd2 = NULL;
  int *foreInd2 = NULL;
  int *isBack2 = NULL;
  int *jInd2 = NULL;
  BurstSegParams *bparams2 = NULL;
  int nSeg2 = 0;

  int nback1, nback2;
  int nfore1, nfore2;

  int *jeId1 = NULL,
    *jeId2 = NULL;

  *nback = 0;
  *nfore = 0;


  if(ProcessFiles(
		  &nback1,
		  &nfore1,
		  files1,
		  nfiles1,
		  &nSeg1,
		  &isBack1,
		  &jInd1,
		  &ntocs1,
		  &tocs1,
		  &bparams1,
		  &events1,
		  &jeId1,
		  &Nevents1,
		  &backInd1,
		  &foreInd1
		  )) {
    fprintf(stderr,"Error in ProcessFiles\n");
    return 1;
  }

  if(ProcessFiles(
		  &nback2,
		  &nfore2,
		  files2,
		  nfiles2,
		  &nSeg2,
		  &isBack2,
		  &jInd2,
		  &ntocs2,
		  &tocs2,
		  &bparams2,
		  &events2,
		  &jeId2,
		  &Nevents2,
		  &backInd2,
		  &foreInd2
		  )) {
    fprintf(stderr,"Error in ProcessFiles\n");
    return 1;
  }

  if(ntocs1 != ntocs2) {
    fprintf(stderr,"1: The two sets of files are incompatible\n");
    return 1;
  }

  /*
  if(nSeg1 != nSeg2) {
    fprintf(stderr,"2: The two sets of files are incompatible\n");
    return 1;
  }
  */

  foreNew1 = (int *)calloc(ntocs1*ntocs2, sizeof(int));
  backNew1 = (int *)calloc(ntocs1*ntocs2, sizeof(int));
  foreParams = (char **)calloc(ntocs1*ntocs2, sizeof(char *));
  infosFore = (char **)calloc(ntocs1*ntocs2, sizeof(char *));
  backParams = (char **)calloc(ntocs1*ntocs2, sizeof(char *));
  infosBack = (char **)calloc(ntocs1*ntocs2, sizeof(char *));

  jInd1Inv = (int **)calloc(ntocs1, sizeof(int *));
  NjInd1Inv = (int *)calloc(ntocs1, sizeof(int));
  for(j1=0;j1<ntocs1;j1++) {
    for(i1=0;i1<nSeg1;i1++) {
      if(j1 == jInd1[i1]) {
	(NjInd1Inv[j1])++;
	jInd1Inv[j1] = (int *)realloc(jInd1Inv[j1], NjInd1Inv[j1] * sizeof(int));
	jInd1Inv[j1][NjInd1Inv[j1]-1] = i1;
      }
    }
  }
  

  jInd2Inv = (int **)calloc(ntocs2, sizeof(int *));
  NjInd2Inv = (int *)calloc(ntocs2, sizeof(int));
  for(j2=0;j2<ntocs2;j2++) {
    for(i2=0;i2<nSeg2;i2++) {
      if(j2 == jInd2[i2]) {
	(NjInd2Inv[j2])++;
	jInd2Inv[j2] = (int *)realloc(jInd2Inv[j2], NjInd2Inv[j2] * sizeof(int));
	jInd2Inv[j2][NjInd2Inv[j2]-1] = i2;
      }
    }
  }

  for(j1=0;j1<ntocs1;j1++) {
    for(j2=0;j2<ntocs2;j2++) {      
      
      int uk2 = NjInd2Inv[j2];
      int *ui2 = jInd2Inv[j2];

      if(mixETGParams ||
	 !paramsdiffIFO(tocs1[j1],tocs2[j2],IFO1,IFO2)) {

	for(k1=0;k1<NjInd1Inv[j1];k1++) {
	  int b1t0s, b1t0ns;
	  char *b1params;

	  i1 = jInd1Inv[j1][k1];

	  b1t0s = bparams1[i1].t0s;
	  b1t0ns = bparams1[i1].t0ns;
	  b1params = bparams1[i1].params;

	  for(k2=0;k2<uk2;k2++) { /* this double loop is SLOW!! */
	    BurstSegParams *bpi2;
 
	    i2 = ui2[k2];
	    bpi2 = bparams2 + i2;

	    if(b1t0s == bpi2->t0s &&
	       b1t0ns == bpi2->t0ns &&
	       !paramsdiffcp(b1params,bpi2->params,0)
	       ) {

		kk = j1*ntocs2+j2;	  
		
		if(isBack1[i1] && isBack2[i2]) {

		  if(!backNew1[j1*ntocs2+j2]) {

		    backNew1[j1*ntocs2+j2] = 1;

		    (*nback)++;

		    backParams[kk] = (char *)malloc(backSize);
		    memcpy(backParams[kk], backParams0, backSize);

		    j1Back = (int *)realloc(j1Back, *nback * sizeof(int));
		    j2Back = (int *)realloc(j2Back, *nback * sizeof(int));
	    
		    j1Back[*nback-1] = j1;
		    j2Back[*nback-1] = j2;

		  }

		  if((*backFun)(infosBack + kk, events1[j1] + jeId1[i1], events2[j2] + jeId2[i2], backParams[kk], bparams1 + i1, bparams2 + i2)) {
		    fprintf(stderr,"Error in backFun\n");
		    return 1;
		  }

		} 

		if(!isBack1[i1] && !isBack2[i2]) {

		  if(!foreNew1[j1*ntocs2+j2]) {

		    foreNew1[j1*ntocs2+j2] = 1;
	    
		    (*nfore)++;

		    foreParams[kk] = (char *)malloc(foreSize);
		    memcpy(foreParams[kk], foreParams0, foreSize);

		    j1Fore = (int *)realloc(j1Fore, *nfore * sizeof(int));
		    j2Fore = (int *)realloc(j2Fore, *nfore * sizeof(int));

		    j1Fore[*nfore-1] = j1;
		    j2Fore[*nfore-1] = j2;

		  } 

		  if((*foreFun)(infosFore + kk, events1[j1] + jeId1[i1], events2[j2] + jeId2[i2], foreParams[kk], bparams1 + i1, bparams2 + i2)) {
		    fprintf(stderr,"Error in foreFun\n");
		    return 1;
		  }

		}
	      }
	    }
	  }
	}
      }
  }


  {
    BurstSegInfo *pbsi = backData;
      
    for(j=0;j<*nback;j++) {
      pbsi->next = (BurstSegInfo *)calloc(1,sizeof(BurstSegInfo));
      pbsi = pbsi->next;
      
      pbsi->params = (char *)calloc(strlen(tocs1[j1Back[j]]) + strlen(tocs2[j2Back[j]]) + 2, sizeof(char));
      sprintf(pbsi->params,"%s,%s",tocs1[j1Back[j]],tocs2[j2Back[j]]);

      pbsi->info = infosBack[j1Back[j]*ntocs2+j2Back[j]];
    }
  }

  {
    BurstSegInfo *pbsi = foreData;
    
    for(j=0;j<*nfore;j++) {
      pbsi->next = (BurstSegInfo *)calloc(1,sizeof(BurstSegInfo));
      pbsi = pbsi->next;
      
      pbsi->params = (char *)calloc(strlen(tocs1[j1Fore[j]]) + strlen(tocs2[j2Fore[j]]) + 2, sizeof(char));
      sprintf(pbsi->params,"%s,%s",tocs1[j1Fore[j]],tocs2[j2Fore[j]]);

      pbsi->info = infosFore[j1Fore[j]*ntocs2+j2Fore[j]];
    }
  }

  for(j1=0;j1<ntocs1;j1++) {
    free(jInd1Inv[j1]);
  }
  free(jInd1Inv);
  free(NjInd1Inv);

  for(j2=0;j2<ntocs2;j2++) {
    free(jInd2Inv[j2]);
  }
  free(jInd2Inv);
  free(NjInd2Inv);

  free(j1Fore);
  free(j1Back);
  free(j2Fore);
  free(j2Back);

  free(jeId1);
  free(jeId2);

  free(foreNew1);
  free(backNew1);

  free(infosFore);
  free(infosBack);

  for(j=0;j<ntocs1*ntocs2;j++) {
    if(backNew1[j]) {
      free(backParams[j]);
    }
    if(foreNew1[j]) {
      free(foreParams[j]);
    }
  }
  free(backParams);
  free(foreParams);

  for(j=0;j<ntocs1;j++) {
    for(i=0;i<Nevents1[j];i++) {
      SnglBurstTableC *eptr = (events1[j]+i)->next, *oeptr;
      while(eptr) {
	oeptr = eptr->next;
	free(eptr);
	eptr = oeptr;
      }
    }
    free(events1[j]);
    free(tocs1[j]);
  }

  free(events1);
  free(Nevents1);
  
  free(tocs1);
  free(foreInd1);
  free(backInd1);

  free(jInd1);
  free(isBack1);

  for(i=0;i<nSeg1;i++) {
    free(bparams1[i].params);
  }
  free(bparams1);


  for(j=0;j<ntocs2;j++) {
    for(i=0;i<Nevents2[j];i++) {
      SnglBurstTableC *eptr = (events2[j]+i)->next, *oeptr;
      while(eptr) {
	oeptr = eptr->next;
	free(eptr);
	eptr = oeptr;
      }
    }
    free(events2[j]);
    free(tocs2[j]);
  }

  free(events2);
  free(Nevents2);
  
  free(tocs2);
  free(foreInd2);
  free(backInd2);

  free(jInd2);
  free(isBack2);

  for(i=0;i<nSeg2;i++) {
    free(bparams2[i].params);
  }
  free(bparams2);
  

  return 0;
}

/******************************************************************/

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
