#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define ILWD_HEAD(fid,fidn) \
  sprintf(fname,"./%s.ilwd",fidn); \
  if(!(fid = fopen(fname,"w"))) { \
    fprintf(stderr,"File error; can't create %s.ilwd\n",fidn); \
    return 1; \
  } \
  fprintf(fid,"<?ilwd?>\n"); \
  fprintf(fid,"\t<ilwd comment='%s' name='%s::sequence' size='2'>\n",fidn,fidn); \
  fprintf(fid,"\t\t<lstring name='real:domain' size='4'>NONE</lstring>\n"); \
  fprintf(fid,"\t\t<real_4 ndim='2' dims='%u,%u' name='data' units='ADC'>",Ntot/Nps,Nps) \

#define ILWD_TAIL(fid) \
  fprintf(fid,"</real_4>\n</ilwd>\n"); \
  fclose(fid)


int main(int argc, char *argv[]) {

  FILE *hamp, *alpha, *delta, *psi, *injTime;

  char fname[256];

  double *h = NULL;
  int nh = 0;
  char *ptr;
  
  int i, j, k, Nps = 0, Ntot = 0;
  double dur;
  double skip = 2.0;

  srand48(time(NULL));

  if(argc != 5) {
    fprintf(stderr,"ShellInjection h N_per_seg N_total dur\n");
    return 1;
  }

  Nps = atoi(argv[2]);
  Ntot = atoi(argv[3]);
  dur = 16384.0*atof(argv[4]);

  ptr = strtok(argv[1],",");
  
  while(ptr) {
    nh++;
    h = (double *)realloc(h, nh * sizeof(double));

    h[nh-1] = atof(ptr);

    ptr = strtok(NULL,",");
  }

  ILWD_HEAD(hamp,"hamp");
  ILWD_HEAD(alpha,"alpha");
  ILWD_HEAD(delta,"delta");
  ILWD_HEAD(psi,"psi");
  ILWD_HEAD(injTime,"injTime");

  for(i=0;i<Ntot/(Nps * nh);i++) {

    for(j=0;j<nh;j++) {

      for(k=0;k<Nps;k++) {

	if(i+1<Ntot/(Nps * nh) ||
	   j+1<nh ||
	   k+1<Nps) {
	  fprintf(hamp,"%.15g,",h[j]);
	  fprintf(alpha,"%.15g,",6.283185307*drand48());
	  fprintf(delta,"%.15g,",-1.570796327+3.141592654*drand48());
	  fprintf(psi,"%.15g,",3.141592654*drand48());
	  fprintf(injTime,"%.15g,",skip*16384.0 + floor((dur-32768.0*skip)*drand48()));
	} else {
	  fprintf(hamp,"%.15g",h[j]);
	  fprintf(alpha,"%.15g",6.283185307*drand48());
	  fprintf(delta,"%.15g",-1.570796327+3.141592654*drand48());
	  fprintf(psi,"%.15g",3.141592654*drand48());
	  fprintf(injTime,"%.15g",skip*16384.0 + floor((dur-32768.0*skip)*drand48()));
	}

      }

    }

  }

  

  ILWD_TAIL(hamp);
  ILWD_TAIL(alpha);
  ILWD_TAIL(delta);
  ILWD_TAIL(psi);
  ILWD_TAIL(injTime);

  return 0;
}
