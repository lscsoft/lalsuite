#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <burstdso.h>


int getFrameCache(char *fQuery, 
		  char *dataserver) {

  char *p0, *p1;
  char *buf;

  /* try to delete file */
  unlink(CACHEFILENAME);


  buf = (char *)calloc(1+strlen(fQuery),sizeof(char));

  /* process line by line */

  p0 = fQuery;
  p1 = strchr(fQuery,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  while(strlen(buf)) {

    char type[256], IFO[256], channel[256], Times[256];
    char *alias = (char *)calloc(256, sizeof(char));
    char *q, *T0, *T1;

    if(buf[0]!='\n') {

#ifdef DEBUGBURST
      fprintf(stderr,"Processing: %s\n", buf);
#endif

      /* Examine line
	 Format: type IFO times channel alias */
      sscanf(buf,"%s\t%s\t%s\t%s\t%s",type,IFO,Times,channel,alias);
    
      q = strchr(Times,'-');
      if(!q) {
	fprintf(stderr,"Times: T0-T1\n");
	return 1;
      } 
      *q=0;
      T0 = Times;
      T1 = q+1;

      { /* get data */
	char tname[] = CACHEFILENAME;
	char cmd[1024];
	char *path;
	char *p;

	/* find the data */
	path = getenv("LSC_DATAGRID_CLIENT_LOCATION");
	if(!path) {
	  fprintf(stderr,"Environment variable LSC_DATAGRID_CLIENT_LOCATION not set\n");
	  return 1;
	}
      
	sprintf(cmd,"source %s/setup.sh; %s/ldg-client/bin/LSCdataFind --server %s --observatory %s --type %s --gps-start-time %s --gps-end-time %s --url-type file --lal-cache >> %s", path, path, dataserver, IFO, type, T0, T1, tname);

	if(system(cmd) == -1) {
	  fprintf(stderr,"system call failed\n");
	  perror("Error");
	  return 1;
	}

      }
    }

    /* update line being processed */
    if(p1) {
      p0 = p1+1;
      p1 = strchr(p0,'\n');
      bzero(buf,strlen(buf));
      if(p1) {
	memcpy(buf,p0,p1-p0+1);
      } else {
	strcpy(buf,p0);
      }
    } else {
      break;
    }

  }

  return 0;
}

int main(int argc, char *argv[]) {

  char *fQuery; /* frame data */

  char *times, *dataserver;

  if(argc<3) {
    fprintf(stderr,"buildcache dataserver framesQuery\n");
    return 1;
  }

  dataserver = argv[1];
  fQuery = argv[2];


  /*****************************************/
  /* parse parameters */
  /*****************************************/
  {
    /* frame query */
    struct stat buf;

    if(stat(fQuery,&buf)) {
      times = (char *)calloc(1 + strlen(fQuery), sizeof(char));
      strcpy(times, fQuery);
    } else {
      /* fQuery is a filename */
      FILE *in;
      char *p;

      if((in = fopen(fQuery,"r"))==NULL) {
	fprintf(stderr,"Can't open %s\n",fQuery);
	return 1;
      }

      p = times = (char *)calloc(1 + buf.st_size, sizeof(char));

      while((p=fgets(p,buf.st_size,in))) {
	p += strlen(p);
      }

      fclose(in);
    }
  }


  if(getFrameCache(times, dataserver)) {
    fprintf(stderr,"ERROR in getFrameCache\n");
    return 1;
  }

  return 0;
}
