#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "datacondAPI/DatacondCaller.h"

#include <burstdso.h>

#include <lalapps/lalapps.h>


int main(int argc, char *argv[]) {

  char nocondor = 0;

  char *fQuery; /* frame data */

  char *algorithms; /* algorithms */

  char *filterParams; /* params */

  char *rFiles = NULL; /* response files */

  float f0=0.0, f1=1.0; /* portion of the job to do */

  int Nsymbols = 0;
  static datacond_symbol_type *symbols = NULL; /* datacond symbols */

  BurstSearchParams bparams; /* burst search parameters */

  char *times, *algo, *params, *fraction, *resp, *dataserver;

  /* Which error handler to use */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "227" ); /* no memory debugging */

  /* look at command line arguments */
  if(argc<5) {
    fprintf(stderr,"burstdso [-nocondor] dataserver framesQuery algorithms filterParams [responseFiles] [fraction]\n");
    fprintf(stderr,"frameQuery, algorithms, filterParams and responseFiles can be files or strings\n");
    fprintf(stderr,"Pass an empty string (\"\") for responseFiles to use fraction without responseFiles\n");
    fprintf(stderr,"fraction: f0-f1\n");
    return 1;
  }

  if(strstr(argv[1],"nocondor")) {
    argc -= 1;
    argv++;
    nocondor=1;
  }

  dataserver = argv[1];
  fQuery = argv[2];
  algo = argv[3];
  params = argv[4];




  /*****************************************/
  /* parse parameters */
  /* If the parameter is a string, it is stored directly. 
     If, however, the parameter is a file, the the content
     of the file is stored in the appropriate variable. */
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


  {
    /* algorithms */
    struct stat buf;

    if(stat(algo,&buf)) {
      algorithms = (char *)calloc(1 + strlen(algo), sizeof(char));
      strcpy(algorithms, algo);
    } else {
      /* algo is a filename */
      FILE *in;
      char *p;

      if((in = fopen(algo,"r"))==NULL) {
	fprintf(stderr,"Can't open %s\n",algo);
	return 1;
      }

      p = algorithms = (char *)calloc(1 + buf.st_size, sizeof(char));

      while((p=fgets(p,buf.st_size,in))) {
	p += strlen(p);
      }

      fclose(in);
    }
  }

  {
    /* filterParams */
    struct stat buf;

    if(stat(params,&buf)) {
      filterParams = (char *)calloc(1 + strlen(params), sizeof(char));
      strcpy(filterParams, params);
    } else {
      /* params is a filename */
      FILE *in;
      char *p;

      if((in = fopen(params,"r"))==NULL) {
	fprintf(stderr,"Can't open %s\n",params);
	return 1;
      }

      p = filterParams = (char *)calloc(1 + buf.st_size, sizeof(char));

      while((p=fgets(p,buf.st_size,in))) {
	p += strlen(p);
      }

      fclose(in);
    }
  }

  if(argc>5)
  {
    /* response Files */
    struct stat buf;

    resp = argv[5];

    if(stat(resp,&buf)) {
      if(strlen(resp)>0) {
	rFiles = (char *)calloc(1 + strlen(resp), sizeof(char));
	strcpy(rFiles, resp);
      }
    } else {
      /* resp is a filename */
      FILE *in;
      char *p;

      if((in = fopen(resp,"r"))==NULL) {
	fprintf(stderr,"Can't open %s\n",resp);
	return 1;
      }

      p = rFiles = (char *)calloc(1 + buf.st_size, sizeof(char));

      while((p=fgets(p,buf.st_size,in))) {
	p += strlen(p);
      }

      fclose(in);
    }
  }

  if(argc>6) {
    fraction = argv[6];
    sscanf(fraction,"%g-%g",&f0,&f1);
  }


#ifdef DEBUGBURST
  fprintf(stderr,"%s\n%s\n%s\n%s\n%g\n%g\n",fQuery,algorithms,filterParams,rFiles,f0,f1);
#endif

  /*****************************************/
  /* prepare output symbols */
  /* For each output statement in the algorithms, we
     have to allocate an entry in the symbols table
     for the stand-alone datacondAPI, so that the data
     can be sent back to the analysis code. In the 
     output statement, the comment field must comtain
     the datatype of the output data in square brackets,
     for instance:
     output(cavfaccplx,_,_,H2:CAL-CAV_FAC,H2 cavity factor [COMPLEX8TimeSeries]);
     Note that this is not necessary if the output type
     is a single precision time series. */
  /*****************************************/
  if(OutputSymbols(algorithms, &Nsymbols, &symbols)) {
    fprintf(stderr,"ERROR: OutputSymbols\n");
    return 1;
  }

  /*****************************************/
  /* acquire data */
  /*****************************************/
  /* loop over lines in fQuery, get data, add to datacond callchain */
  /* Format is:
     frame_type Site [file] GPS_Start-GPS_End channel alias
     file is optional, and is a frame file that overrides LSCdataFind
     channel can be adc(H2:LSC-AS_Q), proc(H2:CAL-OLOOP_FAC), or proc(H2:CAL-RESPONSE!0!7000.0001!), for instance; same format as LDAS
     alias is a name to represent the data in the datacondAPI
  */

  if(nocondor) {

    if(getFrameData(times, dataserver, &Nsymbols, &symbols)) {
      fprintf(stderr,"ERROR: can't get frame data\n");
      return 1;
    }

  } else {
    
    if(getFrameData(times, dataserver, &Nsymbols, &symbols)) {
      fprintf(stderr,"ERROR: can't get frame data\n");
      return 1;
    }

  }

  /*****************************************/
  /* get non-frame data */
  /*****************************************/
  /* loop over lines in rFiles, add to datacond callchain & algorithm */
  /* format is:
     file push/pass datacond_symbol
     pass option not supported yet
     file can be a URL */

  if(rFiles) {
    if(getNonFrameData(rFiles, &algorithms, &Nsymbols, &symbols)) {
      fprintf(stderr,"ERROR: can't get non-frame data\n");
      return 1;
    }
  }

  /*****************************************/
  /* condition data (datacondAPI) */
  /*****************************************/
  if(ConditionData(Nsymbols, symbols, algorithms)) {
    fprintf(stderr,"ERROR: dataconditioning failed\n");
    return 1;
  }

  /*****************************************/
  /* prepare search */
  /*****************************************/
  /* setup parameters */
  if(InitSearch(filterParams, &bparams)) {
    fprintf( stderr, "ERROR: InitSearch\n");
    return 1;
  }
  
  /* run calibration, etc. */
  if(ReConditionData(Nsymbols, symbols, algorithms, &bparams)) {
    fprintf(stderr,"ERROR: dataconditioning failed, level 2\n");
    return 1;
  }

  /*****************************************/
  /* run search */
  /*****************************************/
  if(RunSearch(&bparams, f0, f1)) {
    fprintf( stderr, "ERROR: RunSearch\n");
    return 1;
  }

  /*****************************************/
  /* clean up */
  /*****************************************/
  free(times);
  free(algorithms);
  free(filterParams);

  if(rFiles) free(rFiles);

  symbolsfree(symbols);

  bparamsfree(&bparams);

  return 0;
}
