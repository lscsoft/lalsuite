#include "burstdso.h"
#include <lal/LALStdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <strings.h>

#ifdef DEBUGBURST
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include <lalapps/lalapps.h>

NRCSID( INITSEARCHC, "InitSearch.c" );
RCSID( "InitSearch.c" );

/******** <lalLaTeX file="InitSearchC"> ********
\section{Command Line Arguments}
The \texttt{libldastfclusters.so} requires the following arguments to the wrapperAPI:
\begin{verbatim}

   argv[0] = "-filterparams"
   [1] = channel name; first two letters must be a valid IFO. 
   [2] = total number of data in analysis segment
   [3] = waveform(s) to inject; for each waveform name X, X_p and X_c must be passed from the datacond API.
   [4] = waveform multiplicative amplitude(s)
   [5] = right ascension(s) of waveform (rad, between 0 and 2pi)
   [6] = declination(s) of waveform (rad, between -pi/2 and pi/2)
   [7] = polarization angle(s) (rad)
   [8] = number of independent injections to perform (total)
   [9] = injection time(s) in number of samples. See notes for details.
   [10] = ETG to use (TFCLUSTERS,SLOPE)
   [11] = first ETG parameter
   [11+n] = nth ETG parameter

\end{verbatim}

\section{Notes}
\begin{itemize}
\item Parameters can be numbers, strings, lists, ranges or vectors.
\item Parameters can be a list of values in parantheses, for instance (1,2,3,4) or (H1,H2,L1,V,G). 
\item Numerical values in double parantheses give ranges: ([xlo,xhi,N]) gives N uniformly distributed values between xlo and xhi, ([xlo,xhi,N,r]) gives randomly distributed values (each trial of the simulation correspond to a different realization), ([xlo,xhi,N,l]) gives log distributed values. 
\item Numerical values without dots or "e" or "E" are assumed to be integers, unless they are in square brackets, where they are assumed to be floats (except for random values, which can be integers). 
\item String parameters bounded by underscores (e.g. \_a\_vector\_) are assumed to be names of REAL4 sequences passed by the datacondAPI. 
\item The beginning of the data segment defines the origine of time at the center of the Earth. argv[9] gives the delay between that time and the beginning of the signal time-series, in number of samples, and at the center of the Earth. Given the source position and the interferometer, an additional delay for the wave propagation is added.
\item argv[4-7] do not accept arguments of the form ([xlo,xhi,N,r]) with N>1. When any of these arguments is random and more than one waveform is injected per segment, the same realization of the random argument is used for all injected waveforms. If a list is passed, the list must have the same length for all argv[4-7], and the values of the various arguments will be 'locked' together, i.e. will be used as a single table.
\item If argv[5] or argv[6] are outside the allowed range, the waveform is injected with $F_+ = F_\times = 1$ and no time delay with respect to the earth center.
\item argv[4-7] and argv[9] can be matrices in order to perform coherent injections for many IFOs. They are passed as ilwd files of dimension (number of segments, injections per segment), using double underscores (e.g. \_\_file\_\_). argv[8] must then be equal to the product of the dimensions of the ilwd file.
\item The GW data must be passed with name containing GW\_STRAIN\_DATA .
\item The code uses LALExtractResponse to get the response function. Pass noResponse in the first four arguments to assume a unity response function.
\item The first argument can be set to binOutput while {\tt argv[i] -> argv[i+1]} in order to save the output in a binary format. This leads to up to 2x improvements in speed under LDAS.
\item The first or second argument can alternatively be set to fileOutput:path in order to bypass the ligolwAPI completely and speed up large jobs.
\item The first or second argument can be set to noLALBurstOutput in order to bypass the burst parameter estimation function. The unaltered output of the ETG is then reported.
\item If argv[9] is a list of injection times (esp. random), only those injections for which the signals do not overlap will be performed and reported. The check for overlapping signals is not performed for matrix inputs.
\item The first argument can be set to fileOutput:file in order to bypass the LIGOLWAPI and dump data in binary format to file.
\end{itemize}
********* </lalLaTeX> ********/

#define CBUFSIZE 128

int InitSearch(char *filterParams, 
	       BurstSearchParams *params) {

  /*
    Arguments:
    filterParams: string with filter parameters for the ETG
    params: structure with all parameters needed by the ETG and wrapper
  */

  static LALStatus status;  /* LAL status */

  UINT4 i;  /* counter */
  CHAR cbuf[CBUFSIZE];  /* scratch buffer */
  BurstParameterList *tmpBPL;

  INT4 argc;    /* number of arguments in filterParams */
  CHAR **argv;  /* vector of strings in filterParams */


  /* Note: define DEBUGBURST to have debug info written to disk */
#ifdef DEBUGBURST
 {
   char *uname = getenv("USER");
   char *hname = getenv("HOST");
   char fname[1024];
   pid_t pid = getpid();
   fprintf(stderr,"%s\n",uname);
   /*
   sprintf(fname,"/dso-test/%s%s%uslave.err",hname,uname,pid);
   freopen(fname, "w", stderr);
   fchmod(fileno(stderr),S_IRUSR|S_IRGRP|S_IROTH|S_IWOTH|S_IWGRP|S_IWUSR);
   */
 }

 fprintf(stderr,"Init.......\n");
 fflush(NULL);
#endif


 /*
#ifdef DEBUGBURST
  fprintf(stderr,"Sleeping; my pid is %u\n", getpid());
  usleep(10000000);
  fprintf(stderr,"Waking up\n");
  fflush(NULL);
#endif
 */


 /* initialize params */
 bzero(params, sizeof(BurstSearchParams));


 /* Create argc, argv from filterParams string */
 /* comma separated strings between () parantheses are treated
    as single arguments, and the parantheses are stripped. */
 {
   char *p0 = filterParams, *p1;

   argc = 0;
   argv = NULL;

   if(*p0 == '(') {
     /* list of elements between () parantheses */
     p1 = strchr(p0,')');
     if(!p1) {
       /* no closing ) */
       fprintf(stderr,"Invalid filterParams: %s\n",filterParams);
       return 1;
     }
     p1 = strchr(p1,','); /* p1 is 1st comma after ) */
   } else {
     /* no parantheses, p1 is 1st comma after p0 */
     p1 = strchr(p0,',');
   }

   /* loop over all arguments */
   while(p1) {

     /* allocate memory for new argument */
     argc++;
     argv = (CHAR **)realloc(argv,argc*sizeof(CHAR *));
     argv[argc-1] = (CHAR *)calloc(1+p1-p0,sizeof(CHAR));

     /* copy data */
     if(*p0 != '(') {
       /* no parantheses */
       memcpy(argv[argc-1], p0, p1-p0);
     } else {
       /* parantheses, remove 1st and last character */
       memcpy(argv[argc-1], p0+1, p1-p0-2);
     }

     p0=p1+1; /* move to next argument */

     /* deal with parantheses */
     if(*p0 == '(') {
       p1 = strchr(p0,')');
       if(!p1) {
	 fprintf(stderr, "Invalid filterParams: %s\n",filterParams);
	 return 1;
       }
       p1 = strchr(p1,',');
     } else {
       p1 = strchr(p0,',');
     }
   }

   /* add last argument */
   argc++;
   argv = (CHAR **)realloc(argv,argc*sizeof(CHAR *));
   argv[argc-1] = (CHAR *)calloc(strlen(p0),sizeof(CHAR));

   if(*p0 != '(') {
     /* no parantheses */
     strcpy(argv[argc-1], p0);
   } else {
     /* parantheses, remove 1st and last character */
     char *pl;
     strcpy(argv[argc-1], p0+1);     
     pl  = strrchr(argv[argc-1],')');
     if(pl) {
       *pl = 0;
     } else {
       SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
     }
   }

 }


 /*
  *
  * Sanity checks on the input, output, params
  *
  */

 /* check number of inputs */
  if ( argc < 11){
        SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }
 
  /* check first argument is -filterparams */
  if (strcmp( argv[0], "-filterparams" )){
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS);
  }

  /* is 1st or 2nd argument is binOutput, set flag */
  if(!strcmp(argv[1],"binOutput") ||
     !strcmp(argv[2],"binOutput")) {
    (argv)++;
    (argc)--;
    params->outputMethod = 1;

#ifdef DEBUGBURST
    fprintf(stderr,"binOutput: on\n");
    fflush(NULL);
#endif

  }

#ifdef DEBUGBURST
    fprintf(stderr,"argv[1]:%s\n",argv[1]);
    fflush(NULL);
#endif

    /*********************************************************************/
    /* Parse arguments from filterParams */
    /*********************************************************************/

    /* if 1st argument is fileOutput: */
  if(strstr(argv[1],"fileOutput")) {

    /* extract file name */
    char *p = strchr(argv[1],':')+1;
    strncpy(params->dsoOutput,p,1024);

    /* remove from list */
    (argv)++;
    (argc)--;

    /* update flag */
    params->outputMethod = 2;

#ifdef DEBUGBURST
    fprintf(stderr,"fileOutput: on -> %s\n",params->dsoOutput);
    fflush(NULL);
#endif
  }

  /* same for 2nd argument */
  if(strstr(argv[2],"fileOutput")) {

    char *p = strchr(argv[2],':')+1;

    strncpy(params->dsoOutput,p,1024);

    (argv)++;
    (argc)--;
    params->outputMethod = 2;

#ifdef DEBUGBURST
    fprintf(stderr,"fileOutput: on -> %s\n",params->dsoOutput);
    fflush(NULL);
#endif
  }


  /* if 1st or 2nd argument is noLALBurstOutput, update flag */
  if(!strcmp(argv[1],"noLALBurstOutput") ||
     !strcmp(argv[2],"noLALBurstOutput")) {
    (argv)++;
    (argc)--;
    params->noLALBurstOutput = 1;

#ifdef DEBUGBURST
    fprintf(stderr,"noLALBurstOutput: 1\n");
    fflush(NULL);
#endif

  }


  /* if one of first 4 arguments is noResponse, set flag */
  if(!strcmp(argv[1],"noResponse") ||
     !strcmp(argv[2],"noResponse") ||
     !strcmp(argv[3],"noResponse") ||
     !strcmp(argv[4],"noResponse")) {
    (argv)++;
    (argc)--;
    params->allowNoResponse = 1;

#ifdef DEBUGBURST
    fprintf(stderr,"allowNoResponse: 1\n");
    fflush(NULL);
#endif

  }

  /* 1st argument: channel */
  strncpy(params->channel, argv[1], LIGOMETA_CHANNEL_MAX-1);

  /* 2nd argument: number of data points */
  params->Ndata = atoi(argv[2]);
  if ( params->Ndata <= 0 ) {
    SABORT(INITSEARCHH_ESPOS, INITSEARCHH_MSGESPOS );
  }

  /* 3rd argument: set of waveforms */
  LAL_CALL ( ParseParameter( &status, &(params->waveforms), argv[3], &(params->Nwaveforms) ), &status);

  if(params->waveforms.random) {
    /* waveforms cannot be a random parameter */
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }

  /* 4th argument: waveform amplitudes */
  LAL_CALL ( ParseParameter( &status, &(params->waveform_amplitudes), argv[4], &(params->Nwaveform_amplitude)), &status );

  if(params->Nwaveform_amplitude>1 && params->waveform_amplitudes.random == 1) {
    /* if random, should have only one number */
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }

  /* 5th argument: alpha */
  LAL_CALL ( ParseParameter( &status, &(params->alpha), argv[5], &(params->Nalpha)), &status );

  if(params->Nalpha>1 && params->alpha.random == 1) {
    /* if random, should have only one number */
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }

  /* 6th argument: delta */
  LAL_CALL ( ParseParameter( &status, &(params->delta), argv[6], &(params->Ndelta)), &status );

  if(params->Ndelta>1 && params->delta.random == 1) {
    /* if random, should have only one number */
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }

  /* 7th argument: psi */
  LAL_CALL ( ParseParameter( &status, &(params->psi), argv[7], &(params->Npsi)), &status );

  if(params->Npsi>1 && params->psi.random == 1) {
    /* if random, should have only one number */
    SABORT(INITSEARCHH_EARGS, INITSEARCHH_MSGEARGS ); 
  }

  /* 8th argument: number of injections */
  params->Ninj = atoi(argv[8]);

  /* 9th argument: injection times */
  LAL_CALL ( ParseParameter( &status, &(params->injTimes), argv[9], &(params->Ninj_per_segment)), &status );

  /* 10th argument: ETG */
  strncpy(cbuf,argv[10],CBUFSIZE);

  /* interpret string */
  if(!strcmp(cbuf,"TFCLUSTERS")) {
    params->ETG = TFCLUSTERS;
  } else if(!strcmp(cbuf,"SLOPE")) {
    params->ETG = SLOPE;
  } else {
    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
  }


  /* Create linked list of ETG parameters: */
  tmpBPL = &(params->ETGparams);
  params->Nparams = 0;
  params->NNparams = NULL;

  /* loop over arguments after the 10th argument */
  for(i=11; i<argc; i++) {
    /* allocate memory */
    tmpBPL->next = (BurstParameterList *)LALCalloc(1,sizeof(BurstParameterList));
    SASSERT(tmpBPL->next, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
    tmpBPL = tmpBPL->next;

    (params->Nparams)++;
    params->NNparams = (UINT4 *)LALRealloc(params->NNparams, params->Nparams*sizeof(UINT4));

    /* parse the argument string */
    LAL_CALL ( ParseParameter( &status, &(tmpBPL->param), argv[i], &(params->NNparams[params->Nparams-1]) ), &status );


  }

  /*********************************************************************/

  /* assign proper ETG function */
  switch(params->ETG) {
  case TFCLUSTERS:
    params->ETGfun = LALTFClustersETG;
    break;
  case SLOPE:
    params->ETGfun = LALSlopeETG; 
    break;
  default:
    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
  }

  
  /* no searchMaster */
  params->searchMaster = 0;

  /* sanity check */
  if(params->Nwaveforms==0 && ( params->Nwaveform_amplitude || params->Nalpha || params->Ndelta || params->Ninj || params->Ninj_per_segment)) {
    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
  }

#ifdef DEBUGBURST
      fprintf(stderr,"%u %u %u %u\n",params->Nwaveform_amplitude,params->Nalpha,params->Ndelta,params->Npsi);
      fflush(NULL);
#endif

      /*
  if(params->Nwaveform_amplitude) {
    if(params->Nwaveform_amplitude != params->Nalpha ||
       params->Nwaveform_amplitude != params->Ndelta ||
       params->Nwaveform_amplitude != params->Npsi) {

#ifdef DEBUGBURST
      fprintf(stderr,"Error: %u %u %u %u\n",params->Nwaveform_amplitude,params->Nalpha,params->Ndelta,params->Npsi);
      fflush(NULL);
#endif

      ABORT(status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
    } else {
      params->Nalpha = params->Ndelta = params->Npsi = 1;
    }    
  }
      */
				    

      /* Initialize data structure */
      
  if(params->Ninj_per_segment == 0) {
    params->Ninj_per_segment = 1;
  }

  bzero(&(params->data), sizeof(BurstOutputDataSegment));

  bzero(&(params->oparams), sizeof(BurstOutputParameters));
  params->oparams.data = &(params->data);
  params->oparams.method = 1;


  /* restore pointer to allow clean-up */
  if(params->outputMethod) {
    (argv)--;
    (argc)++;
  }

  if(params->noLALBurstOutput) {
    (argv)--;
    (argc)++;
  }

  if(params->allowNoResponse) {
    (argv)--;
    (argc)++;
  }


  /* clean up */
  {
    int i;
    for(i=0;i<argc;i++) {
      free(argv[i]);
    }
    free(argv);
  }

  return 0;
}




/* The following macro is used to parse the string str and to
   allocate memory for the new parameter, and to add it to the
   link list bp.

   With the __FILE__ notation, bp->random is set to 2, and 
   bp->char_ is set to FILE.

   With the _FILE_ notation, bp->real4vector_ is allocated, 
   and bp->char_ is set to FILE.

   In the absence of underscores, the string is simply parsed:
   only digits: int4_
   only digits and +-.eE: real4_
   otherwise: char_
*/
#define SetParams(bp,str) if(str[0]=='_' && str[strlen(str)-1]=='_') { \
 BOOLEAN mat = 0; \
 if(strlen(str)>=5) { \
  if(str[1]=='_' && str[strlen(str)-2]=='_') { \
   mat = 1; \
   bp->random = 2; \
 }} \
 if(!mat) { \
  bp->real4vector_ = (REAL4Vector *)LALCalloc(1,sizeof(REAL4Vector)); \
  RASSERT( bp->real4vector_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
										 } \
 bp->char_ = (CHAR *)LALMalloc((1+strlen(str))); \
 RASSERT(bp->char_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
 strcpy(bp->char_,str); \
} else {\
for(i=0;i<strlen(str);i++) { \
if(!isdigit((INT4)(str[i])) && str[i]!='e' && str[i]!='E' && str[i]!='+' && str[i]!='-' && str[i]!='.') {break;}} \
if(i<strlen(str)) { \
 bp->char_ = (CHAR *)LALMalloc((1+strlen(str))); \
 RASSERT(bp->char_, status,INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
 strcpy(bp->char_,str); \
} else { \
 for(i=0;i<strlen(str);i++) { \
  if(str[i]=='e' || str[i]=='E' || str[i]=='.') {break;}} \
 if(i<strlen(str)) { \
  bp->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4)); \
  RASSERT(bp->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
  *(bp->real4_) = strtod(str, NULL); \
 } else {\
  bp->int4_ = (INT4 *)LALMalloc(sizeof(INT4)); \
  RASSERT(bp->int4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
  *(bp->int4_) = atoi(str); \
 } \
     } \
	 }

/* This function parse the string argv, and returns the result in par and M */
void ParseParameter( 
		    LALStatus *status,
		    BurstParameter *par,
		    CHAR *argv,
		    UINT4 *M
		    ) {

  /* 
     Arguments:
     status: LAL status pointer
     par: pointer to returned BurstParameter structure
     argv: string to be parsed
     M: number of parameters corresponding to the string (e.g. if [1,2,4], M=4)
  */

  INT4 i,N;
  CHAR *ptr;
  REAL4 x, xlo, xhi;
  BurstParameter *params = par;

  INITSTATUS (status, "ParseParameter", INITSEARCHC);

#ifdef DEBUGBURST
  fprintf(stderr,"argv: %s\n",argv);
#endif

  *M = 0;

  if(!strlen(argv)) {
    /* empty string, do nothing */
    RETURN (status);
  }


  /* first check if we have a range: */
  if(argv[0] == '[') {

    /* check we have closing bracket: */
    if(argv[strlen(argv)-1] != ']') {
      ABORT(status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
    }

    /* use strtok to process the string;
       warning: not thread safe, etc. */
    ptr = strtok(argv,"[,]");

    /* check non-null pointer */
    RASSERT(ptr, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);

    /* set xlo to first element in range, i.e. to lower bound */
    xlo = strtod(ptr, NULL);

    /* second element is upper bound */
    ptr = strtok(NULL,"[,]");
    RASSERT(ptr, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
    xhi = strtod(ptr, NULL);
    
    /* third argument is number of elements */
    ptr = strtok(NULL,"[,]");
    RASSERT(ptr, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
    *M = N = atoi(ptr);

#ifdef DEBUGBURST
  fprintf(stderr,"xlo: %g\n",xlo);
  fprintf(stderr,"xhi: %g\n",xhi);
  fprintf(stderr,"N: %i\n",N);
#endif

  /* check for 4th element */
    ptr = strtok(NULL,"[,]");
    if(!ptr) {

      /* No 4th element */

      /* loop to fill in the values */
      for(i=0;i<N;i++) {
	if(N>1) {
	  x = xlo + (xhi-xlo)*(REAL4)i/(REAL4)(N-1);
	} else {
	  x = xhi;
	}

	/* allocate memory for output (link list) */
	params->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	RASSERT(params->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	params = params->next;

	/* type of each element in link list is real4 */
	params->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4)); 
	RASSERT(params->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); 
	*(params->real4_) = x;
      }
    } else {

      /* Have a 4th element, switch it */

      switch(*ptr) {

      case 'l':
	/* log-spaced interval */
	for(i=0;i<N;i++) {
	  if(N>1) {
	    x = xlo*pow(xhi/xlo,(REAL4)i/(REAL4)(N-1));
	  } else {
	    x = xhi;
	  }

	  /* allocate memory for output (link list) */
	params->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	RASSERT(params->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	params = params->next;

	/* type of each element in link list is real4 */
	params->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4)); 
	RASSERT(params->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); 
	*(params->real4_) = x;
	}
	break;

      case 'r':

	/* random interval */

	params->random = 1; /* set flag */

	/* setup link list; 
	   1st element is lower bound of random interval (xlo),
	   2nd element is upper bound (xhi) */

	params->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	RASSERT(params->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	params = params->next;

	if(fabs((xlo - floor(xlo))/xlo) < 1E-15 ||
	   fabs((xlo - ceil(xlo))/xlo) < 1E-15) {
	  /* xlo is integer type */
	  params->int4_ = (INT4 *)LALMalloc(sizeof(INT4));
	  RASSERT(params->int4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	  *(params->int4_) = (INT4)xlo;
	} else {
	  /* xlo is a REAL4 type */
	  params->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4));
	  RASSERT(params->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	  *(params->real4_) = (REAL4)xlo;
	}

	params->random = 1;

	/* xhi: */
	params->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	RASSERT(params->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	params = params->next;

	if(fabs((xhi - floor(xhi))/xhi) < 1E-15 ||
	   fabs((xhi - ceil(xhi))/xhi) < 1E-15) {
	  /* integer */
	  params->int4_ = (INT4 *)LALMalloc(sizeof(INT4));
	  RASSERT(params->int4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	  *(params->int4_) = (INT4)xhi;
	} else {
	  /* float */
	  params->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4));
	  RASSERT(params->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	  *(params->real4_) = (REAL4)xhi;
	}

	params->random = 1;

	break;
      default:
	RASSERT(ptr, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	break;
      }
    }
  } else {

    /* we don't have a range delimited by square brackets */
    /* NOTE: the wrapper is responsible for stripping the () parantheses */

    CHAR *tptr = argv;

    /* loop over elements separated by commas: */
    while((ptr = strtok(tptr,","))) {
      tptr = NULL;

      /* allocate memory (link list) */
      params->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
      RASSERT(params->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
      params = params->next;

      /* parse parameter using macro */
      SetParams(params,ptr)

      (*M)++;
    }

  }

  /* Normal exit */

  RETURN (status);
}

#undef SetParams
