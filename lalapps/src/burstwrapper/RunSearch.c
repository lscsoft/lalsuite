#include "burstdso.h"

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/Inject.h>
#include <lal/Random.h>
#include <lal/SimulateCoherentGW.h>
#include <strings.h>
#include <unistd.h>
/* #include <BuildDB.h> */
#include <signal.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <time.h>
#include <math.h>
/* #include <LALWrapperConfig.h> */
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "burstdso.h"
#include <lalapps/lalapps.h>

/* largest possible size for TOC */
#define MAXCOMMENT 16384

NRCSID( RUNSEARCHC, "RunSearch.c" );
RCSID( "RunSearch.c" );

/* LAL status used by signal handler */
LALStatus *gStatus;

/* signal handler */
void BurstSigHandler(int);

void BuildOutputFileBinary(LALStatus *status,
			   int fid,
			   EventIDColumn     *burstEventList,
			   CHAR               *tag,
			   UINT4 t0Sec,
			   UINT4 t0NanoSec
			   );

#ifndef WORDS_BIGENDIAN
static void endian_swap(char * pdata, int dsize, int nelements);
#endif

int RunSearch(BurstSearchParams *params,
	      double f0, 
	      double f1) {

  /*
    Applies ETG and PEst functions, write output.
    
    Arguments:
    params: burst search parameters
    f0, f1: fractional range of jobs to do; (f0,f1)=(0,1) does everything, 
    (0,0.5) does first half.
  */

  static LALStatus status;

  EventIDColumn boutput,  /* output from ETG */
    stdOutput;            /* output from PEst */

  /* variables used to parse the burst parameters: */
  BurstParameter bparams;
  BurstParameter *bparamsptr, *bptr, *wave=NULL, *wamp;
  BurstParameterList *ETGparams;
  REAL4 *wbuf = NULL;
  REAL4 wAmp, Alpha, Delta, Psi;
  REAL4Vector *wAmpv, *Alphav, *Deltav, *Psiv;
  BOOLEAN got_wid = 0;
  CHAR comment[MAXCOMMENT];


  UINT4 j, runp, runs;

  /* variables for random numbers */
  RandomParams *rpar = NULL;
  REAL4 arandomnumber;

  /* job handling variables */
  int fid = -1;
  unsigned int i, si, silo, sihi,
    jobsTotal = 0,
    jobsParams = 0;

  /*************************************************************/
  /* Initialize */
  /*************************************************************/
  
  if(params->outputMethod == 2) {

    /* output to a file */

#ifdef DEBUGBURST
    fprintf(stderr,"file to open: %s\n",params->dsoOutput);
    fflush(NULL);
#endif

    /* open with right permission */
    fid = open(params->dsoOutput,O_WRONLY|O_CREAT|O_TRUNC,S_IRUSR|S_IRGRP|S_IROTH|S_IWOTH|S_IWGRP|S_IWUSR);

    /* check error */
    if(fid == -1) {
#ifdef DEBUGBURST
      perror("ERROR from fopen:");
      fflush(NULL);
#endif

      SABORT(INITSEARCHH_EFILE,INITSEARCHH_MSGEFILE);
    }

#ifdef DEBUGBURST
    fprintf(stderr,"MASTER: fid = %i\n",fid);
    fflush(NULL);
#endif
      
  } else {
    /* only support output to file */
    SABORT(INITSEARCHH_EFILE,INITSEARCHH_MSGEFILE);
  }

  
  /* Compute number of jobs */

  /* number of jobs from injections: */
  if(params->MAXinj != 0) {
    jobsTotal = 1 + params->Nwaveforms * params->Ninj / params->Ninj_per_segment;
  } else {
    jobsTotal = 1 + params->Nwaveforms * params->Nwaveform_amplitude * params->Nalpha * params->Ndelta * params->Npsi * params->Ninj / params->Ninj_per_segment;
  }

  if(params->Nwaveforms && !jobsTotal) {
    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
  }

  if(!jobsTotal) {
    jobsTotal = 1;
  }

  /* number of jobs from ETG parameters: */
  jobsParams = 1;
  for(i=0;i<params->Nparams;i++) {
    jobsTotal *= params->NNparams[i];
    jobsParams *= params->NNparams[i];
  }

  /* initialize random number generator */
  LAL_CALL ( LALCreateRandomParams(&status, &rpar, getpid() * clock()), &status );

  /*************************************************************/
  /*************************************************************/
  /* loop over jobs */
  /*************************************************************/
  /*************************************************************/

  /* low and high index of jobs to run */
  silo = (unsigned int)floor(f0 * (double)jobsTotal);
  sihi = (unsigned int)ceil(f1 * (double)jobsTotal);

  /* loop over job Id called si */
  for(si=silo; si<sihi; ++si) {

    /* variables used to construct job parameters */
    UINT4 wid, waid, aid, did, pid, Nid;
    UINT4 *Pid;

    wAmp = Alpha = Delta = Psi = 0.0;
    wAmpv = Alphav = Deltav = Psiv = NULL;

    wid = waid = aid = did = pid = Nid = 0;

    /* initialize */
    bzero(&boutput, sizeof(EventIDColumn));
    bzero(&stdOutput, sizeof(EventIDColumn));
    bzero(&bparams, sizeof(BurstParameter));

    /*************************************************************/
    /* convert job Id si into a set of parameters to run the job */
    /*************************************************************/

    /* injections */
    if(si < jobsParams) {
      Nid = 0;
      runs = si;
      runp = jobsParams;	  
    } else {

      si -= jobsParams;

      if(params->Nwaveforms) {
	got_wid = 1;

	runp = (jobsTotal - jobsParams) / params->Nwaveforms;
	wid = si / runp;

	runs = si - wid * runp;
	runp /= params->Nwaveform_amplitude;
	waid = runs / runp;
	  
	wid++;

	runs -= waid * runp;
	runp /= params->Nalpha;
	aid = runs / runp;

	runs -= aid * runp;
	runp /= params->Ndelta;
	did = runs / runp;

	runs -= did * runp;
	runp /= params->Npsi;
	pid = runs / runp;

	runs -= pid * runp;
	runp /= params->Ninj / params->Ninj_per_segment;
	Nid = runs / runp;
      } else {
	/*
	  Nid = 0;
	  runs = si;
	  runp = jobsTotal;
	*/
	SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
      }
    }

#ifdef DEBUGBURST
    fprintf(stderr,"wid = %u, waid = %u, aid = %u, did = %u, pid = %u, Nid = %u\n",wid,waid,aid,did,pid,Nid);
    fflush(NULL);
#endif


    Pid = (UINT4 *)LALMalloc(params->Nparams * sizeof(UINT4));
    TASSERT(Pid, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

    {
      UINT4 Mid = Nid;

      for(i=0; i<params->Nparams; i++) {
	runs -= Mid * runp;
	runp /= params->NNparams[i];
	Mid = Pid[i] = runs / runp;
      }
    }

    if(wid) {
      wave = &(params->waveforms);
      for(i=0; i<wid; i++) {
	wave = wave->next;
      }

      wamp = params->waveform_amplitudes.next;
      if(wamp) {
	if(wamp->random == 1) {
	  REAL4 lo,hi;
	      
	  if(wamp->real4_) {
	    lo = *(wamp->real4_);
	  } else {
	    lo = (REAL4)(*(wamp->int4_));
	  }
	      
	  wamp = wamp->next;
	    
	  if(wamp->real4_) {
	    hi = *(wamp->real4_);
	  } else {
	    hi = (REAL4)(*(wamp->int4_));
	  }

	  LAL_CALL ( LALUniformDeviate(&status, &arandomnumber, rpar), &status);
	  wAmp = lo + (hi-lo)*arandomnumber;

	} else {
	  if(params->MAXinj) {
	    waid = Nid;
	  }
	  for(i=0; i<waid; i++) {
	    wamp = wamp->next;
	  }
	  if(wamp->real4vector_) {
	    wAmpv = wamp->real4vector_;
	  } else if(wamp->real4_) {
	    wAmp = *(wamp->real4_);
	  } else {
	    wAmp = (REAL4)(*(wamp->int4_));
	  }
	}
      }
	  
      wamp = params->alpha.next;
      if(wamp) {
	if(wamp->random == 1) {
	  REAL4 lo,hi;
	      
	  if(wamp->real4_) {
	    lo = *(wamp->real4_);
	  } else {
	    lo = (REAL4)(*(wamp->int4_));
	  }

	  wamp = wamp->next;
	      
	  if(wamp->real4_) {
	    hi = *(wamp->real4_);
	  } else {
	    hi = (REAL4)(*(wamp->int4_));
	  }

	  LAL_CALL ( LALUniformDeviate(&status, &arandomnumber, rpar), &status);
	  Alpha = lo + (hi-lo)*arandomnumber;

	} else {
	  if(params->MAXinj) {
	    aid = Nid;
	  }
	  for(i=0; i<aid; i++) {
	    wamp = wamp->next;
	  }
	  if(wamp->real4vector_) {
	    Alphav = wamp->real4vector_;
	  } else if(wamp->real4_) {
	    Alpha = *(wamp->real4_);
	  } else {
	    Alpha = (REAL4)(*(wamp->int4_));
	  }
	}
      }
	  
      wamp = params->delta.next;
      if(wamp) {
	if(wamp->random == 1) {
	  REAL4 lo,hi;
	      
	  if(wamp->real4_) {
	    lo = *(wamp->real4_);
	  } else {
	    lo = (REAL4)(*(wamp->int4_));
	  }
	  
	  wamp = wamp->next;
	  
	  if(wamp->real4_) {
	    hi = *(wamp->real4_);
	  } else {
	    hi = (REAL4)(*(wamp->int4_));
	  }
	  
	  LAL_CALL(LALUniformDeviate(&status, &arandomnumber, rpar), &status);
	  Delta = lo + (hi-lo)*arandomnumber;
	      
	} else {
	  if(params->MAXinj) {
	    did = Nid;
	  }
	  for(i=0; i<did; i++) {
	    wamp = wamp->next;
	  }
	  if(wamp->real4vector_) {
	    Deltav = wamp->real4vector_;
	  } else if(wamp->real4_) {
	    Delta = *(wamp->real4_);
	  } else {
	    Delta = (REAL4)(*(wamp->int4_));
	  }
	}
      }

      wamp = params->psi.next;
      if(wamp) {
	if(wamp->random == 1) {
	  REAL4 lo,hi;
	  
	  if(wamp->real4_) {
	    lo = *(wamp->real4_);
	  } else {
	    lo = (REAL4)(*(wamp->int4_));
	  }
	  
	  wamp = wamp->next;
	    
	  if(wamp->real4_) {
	    hi = *(wamp->real4_);
	  } else {
	    hi = (REAL4)(*(wamp->int4_));
	  }

	  LAL_CALL(LALUniformDeviate(&status, &arandomnumber, rpar), &status);
	  Psi = lo + (hi-lo)*arandomnumber;
	      
	} else {
	  if(params->MAXinj) {
	    pid = Nid;
	  }
	  for(i=0; i<pid; i++) {
	    wamp = wamp->next;
	  }
	  if(wamp->real4vector_) {
	    Psiv = wamp->real4vector_;
	  } else if(wamp->real4_) {
	    Psi = *(wamp->real4_);
	  } else {
	    Psi = (REAL4)(*(wamp->int4_));
	  }
	}
      }
    }

#ifdef DEBUGBURST
    fprintf(stderr,"updating comment...\n");
    fflush(NULL);
#endif

    if(wid) {

      sprintf(comment,"%s,",1+wave->char_);

#ifdef DEBUGBURST
      fprintf(stderr,"\t\twAmpv\n");
      fflush(NULL);
#endif

      if(wAmpv) {
	UINT4 kj;
	sprintf(comment,"%s(%g",comment,wAmpv->data[0]);
	for(kj=1;kj<wAmpv->length;kj++) {
	  sprintf(comment,"%s,%g",comment,wAmpv->data[kj]);
	}
	sprintf(comment,"%s),",comment);
      } else {
	sprintf(comment,"%s%g,",comment,wAmp);
      }
	
#ifdef DEBUGBURST
      fprintf(stderr,"\t\tAlphav\n");
      fflush(NULL);
#endif
  
      if(Alphav) {
	UINT4 kj;
	sprintf(comment,"%s(%g",comment,Alphav->data[0]);
	for(kj=1;kj<Alphav->length;kj++) {
	  sprintf(comment,"%s,%g",comment,Alphav->data[kj]);
	}
	sprintf(comment,"%s),",comment);
      } else {
	sprintf(comment,"%s%g,",comment,Alpha);
      }

#ifdef DEBUGBURST
      fprintf(stderr,"\t\tDeltav\n");
      fflush(NULL);
#endif

      if(Deltav) {
	UINT4 kj;
	sprintf(comment,"%s(%g",comment,Deltav->data[0]);
	for(kj=1;kj<Deltav->length;kj++) {
	  sprintf(comment,"%s,%g",comment,Deltav->data[kj]);
	}
	sprintf(comment,"%s),",comment);
      } else {
	sprintf(comment,"%s%g,",comment,Delta);
      }

#ifdef DEBUGBURST
      fprintf(stderr,"\t\tPsiv\n");
      fflush(NULL);
#endif

      if(Psiv) {
	UINT4 kj;
	sprintf(comment,"%s(%g",comment,Psiv->data[0]);
	for(kj=1;kj<Psiv->length;kj++) {
	  sprintf(comment,"%s,%g",comment,Psiv->data[kj]);
	}
	sprintf(comment,"%s)",comment);
      } else {
	sprintf(comment,"%s%g",comment,Psi);
      }
      
    } else {
      sprintf(comment,",,,,,,");
    }

    /************************************************************/
    /* injection code */
    /************************************************************/
    if(wid) {
      LALDetAMResponse F; /* response functions */
      DetTimeAndASource dtS; /* for time delays */
      LALPlaceAndGPS pGPS; /* for time delays */
      LALDetAndSource dAs; /* for responses */
      LALSource source; /* for responses */
      SkyPosition position;  
      LIGOTimeGPS t0;
      REAL8 tDelay;
      REAL4TimeSeries sig;
      BurstParameter *injTime;
      INT4 ndt;
      BOOLEAN F1 = 0;
      REAL4 lo,hi;
      BOOLEAN injRandom;

      INT4 *NDT = NULL;
      INT4 indt, NDTN = 0;

#ifdef DEBUGBURST
      fprintf(stderr,"Beginning Injection...\n");
      fflush(NULL);
#endif

      lo = hi = 0.0;

      /* save a copy of input data */
      wbuf = (REAL4 *)LALMalloc(params->data.data->data->vectorLength * sizeof(REAL4));
      TASSERT(wbuf, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
      memcpy(wbuf, params->data.data->data->data, params->data.data->data->vectorLength * sizeof(REAL4));

      sprintf(comment,"%s,(",comment);

      sig.deltaT = params->data.data->deltaT;
      sig.f0 = params->data.data->f0;
      sig.sampleUnits = params->data.data->sampleUnits;

#ifdef DEBUGBURST
      fprintf(stderr,"**%i\t%i\n",wid,params->wave_pc[wid-1].data->vectorLength);
      fflush(NULL);
#endif 

      sig.data = NULL;
      LAL_CALL ( LALCreateVector(&status, &(sig.data), params->wave_pc[wid-1].data->vectorLength), &status);

      injTime = params->injTimes.next;
      
      wamp = params->injTimes.next;
      injRandom = 0;

      if(wamp->random == 1) {
	    
	injRandom = 1;
	
	if(wamp->real4_) {
	  lo = *(wamp->real4_);
	} else {
	  lo = (REAL4)(*(wamp->int4_));
	}

	wamp = wamp->next;

	if(wamp->real4_) {
	  hi = *(wamp->real4_);
	} else {
	  hi = (REAL4)(*(wamp->int4_));
	}
      } else if(params->MAXinj ==0) {
	for(i=0;i<Nid;i++) {
	  injTime = injTime->next; 
	}
      }


      if(params->channel == strstr(params->channel,"H1")) {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexLHODIFF])); 
      } else if(params->channel == strstr(params->channel,"H2")) {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexLHODIFF])); 
      } else if(params->channel == strstr(params->channel,"L1")) {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexLLODIFF]));
      } else if(params->channel[0] == 'G') {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexGEO600DIFF]));
      } else if(params->channel[0] == 'T') {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexTAMA300DIFF]));
      } else if(params->channel[0] == 'V') {
	pGPS.p_detector = dAs.pDetector = (LALDetector *)(&(lalCachedDetectors[LALDetectorIndexVIRGODIFF]));
      } else {
	SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
      }

#ifdef DEBUGBURST
      fprintf(stderr,"APPLYSEARCH: looping over %i\n",params->Ninj_per_segment);
      fflush(NULL);
#endif

      for(i=0; i<params->Ninj_per_segment; i++) { 

	if(wAmpv) {
	  if(wAmpv->length != params->Ninj_per_segment) {
	    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	  }
	  wAmp = wAmpv->data[i];
	}

	if(Alphav) {
	  if(Alphav->length != params->Ninj_per_segment) {
	    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	  }
	  Alpha = Alphav->data[i];
	}

	if(Deltav) {
	  if(Deltav->length != params->Ninj_per_segment) {
	    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	  }
	  Delta = Deltav->data[i];
	}

	if(Psiv) {
	  if(Psiv->length != params->Ninj_per_segment) {
	    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	  }
	  Psi = Psiv->data[i];
	}

	position.longitude = Alpha;
	position.latitude = Delta;
	position.system = COORDINATESYSTEM_EQUATORIAL;
	dtS.p_source = &position;

	dAs.pSource = &source;
	source.equatorialCoords = position;
	source.orientation = Psi;
	
	/* inject waveform: */

#ifdef DEBUGBURST
	fprintf(stderr,"APPLYSEARCH: before LALUniformDeviate\n");
	fflush(NULL);
#endif

	if(injRandom) {
	  LAL_CALL(LALUniformDeviate(&status, &arandomnumber, rpar), &status);
	  ndt = (INT4)floor(lo + (hi-lo)*arandomnumber);

#ifdef DEBUGBURST
	  fprintf(stderr,"%i\t%g\t%i\n",i,arandomnumber,ndt);
	  fflush(NULL);
#endif
	  
	} else {
	  if(injTime->real4vector_) {
	    if(injTime->real4vector_->length != params->Ninj_per_segment) {
	      SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	    }
	    ndt = (INT4)(injTime->real4vector_->data[i]);
	  } else if(injTime->real4_) {
	    ndt = (INT4)(*(injTime->real4_));
	    injTime = injTime->next;	      
	  } else {
	    ndt = *(injTime->int4_);
	    injTime = injTime->next;
	  }
	}

	/* check no overlap of signals */
	if(params->MAXinj == 0) {
	  for(indt=0;indt<NDTN;indt++) {
	    if(abs(NDT[indt] - ndt) <= sig.data->length) {
	      break;
	    }
	  }
	} else {
	  indt = NDTN;
	}
	    
	if(indt == NDTN) {
	      
	  NDTN++;
	  NDT = (INT4 *)LALRealloc(NDT, NDTN * sizeof(INT4));
	  NDT[indt] = ndt;

	  if(F1) {
	    sprintf(comment,"%s,%i",comment,ndt);
	  } else {
	    sprintf(comment,"%s%i",comment,ndt);
	    F1 = 1;
	  }
	    
	  {
	    REAL8 dt = (REAL8)ndt * sig.deltaT;

	    t0 = params->data.data->epoch;
	    
	    t0.gpsSeconds += (INT4)floor(dt);
	    t0.gpsNanoSeconds += (INT4)floor(1E9 * (dt - floor(dt)));

	    if(t0.gpsNanoSeconds >= 1000000000) {
	      (t0.gpsSeconds)++;
	      (t0.gpsNanoSeconds)-=1000000000;
	    }
	    
	  }

	  dtS.p_det_and_time = &pGPS;
	  pGPS.p_gps = &t0;

#ifdef DEBUGBURST
	  fprintf(stderr,"APPLYSEARCH: before LALTimeDelayFromEarthCenter\n");
	  fflush(NULL);
#endif

	  if(dtS.p_source->longitude < 0.0 ||
	     dtS.p_source->longitude > 6.2832 ||
	     dtS.p_source->latitude < -1.5708 ||
	     dtS.p_source->latitude > 1.5708) {
	    tDelay = 0.0;
	  } else {
	    LAL_CALL(LALTimeDelayFromEarthCenter( &status,&tDelay, &dtS ), &status );
	  }

	  {
	    LALGPSandAcc wi;
	    wi.gps = t0;
	    wi.accuracy = LALLEAPSEC_STRICT;
		
#ifdef DEBUGBURST
	    fprintf(stderr,"APPLYSEARCH: before LALComputeDetAMResponse\n");
	    fflush(NULL);
#endif

	    if(dtS.p_source->longitude < 0.0 ||
	       dtS.p_source->longitude > 6.2832 ||
	       dtS.p_source->latitude < -1.5708 ||
	       dtS.p_source->latitude > 1.5708) {
	      F.plus = F.cross = 1.0;
	    } else {
	      LAL_CALL(LALComputeDetAMResponse(&status,&F,&dAs,&wi), &status);
	    }
	  }

	  {
	    REAL8 dt = tDelay;

	    if(tDelay > 0.0) {
	      t0.gpsSeconds += (INT4)floor(dt);
	      t0.gpsNanoSeconds += (INT4)floor(1E9 * (dt - floor(dt)));
	    } else {
	      dt = fabs(dt);
	      t0.gpsSeconds -= (INT4)floor(dt);
	      t0.gpsNanoSeconds -= (INT4)floor(1E9 * (dt - floor(dt)));
	    }
	    
	    if(t0.gpsNanoSeconds >= 1000000000) {
	      (t0.gpsSeconds)++;
	      (t0.gpsNanoSeconds)-=1000000000;
	    }
	    
	    if(t0.gpsNanoSeconds < 0) {
	      (t0.gpsSeconds)--;
	      (t0.gpsNanoSeconds) = 1000000000 + t0.gpsNanoSeconds;
	    }
		
	  }	    
	      
	  sig.epoch = t0;
	  
	  bzero(sig.data->data, sig.data->length * sizeof(REAL4));
	    
	  for(j=0;j<sig.data->length;j++) {
	    sig.data->data[j] += wAmp * (F.plus * params->wave_pc[wid-1].data->data[j] + F.cross * params->wave_pc[wid-1].data->data[j + sig.data->length]);
	  }
	    

#ifdef DEBUGBURST 
	  {
	    static int dInjFID=0;
	    FILE *InjFID;
	    char InjFName[256];

	    REAL4TimeSeries ats;
	    REAL4Vector av;
	    ats.deltaT = params->data.data->deltaT;
	    ats.f0 = params->data.data->f0;
	    ats.epoch = params->data.data->epoch;
	    ats.data = &av;
	    av.length = params->data.data->data->vectorLength;
	    av.data = (REAL4 *)malloc(av.length * sizeof(REAL4));
	    bzero(av.data, av.length * sizeof(REAL4));

	    LAL_CALL(LALSSInjectTimeSeries(&status, &ats, &sig), &status);

	    sprintf(InjFName,"/dso-test/InjFile.%i",dInjFID);
	    
	    fprintf(stderr,"Dumping injection data to file %s...\n",InjFName);
	    fflush(NULL);

	    if(!(InjFID = fopen(InjFName,"w"))) {
	      fprintf(stderr,"can't write to file %s\n",InjFName);
	      fflush(NULL);
	    } else {
	      REAL4 t;
	      INT4 icnt;

	      fprintf(InjFID,"# Injection: \n");
	      fprintf(InjFID,"# waveform: %s\n",1+wave->char_);
	      fprintf(InjFID,"# longitude = %g\n",dtS.p_source->longitude);
	      fprintf(InjFID,"# latitude = %g\n",dtS.p_source->latitude);
	      fprintf(InjFID,"# orientation = %g\n", source.orientation);
	      fprintf(InjFID,"# Fplus = %g\n",F.plus);
	      fprintf(InjFID,"# Fcross = %g\n",F.cross);
	      fprintf(InjFID,"# t0.sec = %i\n",t0.gpsSeconds);
	      fprintf(InjFID,"# t0.nanosec = %i\n",t0.gpsNanoSeconds);

	      fprintf(InjFID,"# Data: \n");
	      fprintf(InjFID,"# t0.sec = %i\n",ats.epoch.gpsSeconds);
	      fprintf(InjFID,"# t0.nanosec = %i\n",ats.epoch.gpsNanoSeconds);

	      fprintf(InjFID,"# i\ttime\tADC\n");
	      for(t=0.0, icnt=0; icnt<av.length; icnt++, t+=ats.deltaT) {
		fprintf(InjFID,"%i\t%.10g\t%.10g\n",icnt,t,av.data[icnt]);
	      }

	      fclose(InjFID);

	      dInjFID++;
	    }

	    free(av.data);
	  }
#endif

#ifdef DEBUGBURST
	  fprintf(stderr,"APPLYSEARCH: before LALSSInjectTimeSeries\n");
	  fflush(NULL);
#endif

	  {
	    REAL4TimeSeries ats;
	    REAL4Vector av;
	    ats.deltaT = params->data.data->deltaT;
	    ats.f0 = params->data.data->f0;
	    ats.epoch = params->data.data->epoch;
	    ats.data = &av;
	    av.length = params->data.data->data->vectorLength;
	    av.data = params->data.data->data->data;

	    LAL_CALL(LALSSInjectTimeSeries(&status, &ats, &sig),&status);
	  }

#ifdef DEBUGBURST
	  fprintf(stderr,"APPLYSEARCH: after LALSSInjectTimeSeries\n");
	  fflush(NULL);
#endif

	}
      }

      sprintf(comment,"%s)",comment);

      LALFree(NDT);

      LAL_CALL(LALDestroyVector(&status, &(sig.data)), &status);
    }



    /************************************************************/
    /* pack burst parameters */
    /************************************************************/

    bparamsptr = &bparams;
    ETGparams = params->ETGparams.next;
    i=0;
    while(ETGparams) {

      bptr = ETGparams->param.next;
      for(j=0; j<Pid[i]; j++) {
	bptr = bptr->next;
      }

      if(params->ETGparamsLOCK == 0) {

	bparamsptr->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	TASSERT(bparamsptr->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	bparamsptr = bparamsptr->next;

	if(bptr->random == 1) {
	  SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	}

	if(bptr->char_ && !(bptr->real4vector_)) {
	  sprintf(comment,"%s,%s",comment,bptr->char_);

	  bparamsptr->char_ = (CHAR *)LALCalloc(1+strlen(bptr->char_),sizeof(CHAR));
	  TASSERT(bparamsptr->char_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	  strcpy(bparamsptr->char_, bptr->char_);

	} else {
	  if(bptr->real4_) {
	    sprintf(comment,"%s,%g",comment,*(bptr->real4_));
	    
	    bparamsptr->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4));
	    TASSERT(bparamsptr->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	    *(bparamsptr->real4_) = *(bptr->real4_);

	  } else {
	    if(bptr->int4_) {
	      sprintf(comment,"%s,%i",comment,*(bptr->int4_));

	      bparamsptr->int4_ = (INT4 *)LALMalloc(sizeof(INT4));
	      TASSERT(bparamsptr->int4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	      
	      *(bparamsptr->int4_) = *(bptr->int4_);
		  
	    } else {
	      if(bptr->real4vector_) {
		sprintf(comment,"%s,%s",comment,bptr->char_+1);

		bparamsptr->real4vector_ = NULL;
		LAL_CALL(LALCreateVector(&status, &(bparamsptr->real4vector_), bptr->real4vector_->length), &status);

		memcpy(bparamsptr->real4vector_->data, bptr->real4vector_->data, bptr->real4vector_->length*sizeof(REAL4));

	      } else {
		sprintf(comment,"%s,$",comment);
	      }
	    }
	  }
	}
      } else {
	/* have ETGparamsLOCK ! */
	BurstParameter *bptr1;

#ifdef DEBUGBURST
	fprintf(stderr,"SLAVE ETGParams: entering\n");
	fflush(NULL);
#endif 

	if(!(bptr->char_)) {
	  SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	}

	memcpy(&bptr1, bptr->char_, sizeof(BurstParameter *));

	bptr1 = bptr1->next;
	
	while(bptr1) {

	  bparamsptr->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	  TASSERT(bparamsptr->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	  bparamsptr = bparamsptr->next;

	  if(bptr1->random == 1) {
	    SABORT(INITSEARCHH_EIN, INITSEARCHH_MSGEIN);
	  }

	  if(bptr1->char_ && !(bptr1->real4vector_)) {

#ifdef DEBUGBURST
	    fprintf(stderr,"SLAVE ETGParams: str: %s\n",bptr1->char_);
	    fflush(NULL);
#endif 

	    sprintf(comment,"%s,%s",comment,bptr1->char_);
	    
	    bparamsptr->char_ = (CHAR *)LALCalloc(1+strlen(bptr1->char_),sizeof(CHAR));
	    TASSERT(bparamsptr->char_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	    strcpy(bparamsptr->char_, bptr1->char_);

	  } else {
	    if(bptr1->real4_) {

#ifdef DEBUGBURST
	      fprintf(stderr,"SLAVE ETGParams: real4: %g\n",*(bptr1->real4_));
	      fflush(NULL);
#endif 

	      sprintf(comment,"%s,%g",comment,*(bptr1->real4_));
		
	      bparamsptr->real4_ = (REAL4 *)LALMalloc(sizeof(REAL4));
	      TASSERT(bparamsptr->real4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	      *(bparamsptr->real4_) = *(bptr1->real4_);

	    } else {
	      if(bptr1->int4_) {

#ifdef DEBUGBURST
		fprintf(stderr,"SLAVE ETGParams: int4: %i\n",*(bptr1->int4_));
		fflush(NULL);
#endif 

		sprintf(comment,"%s,%i",comment,*(bptr1->int4_));

		bparamsptr->int4_ = (INT4 *)LALMalloc(sizeof(INT4));
		TASSERT(bparamsptr->int4_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

		*(bparamsptr->int4_) = *(bptr1->int4_);
		  
	      } else {
		if(bptr1->real4vector_) {
		  sprintf(comment,"%s,%s",comment,bptr1->char_+1);

		  bparamsptr->real4vector_ = NULL;
		  LAL_CALL(LALCreateVector(&status, &(bparamsptr->real4vector_), bptr1->real4vector_->length), &status);

		  memcpy(bparamsptr->real4vector_->data, bptr1->real4vector_->data, bptr1->real4vector_->length*sizeof(REAL4));

		} else {
		  sprintf(comment,"%s,$",comment);
		}
	      }
	    }
	  } 
	  
	  bptr1 = bptr1->next;
	}
      }
      
      ETGparams = ETGparams->next;
      i++;
    }

    bparamsptr->next = NULL;

#ifdef DEBUGBURST
    fprintf(stderr,"%s\n",comment);
    fflush(NULL);
#endif

    /************************************************************/
    /* run ETG */
    /************************************************************/

    {

      void (*oldSegFault)(int), (*oldSigBus)(int);

      oldSegFault = signal(SIGSEGV, BurstSigHandler);
      oldSigBus = signal(SIGBUS, BurstSigHandler);

#ifdef DEBUGBURST
      fprintf(stderr,"Launching ETG function...\n");
      fflush(NULL);

      {
	clock_t tic = clock(),
	  toc;
	REAL4 *tdata = (REAL4 *)LALMalloc(params->data.data->data->vectorLength * sizeof(REAL4));
	
	memcpy(tdata, params->data.data->data->data, params->data.data->data->vectorLength * sizeof(REAL4));

#endif

	LAL_CALL(params->ETGfun(&status,
				&boutput,
				params->data.data,
				&bparams), &status);

#ifdef DEBUGBURST
	toc = clock();

	fprintf(stderr,"Time in ETG function: %g s\n",(REAL8)(toc-tic)/(REAL8)CLOCKS_PER_SEC);

	/*
	  {
	  UINT4 nbursts = 0;
	  EventIDColumn *optr = boutput.next;
	  
	  fprintf(stderr,"*********************************************\n");
	  while(optr) {
	  nbursts++;
	  fprintf(stderr,"%s %s\t%s\t%u\t%u\t%g\t%g\t%g\n",optr->snglBurstTable->ifo, optr->snglBurstTable->search, optr->snglBurstTable->channel, optr->snglBurstTable->start_time.gpsSeconds, optr->snglBurstTable->start_time.gpsNanoSeconds, optr->snglBurstTable->duration, optr->snglBurstTable->central_freq, optr->snglBurstTable->bandwidth);
	  optr = optr->next;
	  }
	  
	  fprintf(stderr,"%u bursts\n",nbursts);
	  fprintf(stderr,"*********************************************\n");
	  }
	*/

	fflush(NULL);
#endif

	signal(SIGSEGV, oldSegFault);
	signal(SIGBUS, oldSigBus);

#ifdef DEBUGBURST

	if(memcmp(tdata, params->data.data->data->data, params->data.data->data->vectorLength * sizeof(REAL4)) != 0) {
	  SABORT(INITSEARCHH_EETGDATA, INITSEARCHH_MSGEETGDATA);
	}

	LALFree(tdata);
      }

      fprintf(stderr,"ETG function is done...\n");
      fflush(NULL);
#endif

    }

    /************************************************************/
    /* run standardized output function */
    /************************************************************/

    if(!(params->noLALBurstOutput)) {
#ifdef DEBUGBURST
      {
	clock_t tic = clock(),
	  toc;
#endif
	
	LAL_CALL(LALBurstOutput(&status, &stdOutput, &boutput, &(params->oparams)), &status);

#ifdef DEBUGBURST
	toc = clock();
	fprintf(stderr,"Time in LALBurstOutput function: %g s\n",(REAL8)(toc-tic)/(REAL8)CLOCKS_PER_SEC);
      }
#endif
    } else {
	  
      stdOutput.next = boutput.next;
      
    }

    /* remove injected waveform */
    if(wid) {
      memcpy(params->data.data->data->data, wbuf, params->data.data->data->vectorLength * sizeof(REAL4));
      LALFree(wbuf);

      if(!(params->noLALBurstOutput)) {
	LAL_CALL(LALBurstOutput(&status, NULL, NULL, NULL), &status);
      }
    }

    /************************************************************/
    /* prepare output */
    /************************************************************/
    switch(params->outputMethod) {

      /*
	case 0:
	  BuildOutput(status->statusPtr, output, &stdOutput, comment, params->data.data->epoch.gpsSeconds, params->data.data->epoch.gpsNanoSeconds, params->data.data->deltaT * params->data.data->data->vectorLength, numSlaves, zerobid);
	  CHECKSTATUSPTR (status);
	  break;
      */
      /*
	case 1:
	  BuildOutputBinary(status->statusPtr, output, &stdOutput, comment, params->data.data->epoch.gpsSeconds, params->data.data->epoch.gpsNanoSeconds, zerobid);
	  CHECKSTATUSPTR (status);
	  break;
      */
	case 2:
	  LAL_CALL(
		   BuildOutputFileBinary(&status, fid, &stdOutput, comment, params->data.data->epoch.gpsSeconds, params->data.data->epoch.gpsNanoSeconds),
		   &status);
	  break;

	default:
	  SABORT(INITSEARCHH_EIN,INITSEARCHH_MSGEIN);
    }

    /************************************************************/
    /* clean up */
    /************************************************************/
    LALFree(Pid);

    {
      EventIDColumn *optr;

      optr = boutput.next;
      while(optr) {
	EventIDColumn *optr2 = optr->next;

	if(optr->snglTransdataTable) {
	  if(optr->snglTransdataTable->trans_data) {
	    LALFree(optr->snglTransdataTable->trans_data);
	  }
	  LALFree(optr->snglTransdataTable);
	}

	if(optr->snglBurstTable) {
	  LALFree(optr->snglBurstTable);
	}

	LALFree(optr);
	optr = optr2;
      }

      if(!(params->noLALBurstOutput)) {
	optr = stdOutput.next;
	while(optr) {
	  EventIDColumn *optr2 = optr->next;

	  /* no need to free trans_data since pointer copy of boutput */
	  
	  if(optr->snglBurstTable) {
	    LALFree(optr->snglBurstTable);
	  }

	  LALFree(optr);
	  optr = optr2;
	}
      }
    }

    bparamsptr = bparams.next;
    while(bparamsptr) {
      BurstParameter *bptr = bparamsptr->next;

      if(bparamsptr->char_) { 
	LALFree(bparamsptr->char_); 
      } 
      if(bparamsptr->int4_) { 
	LALFree(bparamsptr->int4_); 
      } 
      if(bparamsptr->real4_) { 
	LALFree(bparamsptr->real4_); 
      } 
      if(bparamsptr->real4vector_) { 
	LAL_CALL(LALDestroyVector(&status, &(bparamsptr->real4vector_)), &status); 
      }

      LALFree(bparamsptr);
      bparamsptr = bptr;
    }

  }

  if(!got_wid && !(params->noLALBurstOutput)) {
    LAL_CALL(LALBurstOutput(&status, NULL, NULL, NULL), &status);
  }

  LAL_CALL(LALDestroyRandomParams(&status, &rpar), &status);


  return 0;
}



/***********************************************************************/
/***********************************************************************/


void BuildOutputFileBinary(LALStatus *status,
			   int fid,
			   EventIDColumn     *burstEventList,
			   CHAR               *tag,
			   UINT4 t0Sec,
			   UINT4 t0NanoSec
			   )
{

  /*
    Writes data to disk

    Arguments:
    status: LAL status pointer
    fid: file identifier
    burstEventList: data from ETG
    tag: string describing injection and ETG parameters used (TOC)
    t0Sec: GPS start of analyzed data
    t0NanoSec: nanosec portion of GPS start
  */

#ifndef WORDS_BIGENDIAN
  static void endian_swap(char * pdata, int dsize, int nelements);
#endif

  /*
  FILE *fid;
  */

  /*
  struct flock lock;
  */

  SnglBurstTable *sbtptr;
  EventIDColumn *sbtptrG;

  UCHAR *data;
  UINT4 i,numOutput,datalength, size1, datalength_;

  size_t off;

  INITSTATUS (status, "BuildOutput", RUNSEARCHC);
  ATTATCHSTATUSPTR (status);

  sbtptrG = burstEventList->next;
  numOutput = 0;
  while(sbtptrG) {
    numOutput++;
    sbtptrG = sbtptrG->next;
  }

#ifdef DEBUGBURST
  fprintf(stderr,"Got %u bursts\n",numOutput);
  fflush(NULL);
#endif

  size1 = (LIGOMETA_IFO_MAX + LIGOMETA_SEARCH_MAX + LIGOMETA_CHANNEL_MAX) * sizeof(CHAR) + 2*sizeof(INT4) + 6 * sizeof(REAL4);
  datalength_ = datalength = numOutput * size1 + sizeof(UINT4) + strlen(tag) + 1 + 3*sizeof(UINT4);
  data = (UCHAR *)LALMalloc( datalength );
  if(!data) {
    ABORT(status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
  }

  memcpy(data,&t0Sec,sizeof(UINT4));
  memcpy(data+sizeof(UINT4),&t0NanoSec,sizeof(UINT4));
  memcpy(data+2*sizeof(UINT4),&datalength_,sizeof(UINT4));

#ifndef WORDS_BIGENDIAN
  endian_swap((char *)(data), sizeof(UINT4), 3);
#endif

  strcpy(data+3*sizeof(UINT4),tag);

#ifndef WORDS_BIGENDIAN
  endian_swap((char *)(data+3*sizeof(UINT4)), sizeof(UCHAR), 1+strlen(tag));
#endif

  memcpy(data + strlen(tag) + 1 + 3*sizeof(UINT4), &numOutput, sizeof(UINT4));
#ifndef WORDS_BIGENDIAN
  endian_swap((char *)(data + strlen(tag) + 1 + 3*sizeof(UINT4)), sizeof(UINT4), 1);
#endif

  off = sizeof(UINT4) + strlen(tag) + 1 + 3*sizeof(UINT4);

  i = 0;
  sbtptrG = burstEventList->next;

  while(sbtptrG) {
    UCHAR *p1 = data + i * size1 + off,
      *pst, *pr4;

    if(sbtptrG->snglTransdataTable) {
      ABORT(status,INITSEARCHH_EUI,INITSEARCHH_MSGEUI);
    }

    sbtptr = sbtptrG->snglBurstTable;

    memcpy(p1, sbtptr->ifo,LIGOMETA_IFO_MAX * sizeof(CHAR));

    p1 += LIGOMETA_IFO_MAX * sizeof(CHAR);
    memcpy(p1, sbtptr->search,LIGOMETA_SEARCH_MAX * sizeof(CHAR));

    p1 += LIGOMETA_SEARCH_MAX * sizeof(CHAR);
    memcpy(p1, sbtptr->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR));

    pst = p1 + LIGOMETA_CHANNEL_MAX * sizeof(CHAR);
    memcpy(pst, &(sbtptr->start_time.gpsSeconds), sizeof(INT4));
    memcpy(pst + sizeof(INT4), &(sbtptr->start_time.gpsNanoSeconds), sizeof(INT4));

    pr4 = pst + 2*sizeof(INT4);
    memcpy(pr4, &(sbtptr->duration), sizeof(REAL4));

    p1 = pr4 + sizeof(REAL4);
    memcpy(p1, &(sbtptr->central_freq), sizeof(REAL4));

    p1 += sizeof(REAL4);
    memcpy(p1, &(sbtptr->bandwidth), sizeof(REAL4));

    p1 += sizeof(REAL4);
    memcpy(p1, &(sbtptr->amplitude), sizeof(REAL4));
    
    p1 += sizeof(REAL4);
    memcpy(p1, &(sbtptr->snr), sizeof(REAL4));

    p1 += sizeof(REAL4);
    memcpy(p1, &(sbtptr->confidence), sizeof(REAL4));

#ifndef WORDS_BIGENDIAN
    endian_swap((char *)(pst), sizeof(INT4), 2);
    endian_swap((char *)(pr4), sizeof(REAL4), 6);
#endif

    sbtptrG = sbtptrG->next;
    i++;
  }

  
  if(write(fid,data,datalength) < datalength) {
#ifdef DEBUGBURST
    perror("ERROR from fwrite:");
    fflush(NULL);
#endif
	    
    ABORT(status,INITSEARCHH_EFILE,INITSEARCHH_MSGEFILE);
  }

 LALFree(data);

 DETATCHSTATUSPTR (status);
 RETURN (status);
}


/********************************************************************/
/********************************************************************/


void BurstSigHandler(int s) {

  /* replace standard signal handler */

  switch(s) {
  case SIGSEGV:
    fprintf(stderr,"Got SIGSEGV from ETG function!\n");
    break;
  case SIGBUS:
    fprintf(stderr,"Got SIGBUS from ETG function!\n");
    break;
  default:
    fprintf(stderr,"Got unknown signal %i from ETG function\n",s);
    break;
  }

  fflush(NULL);

#ifdef DEBUGBURST
  exit(-1);
#endif

  ABORT(gStatus,INITSEARCHH_EETG,INITSEARCHH_MSGEETG);

}


/********************************************************************/
/********************************************************************/


#ifndef WORDS_BIGENDIAN
static void endian_swap(char * pdata, int dsize, int nelements)

{

  /* convert between little and big endian */

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


/******************************************************************/
/******************************************************************/


int symbolsfree(datacond_symbol_type *symbols) {

  /* unimplemented clean-up function */

  return 0;
}

int bparamsfree(BurstSearchParams *bparams) {

  /* unimplemented clean-up function */

  return 0;
}
