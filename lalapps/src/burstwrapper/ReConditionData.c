#include <burstdso.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/Random.h>
#include <lal/LALMalloc.h>

#include <lal/Calibration.h>
#include "datacondAPI/DatacondCaller.h"
#include <strings.h>

#include <lalapps/lalapps.h>

#include <Translation.h>

#define RESPONSE_REGEX "CAL-RESPONSE"
#define CAV_GAIN_REGEX "CAL-CAV_GAIN"
#define CAV_FAC_REGEX "CAL-CAV_FAC"
#define OLOOP_FAC_REGEX "CAL-OLOOP_FAC"

NRCSID( RECONDITIONDATAC , "ReConditionData.c" );
RCSID( "ReConditionData.c" );

/* units for calibration */
const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

int ReConditionData(int Nsymbols,
		    datacond_symbol_type *symbols,
		    char *algorithms,
		    BurstSearchParams *params)
{

  /* 
     This function unpack the data, prepares calibration.

     Arguments:
     Nsumbols: number of datacondAPI symbols
     symbols: list of datacondAPI symbols (i.e., data inputs)
     algorithms: not used
     params: burst search parameters
  */

  static LALStatus             status;

  BurstParameter *ptr;
  BurstParameterList *ptrl;
  UINT4 si, Nmat;
  CHAR *cptr;
  CHAR ccptr[8192];

  { 
    INT4 i = 0;
    INT4 ki; 

    REAL4Vector *ptr_p, *ptr_c;

    /************************************************************/
    /* Handle __FILE__ case */
    /* NOTE: this is EXPERIMENTAL code for the stand-alone system */	
    /************************************************************/

    /* look at 1st ETG parameter */
    ptr = params->ETGparams.next->param.next;
    if(ptr) {
      if(ptr->random == 2) {

	/* flag random==2 means that input is a matrix stored
	   in a ilwd file */

	params->ETGparamsLOCK = 1; /* lock ETG parameters */

	/* remove underscores to get symbol name */
	cptr = strtok(ptr->char_,"_"); 

	/* loop over datacond symbols, try to find a match
	   to the symbol */
	for(si=0;si<Nsymbols;si++) {
	  if(strstr(symbols[si].s_symbol_name, cptr)) { 

	    /* we have a matching name */

	    CHARVectorSequence *cvs;

	    CHAR *cbuf;

	    params->NNparams[0] = 0;

	    LALFree(ptr->char_);
	    LALFree(ptr);

	    ptr = &(params->ETGparams.next->param);
	    bzero(ptr, sizeof(BurstParameter));

	    /* process ilwd file, create burstparams lists */
	    if(symbols[si].s_translator != TranslateCHARVectorSequence ||
	       symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT) { 
	      SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN ); 
	    } 

	    cvs = (CHARVectorSequence *)(symbols[si].s_user_data);

	    if(cvs->length != 2) { 
	      SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN ); 
	    } 

	    cbuf = (CHAR *)calloc(1+cvs->vectorLength, sizeof(CHAR));
	    TASSERT(cbuf, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	    for(ki=0;ki<cvs->length;ki++) { 
	      BurstParameter *tmp;

#ifdef DEBUGBURST
	      fprintf(stderr,"ETGParams: entering loop\n");
	      fflush(NULL);
#endif 

	      (params->NNparams[0])++;

	      tmp = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	      TASSERT(tmp, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	      ptr->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter));
	      TASSERT(ptr->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
	      ptr = ptr->next;
	      ptr->char_ = (CHAR *)LALCalloc(1, sizeof(BurstParameter *));
	      TASSERT(ptr->char_, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	      memcpy(ptr->char_, &tmp, sizeof(BurstParameter *));

	      bzero(cbuf, cvs->vectorLength * sizeof(CHAR));
	      memcpy(cbuf, cvs->data + cvs->vectorLength * ki, cvs->vectorLength * sizeof(CHAR));


#ifdef DEBUGBURST
	      fprintf(stderr,"ETGParams: processing string\n%s\n",cbuf);
	      fflush(NULL);
#endif 

	      /* process string */
	      {
		CHAR *arg, *s3;
		UINT4 Npa;

		arg = strtok_r(cbuf,",",&s3);

		while(arg) {
#ifdef DEBUGBURST
		  fprintf(stderr,"ETGParams: start ParseParameter for %s\n",arg);
		  fflush(NULL);
#endif 
		  Npa = 0;

		  LAL_CALL ( ParseParameter(&status, tmp, arg, &Npa), &status);


		  tmp = tmp->next;

#ifdef DEBUGBURST
		  fprintf(stderr,"ETGParams: done ParseParameter (%u)\n",Npa);
		  fprintf(stderr,"ETGParams: %s\n%s\n",cbuf,arg);
		  fflush(NULL);
#endif 

		  arg = strtok_r(NULL,",",&s3);

#ifdef DEBUGBURST
		  fprintf(stderr,"ETGParams: %s\n%s\n",cbuf,arg);
		  fflush(NULL);
#endif 
		}
	      }
	    }

	    free(cbuf);

	    break;

	  }

	}

	if(si == Nsymbols) { 
	  SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL); 
	} 

      }
    }

#ifdef DEBUGBURST
    fprintf(stderr,"ETGParams: done parsing\n");
    fflush(NULL);
#endif 

#define DOMATRIX(par_) \
    Nmat = 0; \
    ptr = params->par_.next; \
    if(ptr) { \
      REAL4VectorSequence *rvs; \
      if(ptr->random == 2) { \
	cptr = strtok(ptr->char_,"_"); \
\
	for(si=0;si<Nsymbols;si++) { \
	  if(strstr(symbols[si].s_symbol_name, cptr)) { \
	    \
	    size_t off; \
\
	    LALFree(ptr->char_); \
	    LALFree(ptr); \
	    ptr = &(params->par_); \
	    \
	    if(symbols[si].s_translator != TranslateREAL4Sequence || symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT ) { \
	      SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN ); \
	    } \
            rvs = (REAL4VectorSequence *)(symbols[si].s_user_data); \
	    if(rvs->length != 2) { \
	      SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN ); \
	    } \
\
	    for(ki=0, off=0;ki<rvs->length;ki++) { \
	      \
	      ptr->next = (BurstParameter *)LALCalloc(1,sizeof(BurstParameter)); \
	      TASSERT(ptr->next, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC); \
	      ptr = ptr->next; \
              Nmat++; \
	      LAL_CALL ( LALCreateVector(&status, &(ptr->real4vector_), rvs->vectorLength), &status); \
\
	      memcpy(ptr->real4vector_->data, rvs->data + off, rvs->vectorLength * sizeof(REAL4)); \
\
	      off += rvs->vectorLength; \
              if(params->MAXinj == 0) {params->MAXinj = rvs->vectorLength;} \
\
	  } \
         params->Ninj_per_segment = params->MAXinj; \
         break; \
	} \
}\
	if(si == Nsymbols) { \
	  SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL); \
	} \
      \
      } \
    } \

    DOMATRIX(waveform_amplitudes)
      /*
    if(params->MAXinj) {
      params->Nwaveform_amplitude = Nmat;
    }
      */
    DOMATRIX(alpha)
    DOMATRIX(delta)
    DOMATRIX(psi)
    DOMATRIX(injTimes)

#undef DOMATRIX

    if(params->MAXinj) {
      params->Nwaveform_amplitude = params->Nalpha = params->Ndelta = params->Npsi = 1;
    }



    /************************************************************/
    /* Handle injected waveforms */
    /************************************************************/
    ptr = params->waveforms.next;

#ifdef DEBUGBURST
for(si=0;si<Nsymbols;si++) {
      fprintf(stderr,"**%i:%s\n",si,symbols[si].s_symbol_name);
      fflush(NULL);
}
#endif 

    if(params->Nwaveforms) {
      params->wave_pc = (REAL4TimeVectorSeries *)LALCalloc(params->Nwaveforms, sizeof(REAL4TimeVectorSeries));
      TASSERT(params->wave_pc, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);
    }

    while(ptr) {

      cptr = strtok(ptr->char_,"_");

      strcpy(ccptr,cptr);
      strcat(ccptr,"_p");

      for(si=0;si<Nsymbols;si++) {
	if(strstr(symbols[si].s_symbol_name, ccptr)) {

	  REAL4Vector *rv;

#ifdef DEBUGBURST
	  fprintf(stderr,"Using waveforms %s for %s\n",symbols[si].s_symbol_name,ccptr);
	  fflush(NULL);
#endif
	  ptr_p = NULL;

	  if(symbols[si].s_translator != TranslateREAL4Sequence || symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT) {
	    SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	  }

	  rv = (REAL4Vector *)(symbols[si].s_user_data);
	  
	  LAL_CALL ( LALCreateVector(&status, &ptr_p, rv->length), &status );


	  
	  memcpy(ptr_p->data, rv->data,ptr_p->length * sizeof(REAL4)); 
	  
	  break;
	}
      }
      
      if(si == Nsymbols) {
	SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
      }
      

     

      cptr = strtok(ptr->char_,"_");
      strcpy(ccptr,cptr);
      strcat(ccptr,"_c");

      for(si=0;si<Nsymbols;si++) {
	if(strstr(symbols[si].s_symbol_name, ccptr)) {

	  REAL4Vector *rv;

#ifdef DEBUGBURST
	  fprintf(stderr,"Using waveforms %s for %s\n",symbols[si].s_symbol_name,ccptr);
	  fflush(NULL);
#endif

	  if(symbols[si].s_translator != TranslateREAL4Sequence || symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT) {
	    SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	  }

	  rv = (REAL4Vector *)(symbols[si].s_user_data);
	  
	  ptr_c = NULL;
	  LAL_CALL ( LALCreateVector(&status, &ptr_c, rv->length), &status );


	  memcpy(ptr_c->data, rv->data, ptr_c->length * sizeof(REAL4)); 
	  
	  break;
	}
      }
      
      if(si == Nsymbols) {
	SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
      }


      /* pack into TimeVectorSeries */
      TASSERT(ptr_c->length == ptr_p->length, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN );

      params->wave_pc[i].data = (REAL4VectorSequence *)LALMalloc(sizeof(REAL4VectorSequence));
      TASSERT(params->wave_pc[i].data, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN );

      params->wave_pc[i].data->length = 2;
      params->wave_pc[i].data->vectorLength = ptr_c->length;

#ifdef DEBUGBURST
      fprintf(stderr,"**%i\t%i\n",i,params->wave_pc[i].data->vectorLength);
      fflush(NULL);
#endif 

      params->wave_pc[i].data->data = (REAL4 *)LALMalloc(2 * ptr_c->length * sizeof(REAL4));
      TASSERT(params->wave_pc[i].data->data, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN );

      memcpy(params->wave_pc[i].data->data, ptr_p->data, ptr_p->length * sizeof(REAL4)); 

      memcpy(params->wave_pc[i].data->data + ptr_p->length, ptr_c->data, ptr_c->length * sizeof(REAL4)); 
      
      ptr = ptr->next;
      i++;

      LAL_CALL ( LALDestroyVector(&status, &ptr_p), &status );


      LAL_CALL ( LALDestroyVector(&status, &ptr_c), &status );


    }
  

    /************************************************************/
    /* handle _VECT_ from ETG parameters */
    /************************************************************/

    ptrl = params->ETGparams.next;
    while(ptrl) {
      ptr = ptrl->param.next;
      while(ptr) {
	
	if(ptr->real4vector_) {
	  
	  cptr = strtok(ptr->char_,"_");
	  
	  for(si=0;si<Nsymbols;si++) {
	    if(strstr(symbols[si].s_symbol_name, cptr)) {

	      REAL4Vector *rv;

	      if(symbols[si].s_translator != TranslateREAL4Sequence || symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT ) {
		SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	      }
	      
	      rv = (REAL4Vector *)(symbols[si].s_user_data);

	      ptr->real4vector_->length = rv->length;
	      ptr->real4vector_->data = (REAL4 *)LALMalloc(ptr->real4vector_->length * sizeof(REAL4));
	      TASSERT(ptr->real4vector_->data, status, INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	      
	      memcpy(ptr->real4vector_->data, rv->data, ptr->real4vector_->length * sizeof(REAL4)); 
	      
	      break;
	    }
	  }

	  if(si == Nsymbols) {
	    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
	  }

	}

	ptr = ptr->next;
      }
      ptrl = ptrl->next;
    }


    /************************************************************/
    /* get GW data */
    /************************************************************/

    /* loop over datacond symbols */
    for(si=0;si<Nsymbols;si++) {

      /* match name */
      if(strstr(symbols[si].s_symbol_name, GW_STRAIN_DATA)) {
	
	REAL4TimeSeries *r4ts;

	/* make sure it's an output symbol, and it's REAL4 */
	if(symbols[si].s_direction != DATACOND_SYMBOL_OUTPUT ||
	   symbols[si].s_translator != TranslateREAL4TimeSeries) {
	  SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	}

	/* pointer to data */
	r4ts = (REAL4TimeSeries *)(symbols[si].s_user_data);

	/* check length */
	if(r4ts->data->length != params->Ndata) {
	  SABORT (  INITSEARCHH_EIN, INITSEARCHH_MSGEIN );
	}

	/* allocate memory */
	params->data.data = (REAL4TimeVectorSeries *)LALCalloc(1, sizeof(REAL4TimeVectorSeries));
	TASSERT(params->data.data, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

	params->data.data->data = NULL;

	{
	  CreateVectorSequenceIn vsin;
	  
	  vsin.length = 1;
	  vsin.vectorLength = r4ts->data->length;

	  LAL_CALL ( LALCreateVectorSequence(&status, &(params->data.data->data), &vsin), &status );

	}

	/* copy data */
	memcpy(params->data.data->data->data, r4ts->data->data,params->data.data->data->vectorLength * sizeof(REAL4)); 

	/* set time */
	params->data.data->epoch = r4ts->epoch;

	/* set sampling frequency */
	params->data.data->deltaT = r4ts->deltaT;

      break;

      }

    }

    if(si == Nsymbols) {
      /* can't find the data! */
      SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
    }





    /************************************************************/
    /* set up calibration */
    /* code stolen from inspiral dso */
    /************************************************************/

    {
      REAL4TimeSeries              ats, *chanPtr = &ats;
      REAL4Vector av;
      COMPLEX8FrequencySeries      *respPtr;
      CHAR                          ifo[3];

      chanPtr->f0 =  params->data.data->f0;
      chanPtr->deltaT =  params->data.data->deltaT;
      chanPtr->epoch =  params->data.data->epoch;
      chanPtr->f0 =  params->data.data->f0;
      chanPtr->data = &av;
      chanPtr->data->length = params->data.data->data->vectorLength;
      chanPtr->data->data = params->data.data->data->data;

      respPtr = params->data.resp = (COMPLEX8FrequencySeries *)LALCalloc(1, sizeof(COMPLEX8FrequencySeries));
      TASSERT(respPtr, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

      strncpy(ifo,params->channel,2);
      ifo[2] = 0;

#ifdef DEBUGBURST
      fprintf(stderr,"ifo = %s\n",ifo);
      fflush(NULL);
#endif 

      /* set the start epoch for the response function from the data */
      respPtr->epoch.gpsSeconds = chanPtr->epoch.gpsSeconds;
      respPtr->epoch.gpsNanoSeconds = chanPtr->epoch.gpsNanoSeconds;

      /* set the frequency step size of the response from the psd */
      respPtr->f0 = 0.0;
      respPtr->deltaF = 1.0 / 
        ( ((REAL8) chanPtr->data->length) * chanPtr->deltaT );

      /* set the units of the desired response function */
      respPtr->sampleUnits = strainPerCount;

      /* allocate memory for the vector in the response function */
      respPtr->data = (COMPLEX8Vector *)LALCalloc(1, sizeof(COMPLEX8Vector));
      TASSERT(respPtr->data, status, INITSEARCHH_EALOC, INITSEARCHH_MSGEALOC);

      respPtr->data->length = chanPtr->data->length / 2 + 1;
      respPtr->data->data = (COMPLEX8 *) 
        LALCalloc( respPtr->data->length, sizeof(COMPLEX8) );

      /* extract it */
      {
	if ( ExtractResponse( respPtr, Nsymbols, symbols, ifo ) ) {
	  if(params->allowNoResponse) {
	    /* No response */
	    LAL_CALL ( LALCDestroyVector(&status, &(params->data.resp->data)), &status );

	    LALFree(params->data.resp);
	    params->data.resp = NULL;
	  } else {
	    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
	  }
	}
      }

/*
#ifdef DEBUGBURST
	{
		int ind;

      		fprintf(stderr,"Response function:\n");

		for(ind=0;ind<respPtr->data->length;ind++) {
			fprintf(stderr,"%i\t%g\t%g\n",ind,respPtr->data->data[ind].re,respPtr->data->data[ind].im);
		}

      		fflush(NULL);
	}
#endif 
*/

    }

  }

  return 0;

}


#include <sys/types.h>
#include <regex.h>

#define FindSequence(x) regcomp(&reg, x, REG_NOSUB); \
    for(i=0;i<Nsymbols;i++) if(!regexec(&reg, symbols[i].s_symbol_name, 0, NULL, 0)) break; \
    regfree(&reg); \
    if(i==Nsymbols) {fprintf(stderr,"Can't find %s\n",regex); \
     return 1;} else (void)(0)
    


int ExtractResponse( COMPLEX8FrequencySeries *output,
		     int Nsymbols, 
		     datacond_symbol_type *symbols, 
		     char *ifo ) {

  static LALStatus status;

  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
  CHAR regex[64];
  CHAR regexhead[64];
  COMPLEX8Vector abvec;
  COMPLEX8Vector avec;

  COMPLEX8FrequencySeries *R0;
  COMPLEX8FrequencySeries *C0;
  COMPLEX8TimeSeries *ab;
  COMPLEX8TimeSeries *a;

  CalibrationFunctions    calfuncs;
  CalibrationUpdateParams calfacts;

  int i;
  regex_t reg;

  strncpy( regexhead, ifo ? ifo : "", sizeof( regexhead ) - 1 );
  strncat( regexhead, ifo ? ".*:" : "",
      sizeof( regexhead ) - strlen( regexhead ) - 1 );

  strncpy( regex, regexhead, sizeof( regex ) - 1 );
  strncat( regex, RESPONSE_REGEX, sizeof( regex ) - strlen( regex ) - 1 );
  FindSequence(regex);
  if(symbols[i].s_direction != DATACOND_SYMBOL_OUTPUT || 
     symbols[i].s_translator != TranslateCOMPLEX8FrequencySeries) {
    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }

  R0 = (COMPLEX8FrequencySeries *)(symbols[i].s_user_data);


  strncpy( regex, regexhead, sizeof( regex ) - 1 );
  strncat( regex, CAV_GAIN_REGEX, sizeof( regex ) - strlen( regex ) - 1 );
  FindSequence(regex);
  if(symbols[i].s_direction != DATACOND_SYMBOL_OUTPUT || 
     symbols[i].s_translator != TranslateCOMPLEX8FrequencySeries) {
    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }

  C0 = (COMPLEX8FrequencySeries *)(symbols[i].s_user_data);


  strncpy( regex, regexhead, sizeof( regex ) - 1 );
  strncat( regex, OLOOP_FAC_REGEX, sizeof( regex ) - strlen( regex ) - 1 );
  FindSequence(regex);
  if(symbols[i].s_direction != DATACOND_SYMBOL_OUTPUT || 
     symbols[i].s_translator != TranslateCOMPLEX8TimeSeries) {
    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }
  ab = (COMPLEX8TimeSeries *)(symbols[i].s_user_data);

  strncpy( regex, regexhead, sizeof( regex ) - 1 );
  strncat( regex, CAV_FAC_REGEX, sizeof( regex ) - strlen( regex ) - 1 );
  FindSequence(regex);
  if(symbols[i].s_direction != DATACOND_SYMBOL_OUTPUT || 
     symbols[i].s_translator != TranslateCOMPLEX8TimeSeries) {
    SABORT( INITSEARCHH_ENULL, INITSEARCHH_MSGENULL);
  }
  a = (COMPLEX8TimeSeries *)(symbols[i].s_user_data);


  calfuncs.responseFunction = R0;
  calfuncs.sensingFunction  = C0;
  calfacts.openLoopFactor   = ab;
  calfacts.sensingFactor    = a;

  calfacts.epoch = output->epoch; /* specify correct epoch for update */
  calfacts.duration.gpsSeconds = 0;     /* keep behaviour the same as  */
  calfacts.duration.gpsNanoSeconds = 0; /* before steve's modification */
                                        /* to LALUpdateCalibration()   */

  /* should be able to update into the same functions... */
  calfuncs.responseFunction->sampleUnits = strainPerCount;
  LAL_CALL ( LALUpdateCalibration( &status, &calfuncs, &calfuncs, &calfacts ), &status );

  if(status.statusPtr)
  {
    if ( status.statusPtr->statusCode == CALIBRATIONH_EZERO )
    {
      LAL_CALL( LALResponseConvert( status.statusPtr, output,
            calfuncs.responseFunction ), status.statusPtr );
    }
  }


  /* now convert response to get output, hardwire units */
  LAL_CALL ( LALResponseConvert( &status, output, calfuncs.responseFunction ), &status );

  return 0;

}
