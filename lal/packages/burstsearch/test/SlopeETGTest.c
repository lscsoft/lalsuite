#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/StdBurstSearch.h>


NRCSID (MAIN, "$Id$");

#define CHKST if(status.statusCode != 0) {REPORTSTATUS (&status); return -1;}

#define SetStringParameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->char_ = (CHAR *)LALCalloc(1+strlen(channel), sizeof(CHAR)); \
  strcpy(par_->char_, str_)

#define SetINT4Parameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->int4_ = (INT4 *)LALCalloc(1, sizeof(INT4)); \
  *(par_->int4_) = str_

#define SetREAL4Parameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->real4_ = (REAL4 *)LALCalloc(1, sizeof(REAL4)); \
  *(par_->real4_) = str_

#define SetREAL4VectorParameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->real4vector_ = NULL; \
  LALCreateVector(&status, &(par_->real4vector_), str_->length); \
  CHKST \
  memcpy(par_->real4vector_->data, str_->data, str_->length * sizeof(REAL4))

INT4 lalDebugLevel = 2; /* set to 2 for full debug */


INT4 main(INT4 argc, CHAR *argv[]) {

  static LALStatus status;
  BurstParameter params, *par;
  EventIDColumn output, *optr;
  REAL4TimeSeries input;
  REAL4Vector *data = NULL;
  RandomParams *rpar = NULL;

  REAL4Vector *coef = NULL;

  CHAR channel[10] = "H1:test";
  UINT4 N = 10*16384; /* length of data */
  UINT4 i;

  /* create ramp */
  LALCreateVector(&status, &coef, 10);
  CHKST
    
  for(i=0;i<coef->length;i++) {
    coef->data[i] = (12.0*16384.0/990.0)*((REAL4)i - 4.5);
  }

  /* set parameters */
  par = &params;
  
  SetStringParameter(par,channel); /* channel name */
  SetREAL4Parameter(par, 7000.0);     /* threshold */
  SetREAL4VectorParameter(par,coef); /* filter coefficients */
  SetINT4Parameter(par, 49);        /* number of points in clustering */
  SetREAL4Parameter(par, 800.0);        /* min frequency */
  SetREAL4Parameter(par, 1200.0);     /* max frequency */

  /* create input */
  LALCreateRandomParams(&status, &rpar, 0);
  CHKST

  LALCreateVector(&status, &data, N);
  CHKST

  input.deltaT = 1.0/16384.0;
  input.epoch.gpsSeconds = 700000000;
  input.epoch.gpsNanoSeconds = 0;
  input.data = data;

  LALNormalDeviates(&status, data, rpar);
  CHKST

  LALDestroyRandomParams(&status, &rpar);
  CHKST

  /* run ETG */
  LALSlopeETG(&status, &output, &input, &params);
  CHKST

  /* report output */
  optr = output.next;
  while(optr) {

      printf("%s %s\t%s\t%u\t%u\t%g\t%g\t%g\n",optr->snglBurstTable->ifo, optr->snglBurstTable->search, optr->snglBurstTable->channel, optr->snglBurstTable->start_time.gpsSeconds, optr->snglBurstTable->start_time.gpsNanoSeconds, optr->snglBurstTable->duration, optr->snglBurstTable->central_freq, optr->snglBurstTable->bandwidth);
    
    optr = optr->next;

  }


  /* clean up */
  LALDestroyVector(&status, &coef);
  CHKST

  LALDestroyVector(&status, &data);
  CHKST


  optr = output.next;
  while(optr) {
    EventIDColumn *optrn = optr->next;

    if(optr->snglTransdataTable) {
      if(optr->snglTransdataTable->trans_data) {
	LALFree(optr->snglTransdataTable->trans_data);
      }
      LALFree(optr->snglTransdataTable);
    }
    
    LALFree(optr->snglBurstTable);
    LALFree(optr);

    optr = optrn;

  }


  par = params.next;
  while(par) {
    BurstParameter *tpar = par->next;

    if(par->char_) {
      LALFree(par->char_);
    }

    if(par->int4_) {
      LALFree(par->int4_);
    }

    if(par->real4_) {
      LALFree(par->real4_);
    }

    if(par->real4vector_) {
      if(par->real4vector_->data) {
	LALDestroyVector(&status, &(par->real4vector_));
	CHKST
      }
    }

    LALFree(par);
    par = tpar;
  }

  return 0;

}
