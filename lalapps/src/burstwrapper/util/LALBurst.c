/* [start_time_s,start_time_ns,duration,central_freq,bandwidth] = LALBurst(n,0,1/16384,'tfclusters','test','string',0.08,'float',1,'int',0,'int',0,'int',0.125,'float',8,'float',8184,'float',0.5,'float',5,'int',0,'int',0,'int',0,'int',0,'int',0,'int',0,'int',2,'int',3,'int',4,'int',4,'int'); */

#include "mex.h"
#include <math.h>
#include <strings.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/StdBurstSearch.h>

#define REPORTSTATUSB( statusptr )                                             \
  do                                                                          \
  {                                                                           \
    LALStatus *ptr_;                                                          \
    for ( ptr_ = (statusptr); ptr_; ptr_ = ptr_->statusPtr )                  \
    {                                                                         \
      sprintf(ebuf, " | Level %i: %s | %s", ptr_->level, ptr_->Id, ebuf );             \
      if ( ptr_->statusCode )                                                 \
      {                                                                       \
        LALPrintError( "\tStatus code %i: %s\n", ptr_->statusCode,            \
                       ptr_->statusDescription );                             \
      }                                                                       \
      else                                                                    \
      {                                                                       \
        sprintf(ebuf, "\tStatus code 0: Nominal\n" );                        \
      }                                                                       \
      sprintf(ebuf, "\tfunction %s, file %s, line %i | %s",                     \
                     ptr_->function, ptr_->file, ptr_->line, ebuf );                \
    }                                                                         \
  } while ( 0 )

#define CHKST if(status.statusCode != 0) {sprintf(ebuf,"LAL Error!\n"); REPORTSTATUSB(&status); mexErrMsgTxt(ebuf); return;}

#define GetString(prhs_,ETG_) cptr = mxGetChars(prhs_); \
  Mstr = mxGetN(prhs_); \
  for(i=0;i<Mstr;i++) { \
    ETG_[i] = (char)cptr[i]; \
  } \
  ETG_[i] = 0;

int lalDebugLevel = 0;

void  mexFunction( int nlhs,       mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{ 
  LALStatus status;
  char ebuf[2048];

  unsigned int i, j, N, Ndata, Mstr;
  double *data;
  double *data_t0;
  double *data_dt;
  char ETG[64];

  mxChar *cptr;

  EventIDColumn output, *ptr;
  REAL4TimeVectorSeries input;
  REAL4VectorSequence vser;
  BurstParameter params, *par;

  double *out0, *out1, *out2, *out3, *out4;

  bzero(&status, sizeof(LALStatus));

  if(nlhs != 5 || nrhs < 5) {
    mexErrMsgTxt("[start_time_s,start_time_ns,duration,central_freq,bandwidth] = LALBurst(data,data_start_time,data_dt,ETG,params,param_type,...)\n");
    return;
  }

  Ndata = mxGetN(prhs[0]);
  if(mxGetM(prhs[0]) > Ndata) {
    Ndata = mxGetM(prhs[0]);
  }

  data = mxGetPr(prhs[0]);
  data_t0 = mxGetPr(prhs[1]);
  data_dt = mxGetPr(prhs[2]);

  GetString(prhs[3],ETG)

  bzero(&params,sizeof(BurstParameter));
  par = &params;

  for(j=4;j<nrhs-1;j+=2) {
    char buf[64];
    par->next = (BurstParameter *)calloc(1,sizeof(BurstParameter));
    par = par->next;
    GetString(prhs[j+1],buf)
    if(!strcmp(buf,"string")) {
      char pbuf[256];
      GetString(prhs[j],pbuf)
      par->char_ = (CHAR *)calloc(1+strlen(pbuf),sizeof(CHAR));
      strcpy(par->char_,pbuf);
    } else {
      if(!strcmp(buf,"int")) {
	par->int4_ = (INT4 *)malloc(sizeof(INT4));
	*(par->int4_) = (INT4)floor(*(mxGetPr(prhs[j])));
      } else {
	if(!strcmp(buf,"float")) {
	  par->real4_ = (REAL4 *)malloc(sizeof(REAL4));
	  *(par->real4_) = (REAL4)(*(mxGetPr(prhs[j])));
	} else {
	  if(!strcmp(buf,"vector")) {
	    unsigned int k, nv = mxGetN(prhs[j]);
	    double *d = mxGetPr(prhs[j]);
	    if(mxGetM(prhs[j])>nv) {
	      nv = mxGetM(prhs[j]);
	    }
	    par->real4vector_ = NULL;
	    LALCreateVector(&status,&(par->real4vector_),nv);
	    CHKST
	    for(k=0;k<nv;k++) {
	      par->real4vector_->data[k] = d[k];
	    }
	  } else {
	    char bu[256];
	    sprintf(bu,"Unknown parameter type: %s\n",buf);
	    mexErrMsgTxt(bu);
	    return;
	  }
	}
      }
    }
  }

  input.epoch.gpsSeconds = (INT4)floor(*data_t0);
  input.epoch.gpsNanoSeconds = (INT4)floor(1e9*(*data_t0 - floor(*data_t0)));
  input.deltaT = *data_dt;
  input.data = &vser;

  vser.length = 1;
  vser.vectorLength = Ndata;
  vser.data = (REAL4 *)LALMalloc(Ndata * sizeof(REAL4));
  if(!(vser.data)) {
    mexErrMsgTxt("Memory Error!");
    return;
  }

  for(i=0;i<Ndata;i++) {
    input.data->data[i] = data[i];
  }


  if(!strcmp(ETG,"tfclusters")) {
    LALTFClustersETG(&status, &output, &input, &params);
    CHKST
  } else {
    if(!strcmp(ETG,"slope")) {
      LALSlopeETG(&status, &output, &input, &params);
      CHKST
    } else {
      if(!strcmp(ETG,"power")) {
	LALPowerETG(&status, &output, &input, &params);
	CHKST
      } else {
	mexErrMsgTxt("Unknown ETG!\n");
	return;
      }
    }
  }

  N = 0;
  ptr = output.next;
  while(ptr) {
    ptr = ptr->next;
    N++;
  }

  plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(N,1,mxREAL);

  out0 = mxGetPr(plhs[0]);
  out1 = mxGetPr(plhs[1]);
  out2 = mxGetPr(plhs[2]);
  out3 = mxGetPr(plhs[3]);
  out4 = mxGetPr(plhs[4]);

  ptr = output.next;

  i=0;
  while(ptr) {
    EventIDColumn *optr = ptr;

    out0[i] = ptr->snglBurstTable->start_time.gpsSeconds;
    out1[i] = ptr->snglBurstTable->start_time.gpsNanoSeconds;
    out2[i] = ptr->snglBurstTable->duration;
    out3[i] = ptr->snglBurstTable->central_freq;
    out4[i] = ptr->snglBurstTable->bandwidth;

    ptr = ptr->next;
    i++;

    LALFree(optr->snglBurstTable);
    LALFree(optr);
  }

  LALFree(vser.data);

  par = params.next;
  while(par) {
    BurstParameter *np;

    if(par->char_) {
      free(par->char_);
    }

    if(par->int4_) {
      free(par->int4_);
    }

    if(par->real4_) {
      free(par->real4_);
    }

    if(par->real4vector_) {
      LALDestroyVector(&status,&(par->real4vector_));
      CHKST
    }

    np = par->next;
    free(par);
    par = np;
  }

  return;
}
