#include "mex.h"
#include <math.h>
#include <strings.h>
#include <lal/LALStdlib.h>
#include <lal/StdBurstSearch.h>
#include <lal/AVFactories.h>

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

int lalDebugLevel = 0;

void  mexFunction( int nlhs,       mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{ 
  LALStatus status;
  char ebuf[2048];

  unsigned int i, N, Ndata, Nresponse;
  double *data;
  double *data_t0;
  double *data_dt;
  double *response_re;
  double *response_im;
  double *response_f0;
  double *response_df;
  double *start_time_s;
  double *start_time_ns;
  double *duration;
  double *central_freq;
  double *bandwidth;

  BurstOutputParameters oparams;
  EventIDColumn output, input, *ptr;
  BurstOutputDataSegment Mdata;
  REAL4TimeSeries inputdata;
  REAL4TimeVectorSeries inputdatav;
  COMPLEX8FrequencySeries resp;

  double *out0, *out1, *out2, *out3, *out4, *out5, *out6, *out7;

  REAL4VectorSequence vsdata;

  bzero(&status, sizeof(LALStatus));

  if(nlhs != 8 || nrhs != 12) {
    mexErrMsgTxt("[start_time_s,start_time_ns,duration,central_freq,bandwidth,amplitude,snr,confidence] = LALBurstOutput(data,data_start_time,data_dt,response_re,response_im,response_f0,response_df,start_time_s,start_time_ns,duration,central_freq,bandwidth)\n");
    return;
  }

  Ndata = mxGetN(prhs[0]);
  if(mxGetM(prhs[0]) > Ndata) {
    Ndata = mxGetM(prhs[0]);
  }

  data = mxGetPr(prhs[0]);
  data_t0 = mxGetPr(prhs[1]);
  data_dt = mxGetPr(prhs[2]);

  Nresponse = mxGetN(prhs[3]);
  if(mxGetM(prhs[3]) > Nresponse) {
    Nresponse = mxGetM(prhs[3]);
  }

  response_re = mxGetPr(prhs[3]);
  response_im = mxGetPr(prhs[4]);
  response_f0 = mxGetPr(prhs[5]);
  response_df = mxGetPr(prhs[6]);

  N = mxGetN(prhs[7]);
  if(mxGetM(prhs[7]) > N) {
    N = mxGetM(prhs[7]);
  }

  start_time_s = mxGetPr(prhs[7]);
  start_time_ns = mxGetPr(prhs[8]);
  duration = mxGetPr(prhs[9]);
  central_freq = mxGetPr(prhs[10]);
  bandwidth = mxGetPr(prhs[11]);

  bzero(&input, sizeof(EventIDColumn));
  bzero(&output, sizeof(EventIDColumn));

  ptr = &input;
  for(i=0;i<N;i++) {
    ptr->next = (EventIDColumn *)calloc(1,sizeof(EventIDColumn));
    ptr = ptr->next;

    ptr->snglBurstTable = (SnglBurstTable *)calloc(1,sizeof(SnglBurstTable));

    ptr->snglBurstTable->start_time.gpsSeconds = (INT4)start_time_s[i];
    ptr->snglBurstTable->start_time.gpsNanoSeconds = (INT4)start_time_ns[i];
    ptr->snglBurstTable->duration = duration[i];
    ptr->snglBurstTable->central_freq = central_freq[i];
    ptr->snglBurstTable->bandwidth = bandwidth[i];
  }

  resp.data = NULL;
  resp.f0 = *response_f0;
  resp.deltaF = *response_df;

  LALCCreateVector(&status, &(resp.data), Nresponse);
  CHKST

  for(i=0;i<Nresponse;i++) {
    resp.data->data[i].re = response_re[i];
    resp.data->data[i].im = response_im[i];
  }

  inputdata.data = NULL;
  inputdata.epoch.gpsSeconds = (INT4)floor(*data_t0);
  inputdata.epoch.gpsNanoSeconds = (INT4)floor(1e9*(*data_t0 - floor(*data_t0)));
  inputdata.deltaT = *data_dt;
  LALCreateVector(&status,&(inputdata.data),Ndata);
  CHKST

  for(i=0;i<Ndata;i++) {
    inputdata.data->data[i] = data[i];
  }

  inputdatav.epoch = inputdata.epoch;
  inputdatav.deltaT = inputdata.deltaT;
  inputdatav.data = &vsdata;
  vsdata.length = 1;
  vsdata.vectorLength = inputdata.data->length;
  vsdata.data = inputdata.data->data;

  bzero(&Mdata, sizeof(BurstOutputDataSegment));
  Mdata.resp = &resp;
  Mdata.data = &inputdatav;
  oparams.data = &Mdata;
  oparams.method = 1;

  LALBurstOutput(&status, &output, &input, &oparams);
  CHKST

  LALBurstOutput(&status, NULL, NULL, NULL);
  CHKST

  plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[3] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[4] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[5] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[6] = mxCreateDoubleMatrix(N,1,mxREAL);
  plhs[7] = mxCreateDoubleMatrix(N,1,mxREAL);

  out0 = mxGetPr(plhs[0]);
  out1 = mxGetPr(plhs[1]);
  out2 = mxGetPr(plhs[2]);
  out3 = mxGetPr(plhs[3]);
  out4 = mxGetPr(plhs[4]);
  out5 = mxGetPr(plhs[5]);
  out6 = mxGetPr(plhs[6]);
  out7 = mxGetPr(plhs[7]);

  ptr = output.next;
  
  i=0;
  while(ptr) {
    EventIDColumn *optr = ptr;

    out0[i] = ptr->snglBurstTable->start_time.gpsSeconds;
    out1[i] = ptr->snglBurstTable->start_time.gpsNanoSeconds;
    out2[i] = ptr->snglBurstTable->duration;
    out3[i] = ptr->snglBurstTable->central_freq;
    out4[i] = ptr->snglBurstTable->bandwidth;
    out5[i] = ptr->snglBurstTable->amplitude;
    out6[i] = ptr->snglBurstTable->snr;
    out7[i] = ptr->snglBurstTable->confidence;

    ptr = ptr->next;
    i++;

    LALFree(optr->snglBurstTable);
    LALFree(optr);
  }

  LALDestroyVector(&status, &(inputdata.data));
  CHKST
  LALCDestroyVector(&status, &(resp.data));
  CHKST

  ptr = input.next;
  while(ptr) {
    EventIDColumn *optr = ptr;
    ptr = ptr->next;
    LALFree(optr->snglBurstTable);
    LALFree(optr);
  }

  return;
}
