/*
gcc -c ReadBurstBin.c -O
gcc -c  -I/home/jsylvest/pro/include -I/apps/matlab/extern/include -I/apps/matlab/simulink/include -DMATLAB_MEX_FILE -fPIC  -O -DNDEBUG reab_m.c
gcc -c  -I/home/jsylvest/pro/include -I/apps/matlab/extern/include -I/apps/matlab/simulink/include -DMATLAB_MEX_FILE -fPIC  -O -DNDEBUG /apps/matlab/extern/src/mexversion.c
gcc -O -shared -o reab.mexsol reab_m.o ReadBurstBin.o mexversion.o -L/home/jsylvest/pro/lib -llal -L/apps/matlab/bin/sol2 -lmx -lmex -lmatlb -lmat -lmwservices -lut -lm -lm
cp reab.mexsol ~/matlab

*/


#include "mex.h"

#define LIGOMETA_IFO_MAX 3
#define LIGOMETA_SEARCH_MAX 25
#define LIGOMETA_CHANNEL_MAX 65

int ReadBurstBin(
		 char *file,
		 unsigned int *N,
		 char ***ifo,
		 char ***search,
		 char ***channel,
		 unsigned int **start_time_s,
		 unsigned int **start_time_ns,
		 float **duration,
		 float **central_freq,
		 float **bandwidth,
		 float **amplitude,
		 float **snr,
		 float **confidence
		 );

void  mexFunction( int nlhs,       mxArray *plhs[],
                   int nrhs, const mxArray *prhs[])
{ 
  mxChar *mfile;
  char file[2048];
  unsigned int i, j, N;
  char **ifo;
  char **search;
  char **channel;
  unsigned int *start_time_s;
  unsigned int *start_time_ns;
  float *duration;
  float *central_freq;
  float *bandwidth;
  float *amplitude;
  float *snr;
  float *confidence;

  double *dat;
  mxChar *dac;
  int dims[2];

  if(nrhs != 1 || nlhs != 11) {
    mexErrMsgTxt("[ifo,search,channel,start_time_s,start_time_ns,duration,central_freq,bandwidth,amplitude,snr,confidence] = reab(file)\n");
    return;
  }

  mfile = mxGetChars(prhs[0]);
  for(i=0;i<mxGetN(prhs[0]);i++) {
    file[i] = mfile[i];
  }
  file[i] = 0;
  
  if(ReadBurstBin(file, &N, &ifo, &search, &channel, &start_time_s, &start_time_ns, &duration, &central_freq, &bandwidth, &amplitude, &snr, &confidence)) {
    mexErrMsgTxt("Error processing file\n");
    return;
  }

  dims[0] = N;
  dims[1] = LIGOMETA_IFO_MAX;
  plhs[0] = mxCreateCharArray(2,dims);
  dac = mxGetChars(plhs[0]);
  for(i=0;i<N;i++) {
    for(j=0;j<dims[1];j++) {
      dac[j*dims[0]+i] = (mxChar)ifo[i][j];
    }
  }
  
  dims[0] = N;
  dims[1] = LIGOMETA_SEARCH_MAX;
  plhs[1] = mxCreateCharArray(2,dims);
  dac = mxGetChars(plhs[1]);
  for(i=0;i<N;i++) {
    for(j=0;j<dims[1];j++) {
      dac[j*dims[0]+i] = (mxChar)search[i][j];
    }
  }
  
  dims[0] = N;
  dims[1] = LIGOMETA_CHANNEL_MAX;
  plhs[2] = mxCreateCharArray(2,dims);
  dac = mxGetChars(plhs[2]);
  for(i=0;i<N;i++) {
    for(j=0;j<dims[1];j++) {
      dac[j*dims[0]+i] = (mxChar)channel[i][j];
    }
  }

  plhs[3] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[3]);
  for(i=0;i<N;i++) {
    dat[i] = (double)start_time_s[i];
  }

  plhs[4] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[4]);
  for(i=0;i<N;i++) {
    dat[i] = (double)start_time_ns[i];
  }

  plhs[5] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[5]);
  for(i=0;i<N;i++) {
    dat[i] = (double)duration[i];
  }
  
  plhs[6] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[6]);
  for(i=0;i<N;i++) {
    dat[i] = (double)central_freq[i];
  }

  plhs[7] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[7]);
  for(i=0;i<N;i++) {
    dat[i] = (double)bandwidth[i];
  }

  plhs[8] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[8]);
  for(i=0;i<N;i++) {
    dat[i] = (double)amplitude[i];
  }

  plhs[9] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[9]);
  for(i=0;i<N;i++) {
    dat[i] = (double)snr[i];
  }

  plhs[10] = mxCreateDoubleMatrix(N,1,mxREAL);
  dat = mxGetPr(plhs[10]);
  for(i=0;i<N;i++) {
    dat[i] = (double)confidence[i];
  }

  for(i=0;i<N;i++) {
    free(ifo[i]);
    free(channel[i]);
    free(search[i]);
  }
  free(ifo);
  free(channel);
  free(search);
  free(start_time_s);
  free(start_time_ns);
  free(duration);
  free(central_freq);
  free(bandwidth);
  free(amplitude);
  free(snr);
  free(confidence);

  return;
}
