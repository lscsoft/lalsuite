#include <stdio.h>


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


int main(int argc, char *argv[]) {

  char *file;
  unsigned int i;

  unsigned int N;
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

  if(argc != 2) {
    fprintf(stderr,"reab file\n");
    return 1;
  }

  file = argv[1];

  if(ReadBurstBin(file, &N, &ifo, &search, &channel, &start_time_s, &start_time_ns, &duration, &central_freq, &bandwidth, &amplitude, &snr, &confidence)) {
    fprintf(stderr,"Error processing file %s\n",file);
    return 1;
  }

  for(i=0;i<N;i++) {

    printf("%.3u %.3s %.25s %.65s %.9i %.9i %.9g %.4g %.4g %.6g %.6g %.6g\n",i,ifo[i],search[i],channel[i],start_time_s[i],start_time_ns[i],duration[i],central_freq[i],bandwidth[i],amplitude[i],snr[i],confidence[i]);

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

  return 0;
}
