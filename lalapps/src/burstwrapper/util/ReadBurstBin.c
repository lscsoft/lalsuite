#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIGOMETA_IFO_MAX 3
#define LIGOMETA_SEARCH_MAX 25
#define LIGOMETA_CHANNEL_MAX 65

#ifdef linux
  static void endian_swap(char * pdata, int dsize, int nelements);
#endif

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
		 )
{

  FILE *in;
  unsigned int i, length;
  char *buffer;
  size_t size1 = (LIGOMETA_IFO_MAX + LIGOMETA_SEARCH_MAX + LIGOMETA_CHANNEL_MAX) * sizeof(char) + 2*sizeof(int) + 6 * sizeof(float);

  *N = 0;
  *ifo = *search = *channel = NULL;
  *start_time_s = *start_time_ns = NULL;
  *duration = *central_freq = *bandwidth = *amplitude = *snr = *confidence = NULL;

  if(!(in = fopen(file,"r"))) {
    fprintf(stderr,"Can't open file %s\n",file);
    return 1;
  }

  if(!fread(N, sizeof(unsigned int), 1, in)) {
    fclose(in);
    fprintf(stderr,"Can't read length from file %s\n",file);
    return 1;
  }

#ifdef linux
  endian_swap((char *)N, sizeof(unsigned int), 1);
#endif

  if(!(*N)) {
    fclose(in);
    return 0;
  }

  length = *N * size1;

  if(!(*ifo = (char **)calloc(*N, sizeof(char *)))) {
    fprintf(stderr,"1 Malloc error: %i\n",*N);
    fclose(in);
    return 1;
  }

  if(!(*search = (char **)calloc(*N, sizeof(char *)))) {
    fprintf(stderr,"2 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*channel = (char **)calloc(*N, sizeof(char *)))) {
    fprintf(stderr,"3 Malloc error\n");
    fclose(in);
    return 1;
  }


  for(i=0;i<*N;i++) {
      if(!((*ifo)[i] = (char *)calloc(LIGOMETA_IFO_MAX, sizeof(char)))) {
	fprintf(stderr,"4 Malloc error\n");
	fclose(in);
	return 1;
      }

      if(!((*search)[i] = (char *)calloc(LIGOMETA_SEARCH_MAX, sizeof(char)))) {
	fprintf(stderr,"5 Malloc error\n");
	fclose(in);
	return 1;
      }

      if(!((*channel)[i] = (char *)calloc(LIGOMETA_CHANNEL_MAX, sizeof(char)))) {
	fprintf(stderr,"6 Malloc error\n");
	fclose(in);
	return 1;
      }

  }


  if(!(*start_time_s = (unsigned int *)malloc(*N * sizeof(unsigned int)))) {
    fprintf(stderr,"7 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*start_time_ns = (unsigned int *)malloc(*N * sizeof(unsigned int)))) {
    fprintf(stderr,"8 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*duration = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"9 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*central_freq = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"10 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*bandwidth = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"11 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*amplitude = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"12 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*snr = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"13 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(*confidence = (float *)malloc(*N * sizeof(float)))) {
    fprintf(stderr,"14 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!(buffer=(char *)malloc(length * sizeof(char)))) {
    fprintf(stderr,"15 Malloc error\n");
    fclose(in);
    return 1;
  }

  if(!fread(buffer, sizeof(char), length, in)) {
    fclose(in);
    fprintf(stderr,"Can't read data from file %s\n",file);
    return 1;
  }


  for(i=0;i<*N;i++) {
    char *p1 = buffer + i*size1, *pst, *pr4;

    memcpy((*ifo)[i], p1, LIGOMETA_IFO_MAX * sizeof(char));

    p1 += LIGOMETA_IFO_MAX * sizeof(char);
    memcpy((*search)[i], p1,LIGOMETA_SEARCH_MAX * sizeof(char));

    p1 += LIGOMETA_SEARCH_MAX * sizeof(char);
    memcpy((*channel)[i], p1, LIGOMETA_CHANNEL_MAX * sizeof(char));

    pst = p1 + LIGOMETA_CHANNEL_MAX * sizeof(char);

#ifdef linux
    endian_swap((char *)(pst), sizeof(int), 2);
#endif

    memcpy(*start_time_s + i, pst, sizeof(int));
    memcpy(*start_time_ns + i, pst + sizeof(int), sizeof(int));

    pr4 = pst + 2*sizeof(int);

#ifdef linux
    endian_swap((char *)(pr4), sizeof(float), 6);
#endif

    memcpy(*duration + i, pr4, sizeof(float));

    p1 = pr4 + sizeof(float);
    memcpy(*central_freq + i, p1, sizeof(float));

    p1 += sizeof(float);
    memcpy(*bandwidth + i, p1, sizeof(float));

    p1 += sizeof(float);
    memcpy(*amplitude + i, p1, sizeof(float));
    
    p1 += sizeof(float);
    memcpy(*snr + i, p1, sizeof(float));

    p1 += sizeof(float);
    memcpy(*confidence + i, p1, sizeof(float));

  }

  fclose(in);

  free(buffer);

  return 0;
}


#ifdef linux
static void endian_swap(char * pdata, int dsize, int nelements)

{

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
