#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/StdBurstSearch.h>
#include <lal/AVFactories.h>

NRCSID (MAIN, "$Id$");

#define CHKST if(status.statusCode != 0) {REPORTSTATUS (&status); return -1;}

#define SetStringParameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->char_ = (CHAR *)LALCalloc(1+strlen(str_), sizeof(CHAR)); \
  strcpy(par_->char_, str_)

#ifndef WORDS_BIGENDIAN
  static void endian_swap(char * pdata, int dsize, int nelements);
#endif

INT4 lalDebugLevel = 2; /* set to 2 for full debug */


INT4 main(INT4 argc, CHAR *argv[]) {

  static LALStatus status;
  BurstParameter params, *par;
  EventIDColumn output, *optr;
  REAL4TimeVectorSeries input;
  REAL4VectorSequence vseq;
  REAL4Vector *data = NULL;
  RandomParams *rpar = NULL;

  CHAR channel[10] = "H1:test";
  UINT4 N = 10*16384; /* length of data */
  UINT4 i;

  /* set parameters */
  par = &params;
  
  SetStringParameter(par,"-filterparams");
  SetStringParameter(par,"163840");    /* number of points per segment */
  SetStringParameter(par,"1");         /* number of segments */
  SetStringParameter(par,"0");         /* overlap of each segment */
  SetStringParameter(par,"3");         /* ? */
  SetStringParameter(par,"2");         /* ? */
  SetStringParameter(par,"2");         /* ? */
  SetStringParameter(par,"100.0");     /* low frequency bound of search */
  SetStringParameter(par,"1.0");       /* frequency steps */
  SetStringParameter(par,"1024.0");    /* bandwidth above low freq */
  SetStringParameter(par,"2.0");       /* ? */
  SetStringParameter(par,"0.5");       /* ? */
  SetStringParameter(par,"1");         /* segment analyzed at one time */
  SetStringParameter(par,"1.0e-36");   /* dso threshold */
  SetStringParameter(par,"10");        /* events to communicate to master */
  SetStringParameter(par,channel);     /* channel name */
  SetStringParameter(par,"0");         /* ? */
  SetStringParameter(par,"useMedian"); /* ? */
  SetStringParameter(par,"2");          /* ? */

  /* create input */
  LALCreateRandomParams(&status, &rpar, 0);
  CHKST

  LALCreateVector(&status, &data, N);
  CHKST

  input.deltaT = 1.0/16384.0;
  input.epoch.gpsSeconds = 700000000;
  input.epoch.gpsNanoSeconds = 0;
  input.data = &vseq;

  vseq.length = 1;
  vseq.vectorLength = N;
  vseq.data = data->data;

  LALNormalDeviates(&status, data, rpar);
  CHKST

  LALDestroyRandomParams(&status, &rpar);
  CHKST

  /* run ETG */
  LALPowerETG(&status, &output, &input, &params);
  CHKST

  /* report output */
  optr = output.next;
  while(optr) {

    if(optr->snglTransdataTable) {
      UINT4 nclusters;

      printf("%s %s\t%s\t%u\t%u\t%g\t%g\t%g\t1\n",optr->snglBurstTable->ifo, optr->snglBurstTable->search, optr->snglBurstTable->channel, optr->snglBurstTable->start_time.gpsSeconds, optr->snglBurstTable->start_time.gpsNanoSeconds, optr->snglBurstTable->duration, optr->snglBurstTable->central_freq, optr->snglBurstTable->bandwidth);
      printf("***************************************************\n");
      printf("%s %s\t%i\t%i\t%g\t%g\t%s\t%i\t%g\t%g\t%s\t%s\t%s\t%i\n", optr->snglTransdataTable->ifo, optr->snglTransdataTable->name, optr->snglTransdataTable->dimensions, optr->snglTransdataTable->x_bins, optr->snglTransdataTable->x_start, optr->snglTransdataTable->x_end, optr->snglTransdataTable->x_units, optr->snglTransdataTable->y_bins, optr->snglTransdataTable->y_start, optr->snglTransdataTable->y_end, optr->snglTransdataTable->y_units, optr->snglTransdataTable->data_type, optr->snglTransdataTable->data_units, optr->snglTransdataTable->transdata_length);

      memcpy(&nclusters, optr->snglTransdataTable->trans_data, sizeof(UINT4));

#ifndef WORDS_BIGENDIAN
      endian_swap((char *)(&nclusters), sizeof(UINT4), 1);
#endif

      for(i=0;i<nclusters;i++) {

	UINT4 ti,fi;
	REAL8 P;

	memcpy(&ti, optr->snglTransdataTable->trans_data+sizeof(UINT4)+i*(2*sizeof(UINT4)+sizeof(REAL8)), sizeof(UINT4));
	memcpy(&fi, optr->snglTransdataTable->trans_data+sizeof(UINT4)+sizeof(UINT4)+i*(2*sizeof(UINT4)+sizeof(REAL8)), sizeof(UINT4));
	memcpy(&P, optr->snglTransdataTable->trans_data+2*sizeof(UINT4)+sizeof(UINT4)+i*(2*sizeof(UINT4)+sizeof(REAL8)), sizeof(REAL8));

#ifndef WORDS_BIGENDIAN
      endian_swap((char *)(&ti), sizeof(UINT4), 1);
      endian_swap((char *)(&fi), sizeof(UINT4), 1);
      endian_swap((char *)(&P), sizeof(REAL8), 1);
#endif

	printf("%u\t%u\t%g\n", ti, fi, P);

      }
      printf("***************************************************\n");
    } else {
      printf("%s %s\t%s\t%u\t%u\t%g\t%g\t%g\t0\n",optr->snglBurstTable->ifo, optr->snglBurstTable->search, optr->snglBurstTable->channel, optr->snglBurstTable->start_time.gpsSeconds, optr->snglBurstTable->start_time.gpsNanoSeconds, optr->snglBurstTable->duration, optr->snglBurstTable->central_freq, optr->snglBurstTable->bandwidth);
    }

    optr = optr->next;

  }


  /* clean up */
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


#ifndef WORDS_BIGENDIAN
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

