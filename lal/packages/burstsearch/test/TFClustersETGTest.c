/*
*  Copyright (C) 2007 Jolien Creighton, Julien Sylvestre
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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

#define SetINT4Parameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->int4_ = (INT4 *)LALCalloc(1, sizeof(INT4)); \
  *(par_->int4_) = str_

#define SetREAL4Parameter(par_, str_) par_->next = (BurstParameter *)LALCalloc(1, sizeof(BurstParameter)); \
  par_ = par_->next; \
  par_->real4_ = (REAL4 *)LALCalloc(1, sizeof(REAL4)); \
  *(par_->real4_) = str_

#ifndef WORDS_BIGENDIAN
  static void endian_swap(char * pdata, int dsize, int nelements);
#endif

INT4 lalDebugLevel = 0; /* set to 2 for full debug, painfully slow with tfclusters */


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

  SetStringParameter(par,channel); /* channel name */
  SetREAL4Parameter(par, 0.08);    /* black pixel probability */
  SetINT4Parameter(par, 1);        /* windowing? */
  SetINT4Parameter(par, 0);        /* threshold method (1=-log(p), 0=rank) */
  SetINT4Parameter(par, 1);        /* set to 1 to report tf information */
  SetREAL4Parameter(par, 0.125);   /* time resolution */
  SetREAL4Parameter(par, 8);        /* min frequency */
  SetREAL4Parameter(par, 8184);     /* max frequency */
  SetREAL4Parameter(par, 0.5);     /* alpha */
  SetINT4Parameter(par, 5);        /* sigma */
  SetINT4Parameter(par, 0);        /* delta(1,1) */
  SetINT4Parameter(par, 0);        /* delta(2,1) */
  SetINT4Parameter(par, 0);        /* delta(3,1) */
  SetINT4Parameter(par, 0);        /* delta(4,1) */
  SetINT4Parameter(par, 0);        /* delta(2,2) */
  SetINT4Parameter(par, 0);        /* delta(3,2) */
  SetINT4Parameter(par, 2);        /* delta(4,2) */
  SetINT4Parameter(par, 3);        /* delta(3,3) */
  SetINT4Parameter(par, 4);        /* delta(4,3) */
  SetINT4Parameter(par, 4);       /* delta(4,4) */

  /* create input */
  LALCreateRandomParams(&status, &rpar, 0);
  CHKST

  LALCreateVector(&status, &data, 2*N);
  CHKST

  input.deltaT = 1.0/16384.0;
  input.epoch.gpsSeconds = 700000000;
  input.epoch.gpsNanoSeconds = 0;
  input.data = &vseq;

  vseq.length = 2;
  vseq.vectorLength = N;
  vseq.data = data->data;

  LALNormalDeviates(&status, data, rpar);
  CHKST

  LALDestroyRandomParams(&status, &rpar);
  CHKST

  /* run ETG */
  LALTFClustersETG(&status, &output, &input, &params);
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

