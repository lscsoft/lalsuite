/*
*  Copyright (C) 2007 Jolien Creighton, Julien Sylvestre, Kipp Cannon
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

/*-----------------------------------------------------------------------
 *
 * File Name: StdBurstSearch.h
 *
 * Author: Julien Sylvestre
 *
 * Revision: $$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="StdBurstSearchHV">
Author: Julien Sylvestre
$$
</lalVerbatim>
<lalLaTeX file="StdBurstSearchH">
\section{Header \texttt{StdBurstSearch.h}}
\label{s:StdBurstSearch.h}

Provides a standard environment for burst event trigger generators.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/StdBurstSearch.h>
\end{verbatim}

\subsection*{Error Conditions}
\input{STDBURSTSEARCHHErrTab}

</lalLaTeX>
#endif

#ifndef _STDBURSTSEARCH_H
#define _STDBURSTSEARCH_H

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

#include <lal/LALStdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/IIRFilter.h>
#include <lal/RealFFT.h>


NRCSID( STDBURSTSEARCHH, "$Id$" );

#include <lal/LIGOMetadataTables.h>

/******************************** <lalErrTable file="STDBURSTSEARCHHErrTab"> */
#define STDBURSTSEARCHH_ENULLP       1
#define STDBURSTSEARCHH_EMEM 2
#define STDBURSTSEARCHH_EUINP 3
#define STDBURSTSEARCHH_EIN 4
#define STDBURSTSEARCHH_ENULLPI       5

#define STDBURSTSEARCHH_MSGENULLP "Null pointer"
#define STDBURSTSEARCHH_MSGEMEM "memory error"
#define STDBURSTSEARCHH_MSGEUINP "unimplemented feature requested"
#define STDBURSTSEARCHH_MSGEIN "invalid input"
#define STDBURSTSEARCHH_MSGENULLPI "Parameter is not of the right type, or null pointer"
/*************************************************** </lalErrTable> */
/* <lalVerbatim file="StdBurstSearchH"> */

#define STDBURSTSEARCHSKIP_STARTTIME 1
#define STDBURSTSEARCHSKIP_CENTRALFREQ 2
#define STDBURSTSEARCHSKIP_DURATION 4
#define STDBURSTSEARCHSKIP_BANDWIDTH 8
#define STDBURSTSEARCHSKIP_AMPLITUDE 16
#define STDBURSTSEARCHSKIP_CONFIDENCE 32
#define STDBURSTSEARCHSKIP_SNR 128


typedef struct tagBurstParameter {
  struct tagBurstParameter *next;
  CHAR *char_;
  INT4 *int4_;
  REAL4 *real4_;
  REAL4Vector *real4vector_;

  INT4 random;
} BurstParameter;
/* </lalVerbatim> */

/* <lalVerbatim file="StdBurstSearchH"> */
void
LALTFClustersETG(
		 LALStatus *status,
		 EventIDColumn *output,
		 REAL4TimeVectorSeries *input,
		 BurstParameter *params);
/* </lalVerbatim> */

/* <lalVerbatim file="StdBurstSearchH"> */
void
LALSlopeETG(
		 LALStatus *status,
		 EventIDColumn *output,
		 REAL4TimeVectorSeries *input,
		 BurstParameter *params
		 );
/* </lalVerbatim> */

/* <lalVerbatim file="StdBurstSearchH"> */
typedef struct
tagBurstOutputDataSegment
{
  REAL4TimeVectorSeries              *data;
  REAL4FrequencySeries         *spec;
  COMPLEX8FrequencySeries      *resp;
  REAL4IIRFilter               *preprocessing_filter;
}
BurstOutputDataSegment;
/* </lalVerbatim> */


/* <lalVerbatim file="StdBurstSearchH"> */
typedef struct tagBurstOutputParameters {

  BurstOutputDataSegment *data;    /* input data and metadata */

  INT4 method; /* 0 for plain copy; 1 for standard */

  UINT4 skip; /* estimation to skip */

} BurstOutputParameters;
/* </lalVerbatim> */


typedef struct tagBurstOutputSpecStat {

  struct tagBurstOutputSpecStat *next;

  REAL8 duration;

  UINT4 nTime;

  UINT4 nFreq;

  REAL4 *wwin;

  REAL8 norm;

  RealFFTPlan *pfwd;

  COMPLEX8FrequencySeries *resp;

  REAL8 *P0;

  REAL8 *Q;

} BurstOutputSpecStat;

typedef struct tagRiceLikelihoodParams {

  REAL8 P;

  REAL8 P0;

  REAL8 Q;

} RiceLikelihoodParams;

void
LALRiceLikelihood(
		  LALStatus *status,
		  REAL8 *llik,
		  REAL8 P,
		  void *params
		  );


/* <lalVerbatim file="StdBurstSearchH"> */
void
LALBurstOutput(
	       LALStatus *status,
	       EventIDColumn *output,  /* output linked list of events */
	       EventIDColumn *input,   /* linked list of events from ETG */
	       BurstOutputParameters *params
);
/* </lalVerbatim> */


#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _STDBURSTSEARCH_H */
