/*
*  Copyright (C) 2007 Badri Krishnan, Reinhard Prix, Alicia Sintes Olives
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
 * File Name: SFTbin.c
 *
 * Authors: Sintes, A.M.,  Krishnan, B. & inspired from Siemens, X.
 *
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified by Krishnan on Feb 22, 2004
 *
 *-----------------------------------------------------------------------
 */

/** OBSOLETE -- Use LAL functions from SFTfileIO.c instead */

/**
\author Sintes, A.M., Krishnan, B.

\heading{\ref SFTbin.c}
Routines for reading SFT binary files

\heading{Prototypes}
<tt>ReadSFTbinHeader1()</tt>
<tt>ReadCOMPLEX8SFTbinData1()</tt>
<tt>ReadCOMPLEX16SFTbinData1()</tt>
<tt>COMPLEX8SFT2Periodogram1()</tt>
<tt>COMPLEX16SFT2Periodogram1()</tt>

\heading{Description}

the output of the periodogram should be properly normalized !!!

\heading{Uses}
\code
LALHO()
\endcode
*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "./SFTbin.h"

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void COMPLEX8SFT2Periodogram1 (LALStatus  *status,
                   REAL8Periodogram1     *peri,
		   COMPLEX8SFTData1      *sft)
{ 

  INT4       length;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (peri,  status, SFTBINH_ENULL, SFTBINH_MSGENULL);

  peri->epoch.gpsSeconds = sft->epoch.gpsSeconds;
  peri->epoch.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  peri->timeBase = sft->timeBase;
  peri->fminBinIndex = sft->fminBinIndex;
  
  length = sft->length;
  ASSERT (length==peri->length,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);

  if (length > 0){
    REAL8      *out;
    COMPLEX8   *in1;
    REAL8      re,im, factor;
    INT4      n;

    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (peri->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    out = peri->data;
    in1 = sft->data;
    n= length;

    /* if data was properly normalized */
    factor = 1./peri->timeBase;
    
    /* if data is not normalized... this factor needs to be clarified  */
    /* note inconsistency with the previous function and quantity *factor*  there */
    
    /* if bin zero is included should be treated properly because of factor 2 */
    
    while (n-- >0){
      re = crealf(*in1);
      im = cimagf(*in1);
      *out = (re*re + im*im) *factor; /* factor 2 still missing if one-sided*/
      ++out;
      ++in1;
    }
  }
  
  

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
 



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void SFT2Periodogram (LALStatus  *status,
                   REAL8Periodogram1     *peri,
		   SFTtype      *sft)
{ 

  INT4       length;
  REAL8      f0, deltaF;
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (peri,  status, SFTBINH_ENULL, SFTBINH_MSGENULL);

  f0 = sft->f0;
  deltaF = sft->deltaF;  
  length = sft->data->length;

  peri->epoch.gpsSeconds = sft->epoch.gpsSeconds;
  peri->epoch.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  peri->timeBase = 1.0/deltaF;
  peri->fminBinIndex = floor(f0/deltaF + 0.5);
  
  ASSERT (length==peri->length,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);

  if (length > 0){
    REAL8      *out;
    COMPLEX8   *in1;
    REAL8      re,im, factor;
    INT4      n;

    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (peri->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    out = peri->data;
    in1 = sft->data->data;
    n= length;

    /* if data was properly normalized */
    factor = 1./peri->timeBase;
    
    /* if data is not normalized... this factor needs to be clarified  */
    /* note inconsistency with the previous function and quantity *factor*  there */
    
    /* if bin zero is included should be treated properly because of factor 2 */
    
    while (n-- >0){
      re = crealf(*in1);
      im = cimagf(*in1);
      *out = (re*re + im*im) *factor; /* factor 2 still missing if one-sided*/
      ++out;
      ++in1;
    }
  }
  
  

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
 








