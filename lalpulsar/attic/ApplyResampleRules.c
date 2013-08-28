/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

/**
 * \author Creighton, T. D.
 * \ingroup Resample_h
 * \brief Resamples a time series according to a set of resampling rules.
 *
 * This function sets <tt>output->deltaT</tt> and fills <tt>output->data</tt>
 * with data from <tt>*input</tt>, using the resampling rules specified in
 * <tt>*rules</tt>.  If the timespan required to fill <tt>output->data</tt>
 * is not a subset of the timespan covered by <tt>*input</tt> or
 * <tt>*rules</tt>, the data at the nonintersecting times are set to zero.
 *
 * \heading{Algorithm}
 *
 * At present this routine is just a stub.  It does not apply or even
 * check <tt>*rules</tt>, and instead simply makes <tt>*output</tt>
 * equivalent to (a subset of) <tt>*input</tt>.
 */
void
LALApplyResampleRules( LALStatus       *stat,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *input,
		       ResampleRules   *rules )
{
  INT4 nStart, nStop; /* output domain for which we can get data */

  INITSTATUS(stat);

  /* Check that the inputs all exist. */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(output->data->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(input->data->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Set the output sampling time. */
  output->deltaT=input->deltaT;

  /* Find the difference between the input and output start times, in
     samples. */
  nStart=(INT4)((input->epoch.gpsSeconds-output->epoch.gpsSeconds)
		/output->deltaT);
  nStart+=(INT4)((input->epoch.gpsNanoSeconds
		  -output->epoch.gpsNanoSeconds)
		 /(1e9*output->deltaT));
  if(nStart>(INT4)(output->data->length))
    nStart=output->data->length;

  /* Ditto for stop times. */
  nStop=nStart+input->data->length; /* since deltaT's are equal */
  if(nStop>(INT4)(output->data->length))
    nStop=output->data->length;

  if(nStart>0)
    memset(output->data->data,0,nStart*sizeof(REAL4));
  if(nStop>nStart)
    memcpy(output->data->data+nStart,input->data->data,
	   (nStop-nStart)*sizeof(REAL4));
  if((INT4)(output->data->length)>nStop)
    memset(output->data->data+nStop,0,
	   (output->data->length-nStop)*sizeof(REAL4));

  (void)rules;

  /* That's all for the current stub. */
  RETURN(stat);
}
