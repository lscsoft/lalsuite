/*
 * Copyright (C) 2013 Evan Ochsner, C. Pankow
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

#include <lal/LALSimInspiralSphHarmSeries.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Prepend a node to a linked list of SphHarmTimeSeries, or create a new head
 */
SphHarmTimeSeries* XLALSphHarmTimeSeriesAddMode(
            SphHarmTimeSeries *appended, /**< Linked list to be prepended */
            const COMPLEX16TimeSeries* inmode, /**< Time series of h_lm mode being prepended */
            UINT4 l, /**< l index of h_lm mode being prepended */
            INT4 m /**< m index of h_lm mode being prepended */
            )
{
    SphHarmTimeSeries* ts;

	// Check if the node with this l, m already exists
	ts = appended;
	while( ts ){
		if( l == ts->l && m == ts->m ){
			break;
		}
		ts = ts->next;
	}

	if( ts ){
		XLALDestroyCOMPLEX16TimeSeries( ts->mode );
		ts->mode = XLALCutCOMPLEX16TimeSeries( inmode, 0, inmode->data->length);
		return appended;
	} else {
    	ts = XLALMalloc( sizeof(SphHarmTimeSeries) );
	}

    ts->l = l;
    ts->m = m;
	// Cut returns a new series using a slice of the original. I ask it to
	// return a new one for the full data length --- essentially a duplication
	if( inmode ){
		ts->mode = XLALCutCOMPLEX16TimeSeries( inmode, 0, inmode->data->length);
	} else {
		ts->mode = NULL;
	}

    if( appended ){
        ts->next = appended;
		ts->tdata = appended->tdata;
    } else {
        ts->next = NULL;
		ts->tdata = NULL;
    }

    return ts;
}

/**
 * Set the tdata member for *all* nodes in the list.
 */
void XLALSphHarmTimeSeriesSetTData(
            SphHarmTimeSeries *ts, /**< Linked list to be prepended */
            REAL8Sequence* tdata /**< series of time data*/
            )
{
	while( ts ){
		ts->tdata = tdata;
		ts = ts->next;
	}
}

/**
 * Get the tdata member for nodes in the list.
 */
REAL8Sequence* XLALSphHarmTimeSeriesGetTData(
            SphHarmTimeSeries *ts /**< Get tdata from this list */
            )
{
	if( ts ){
		return ts->tdata;
	}
	return NULL;
}

/** Delete list from current pointer to the end of the list */
void XLALDestroySphHarmTimeSeries(
            SphHarmTimeSeries* ts /**< Head of linked list to destroy */
            )
{
    SphHarmTimeSeries* pop;
    while( (pop = ts) ){
		if( pop->mode ){
        	XLALDestroyCOMPLEX16TimeSeries( pop->mode );
		}
		// The tdata pointer is shared so we delete on the last node
		if( pop->next == NULL && pop->tdata ){
        	XLALDestroyREAL8Sequence( pop->tdata );
		}
        ts = pop->next;
        XLALFree( pop );
    }
}

/**
 * Get the time series of a waveform's (l,m) spherical harmonic mode from a
 * SphHarmTimeSeries linked list. Returns a pointer to its COMPLEX16TimeSeries
 */
COMPLEX16TimeSeries* XLALSphHarmTimeSeriesGetMode(
            SphHarmTimeSeries *ts, /** linked list to extract mode from */
            UINT4 l, /**< l index of h_lm mode to get */
            INT4 m /**< m index of h_lm mode to get */
            )
{
    if( !ts ) return NULL;

    SphHarmTimeSeries *itr = ts;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr->mode;
}

/**
 * Get the largest l index of any mode in the SphHarmTimeSeries linked list
 */
UINT4 XLALSphHarmTimeSeriesGetMaxL( SphHarmTimeSeries* ts ){
    SphHarmTimeSeries *itr = ts;
	UINT4 maxl=0;

    while( itr ){
		maxl = itr->l > maxl ? itr->l : maxl;
		itr = itr ->next;
    }
    return maxl;
}

/**
 * For every (l,m) node in the SphHarmTimeSeries linked list,
 * call XLALResizeCOMPLEX16TimeSeries(ts->mode, first, length)
 *
 * The TimeSeries of each (l,m) mode will have the given length,
 * and its contents will consist of that part of the original time series
 * that started at sample first. If first is negative, then the new time
 * series is padded at the start by that many samples. The time series' epoch
 * is adjusted appropriately.
 */
SphHarmTimeSeries *XLALResizeSphHarmTimeSeries(
        SphHarmTimeSeries *ts, /**< SphHarmTimeSeries to be resized */
        int first, /**< index of first time sample to be copied over */
        size_t length /**< length to resize all COMPLEX16TimeSeries to */
        )
{
    SphHarmTimeSeries *this = ts;
    while( this ) {
        this->mode = XLALResizeCOMPLEX16TimeSeries(this->mode, first, length);
        this = this->next;
    }

    return ts;
}

/**
 * Create a SphHarmFrequencySeries from a SphHarmTimeSeries
 * by performing an FFT on each mode in the SphHarmTimeSeries.
 */
SphHarmFrequencySeries *XLALSphHarmFrequencySeriesFromSphHarmTimeSeries(
        SphHarmTimeSeries *hlms_TD /**< SphHarmTimeSeries to be FFT'd */
        )
{
    UINT4 l, Lmax, length;
    int m;
    COMPLEX16TimeSeries *ht;
    COMPLEX16FrequencySeries *hf;
    SphHarmFrequencySeries *hlms_FD = NULL;
    REAL8 deltaF;
    if( !hlms_TD ) // Check head of linked list is valid
        XLAL_ERROR_NULL(XLAL_EINVAL);

    Lmax = XLALSphHarmTimeSeriesGetMaxL(hlms_TD);
    length = hlms_TD->mode->data->length; // N.B. Assuming all hlms same length
    deltaF = 1./hlms_TD->mode->deltaT/length;
    COMPLEX16FFTPlan *fwdplan = XLALCreateForwardCOMPLEX16FFTPlan(length, 0);
    hf = XLALCreateCOMPLEX16FrequencySeries( "FD Mode", &hlms_TD->mode->epoch,
            0., deltaF, &lalHertzUnit, length);
    // Loop over TD modes, FFT, add to SphHarmFrequencySeries
    for(l = 2; l <= Lmax; l++) {
        for(m = -l; m <= (int) l; m++) {
            ht = XLALSphHarmTimeSeriesGetMode(hlms_TD, l, m);
            if( ht ) {
                XLALCOMPLEX16TimeFreqFFT(hf, ht, fwdplan);
                hlms_FD = XLALSphHarmFrequencySeriesAddMode(hlms_FD, hf, l, m);
            }
        }
    }

    return hlms_FD;

}

SphHarmTimeSeries *XLALSphHarmTimeSeriesFromSphHarmFrequencySeriesDataAndPSD(
                                                                             SphHarmFrequencySeries *hlms, 
                                                                             COMPLEX16FrequencySeries* data,
                                                                             COMPLEX16FrequencySeries* psd
        )
{
  UINT4 l, Lmax, length,i;
    int m;
    COMPLEX16TimeSeries *rhoT;
    COMPLEX16FrequencySeries *hf,*hfBuffer;
    SphHarmTimeSeries *rhoTlm = NULL;
    REAL8 deltaF;
    COMPLEX16 wt;
    if( !hlms ) // Check head of linked list is valid
        XLAL_ERROR_NULL(XLAL_EINVAL);

    Lmax = XLALSphHarmFrequencySeriesGetMaxL(hlms);
    length = hlms->mode->data->length; // N.B. Assuming all hlms same length
    deltaF = hlms->mode->deltaF;
    COMPLEX16FFTPlan *revplan = XLALCreateReverseCOMPLEX16FFTPlan(length, 0);
    // Output working buffer : should be copied
    rhoT = XLALCreateCOMPLEX16TimeSeries( "rhoTD", &hlms->mode->epoch,
            0., deltaF, &lalDimensionlessUnit, length);
    hfBuffer = XLALCreateCOMPLEX16FrequencySeries( "FD Mode", &hlms->mode->epoch,
            0., deltaF, &lalHertzUnit, length);
    // Loop over TD modes, FFT, add to SphHarmFrequencySeries
    for(l = 2; l <= Lmax; l++) {
        for(m = -l; m <= (int) l; m++) {
            hf = XLALSphHarmFrequencySeriesGetMode(hlms, l, m);
            if( hf ) {
              hfBuffer->epoch = hf->epoch;
              hfBuffer->deltaF = hf->deltaF;
              for (i =0; i<length; i++) {
                if (psd->data->data[i]) {
                  wt = 1./psd->data->data[i];
                } else {
                  wt = 0;
                    }
                hfBuffer->data->data[i] = conj(hf->data->data[i])* data->data->data[i]*wt;
              }
              XLALCOMPLEX16FreqTimeFFT(rhoT, hfBuffer, revplan);
              rhoTlm = XLALSphHarmTimeSeriesAddMode(rhoTlm, rhoT, l, m);
            }
        }
    }
    return rhoTlm;
}


/**
 * Prepend a node to a linked list of SphHarmFrequencySeries, or create a new head
 */
SphHarmFrequencySeries* XLALSphHarmFrequencySeriesAddMode(
            SphHarmFrequencySeries *appended, /**< Linked list to be prepended */
            const COMPLEX16FrequencySeries* inmode, /**< Time series of h_lm mode being prepended */
            UINT4 l, /**< l index of h_lm mode being prepended */
            INT4 m /**< m index of h_lm mode being prepended */
            )
{
    SphHarmFrequencySeries* ts;

	// Check if the node with this l, m already exists
	ts = appended;
	while( ts ){
		if( l == ts->l && m == ts->m ){
			break;
		}
		ts = ts->next;
	}

	if( ts ){
		XLALDestroyCOMPLEX16FrequencySeries( ts->mode );
		ts->mode = XLALCutCOMPLEX16FrequencySeries( inmode, 0, inmode->data->length);
		return appended;
	} else {
    	ts = XLALMalloc( sizeof(SphHarmFrequencySeries) );
	}

    ts->l = l;
    ts->m = m;
	// Cut returns a new series using a slice of the original. I ask it to
	// return a new one for the full data length --- essentially a duplication
	if( inmode ){
		ts->mode = XLALCutCOMPLEX16FrequencySeries( inmode, 0, inmode->data->length);
	} else {
		ts->mode = NULL;
	}

    if( appended ){
        ts->next = appended;
		ts->fdata = appended->fdata;
    } else {
        ts->next = NULL;
		ts->fdata = NULL;
    }

    return ts;
}

/**
 * Set the tdata member for *all* nodes in the list.
 */
void XLALSphHarmFrequencySeriesSetFData(
            SphHarmFrequencySeries *ts, /**< Linked list to be prepended */
            REAL8Sequence* fdata /**< series of frequency data*/
            )
{
	while( ts ){
		ts->fdata = fdata;
		ts = ts->next;
	}
}

/**
 * Get the fdata member.
 */
REAL8Sequence* XLALSphHarmFrequencySeriesGetFData(
            SphHarmFrequencySeries *ts /**< Get tdata from this list */
            )
{
	if( ts ){
		return ts->fdata;
	}
	return NULL;
}

/** Delete list from current pointer to the end of the list */
void XLALDestroySphHarmFrequencySeries(
            SphHarmFrequencySeries* ts /**< Head of linked list to destroy */
            )
{
    SphHarmFrequencySeries* pop;
    while( (pop = ts) ){
		if( pop->mode ){
        	XLALDestroyCOMPLEX16FrequencySeries( pop->mode );
		}
		// The fdata pointer is shared so we delete on the last node
		if( pop->next == NULL && pop->fdata ){
        	XLALDestroyREAL8Sequence( pop->fdata );
		}
        ts = pop->next;
        XLALFree( pop );
    }
}

/**
 * Get the time series of a waveform's (l,m) spherical harmonic mode from a
 * SphHarmFrequencySeries linked list. Returns a pointer to its COMPLEX16FrequencySeries
 */
COMPLEX16FrequencySeries* XLALSphHarmFrequencySeriesGetMode(
            SphHarmFrequencySeries *ts, /** linked list to extract mode from */
            UINT4 l, /**< l index of h_lm mode to get */
            INT4 m /**< m index of h_lm mode to get */
            )
{
    if( !ts ) return NULL;

    SphHarmFrequencySeries *itr = ts;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr->mode;
}

/**
 * Get the largest l index of any mode in the SphHarmFrequencySeries linked list
 */
UINT4 XLALSphHarmFrequencySeriesGetMaxL( SphHarmFrequencySeries* ts ){
    SphHarmFrequencySeries *itr = ts;
	UINT4 maxl=0;

    while( itr ){
		maxl = itr->l > maxl ? itr->l : maxl;
		itr = itr ->next;
    }
    return maxl;
}
