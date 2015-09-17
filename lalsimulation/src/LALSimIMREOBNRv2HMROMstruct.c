/*
 *  Copyright (C) 2014 Sylvain Marsat
 *  Reduced Order Model for EOBNRv2HM
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

/**
 * \author Sylvain Marsat
 *
 * \file
 *
 * \brief C code for structures EOBNRv2HM reduced order model (non-spinning version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details on the reduced order method.
 * See arXiv:1106.1021 for the EOBNRv2HM model.
 *
 * Borrows from the SEOBNR ROM LAL code written by Michael Puerrer and John Veitch.
 *
 * The binary data files are available at [TBD].
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * Parameter ranges:
 *   q = 1-6
 *   No spin
 *   Mtot >= 20Msun for fstart=9Hz
 *
 */


#define _XOPEN_SOURCE 500

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <stdbool.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>

#include "LALSimIMREOBNRv2HMROMstruct.h"
#include "LALSimIMREOBNRv2HMROM.h"


/************************************************************************/
/********************* Functions for list structures ********************/

/***************** Functions for the SplineList structure ****************/

/* Prepend a node to a linked list of splines, or create a new head */
SplineList* SplineList_AddElementNoCopy(
	   SplineList* appended,  /* List structure to prepend to */
	   gsl_spline* spline,  /* spline to contain */
           gsl_interp_accel* accel,  /* accelerator to contain */
	   UINT4 i /* index in the list */)
{
    SplineList* splinelist;
    /* Check if the node with this index already exists */
    splinelist = appended;
    while( splinelist ){
      if( i == splinelist->i ){
	break;
      }
      splinelist = splinelist->next;
    }
    if( splinelist ){ /* We don't allow for the case where the index already exists*/
      XLALPrintError("Error: Tried to add an already existing index to a SplineList");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input spline, which therefore can't be
                 used anywhere else; this will be acceptable as these operations will only be done
                 when initializing the data */
      splinelist = XLALMalloc( sizeof(SplineList) );
    }
    splinelist->i = i;
    if( spline ){
      splinelist->spline = spline;
    } else {
      splinelist->spline = NULL;
    }
    if( accel ){
      splinelist->accel = accel;
    } else {
      splinelist->accel = NULL;
    }
    if( appended ){
      splinelist->next = appended;
    } else {
        splinelist->next = NULL;
    }
    return splinelist;
}
/* Get the element of a SplineList with a given index */
SplineList* SplineList_GetElement(
	      SplineList* const splinelist,  /* List structure to get element from */
              const UINT4 i ) /* Index looked for */
{
    if( !splinelist ) return NULL;

    SplineList* itr = splinelist;
    while( itr->i != i ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a SplineList */
}
/* Delete list from given pointer to the end of the list */
void SplineList_Destroy( SplineList* splinelist ) /* Head of linked list to destroy */
{
  SplineList* pop;
  while( (pop = splinelist) ){
    if( pop->spline ){ /* Internal spline and accelerator are freed */
      gsl_spline_free( pop->spline );
    }
    if( pop->accel ){
      gsl_interp_accel_free( pop->accel );
    }
    /* Notice that the index i is not freed, like in SphHarmTimeSeries struct indices l and m */
    splinelist = pop->next;
    XLALFree( pop );
  }
}

/***************** Functions for the EOBNRHMROMdata structure ****************/
ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata* data,  /* data to contain */
	   const UINT4 l, /* major mode number */
	   const INT4 m  /* minor mode number */)
{
    ListmodesEOBNRHMROMdata* list;
    /* Check if the node with this mode already exists */
    list = appended;
    while( list ){
      if( l == list->l && m == list->m ){
	break;
      }
      list = list->next;
    }
    if( list ){ /* We don't allow for the case where the mode already exists in the list*/
      XLALPrintError("Error: Tried to add an already existing mode to a ListmodesEOBNRHMROMdata ");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input interpolated data, which therefore can't be
		used anywhere else; this will be acceptable as these operations will only be done
		when interpolating the initialization data */
      list = XLALMalloc( sizeof(ListmodesEOBNRHMROMdata) );
    }
    list->l = l;
    list->m = m;
    if( data ){
      list->data = data;
    } else {
      list->data = NULL;
    }
    if( appended ){
      list->next = appended;
    } else {
        list->next = NULL;
    }
    return list;
}
/* Get the element of a ListmodesEOBNRHMROMdata with a given index */
ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_GetMode(
	   ListmodesEOBNRHMROMdata* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */ )
{
    if( !list ) return NULL;

    ListmodesEOBNRHMROMdata *itr = list;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a ListmodesEOBNRHMROMdata */
}
void ListmodesEOBNRHMROMdata_Destroy(
	   ListmodesEOBNRHMROMdata* list  /* List structure to destroy; notice that the data is destroyed too */
)
{
  ListmodesEOBNRHMROMdata* pop;
  while( (pop = list) ){
    if( pop->data ){ /* Destroying the EOBNRHMROMdata data */
      EOBNRHMROMdata_Cleanup( pop->data );
    }
    /* Notice that the mode indices l and m are not freed, like in SphHarmTimeSeries struct indices l and m */
    list = pop->next;
    XLALFree( pop );
  }
}

/***************** Functions for the EOBNRHMROMdata_interp structure ****************/
ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata_interp* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata_interp* data_interp,  /* data to contain */
	   UINT4 l, /* major mode number */
	   INT4 m  /* minor mode number */)
{
    ListmodesEOBNRHMROMdata_interp* list;
    /* Check if the node with this mode already exists */
    list = appended;
    while( list ){
      if( l == list->l && m == list->m ){
	break;
      }
      list = list->next;
    }
    if( list ){ /* We don't allow for the case where the mode already exists in the list*/
      XLALPrintError("Error: Tried to add an already existing mode to a ListmodesEOBNRHMROMdata_interp ");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input interpolated data, which therefore can't be
		used anywhere else; this will be acceptable as these operations will only be done
		when interpolating the initialization data */
      list = XLALMalloc( sizeof(ListmodesEOBNRHMROMdata_interp) );
    }
    list->l = l;
    list->m = m;
    if( data_interp ){
      list->data_interp = data_interp;
    } else {
      list->data_interp = NULL;
    }
    if( appended ){
      list->next = appended;
    } else {
        list->next = NULL;
    }
    return list;
}
/* Get the element of a ListmodesEOBNRHMROMdata with a given index */
ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_GetMode(
	   ListmodesEOBNRHMROMdata_interp* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */ )
{
    if( !list ) return NULL;

    ListmodesEOBNRHMROMdata_interp *itr = list;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a ListmodesEOBNRHMROMdata_interp */
}
void ListmodesEOBNRHMROMdata_interp_Destroy(
	   ListmodesEOBNRHMROMdata_interp* list  /* List structure to destroy; notice that the data is destroyed too */
)
{
  ListmodesEOBNRHMROMdata_interp* pop;
  while( (pop = list) ){
    if( pop->data_interp ){ /* Destroying the EOBNRHMROMdata_interp data */
      EOBNRHMROMdata_interp_Cleanup( pop->data_interp );
    }
    /* Notice that the mode indices l and m are not freed, like in SphHarmTimeSeries struct indices l and m */
    list = pop->next;
    XLALFree( pop );
  }
}

/***************** Functions for the EOBNRHMROMdata_coeff structure ****************/
ListmodesEOBNRHMROMdata_coeff* ListmodesEOBNRHMROMdata_coeff_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata_coeff* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata_coeff* data_coeff,  /* data to contain */
	   UINT4 l, /* major mode number */
	   INT4 m  /* minor mode number */)
{
    ListmodesEOBNRHMROMdata_coeff* list;
    /* Check if the node with this mode already exists */
    list = appended;
    while( list ){
      if( l == list->l && m == list->m ){
	break;
      }
      list = list->next;
    }
    if( list ){ /* We don't allow for the case where the mode already exists in the list*/
      XLALPrintError("Error: Tried to add an already existing mode to a ListmodesEOBNRHMROMdata_coeff ");
      return(NULL);
    } else { /* In that case, we do NOT COPY the input interpolated data, which therefore can't be
		used anywhere else; this will be acceptable as these operations will only be done
		when interpolating the initialization data */
      list = XLALMalloc( sizeof(ListmodesEOBNRHMROMdata_coeff) );
    }
    list->l = l;
    list->m = m;
    if( data_coeff ){
      list->data_coeff = data_coeff;
    } else {
      list->data_coeff = NULL;
    }
    if( appended ){
      list->next = appended;
    } else {
        list->next = NULL;
    }
    return list;
}
/* Get the element of a ListmodesEOBNRHMROMdata_coeff with a given index */
ListmodesEOBNRHMROMdata_coeff* ListmodesEOBNRHMROMdata_coeff_GetMode(
	   ListmodesEOBNRHMROMdata_coeff* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */ )
{
    if( !list ) return NULL;

    ListmodesEOBNRHMROMdata_coeff *itr = list;
    while( itr->l != l || itr->m != m ){
        itr = itr->next;
        if( !itr ) return NULL;
    }
    return itr; /* The element returned is itself a pointer to a ListmodesEOBNRHMROMdata_coeff */
}
void ListmodesEOBNRHMROMdata_coeff_Destroy(
	   ListmodesEOBNRHMROMdata_coeff* list  /* List structure to destroy; notice that the data is destroyed too */
)
{
  ListmodesEOBNRHMROMdata_coeff* pop;
  while( (pop = list) ){
    if( pop->data_coeff ){ /* Destroying the EOBNRHMROMdata_coeff data */
      EOBNRHMROMdata_coeff_Cleanup( pop->data_coeff );
    }
    /* Notice that the mode indices l and m are not freed, like in SphHarmTimeSeries struct indices l and m */
    list = pop->next;
    XLALFree( pop );
  }
}

/***************** Other structure functions ****************/

/* Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structures */
INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym) {
  /* Deleted the definition of the string 'func': usage ? */
  COMPLEX16 Y;
  UINT4 j;
  COMPLEX16 hlmtildevalue;

  /* Checks LAL_CHECK_VALID_SERIES and LAL_CHECK_CONSISTENT_TIME_SERIES removed
   * - they do not seem available for frequency series ? */

  INT4 minus1l; /* (-1)^l */
  if ( l%2 ) minus1l = -1;
  else minus1l = 1;
  if ( sym ) { /* equatorial symmetry: add in -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
    COMPLEX16 factorp = 1./2*(Y + minus1l*Ymstar);
    COMPLEX16 factorc = I/2*(Y - minus1l*Ymstar);
    COMPLEX16* datap = hptilde->data->data;
    COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      datap[j] += factorp*hlmtildevalue;
      datac[j] += factorc*hlmtildevalue;
    }
  }
  else { /* not adding in the -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 factorp = 1./2*Y;
    COMPLEX16 factorc = I/2*Y;
    COMPLEX16* datap = hptilde->data->data;
    COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      datap[j] += factorp*hlmtildevalue;
      datac[j] += factorc*hlmtildevalue;
    }
  }

  return 0;
}
