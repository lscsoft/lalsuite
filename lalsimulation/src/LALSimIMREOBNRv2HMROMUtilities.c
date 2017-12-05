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
 *   q = 1-11.98
 *   No spin
 *   Mtot >= 8Msun for fstart=10Hz
 *
 */


/*****************************************************************************/
/**************************** General definitions ****************************/

#define nk_amp 10  /* number of SVD-modes == number of basis functions for amplitude */
#define nk_phi 20  /* number of SVD-modes == number of basis functions for phase */

/* Contrarily to SEOBNR, the frequencies used by the SVD for each mode are to be loaded as data in the mode-by-mode loop (like the coefficients and reduced basis) */
/* Define the number of points in frequency used by the SVD, identical for all modes */
#define nbfreq 300
/* Define the number of training waveforms used by the SVD, identical for all modes */
#define nbwf 301


/***************************************************/
/*************** Type definitions ******************/

typedef struct tagSplineList {
    gsl_spline*            spline; /* The gsl spline */
    gsl_interp_accel*      accel; /* The gsl accelerator */
    UINT4                  i; /* Index in the list  */
    struct tagSplineList*  next; /* Pointer to the next list element */
} SplineList;

typedef struct tagEOBNRHMROMdata
{
  gsl_vector* q;
  gsl_vector* freq;
  gsl_matrix* Camp;
  gsl_matrix* Cphi;
  gsl_matrix* Bamp;
  gsl_matrix* Bphi;
  gsl_vector* shifttime;
  gsl_vector* shiftphase;
} EOBNRHMROMdata;

typedef struct tagEOBNRHMROMdata_interp
{
  SplineList* Camp_interp; /* List of splines for the amp coefficients - SplineList, index of reduced basis */
  SplineList* Cphi_interp; /* List of splines for the phase coefficients - SplineList, index of reduced basis */
  SplineList* shifttime_interp; /* interpolated shift in time - SplineList with one element */
  SplineList* shiftphase_interp; /* interpolated shift in phase - SplineList with one element */
} EOBNRHMROMdata_interp;

typedef struct tagEOBNRHMROMdata_coeff
{
  gsl_vector* Camp_coeff;
  gsl_vector* Cphi_coeff;
  double*     shifttime_coeff;
  double*     shiftphase_coeff;
} EOBNRHMROMdata_coeff;

typedef struct tagListmodesEOBNRHMROMdata
{
    EOBNRHMROMdata*                     data; /* The ROM data. */
    UINT4                               l; /* Node mode l  */
    INT4                                m; /* Node submode m  */
    struct tagListmodesEOBNRHMROMdata*  next; /* next pointer */
} ListmodesEOBNRHMROMdata;

typedef struct tagListmodesEOBNRHMROMdata_interp
{
    EOBNRHMROMdata_interp*                     data_interp; /* The splines built from the coefficients. */
    UINT4                                      l; /* Node mode l  */
    INT4                                       m; /* Node submode m  */
    struct tagListmodesEOBNRHMROMdata_interp*  next; /* next pointer */
} ListmodesEOBNRHMROMdata_interp;


/****************************************************************************/
/******************************** Prototypes ********************************/

/* Function to read data from files */
static INT4 Read_Data_Mode(const char dir[], const INT4 mode[2], EOBNRHMROMdata *data);

/* Functions to initialize and cleanup data structures */
static void EOBNRHMROMdata_Init(EOBNRHMROMdata **data);
static void EOBNRHMROMdata_interp_Init(EOBNRHMROMdata_interp **data_interp);
static void EOBNRHMROMdata_coeff_Init(EOBNRHMROMdata_coeff **data_coeff);

static void EOBNRHMROMdata_Cleanup(EOBNRHMROMdata *data);
static void EOBNRHMROMdata_interp_Cleanup(EOBNRHMROMdata_interp *data_interp);
static void EOBNRHMROMdata_coeff_Cleanup(EOBNRHMROMdata_coeff *data_coeff);

/* Function to add modes for frequency-domain structures */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym);

/* Functions associated to list manipulations */
static SplineList* SplineList_AddElementNoCopy(
	   SplineList* appended,  /* List structure to prepend to */
	   gsl_spline* spline,  /* spline to contain */
           gsl_interp_accel* accel,  /* accelerator to contain */
	   UINT4 i /* index in the list */
);
static SplineList* SplineList_GetElement(
	   SplineList* const splinelist,  /* List structure to get a particular mode from */
	   const UINT4 i /* index in the list */
);
static void SplineList_Destroy(
	   SplineList* list  /* List structure to destroy; notice that the content is destroyed too */
);
static ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata* indata,  /* data to contain */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
static ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_GetMode(
	   ListmodesEOBNRHMROMdata* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
/* Note: we do not add a ListmodesEOBNRHMROMdata_Destroy function,
 * as the only ListmodesEOBNRHMROMdata will be persistent and never destroyed until the program ends */
static ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_AddModeNoCopy(
	   ListmodesEOBNRHMROMdata_interp* appended,  /* List structure to prepend to */
	   EOBNRHMROMdata_interp* data,  /* data to contain */
	   UINT4 l, /* major mode number */
	   INT4 m  /* minor mode number */
);
static ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_GetMode(
	   ListmodesEOBNRHMROMdata_interp* const list,  /* List structure to get a particular mode from */
	   UINT4 l, /*< major mode number */
	   INT4 m  /*< minor mode number */
);
/* Note: we do not add a ListmodesEOBNRHMROMdata_interp_Destroy function,
 * as the only ListmodesEOBNRHMROMdata_interp will be persistent and never destroyed until the program ends */


/************************************************************************************/
/******************************** Internal functions ********************************/

/* Read binary ROM data for frequency vectors, coefficients matrices, basis functions matrices, and shiftvectors in time and phase */
static INT4 Read_Data_Mode(const char dir[], const INT4 mode[2], EOBNRHMROMdata *data) {
  /* Load binary data for amplitude and phase spline coefficients as computed in Mathematica */
  INT4 ret = XLAL_SUCCESS;
  size_t size = strlen(dir) + 64;
  char *file_q = XLALMalloc(size);
  char *file_freq = XLALMalloc(size);
  char *file_Camp = XLALMalloc(size);
  char *file_Cphi = XLALMalloc(size);
  char *file_Bamp = XLALMalloc(size);
  char *file_Bphi = XLALMalloc(size);
  char *file_shifttime = XLALMalloc(size);
  char *file_shiftphase = XLALMalloc(size);
  snprintf(file_q, size, "%s", "EOBNRv2HMROM_q.dat"); /* The q vector is the same for all modes */
  snprintf(file_freq, size, "%s%d%d%s", "EOBNRv2HMROM_freq_", mode[0], mode[1], ".dat");
  snprintf(file_Camp, size, "%s%d%d%s", "EOBNRv2HMROM_Camp_", mode[0], mode[1], ".dat");
  snprintf(file_Cphi, size, "%s%d%d%s", "EOBNRv2HMROM_Cphi_", mode[0], mode[1], ".dat");
  snprintf(file_Bamp, size, "%s%d%d%s", "EOBNRv2HMROM_Bamp_", mode[0], mode[1], ".dat");
  snprintf(file_Bphi, size, "%s%d%d%s", "EOBNRv2HMROM_Bphi_", mode[0], mode[1], ".dat");
  snprintf(file_shifttime, size, "%s%d%d%s", "EOBNRv2HMROM_shifttime_", mode[0], mode[1], ".dat");
  snprintf(file_shiftphase, size, "%s%d%d%s", "EOBNRv2HMROM_shiftphase_", mode[0], mode[1], ".dat");
  ret |= read_vector(dir, file_q, data->q);
  ret |= read_vector(dir, file_freq, data->freq);
  ret |= read_matrix(dir, file_Camp, data->Camp);
  ret |= read_matrix(dir, file_Cphi, data->Cphi);
  ret |= read_matrix(dir, file_Bamp, data->Bamp);
  ret |= read_matrix(dir, file_Bphi, data->Bphi);
  ret |= read_vector(dir, file_shifttime, data->shifttime);
  ret |= read_vector(dir, file_shiftphase, data->shiftphase);
  XLALFree(file_q);
  XLALFree(file_freq);
  XLALFree(file_Camp);
  XLALFree(file_Cphi);
  XLALFree(file_Bamp);
  XLALFree(file_Bphi);
  XLALFree(file_shifttime);
  XLALFree(file_shiftphase);
  return(ret);
}


/********************* Functions to initialize and cleanup data structures ********************/
static void EOBNRHMROMdata_Init(EOBNRHMROMdata **data) {
  if(!data) exit(1);
  /* Create storage for structures */
  if(!*data) *data=XLALMalloc(sizeof(EOBNRHMROMdata));
  else
  {
    EOBNRHMROMdata_Cleanup(*data);
  }
  (*data)->q = gsl_vector_alloc(nbwf);
  (*data)->freq = gsl_vector_alloc(nbfreq);
  (*data)->Camp = gsl_matrix_alloc(nk_amp,nbwf);
  (*data)->Cphi = gsl_matrix_alloc(nk_phi,nbwf);
  (*data)->Bamp = gsl_matrix_alloc(nbfreq,nk_amp);
  (*data)->Bphi = gsl_matrix_alloc(nbfreq,nk_phi);
  (*data)->shifttime = gsl_vector_alloc(nbwf);
  (*data)->shiftphase = gsl_vector_alloc(nbwf);
}
static void EOBNRHMROMdata_interp_Init(EOBNRHMROMdata_interp **data_interp) {
  if(!data_interp) exit(1);
  /* Create storage for structures */
  if(!*data_interp) *data_interp=XLALMalloc(sizeof(EOBNRHMROMdata_interp));
  else
  {
    EOBNRHMROMdata_interp_Cleanup(*data_interp);
  }
  (*data_interp)->Camp_interp = NULL;
  (*data_interp)->Cphi_interp = NULL;
  (*data_interp)->shifttime_interp = NULL;
  (*data_interp)->shiftphase_interp = NULL;
}
static void EOBNRHMROMdata_coeff_Init(EOBNRHMROMdata_coeff **data_coeff) {
  if(!data_coeff) exit(1);
  /* Create storage for structures */
  if(!*data_coeff) *data_coeff=XLALMalloc(sizeof(EOBNRHMROMdata_coeff));
  else
  {
    EOBNRHMROMdata_coeff_Cleanup(*data_coeff);
  }
  (*data_coeff)->Camp_coeff = gsl_vector_alloc(nk_amp);
  (*data_coeff)->Cphi_coeff = gsl_vector_alloc(nk_phi);
  (*data_coeff)->shifttime_coeff = XLALMalloc(sizeof(double));
  (*data_coeff)->shiftphase_coeff = XLALMalloc(sizeof(double));
}
static void EOBNRHMROMdata_Cleanup(EOBNRHMROMdata *data /* data to destroy */) {
  if(data->q) gsl_vector_free(data->q);
  if(data->freq) gsl_vector_free(data->freq);
  if(data->Camp) gsl_matrix_free(data->Camp);
  if(data->Cphi) gsl_matrix_free(data->Cphi);
  if(data->Bamp) gsl_matrix_free(data->Bamp);
  if(data->Bphi) gsl_matrix_free(data->Bphi);
  if(data->shifttime) gsl_vector_free(data->shifttime);
  if(data->shiftphase) gsl_vector_free(data->shiftphase);
  XLALFree(data);
}
static void EOBNRHMROMdata_coeff_Cleanup(EOBNRHMROMdata_coeff *data_coeff) {
  if(data_coeff->Camp_coeff) gsl_vector_free(data_coeff->Camp_coeff);
  if(data_coeff->Cphi_coeff) gsl_vector_free(data_coeff->Cphi_coeff);
  if(data_coeff->shifttime_coeff) free(data_coeff->shifttime_coeff);
  if(data_coeff->shiftphase_coeff) free(data_coeff->shiftphase_coeff);
  XLALFree(data_coeff);
}
static void EOBNRHMROMdata_interp_Cleanup(EOBNRHMROMdata_interp *data_interp) {
  if(data_interp->Camp_interp) SplineList_Destroy(data_interp->Camp_interp);
  if(data_interp->Cphi_interp) SplineList_Destroy(data_interp->Cphi_interp);
  if(data_interp->shifttime_interp) SplineList_Destroy(data_interp->shifttime_interp);
  if(data_interp->shiftphase_interp) SplineList_Destroy(data_interp->shiftphase_interp);
  XLALFree(data_interp);
}


/********************* Function to add modes for frequency-domain structures ********************/

/* Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structures */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym) {
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

/********************* Functions for list structures ********************/

/***************** Functions for the SplineList structure ****************/

/* Prepend a node to a linked list of splines, or create a new head */
static SplineList* SplineList_AddElementNoCopy(
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
static SplineList* SplineList_GetElement(
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
static void SplineList_Destroy( SplineList* splinelist ) /* Head of linked list to destroy */
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

/***************** Functions for the ListmodesEOBNRHMROMdata structure ****************/
static ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_AddModeNoCopy(
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
static ListmodesEOBNRHMROMdata* ListmodesEOBNRHMROMdata_GetMode(
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
/* Note: we do not add a ListmodesEOBNRHMROMdata_Destroy function,
 * as the only ListmodesEOBNRHMROMdata will be persistent and never destroyed until the program ends */


/***************** Functions for the ListmodesEOBNRHMROMdata_interp structure ****************/
static ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_AddModeNoCopy(
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
static ListmodesEOBNRHMROMdata_interp* ListmodesEOBNRHMROMdata_interp_GetMode(
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
/* Note: we do not add a ListmodesEOBNRHMROMdata_interp_Destroy function,
 * as the only ListmodesEOBNRHMROMdata_interp will be persistent and never destroyed until the program ends */
