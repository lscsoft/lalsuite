/*----------------------------------------------------------------------- 
 * 
 * File Name: TrackSearch.c 
 * 
 * Author: R. Balasubramanian
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * TrackSearch.c
 * 
 * SYNOPSIS 
 * (void) LALSignalTrackSearch()
 * 
 * DESCRIPTION 
 * This file contains the routines to detect curves on the time frequency plane. 
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 * TS_MSGNON_NULL_POINTER  "Allocation of memory to Non null pointer" 
 * TS_MSG_ALLOC            "flag for allocating storage space does not contain a valid value (0,1,2)"
 * TS_MSG_UNINITIALIZED_MASKS "Gaussian masks appear not to have been initialized"
 * TS_MSG_MAP_DIMENSION     "Dimensions of Map have to be each greater than unity"
 * TS_MSG_ILLEGAL_PARAMS   "Illegal parameter values ( high and low should be positive and high>low)"
 * TS_MSG_LINE_START "The number of line start points is not consistent with previous module" 
 * MSG_TS_ARRAY_OVERFLOW  "Array bounds can be crossed "
 * MSG_TS_TOO_MANY_CURVES " Too many curves found in map "
 * MSG_TS_OCTANT "octant cannot be greater than 4 "
 * MSG_TS_MASKSIZE "maskSize too small"
 
 * CALLS
 * (list of LLAL, LDAS, other non-system functions/procedures called. 
 * NONE
 *
 * NOTES
 * (Other notes)
 * 
 *-----------------------------------------------------------------------
 */


#include <lal/TrackSearch.h>

NRCSID (TRACKSEARCHC, "$Id$");


/* Local Structure definitions */
typedef struct tagLinePoints {
  INT4 pValue; /* value of the Line Point: 2 if a line Start Point else 1 if an ordinary line point
		  A line start point is not necessarily a point where a curve begins but just the point at 
		  which we begin to search for a curve*/
  INT4 r; /* row of the line point*/
  INT4 c; /* column of the line Point */
  REAL4 angle; /* angle at the line Point */
  REAL4 pRow; /* the subpixel position of the line */
  REAL4 pCol; /* the subpixel position of the line */
  REAL4 eigen; /* the principal eigen value of the Hessian matrix  */
  CHAR flag; /* to flag the line start points if they are already processed */
}
LinePoints;

/* A structure to use in searching for line points in the neighbourhood */
typedef struct tagOffset {
  INT4 x;
  INT4 y;
} Offset;

/*
 * Local Function prototypes; not visible outside;
 */
static void InitializeStorage(LALStatus *, TrackSearchStore *, TrackSearchParams *);
static void DestroyStorage(LALStatus *, TrackSearchStore *);
static void ComputeConvolutions(LALStatus *,  TrackSearchStore *, const TimeFreqRep *, TrackSearchParams *);
static void ComputeLinePoints(LALStatus *, TrackSearchStore * , TrackSearchParams *);
static void ConnectLinePoints(LALStatus *, TrackSearchOut *, TrackSearchParams *, const TimeFreqRep *);
static REAL4 GetAngle(REAL4 , REAL4 );
static REAL4 Gauss0(INT4,REAL4);
static REAL4 Gauss1(INT4,REAL4);
static REAL4 Gauss2(INT4,REAL4);
static INT4 CompareLinePoints(const void *,const void *);
/*
  Main Routine to detect curves. 
  Convolves the image with the derivatives of Gaussian kernels, identifies line points
  and then connects them into curves
*/
void 
LALSignalTrackSearch(LALStatus *status,
                  TrackSearchOut *out,
		     const TimeFreqRep *tfMap, /* type defined in TimeFreq.h */
                  TrackSearchParams *params)
{
    
  /* Initialize status structure   */

  INITSTATUS(status,"LALSignalTrackSearch",TRACKSEARCHC);
  ATTATCHSTATUSPTR (status);
  
  /* Check the the arguments are not null pointers */
  ASSERT(tfMap     != NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  ASSERT(tfMap->map!= NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);  
  ASSERT(params    != NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);  
  ASSERT(out       != NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  
  /* check if the dimensions of the input image are ok */
  ASSERT(((params->width > 8)&&(params->height>8)),status, TS_MAP_DIMENSION,TS_MSG_MAP_DIMENSION);
        /* check if the alloc flag has only allowed values */
  /* this flag is used to determine whether memory is to be allocated 
     for the temporary storage space. (allowed values are 0 for no allocation, 1 for allocation
     and 2 for deallocation*/
  ASSERT(((params->allocFlag>=0)&&(params->allocFlag<=2)), status, TS_ALLOC, TS_MSG_ALLOC);
  /* If the alloc flag==2 then free storage space and return */
  if(params->allocFlag==2){
    DestroyStorage(status->statusPtr, &(out->store));
    CHECKSTATUSPTR(status);
    DETATCHSTATUSPTR(status);
    RETURN(status);    
  }
  /* if the allocFlag==1 then allocate memory and **SET*** allocFlag to 0; in case the caller is forgetful*/
  if(params->allocFlag==1){
    InitializeStorage(status->statusPtr, &(out->store), params);
    CHECKSTATUSPTR(status);
    params->allocFlag=0;
  }
  /* check if the values for the threshold on the second derivative are OK */
  ASSERT(((params->low>0)&&(params->high>params->low)), status, TS_ILLEGAL_PARAMS, TS_MSG_ILLEGAL_PARAMS);
  /* check if sigma is greater than unity. Otherwise the Gaussian masks will be too small to be useful. */
  ASSERT(params->sigma>=1,status,TS_ILLEGAL_PARAMS, TS_MSG_ILLEGAL_PARAMS);    
  /* Now we are ready to begin our work */
  /* Compute Convolutions of the image with the Gaussian derivative Kernels */
  ComputeConvolutions(status->statusPtr, &(out->store), tfMap, params);
  CHECKSTATUSPTR(status);
  /* Compute Possible Line Points */
  ComputeLinePoints(status->statusPtr,  &(out->store), params);
  CHECKSTATUSPTR(status);
  /* initialize the number of curves found to be zero */
  out->numberOfCurves=0;
  /* Connect the Points into curves  only if there are possible line starting points and 
     the number of line points is more than a given threshold */  
  if((out->store).numLStartPoints&&((out->store).numLPoints>=LENGTH_THRESHOLD)){
    ConnectLinePoints(status->statusPtr, out, params, tfMap); /**/
    CHECKSTATUSPTR(status);
  }
  /* Return to Calling program */
  DETATCHSTATUSPTR(status);
    RETURN(status);
}

/* Functions not accessible to the user are defined below */
/* Routine to allocate and initialise the Gaussian mask arrays */

static void 
InitializeStorage( LALStatus * status, 
                   TrackSearchStore *store,
                   TrackSearchParams *params
                   )
{
  INT4 i,j; /* temporary loop variables */
  INT4 maskSize; /* size of the mask arrays (gaussMask.k[i]) */
  
  /* Initialize status structure   */
  INITSTATUS(status,"InitializeStorage",TRACKSEARCHC);
  
  /* check for null pointers */
  ASSERT (store !=NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  ASSERT (params!=NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  
  /* height and width of the imap */
  store->height=params->height;
  store->width=params->width;
  
  /* allocate space for isLine, a char Map to flag possible line points*/
  store->isLine = (CHAR **)LALMalloc(store->height*sizeof(CHAR *));
  for(i=0;i<store->height;i++) 
    *(store->isLine + i ) = (CHAR *) LALMalloc(store->width*sizeof(CHAR));
  
  /* Arrays to store convolved images.. 5 in number for the 2 first derivatives
     and 3 second derivatives */
  for(i=0;i<5;i++){
    store->k[i] = (REAL4 **)LALMalloc(store->height*sizeof(REAL4 *));
    for(j=0;j<store->height;j++) 
      *(store->k[i]+j) = (REAL4 *)  LALMalloc(store->width*sizeof(REAL4));    
  }
  /* Array to store the Eigen vectors and line positions */
  store->eigenVec = LALMalloc(store->height * store->width * 4 * sizeof(REAL4));
  
  /* Array to store temporary images */
  store->imageTemp = (REAL4 **)LALMalloc(store->height*sizeof(REAL4 *));
  for(i=0;i<store->height;i++) 
    *(store->imageTemp + i ) = (REAL4 *) LALMalloc(store->width*sizeof(REAL4));
   
  /* Allocate memory for the Gaussian masks and compute the size of the arrays */
  maskSize = ceil(params->sigma * MASK_SIZE);
  ASSERT(maskSize>2,status,TS_MASKSIZE,MSG_TS_MASKSIZE);
  /* 3 arrays; gaussian function, its first and second derivative */
  for(i=0;i<3;i++)
    store->gaussMask[i] = LALMalloc(sizeof(REAL8)*(2*maskSize+1));
  for(i=-maskSize;i<=maskSize;i++){
    j = i + maskSize;
    /* Averaged Gaussian function */
    *(store->gaussMask[0]+j) = Gauss0(i,params->sigma);
    /* Averaged first derivative of Gaussian function */
    *(store->gaussMask[1]+j) = -1.0 * Gauss1(i,params->sigma);
          /* Averaged second derivative of Gaussian function */
    *(store->gaussMask[2]+j) = Gauss2(i,params->sigma);
  }  
  RETURN(status);
}


static void 
DestroyStorage (LALStatus *status,
                TrackSearchStore *store
                )
{
  int i,j; /* temporary loop variables */
  
  
  /* Initialize status structure   */
    INITSTATUS(status,"DestroyStorage",TRACKSEARCHC);
  
  /* check for null pointers */
  ASSERT (store !=NULL, status, TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  /* Free all pointers */
  for(i=0;i<store->height;i++)
    LALFree(*(store->isLine + i));
  LALFree(store->isLine);
    
  for(i=0;i<=4;i++){
    for(j=0;j<store->height;j++){
      LALFree(*(store->k[i] + j));
    }
    LALFree(store->k[i]);
  }
  
  LALFree(store->eigenVec);
  
  for(i=0;i<store->height;i++)  
    LALFree(*(store->imageTemp + i));
  LALFree(store->imageTemp);  
  
  for(i=0;i<=2;i++)
    LALFree(*(store->gaussMask+i));  
  /* set the storage structure to 0 */
  memset((void *)store,sizeof(*store),0);
  RETURN(status);
}


static void 
ComputeConvolutions (LALStatus *status, 
                     TrackSearchStore *store, 
                     const TimeFreqRep *tfMap, 
                     TrackSearchParams *params)
{
  REAL8 sum;                /* temporary variable to perform the convolution sum */
  INT4 maskSize;            /* One sided length of the Gaussian mask arrays */
  REAL8 *rowMask=0, *colMask=0; /* Pointers to the row and column masking arrays*/
  INT4 i,j,k;               /* temporary loop variables*/
  REAL4 **convolvedMatrix=0;  /* A pointer used to point at one of the 5 smoothed derivative images k[i]*/
  INT4 derivType;           /* Is one of the 5 possible 1st and 2nd directional derivatives 
			       along the row and column directions*/
  REAL4 check;              /* a variable used to check if the Gaussian mask arrays have been Allocated */
  INT4 heightTemp, widthTemp;/* variables used for reflecting images at its edges*/
  
  /* Initialize status structure   */ 
  INITSTATUS(status,"ComputeConvolutions",TRACKSEARCHC);
  
  /* check for null pointers */
  ASSERT(store->k[0],status,TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  ASSERT(tfMap, status,TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  ASSERT(params,status,TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  ASSERT(store->gaussMask[0],status,TS_NULL_POINTER, TS_MSG_NULL_POINTER);
  
  /* check whether the Gaussian mask arrays have been set */
  check = (*(store->gaussMask[0] + (INT4)(MASK_SIZE*params->sigma)));
  ASSERT(check, status,TS_UNITIALIZED_MASKS, TS_MSG_UNINITIALIZED_MASKS); 
  
  /* Computing Derivatives
     store->k[0] 1st derivative along column direction          (dI/df)
     store->k[1] 1st derivative along row direction             (dI/dt)
     store->k[2] 2nd derivative along column direction          (d^2 I/df^2)
     store->k[3] 2nd derivative along row direction             (d^2 I/dt^2) 
     store->k[4] 2nd derivative along row/column direction      (d^2 I/dtdf)
     */
  
  maskSize=ceil(MASK_SIZE*params->sigma);
  for(derivType=0;derivType<=4;derivType++){
    switch(derivType){
    case 0:
      colMask=store->gaussMask[0] + maskSize; /* offset the pointer to point to the central value*/
      rowMask=store->gaussMask[1] + maskSize; /* offset the pointer to point to the central value*/
      convolvedMatrix = store->k[0];
              break;
    case 1:
              colMask=store->gaussMask[1] + maskSize; /* offset the pointer to point to the central value*/
      rowMask=store->gaussMask[0] + maskSize; /* offset the pointer to point to the central value*/
      convolvedMatrix = store->k[1];
      break;      
    case 2:
      colMask=store->gaussMask[0] + maskSize; /* offset the pointer to point to the central value*/
      rowMask=store->gaussMask[2] + maskSize; /* offset the pointer to point to the central value*/
      convolvedMatrix = store->k[2];
      break;      
    case 3:
      colMask=store->gaussMask[2] + maskSize; /* offset the pointer to point to the central value*/
      rowMask=store->gaussMask[0] + maskSize; /* offset the pointer to point to the central value*/
      convolvedMatrix = store->k[3];
      break;      
    case 4:
      colMask=store->gaussMask[1] + maskSize; /* offset the pointer to point to the central value*/
      rowMask=store->gaussMask[1] + maskSize; /* offset the pointer to point to the central value*/
      convolvedMatrix = store->k[4];
      break;      
    default:/* this would never happen */
      break;
    }
    /* variables used for reflection at edges */
    widthTemp = 2*(params->width - 1);
    heightTemp = 2 * (params->height - 1);
    /* Column convolution */
    for(i=0;i<params->height;i++){
      /* take care of left hand edge reflections */
      for(j=0;j<maskSize;j++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++){
	  if(k>=-j)
	    sum += (*(tfMap->map[i]+j + k)) * colMask[k];
	  else
	    sum += (*(tfMap->map[i]-j - k)) * colMask[k];
	}
	*(store->imageTemp[i] + j) = (REAL4)sum;
      }
      /* the central region */
      for(j=maskSize;j<params->width-maskSize;j++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++)
	  sum += (*(tfMap->map[i]+j + k)) * colMask[k];
	*(store->imageTemp[i] + j) = (REAL4)sum;
      }      
      /* take care of reflection at right hand edge */
      for(j=params->width-maskSize;j<params->width;j++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++){
	  if(k<params->width-j)
	    sum += (*(tfMap->map[i]+j + k)) * colMask[k];
	  else
	    sum += (*(tfMap->map[i]+ widthTemp - j - k)) * colMask[k];
	}
	*(store->imageTemp[i] + j) = (REAL4)sum;
      }
    }
    /* Row Convolution */
    for(j=0;j<params->width;j++){
      /* take care of reflection at lower edge */
      for(i=0;i<maskSize;i++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++){
	  if(k>=-i)
	    sum += (*(store->imageTemp[i+k]+j)) * rowMask[k];
	  else
	    sum += (*(store->imageTemp[-i-k]+j)) * rowMask[k];
	}
	*(convolvedMatrix[i] + j) = (REAL4)sum;
      }      
      /* the central region */
      for(i=maskSize;i<params->height-maskSize;i++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++)
	  sum += (*(store->imageTemp[i+k]+j)) * rowMask[k];
	*(convolvedMatrix[i] + j) = (REAL4)sum;
      }
      /* take care of reflection at upper edge */
      for(i=params->height-maskSize;i<params->height;i++){
	sum=0.0;
	for(k=-maskSize;k<=maskSize;k++){
	  if(k<params->height-i)
	    sum += (*(store->imageTemp[i+k]+j)) * rowMask[k];
	  else
	    sum += (*(store->imageTemp[heightTemp - i - k]+j)) * rowMask[k];
	}
	*(convolvedMatrix[i] + j) = (REAL4)sum;
      }      
    }    
  }
  RETURN(status);
}



static void
ComputeLinePoints (LALStatus *status,
                   TrackSearchStore *store, 
                   TrackSearchParams *params)
{ 
  INT4 i,j; /* temporary loop variables */
  REAL8 cotTheta, cosTheta, sinTheta, temp, eVal1, eVal2; /* variables used in eigen value computation */
  REAL8 k0,k1,k2,k3,k4;/* variables to contain the 5 first and second directional derivatives */
  REAL8 a,b,t; /* temporary variables */ 
  REAL8 px,py;  /* the subpixel line point positions */
  REAL4 *eigenVec; /* a pointer to the array containing the eigen vector and subpixel position */  
  /*FILE* fp;*/ /* temp added by cwt */

  /* Initialize status structure   */ 
  INITSTATUS(status,"ComputeLinePoints",TRACKSEARCHC);
  /* initialize the number of possible line start points and line points */
  store->numLStartPoints=0;
  store->numLPoints=0;
  /* Initialize the eigenVec array */
  eigenVec = store->eigenVec;
  memset(eigenVec,0,params->height*params->width*sizeof(REAL4)*4);
  for(i=0;i<params->height;i++){
    for(j=0;j<params->width;j++){
      /* The five first the second directional derivatives at each point */
      k0 = (REAL8)*((store->k[0])[i]+j); /* first along row */
      k1 = (REAL8)*((store->k[1])[i]+j); /* first along column */
            k2 = (REAL8)*((store->k[2])[i]+j); /* second along row */
      k3 = (REAL8)*((store->k[3])[i]+j); /* second along column */
      k4 = (REAL8)*((store->k[4])[i]+j); /* second along row/column */
      eigenVec = store->eigenVec + 4 * (i * params->width + j);
      
      /* check for a already diagonalised Hessian matrix*/
      if(k4){
	cotTheta=0.5*(k3-k2)/(k4);
	temp = 1.0/(fabs(cotTheta) + sqrt(cotTheta*cotTheta+1));
	if(cotTheta<0) temp=-temp;
	/* computing direction */
	cosTheta=1.0/sqrt(temp*temp + 1.0);
	sinTheta=-temp*cosTheta;
	/* computing the eigen values */
	eVal1 = k2 - temp * k4;
	eVal2 = k3 + temp * k4;
      }
      else{
	/* computing the eigen values */
	eVal1 = k2;
	eVal2 = k3;
	/* computing direction */
	cosTheta=1.0;
	sinTheta=0.0;
      }
      /* Select the Eigen Value with maximum absolute value first; if 
	 absolute values are equal then put the one with negative value first */
      if((fabs(eVal1)>fabs(eVal2)) || ((fabs(eVal1)==fabs(eVal2))&&(eVal1<eVal2))){
	*((store->imageTemp)[i] + j) = (REAL4)eVal1;
	eigenVec[0] = cosTheta;
	eigenVec[1] = sinTheta;
      }
      else{
	*((store->imageTemp)[i] + j) = (REAL4)eVal2;
	eigenVec[0] = -sinTheta;
	eigenVec[1] = cosTheta;	
      }
      /* compute the line positions */
      if(*((store->imageTemp)[i] + j)<0.0){
	/* a and b are second order and first order increments of the image function respectively 
	   along the eigen direction */
	a = k2 * eigenVec[0] * eigenVec[0] + 2.*k4*eigenVec[0]*eigenVec[1] + k3*eigenVec[1]*eigenVec[1];
	b = k0*eigenVec[0] + k1*eigenVec[1];
	if(a!=0){
	  t = -b/a;
	  /* position of maximum relative to current position */
	  px = t*eigenVec[0];
	  py = t*eigenVec[1];	
                        /* check if the maximum of the function lies within this pixel */
	  if((fabs(px) <= PIXEL_BOUNDARY) && (fabs(py) <= PIXEL_BOUNDARY)){
	    if(fabs(*((store->imageTemp)[i] + j))>=params->low){
	      store->numLPoints++; /* variable stores the number of possible line points */
	      /* possible line starting point */
	      if(fabs(*((store->imageTemp)[i] + j))>=params->high){
		store->numLStartPoints++; /* variable stores the number of line Start points */
		*((store->isLine)[i]+j) = 2;
	      }
                            else
			      *((store->isLine)[i]+j) = 1;/* possible line Point */
	      /* store line point position */
	      eigenVec[2]=i+px;
	      eigenVec[3]=j+py;
	    }
	    else{
	      *((store->isLine)[i]+j) = 0;
	    }
	  }
	  else{
	    *((store->isLine)[i]+j) = 0;
	  }
	}
	else{
	  *((store->isLine)[i]+j) = 0;
	}	
      }
      else{
	*((store->isLine)[i]+j) = 0;
      }
    }
  }  
}

static const Offset offset[9]={{0,0},{1,0},{1,1},{0,1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1}}; /*neighbouring pixels*/

static void 
ConnectLinePoints(LALStatus *status,
                  TrackSearchOut *out,
                  TrackSearchParams *params,
		  const TimeFreqRep *TimeFreqMap)
{
  LinePoints *linePoints; /* An array containing the possible line Points */ 
  TrackSearchStore *store; /* a pointer defined for convenience, points to the temporary storage space */
  Curve contour[2];  /* 2 Contours, one for each direction; will be linked to a single contour */ 
  CHAR lineFlag[3]; /* used to flag valid curve extension points */
  CHAR flag; /* used to break out of the loop searching for more points for the curve */
  CHAR directionSwitch; /* used to switch the direction of the line point search if the value of ocant is discontinuous*/
  INT4 **label;/*map whose entries point to the approapriate element of the linePoints array*/
  INT4 **mapContour; /* Each Contour is assigned an iteger value and the points of this map depict the contours*/
  INT4 i,j; /* Temporary loop variables */
  INT4 myindex,nextIndex; /* used as an indices to the lineStart array */
  INT4 processed; /*the number of line Points processed so far */
  INT4 octant,lastOctant; /* to use in identifying which neighbouring pixels to search for line points */
  INT4 currentRow,currentCol,nextRow,nextCol; /* variables pointing to the current and next location*/
  INT4 curveLength;    /* the length of the current curve segment */
  INT4 maxCurveLen=MAX_CURVE_LENGTH;/* the maximum curve lenght allowed */
  INT4 direction; /* the 2 opposite directions to traverse to obtain the complete line */
  INT4 reject[2],check[3],which;/* the list of surrounding points to reject and select */
  REAL4 *eigenVec; /* a pointer defined for convenience; points to a particular element of an array */
  REAL4 metric[3]; /* a metric defined to choose between 3 different line points */
  REAL4 minimum; /* the minimum of the metric[i]*/
  REAL4 differ; /* the difference between 2 angles */
  REAL4 powerHalfContourA; /* The resulting power of half the contour from tfMap*/
  REAL4 powerHalfContourB; /* The resulting power of half the contour from tfMap*/
  
  
  /* Initialize status structure   */ 
  INITSTATUS(status,"ConnectLinePoints",TRACKSEARCHC);  
  /* check if the caller has passed a non NULL pointer for Curves */
  ASSERT(out->curves==NULL,status,TS_NON_NULL_POINTER,TS_MSGNON_NULL_POINTER); 
  /* a pointer to the storage structure */
  store = &(out->store);
  /* Allocate arrays */
  linePoints = LALMalloc(sizeof(LinePoints)*(store->numLPoints));
  label= LALMalloc(sizeof(INT4 *) * params->height);
  mapContour= LALMalloc(sizeof(INT4 *) * params->height);
  for(i=0;i<params->height;i++){
    label[i]=LALMalloc(sizeof(INT4)*params->width);
    mapContour[i]=LALMalloc(sizeof(INT4)*params->width);
    memset(label[i],0,sizeof(INT4)*params->width);
    memset(mapContour[i],0,sizeof(INT4)*params->width);
  }
  /* initialize the contours;*/
  contour[0].row = LALMalloc(sizeof(INT4) * MAX_CURVE_LENGTH);
  contour[0].col = LALMalloc(sizeof(INT4) * MAX_CURVE_LENGTH);
  contour[1].row = LALMalloc(sizeof(INT4) * MAX_CURVE_LENGTH);
  contour[1].col = LALMalloc(sizeof(INT4) * MAX_CURVE_LENGTH);
  /* we do not know whether the curve will have a junction so we initialize to zero(no junction)*/
  contour[0].junction=contour[1].junction=0;
  /* Fill up linePoints structure and set label to 0 */
  myindex=0;
  for(i=0;i<params->height;i++){
    for(j=0;j<params->width;j++){
      label[i][j]=0;
      if(store->isLine[i][j]>0){
	linePoints[myindex].pValue = store->isLine[i][j];
	linePoints[myindex].r = i;
	linePoints[myindex].c = j;
	linePoints[myindex].flag = 1;
	linePoints[myindex].eigen = store->imageTemp[i][j];
	eigenVec = store->eigenVec + 4*(params->width*i + j);	
	/* calculate the angles at all the line points (between 0 and PI)*/
	linePoints[myindex].angle=GetAngle(eigenVec[1],eigenVec[0]);	
	/* store the subpixel line position */
	linePoints[myindex].pRow=eigenVec[2];
	linePoints[myindex].pCol=eigenVec[3];	
	myindex++;
      }	
    }
  }
  /* check if the number of line start points is the same as found in the previous module */
  ASSERT(myindex==store->numLPoints,status,TS_LINE_START,TS_MSG_LINE_START);
  /* sort the line point structures based on their values (2's come first ) 
     i.e. The first store->numLStartPoints points are line starting points and the remaining are ordinary points
     Also among the line start points the sorting order is one of decreasing value of linePoint.eigen*/
  qsort (linePoints,store->numLPoints,sizeof(LinePoints), CompareLinePoints);
  /* For easy access initialize map contour */
  for(i=0;i<store->numLPoints;i++)
    label[linePoints[i].r][linePoints[i].c] = i + 1;
  /* loop while there remain unprocessed line Start Points */
  for(processed=0;processed<store->numLStartPoints;processed++){
    /* skip this point if it has already been processed */
    if(!linePoints[processed].flag)
      continue;
    /* set the flag to 0 so that it is not processed again */
    linePoints[processed].flag=0;
    mapContour[linePoints[processed].r][linePoints[processed].c]=out->numberOfCurves+1;
    /* go of in the 2 directions for every line start point */
    for(direction=0;direction<=1;direction++){
      /* Initialise the contour */
      curveLength=0;
      currentRow=linePoints[processed].r;
      currentCol=linePoints[processed].c;
      contour[direction].row[curveLength]=linePoints[processed].r;
      contour[direction].col[curveLength]=linePoints[processed].c;
      curveLength++;
      flag=1;
      /* identify which direction the line is normal to */
      lastOctant=floor((double)linePoints[processed].angle/(LAL_PI/4.0)) + 1;
      /* check if octant is less than or equal to 4 */
      ASSERT(lastOctant<=4,status,TS_OCTANT,MSG_TS_OCTANT);      
      /*initialize directionSwitch to 1 */
      directionSwitch=1;
      /* loop while there are points to be joined to the current curve */
      while(flag){      
	myindex=label[currentRow][currentCol]-1;
	/* identify which direction the line is normal to */
	octant = floor((double)linePoints[myindex].angle/(LAL_PI/4.0)) + 1;
	/* check if octant is less than or equal to 4 */
	ASSERT(octant<=4,status,TS_OCTANT,MSG_TS_OCTANT);
	/* now we have to find whether the octant has changed discontinuously for 1 to 4 or from 4 to 1 
	   This is to avoid infinite loops; if it has done so then  the idea is to switch the direction*/
	if(((octant==1)&&(lastOctant==4))||((octant==4)&&(lastOctant==1))){
	  if(directionSwitch)
	    directionSwitch=0;
	  else
	    directionSwitch=1;
	}	
	/* set the last octant to octant */
	lastOctant=octant;	
	/* identify the surrounding pixels to reject i.e. reject multiple responses */ 
	reject[0]=octant;
	reject[1]=octant+4;	
	/* select the points to search for continuation */
	if(!direction){
	  if(directionSwitch){
	    check[0]=octant+1;
	    check[1]=octant+2;
	    check[2]=octant+3;
	  }
	  else{
	    check[0]=octant+5;
	    if(check[0]>8) check[0]-=8;
	    check[1]=octant+6;
	    if(check[1]>8) check[1]-=8;
	    check[2]=octant+7;
	    if(check[2]>8) check[2]-=8;	    
	  }
	}
	else{
	  if(directionSwitch){
	    check[0]=octant+5;
	    if(check[0]>8) check[0]-=8;
	    check[1]=octant+6;
	    if(check[1]>8) check[1]-=8;
	    check[2]=octant+7;
	    if(check[2]>8) check[2]-=8;	    	    
	  }
	  else{
	    check[0]=octant+1;
	    check[1]=octant+2;
	    check[2]=octant+3;	    	    
	  }
	}
	/* loop through rejected pixels */
	for(i=0;i<2;i++){
	  nextRow=currentRow+offset[reject[i]].x;
	  nextCol=currentCol+offset[reject[i]].y;
	  /* if pixel position is outside the map then ignore this pixel */
	  if((nextCol>=params->width)||(nextCol<0)||(nextRow>=params->height)||(nextRow<0))
	    continue;
	  /* if this pixel is not marked as a line point ignore */
	  if(!label[nextRow][nextCol])
	    continue;
	  /* check whether the pixel is a multiple response*/
	  nextIndex=label[nextRow][nextCol]-1;
	  differ=(REAL4)fabs((double)linePoints[myindex].angle - linePoints[nextIndex].angle);
	  if(differ>LAL_PI/2.0)
	    differ=LAL_PI - differ;
	  if(differ<MAX_ANGLE_DIFFERENCE){	
	    linePoints[nextIndex].flag=0;
	    /* check if one of the rejected points is the original starting point and if it is not then we
	       can safely put the entry in the label equal to zero */
	    if(linePoints[nextIndex].pValue==2)
	      if(!((nextRow==linePoints[processed].r)&&(nextCol==linePoints[processed].c)))
		label[nextRow][nextCol]=0;	
	  }
	}
	/* now go through the possible neighbours */
	for(i=0;i<3;i++){
	  nextRow=currentRow+offset[check[i]].x;
	  nextCol=currentCol+offset[check[i]].y;
	  /* if pixel position is outside the map then ignore this pixel */
	  if((nextCol>=params->width)||(nextCol<0)||(nextRow>=params->height)||(nextRow<0)){
	    lineFlag[i]=0;
	    continue;
	  }
	  /* if this pixel is not marked as a line point ignore */
	  if(!label[nextRow][nextCol]){
	    lineFlag[i]=0;
	    continue;	  
	  }
	  /* check if the pixel is a valid continuation */
	  nextIndex=label[nextRow][nextCol]-1;
	  differ=(REAL4)fabs((double)linePoints[myindex].angle-linePoints[nextIndex].angle);
	  if(differ>LAL_PI/2.0)
	    differ=LAL_PI - differ;
	  if(differ>MAX_ANGLE_DIFFERENCE){
	    lineFlag[i]=0;
	  }
	  else {
	    lineFlag[i]=1;
	    /* the point already belongs to a curve ignore */
	    if(mapContour[nextRow][nextCol]){
	      if(mapContour[nextRow][nextCol]!=out->numberOfCurves+1){
		which=mapContour[nextRow][nextCol]-1;
		out->curves[which].junction=1;
		contour[direction].junction=1;
	      }
	      else{
		contour[direction].junction=1;
	      }
                            lineFlag[i]=0;
	      continue;	  	    
                        }    
	    /* the metric used to select the most plausible adjacent point
	       d = sqrt((px1-px2)^2 + (py1-py2)^2) + abs(angle1 - angle2)*/
	    metric[i]=sqrt((linePoints[myindex].pRow-linePoints[nextIndex].pRow)*(linePoints[myindex].pRow-linePoints[nextIndex].pRow)+\
			   (linePoints[myindex].pCol-linePoints[nextIndex].pCol)*(linePoints[myindex].pCol-linePoints[nextIndex].pCol)) + differ;
	  }
	}
	/* if none of  the pixels are not a valid continuation, the curve stops here */
	if(!(lineFlag[0]+lineFlag[1]+lineFlag[2])){
	  flag=0;
	  continue;
	}
	/* select which pixel to choose as the continuation */
	/* you surely cant have a metric greater than 100!! */
	minimum=100.0;
	which=0;
	for(i=0;i<3;i++){
	  if(lineFlag[i]){
	    if(metric[i]<minimum){
	      which=i;
	      minimum=metric[i];
	    }
	  }
	}
	/* include the line point in the curve */
	nextRow=currentRow+offset[check[which]].x;
	nextCol=currentCol+offset[check[which]].y;
	nextIndex=label[nextRow][nextCol]-1;
	linePoints[nextIndex].flag=0;
	contour[direction].row[curveLength]=nextRow;
	contour[direction].col[curveLength]=nextCol;
	curveLength++;
	/* check if the maximum curve length is reached */
	ASSERT(curveLength<=maxCurveLen,status,TS_ARRAY_OVERFLOW,MSG_TS_ARRAY_OVERFLOW);
	/* if the line point is also a line start point mark it is processed */
	if(linePoints[nextIndex].pValue==2)
	  linePoints[nextIndex].flag=0;	
	mapContour[nextRow][nextCol]=out->numberOfCurves+1;	
	/* redefine current position and make the new point included as the currect position */
	currentRow=nextRow;
	currentCol=nextCol;	      
      }
      /* end of the Contour reached */
      contour[direction].n=curveLength;
    }
    /* check if the length of the curve is greater than the threshhold*/
    if(contour[0].n+contour[1].n-1 >= LENGTH_THRESHOLD){
      /* record the curve found in the output structure by joining the left and right Contours*/
      out->curves = LALRealloc(out->curves,sizeof(Curve)*(out->numberOfCurves+1));
            ((out->curves)[out->numberOfCurves]).n = contour[0].n+contour[1].n-1;
      out->curves[out->numberOfCurves].row=LALMalloc(sizeof(INT4)* out->curves[out->numberOfCurves].n);
      out->curves[out->numberOfCurves].col=LALMalloc(sizeof(INT4)* out->curves[out->numberOfCurves].n);
      out->curves[out->numberOfCurves].depth=LALMalloc(sizeof(REAL4)* out->curves[out->numberOfCurves].n);
      powerHalfContourA = 0;
      for(i=0;i<contour[0].n;i++){
	out->curves[out->numberOfCurves].row[i] = contour[0].row[contour[0].n - 1 - i];
	out->curves[out->numberOfCurves].col[i] = contour[0].col[contour[0].n - 1 - i];	
	/* INSERT CODE TO PLUCK OUT VALUES FROM TFMAP */
      	out->curves[out->numberOfCurves].depth[i] = 
	  TimeFreqMap->map[contour[0].row[contour[0].n - 1 - i]][contour[0].col[contour[0].n - 1 - i]];
	powerHalfContourA = powerHalfContourA + out->curves[out->numberOfCurves].depth[i];
		}
      powerHalfContourB = 0;
      for(i=1;i<contour[1].n;i++){
	out->curves[out->numberOfCurves].row[i+contour[0].n-1] = contour[1].row[i];     
	out->curves[out->numberOfCurves].col[i+contour[0].n-1] = contour[1].col[i]; 
	out->curves[out->numberOfCurves].depth[i+contour[0].n-1] = 
		  TimeFreqMap->map[contour[1].row[i]][contour[1].col[i]];
        powerHalfContourB = powerHalfContourB + out->curves[out->numberOfCurves].depth[i+contour[0].n-1];
      }
      out->curves[out->numberOfCurves].totalPower = powerHalfContourA + powerHalfContourB;
      if(contour[0].junction+contour[1].junction)
	out->curves[out->numberOfCurves].junction=1;
      /* increment the number of curves */
      out->numberOfCurves++;
      ASSERT(out->numberOfCurves<=MAX_NUMBER_OF_CURVES,status,TS_TOO_MANY_CURVES,MSG_TS_TOO_MANY_CURVES);
    }
  }  
  /* Free up used memory */
  LALFree(linePoints);
  LALFree(contour[0].row);
  LALFree(contour[0].col);
  LALFree(contour[1].row);
  LALFree(contour[1].col);
  for(i=0;i<params->height;i++)
    LALFree(label[i]);
  LALFree(label);
  for(i=0;i<params->height;i++)
    LALFree(mapContour[i]);
  LALFree(mapContour);  
  /* Return to Calling program */
  RETURN(status);  
}


/* Computes the average Gaussian function over the pixel 
   Uses a simple averaging process 
*/
REAL4 
Gauss0(INT4 pixel, REAL4 sigma)
{
    REAL8 limit1,limit2; /* the limits over which to average the function */
    INT4 numberInt=20; /* number of bins; */
    REAL8 interval; /* the bin width */
    REAL8 loop;     /* the loop variable */
    INT4 i;         /* loop variable */
    REAL8 sum;     /* the sum */
    REAL4 result; /* the value to return */
    REAL8 norm;   /* normalizing constant */
    REAL8 temp;   /* temporary variable */
    
    limit1 = (REAL4)pixel - 0.5;
    limit2 = (REAL4)pixel + 0.5;
    interval=(limit2-limit1)/numberInt;
    sum=0.0;
    norm = 1.0/(sigma*sqrt(LAL_PI*2.0));
    for(i=0;i<=numberInt;i++){
        loop = limit1 + i*interval;
        temp = -0.5*((loop*loop)/(sigma*sigma));
        sum += norm * exp(temp);
    }
    sum /= (numberInt+1);
    result = sum;
    return result;
}

/* Computes the average of the derivative of the Gaussian function over the pixel 
   Uses a simple averaging process 
*/
REAL4 
Gauss1(INT4 pixel, REAL4 sigma)
{
    REAL8 limit1,limit2; /* the limits over which to average the function */
    INT4 numberInt=20; /* number of bins; */
    REAL8 interval; /* the bin width */
    REAL8 loop;     /* the loop variable */
    INT4 i;         /* loop variable */
    REAL8 sum;     /* the sum */
    REAL4 result; /* the value to return */
    REAL8 norm;   /* normalizing constant */
    REAL8 temp;   /* temporary variable */
    
    limit1 = (REAL4)pixel - 0.5;
    limit2 = (REAL4)pixel + 0.5;
    interval=(limit2-limit1)/numberInt;
    sum=0.0;
    norm = -1.0/(sigma*sigma*sigma*sqrt(LAL_PI*2.0));
    for(i=0;i<=numberInt;i++){
        loop = limit1 + i*interval;
        temp = -0.5*((loop*loop)/(sigma*sigma));
        sum += norm * loop * exp(temp);
    }
    sum /= (numberInt+1);
    result = sum;
    return result;
}

/* Computes the average of the second derivative of the Gaussian function over the pixel 
   Uses a simple averaging process 
*/
REAL4 
Gauss2(INT4 pixel, REAL4 sigma)
{
    REAL8 limit1,limit2; /* the limits over which to average the function */
    INT4 numberInt=20; /* number of bins;*/
    REAL8 interval; /* the bin width */
    REAL8 loop;     /* the loop variable */
    INT4 i;         /* loop variable */
    REAL8 sum;     /* the sum */
    REAL4 result; /* the value to return */
    REAL8 norm;   /* normalizing constant */
    REAL8 temp;   /* temporary variable */
    
    limit1 = (REAL4)pixel - 0.5;
    limit2 = (REAL4)pixel + 0.5;
    interval=(limit2-limit1)/numberInt;
    sum=0.0;
    norm = 1.0/(sigma*sigma*sigma*sqrt(LAL_PI*2.0));
    for(i=0;i<=numberInt;i++){
        loop = limit1 + i*interval;
        temp = ((loop*loop)/(sigma*sigma));
        sum += norm * exp(-0.5*temp) * (temp - 1.0);
    }
    sum /= (numberInt+1);
    result = sum;
    return result;
}

/* computes the angle given the 2 quadratures to within [0,pi] */
REAL4 
GetAngle(
    REAL4 y,
    REAL4 x)
{
    REAL4 angle;
    REAL4 epsilon = 1.e-4;

    if(x==0.0)
        angle=LAL_PI/2.0;
    else
        angle=atan2(y,x);
    if(angle<0) angle += 2*LAL_PI;
    if(angle>LAL_PI) angle -= LAL_PI;
    if(angle>=(LAL_PI-epsilon))
        angle=LAL_PI-epsilon;
    if(angle<epsilon) angle=epsilon;
    return angle;
}

/* routine called by qsort  
   Compare the elements of the LinePoint array. 
*/
static int
CompareLinePoints(const void *element1, const void *element2)
{
    if( ((const LinePoints *)element2)->pValue > ((const LinePoints *)element1)->pValue)
        return 1;
    if( ((const LinePoints *)element1)->pValue > ((const LinePoints *)element2)->pValue)
        return -1;    
    if(((const LinePoints *)element2)->eigen > ((const LinePoints *)element1)->eigen)
        return 1;
    else
        return -1;      
} 







