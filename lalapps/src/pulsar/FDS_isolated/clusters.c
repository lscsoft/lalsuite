/** \file clusters.c 
 * Set of routines that have been used to calculate outliers and clusters of outliers in the data.*/
/* Author: M. A. Papa - AEI August 2003 */
/* Revision: Y. Itoh - AEI December 2003  */
/*           Commented out "if Nclust==0. NclustPoints[Nclust]=k"  */

/* $Id$ */

/* #include <stdlib.h> */
#include <lal/LALStdlib.h>
#include <math.h>
#include "clusters.h"
#include "rngmed.h"

/*this code estimates the floor of a givendata set by the running median.*/
int EstimateFloor(REAL8Vector *input, INT2 windowSize, REAL8Vector *output){
  /* input : vector (N points) over which the running median code is ran with a  */
  /* blocksize = windowsize */
  /* the output of the running median code has N - blocksize+1 points.  */
  /* the output of this routine has N points. The missing blocksize points */
  /* are set equal to the nearest last value of the output of th erunnning median. */


  UINT4 start,start2;
  UINT4 lpc;
  UINT4 nbins=input->length;
  REAL4 halfWindow;

  REAL8 *dmp;

  INT4 M;


  M = nbins-windowSize+1;
  if(!(dmp= (double *) LALCalloc(M,sizeof(double)))){ 
    printf("dmp Memory allocation failure");
    return 1;
  }
  
  /* Call rngmed */
  if (rngmed(input->data, nbins, windowSize, dmp)) {
    fprintf(stderr,"Problem in rngmed()\n");
    return 1;
  }

  /* start is the index of outdata at which the actual rmgmed begins*/
  halfWindow=windowSize/2.0;
  start=((int)(halfWindow))-1;
  start2=(start+nbins-windowSize);
  
  /*Fill-in RngMdnSp */
  for (lpc=0;lpc<start;lpc++){
    output->data[lpc]=dmp[0];
  }
  for (lpc=start;lpc<start2;lpc++){
    output->data[lpc]=dmp[lpc-start];
  }
  for (lpc=start2;lpc<nbins;lpc++){
    output->data[lpc]=dmp[nbins-windowSize];
  }
  
  LALFree(dmp);

  return 0;
}



/*Given a set of outliers and their parameters, this routine finds clusters defined by certain parameters. The main parameter of the cluster algorithm is Dmax: the maximum distance between points above threshold that belong to the same cluster. The other parameter is "smallblock"*/
int DetectClusters(ClustersInput *input, ClustersParams *clParams, Clusters *output){

  INT4 Dist,Dmax;
  INT2 smallBlock;
  INT4 wings;
  INT2 Nclust;
  UINT4 Ntot,NtotCheck;

  UINT4 *Iclust, *NclustPoints; /* first index of cluster, how many points */

  REAL8 *RDMP1, *RDMP2;

  INT4 rwing,lwing;
  REAL8 RngMdn,max1,max2;

  int k,i,i0,lpc,j,imax1,imax2,shift;

  wings      = clParams->wings;
  smallBlock = clParams->smallBlock;

  Dmax = 2*wings;

  /*  find how many clusters there are (Nclust+1) */

  Dist=0; /*distance in the array between adjacent outliers; if > Dmax then 
	   we have found a new cluster.*/
  Nclust=0; /* number of clusters -1 */
  k=0;
  for (i=0;i<input->outliers->Noutliers-1;i++){
    Dist=input->outliers->outlierIndexes[i+1]-input->outliers->outlierIndexes[i];
    k++;
    if (Dist >= Dmax){
      k=0;
      Nclust=Nclust+1;
    }
  }
  
  output->Nclusters=Nclust+1;

  /*  Now that we know how many clusters, allocate space for Iclust */
  if (!(Iclust = (UINT4 *) LALCalloc(Nclust+1,sizeof(UINT4)))){
    printf("Memory allocation failure in Iclust");
    return 0;
  }
  /*  Now that we know how many clusters, allocate space for NclustPoints */
  if (!(NclustPoints = (UINT4 *) LALCalloc(Nclust+1,sizeof(UINT4)))){
    printf("Memory allocation failure in NClustPoints");
    return 0;
  }  
  /*  Now that we know how many clusters, allocate space 
      for output->NclustPoints */
  if (!(output->NclustPoints = (UINT4 *) LALRealloc(output->NclustPoints,(Nclust+1)*sizeof(UINT4)))){
    printf("Memory allocation failure in output->NclustPoints");
    return 0;
  }  

  /*  Populate Iclust and NclustPoints */
  Dist=0;
  Nclust=0; /* number of clusters - 1 */
  k=0;      /* number of outlier values in a single cluster */
            /* note: in general number of points is NOT this*/
            /* it is only if all outliers are neighbouring. */
  
  i0=0;
  for (i=0;i< input->outliers->Noutliers-1;i++){
    Dist=input->outliers->outlierIndexes[i+1]-input->outliers->outlierIndexes[i];
    k++;
    if (Dist >= Dmax){
      NclustPoints[Nclust]=input->outliers->outlierIndexes[i]-input->outliers->outlierIndexes[i0]+1; /* i0: first outlier of current (Nclust) cluster*/
      i0=i+1; /* now i0 is the first outlier index of next cluster */
      Iclust[Nclust]=input->outliers->outlierIndexes[i+1-k];
      k=0;
      Nclust=Nclust+1;
    }
  }/*end loop on outliers*/
  
  NclustPoints[Nclust]=input->outliers->outlierIndexes[i]-input->outliers->outlierIndexes[i0]+1;
  Iclust[Nclust]=input->outliers->outlierIndexes[i+1-k];

  if (Nclust == 0){
    Iclust[Nclust]=input->outliers->outlierIndexes[Nclust];
  }

  if (Nclust != 0 && k!=0){ 	/*  if the last considered outlier was part of a cluster */
    				/*  and the really last outlier belongs to */
	    			/*  the same cluster */
     k++;
     Iclust[Nclust]=input->outliers->outlierIndexes[i+1-k]; 
     NclustPoints[Nclust]=input->outliers->outlierIndexes[i]-input->outliers->outlierIndexes[i0]+1; 
   } 
  
  if (Nclust != 0 && k==0){ /* if the last outlier is isolated */
    Iclust[Nclust]=input->outliers->outlierIndexes[i];
    NclustPoints[Nclust]=1;
  }


  /*  count how many points per cluster and in total  */
  NtotCheck=0;
  for (lpc=0;lpc<output->Nclusters;lpc++){
    rwing=wings;
    lwing=wings;

    if (lpc == 0)
      lwing=input->outliers->leftwing;
    if (lpc == Nclust)
      rwing=input->outliers->rightwing;

    output->NclustPoints[lpc]=NclustPoints[lpc]+lwing+rwing;
    NtotCheck=NtotCheck+output->NclustPoints[lpc];
  }
  

  /*  Allocate space for output vectors */
  if (!(output->clusters = (REAL8 *) LALRealloc(output->clusters,NtotCheck*sizeof(REAL8)))){
    printf("Memory allocation failure in output->clusters");
    LALFree(Iclust);
    LALFree(NclustPoints);
    return 1;
  }
   if (!(output->Iclust = (UINT4 *) LALRealloc(output->Iclust,NtotCheck*sizeof(UINT4)))){
    printf("Memory allocation failure in output->Iclust");
    return 0;
  }



   /*  compute line profiles  */
   /*  for every cluster:*/
  Ntot=0;
  for (lpc=0;lpc<output->Nclusters;lpc++){
    
    rwing=wings;
    lwing=wings;
    
    if (lpc == 0)
      lwing=input->outliers->leftwing;
    if (lpc == Nclust)
      rwing=input->outliers->rightwing;


    /*The next block is NOT relevant for F stat clusters searches because there smallBlock =1 */
    if (smallBlock > 1){

   /*  k : number of samples to feed into rngmed code.  */
    /*  They are smallBlock-1 more than the output.  */
    /*  The output is output->NclustPoints[lpc].     */

    k= output->NclustPoints[lpc]+smallBlock-1; 


    /*  RDMP2 has the relevant part (i.e. the outliers + wings for  */
    /*  this profile) of the input data. */
    /*  RDMP1 will have the profile */

    if (!(RDMP2 = (REAL8 *) LALCalloc(k,sizeof(REAL8)))){
      printf("RDMP2 memory allocation failure");
      return 0;
    }
    if (!(RDMP1 = (REAL8 *) LALCalloc(output->NclustPoints[lpc],sizeof(REAL8)))){
      printf("RDMP1 memory allocation failure");
      return 0;
    }

      /*  compute max of input data */
      max2=0;
      imax2=0;
      for (i=0;i<k;i++){
	/*       j=Iclust[lpc]-rwing-lwing+i; */
	j=Iclust[lpc]-lwing+i;
	RDMP2[i]=input->outliersInput->data->data[j];
	if (input->outliersInput->data->data[j] > max2){
	  max2 = input->outliersInput->data->data[j];
	  imax2 = i;
	}
      }
      
      rngmed(RDMP2, k, smallBlock, RDMP1);
      
      /*  compute max of output data */
      max1=0;
      for (i=0;i<k-smallBlock+1;i++){
	if (RDMP1[i] > max1){
	  max1 = RDMP1[i];
	  imax1 = i;
	  if (imax1 == output->NclustPoints[lpc]-1)
	    imax1=imax1-1;
	}
      }
      
      /*  put maximum of input vector in place of max of running median */
      /* RDMP1[imax1+1]=RDMP2[imax2]; */
      shift=(UINT4)((smallBlock+0.5)/2.0);
      if (shift > imax2)
	shift=imax2;
      if ((imax2-shift) > output->NclustPoints[lpc]-1)
	shift=imax2-output->NclustPoints[lpc]+1;
      RDMP1[imax2-shift]=RDMP2[imax2];
      
      
      /* write output */
      for (i=0;i<k-smallBlock+1;i++){
	j=Iclust[lpc]+smallBlock/2+i-lwing;
	output->Iclust[Ntot]=j;
	RngMdn=input->outliersInput->data->data[j]/input->outliers->ratio[j];
	output->clusters[Ntot]=RDMP1[i]/RngMdn;
	Ntot++;
      }

      LALFree(RDMP2);
      LALFree(RDMP1);
  
    }


    /*Fill-in the output structures */
    if (smallBlock == 1){
      for (i=0;i< output->NclustPoints[lpc];i++){
	/* for (i=0;i< output->Nclusters;i++){ */
	j=i+Iclust[lpc]-lwing;
	output->Iclust[Ntot]=j;
	output->clusters[Ntot]=input->outliersInput->data->data[j];
	Ntot++;
      }
    }
    
  

  }/* loop over clusters*/

  if (Ntot != NtotCheck){
    fprintf(stderr,"Ntot not equal NtotCheck In DetectClusters\n");
    LALFree(Iclust);
    LALFree(NclustPoints);
    LALFree(output->Iclust);
    LALFree(output->NclustPoints);
    LALFree(output->clusters);
    return 1;
  }
   
  LALFree(Iclust);
  LALFree(NclustPoints);

  return 0;
}





/* Given the input (the data and the index of the first datum) and the parameters (the "floor" of the data used to normalize the data, the threshold on the normalized data, the wings that we would like to leave at the edges and the index of the first datum) this routine determines the outliers (how many, how many wing-bins at the edges, the indexes of the outliers, their value)*/
int ComputeOutliers(OutliersInput *input, OutliersParams *outliersParams, Outliers *outliers){

  INT4 j,imin;
  INT4 leftwing, rightwing,wings;
  /*ITOH*/
  /*  UINT4 ileft,iright; */
  INT4 ileft,iright; 
  /*ITOH*/

  REAL4 thresh;

  UINT4 i,lpc,nbins, NI;

  REAL4 *RDMP;
  UINT4 *IDMP;

  nbins=input->data->length;
  wings=outliersParams->wings;
  imin=input->ifmin;
  thresh = outliersParams->Thr;

  if (!(outliers->ratio=(REAL8 *)LALMalloc(nbins*sizeof(REAL8)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers.ratio\n");
    return 1;
  }

  if (!(RDMP = (REAL4 *) LALCalloc(nbins,sizeof(REAL4)))){
    printf("RDMP memory allocation failure in ComputeOutliers");
    return 0;
  }
  if (!(IDMP = (UINT4 *) LALCalloc(nbins,sizeof(UINT4)))){
    printf("IDMP memory allocation failure in ComputeOutliers");
    return 0;
  }


  /*compute ratio of the data to its floor. If this ratio is > threshold then save the value of the index in IDMP and the corresponding value of the ration in RDMP*/
  i=0;     
  for (lpc=0;lpc<nbins;lpc++){
    outliers->ratio[lpc]=input->data->data[lpc]/outliersParams->Floor->data[lpc];
    /*     fprintf(outfile1,"%f %E %E %E\n",(GV.ifmin+lpc)/GV.tsft,Sp[lpc],RngMdnSp[lpc],R[lpc]); */
    if (outliers->ratio[lpc] > thresh){
      IDMP[i]=lpc;
      RDMP[i]=outliers->ratio[lpc];
      i++;
    }
  }

  /*total number of outliers*/
  NI=i;
  outliers->Noutliers=NI;

  /*if no outliers are found clean up and exit*/
  if (NI == 0){
    LALFree(RDMP);
    LALFree(IDMP);
    return 0;
  }

  if (!(outliers->outlierIndexes=(UINT4 *)LALMalloc(NI*sizeof(UINT4)))){
    fprintf(stderr,"Memory allocation failure for SpOutliers.Indexes\n");
    return 1;
  }

  /* check left border */
  /*  leftwing : how many points on the left wing */
  leftwing = wings;
  ileft = IDMP[0]-wings;
  if (ileft < 0){
    ileft = 0;
    leftwing = IDMP[0];
  }
  
  /*  check right border */
  /*  rightwing : how many points on the right wing */
  rightwing = wings;
  iright = IDMP[NI-1]+ wings;
  if (iright > (nbins-1)){
    iright = nbins-1;
    rightwing = nbins - IDMP[NI-1]-1;
    /* ITOH */
    /* rightwing = nbins - IDMP[NI-1] - 1; */
  }

  /*fill-in outlierIndexes, the wings and exit */
  for (j=0;j<NI;j++){
    outliers->outlierIndexes[j] = IDMP[j];
  }
  
  outliers->rightwing = rightwing;
  outliers->leftwing = leftwing;  

  LALFree(RDMP);
  LALFree(IDMP);

  return 0;
}
