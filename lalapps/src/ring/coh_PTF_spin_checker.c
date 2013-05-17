#include "coh_PTF.h"

#define PROGRAM_NAME "lalapps_coh_PTF_spin_checker"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE   "$Source$"
#define CVS_DATE     "$Date$"

/* This function should be migrated to option.c */
/* warning: returns a pointer to a static variable... not reenterant */
/* only call this routine once to initialize params! */
/* also do not attempt to free this pointer! */
static struct coh_PTF_params *coh_PTF_get_params( int argc, char **argv )
{
  static struct coh_PTF_params params;
  static char programName[] = PROGRAM_NAME;
  static char cvsRevision[] = CVS_REVISION;
  static char cvsSource[]   = CVS_SOURCE;
  static char cvsDate[]     = CVS_DATE;
  coh_PTF_parse_options( &params, argc, argv );
  coh_PTF_params_sanity_check( &params ); /* this also sets various params */
  coh_PTF_params_spin_checker_sanity_check( &params );
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}


int main( int argc, char **argv )
{
  INT4 i,j;
  struct coh_PTF_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *psdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  COMPLEX8FFTPlan         *invPlan   = NULL;
  REAL4TimeSeries         *channel[LAL_NUM_IFO+1];
  REAL4FrequencySeries    *invspec[LAL_NUM_IFO+1];
  RingDataSegments        *segments[LAL_NUM_IFO+1];
  INT4                    numTmplts = 0;
  INT4  startTemplate     = -1;           /* index of first template      */
  INT4  stopTemplate      = -1;           /* index of last template       */
  INT4 numSegments        = 0;
  InspiralTemplate        *PTFtemplate = NULL;
  InspiralTemplate        *PTFbankhead = NULL;
  SnglInspiralTable       *PTFSpinTmplt = NULL;
  SnglInspiralTable       *PTFSpinTmpltHead = NULL;
  SnglInspiralTable       *PTFNoSpinTmplt = NULL;
  SnglInspiralTable       *PTFNoSpinTmpltHead = NULL;
  SnglInspiralTable       *PTFLastTmplt = NULL;
  FindChirpTemplate       *fcTmplt     = NULL;
  FindChirpTmpltParams    *fcTmpltParams      = NULL;
  FindChirpInitParams     *fcInitParams = NULL;
  UINT4                   numPoints,ifoNumber,spinTemplate;
  REAL8Array              *PTFM[LAL_NUM_IFO+1];
  REAL8Array              *PTFN[LAL_NUM_IFO+1];
  COMPLEX8VectorSequence  *PTFqVec[LAL_NUM_IFO+1];
  time_t                  startTime;
  LALDetector             *detectors[LAL_NUM_IFO+1];
  REAL4                   *timeOffsets;
  REAL4                   *Fplus;
  REAL8                   FplusTmp;
  REAL4                   *Fcross;
  REAL8                   FcrossTmp;
  REAL8                   detLoc[3];
  REAL4TimeSeries         UNUSED *pValues[10];
  REAL4TimeSeries         UNUSED *gammaBeta[2];
  LIGOTimeGPS             segStartTime;
  MultiInspiralTable      *eventList = NULL;
  UINT4                   numDetectors = 0;
  UINT4                   singleDetector = 0;
  char                    bankFileName[256];
  
  startTime = time(NULL);

  /* set error handlers to abort on error */
  set_abrt_on_error();

  /* options are parsed and debug level is set here... */
  

  /* no lal mallocs before this! */
  params = coh_PTF_get_params( argc, argv );

  /* create process params */
  procpar = create_process_params( argc, argv, PROGRAM_NAME );

  verbose("Read input params %ld \n", time(NULL)-startTime);

  /* create forward and reverse fft plans */
  fwdplan = coh_PTF_get_fft_fwdplan( params );
  psdplan = coh_PTF_get_fft_psdplan( params );
  revplan = coh_PTF_get_fft_revplan( params );

  verbose("Made fft plans %ld \n", time(NULL)-startTime);

  /* Determine if we are analyzing single or multiple ifo data */

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      numDetectors++;
    }
  }
  /* NULL out the pValues pointer array */
  for ( i = 0 ; i < 10 ; i++ )
  {
    pValues[i] = NULL;
  }   
  gammaBeta[0] = NULL;
  gammaBeta[1] = NULL;

  if (numDetectors == 0 )
  {
    fprintf(stderr,"You have not specified any detectors to analyse");
    return 1;
  }
  else if (numDetectors == 1 )
  {
    fprintf(stdout,"You have only specified one detector, why are you using the coherent code? \n");
    singleDetector = 1;
  }

  /* Convert the file names */
  if ( params->bankFile )
    strncpy(bankFileName,params->bankFile,sizeof(bankFileName)-1);

  REAL4                    *timeSlideVectors;
  timeSlideVectors=LALCalloc(1, (LAL_NUM_IFO+1)*
                                params->numOverlapSegments*sizeof(REAL4));
  memset(timeSlideVectors, 0,
         (LAL_NUM_IFO+1) * params->numOverlapSegments * sizeof(REAL4));

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    /* Initialize some of the structures */
    channel[ifoNumber] = NULL;
    invspec[ifoNumber] = NULL;
    segments[ifoNumber] = NULL;
    PTFM[ifoNumber] = NULL;
    PTFN[ifoNumber] = NULL;
    PTFqVec[ifoNumber] = NULL;
    if ( params->haveTrig[ifoNumber] )
    {
      /* Read in data from the various ifos */
      channel[ifoNumber] = coh_PTF_get_data(params,params->channel[ifoNumber],\
                               params->dataCache[ifoNumber],ifoNumber );
      coh_PTF_rescale_data (channel[ifoNumber],1E20);

      /* compute the spectrum */
      invspec[ifoNumber] = coh_PTF_get_invspec( channel[ifoNumber], fwdplan,\
                               revplan, psdplan, params );

      /* create the segments */
      segments[ifoNumber] = coh_PTF_get_segments( channel[ifoNumber],\
           invspec[ifoNumber], fwdplan, ifoNumber, timeSlideVectors, params );
      
      numSegments = segments[ifoNumber]->numSgmnt;

      verbose("Created segments for one ifo %ld \n", time(NULL)-startTime);
    }
  }

  /* Determine time delays and response functions */ 

  timeOffsets = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL4 ));
  Fplus = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL4 ));
  Fcross = LALCalloc(1, numSegments*LAL_NUM_IFO*sizeof( REAL4 ));
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    detectors[ifoNumber] = LALCalloc( 1, sizeof( *detectors[ifoNumber] ));
    XLALReturnDetector(detectors[ifoNumber] ,ifoNumber);
    for ( i = 0; i < 3; i++ )
    {
      detLoc[i] = (double) detectors[ifoNumber]->location[i];
    }
    for ( j = 0; j < numSegments; ++j )
    {
      /* Despite being called segStartTime we use the time at the middle 
      * of a segment */
      segStartTime = params->startTime;
      
      /*XLALGPSAdd(&segStartTime,(j+1)*params->segmentDuration/2.0);*/
      XLALGPSAdd(&segStartTime,8.5*params->segmentDuration/2.0);
      /*XLALGPSMultiply(&segStartTime,0.);
      XLALGPSAdd(&segStartTime,874610713.072549154);*/
      timeOffsets[j*LAL_NUM_IFO+ifoNumber] = (REAL4)
          XLALTimeDelayFromEarthCenter(detLoc,params->rightAscension,
          params->declination,&segStartTime);
      XLALComputeDetAMResponse(&FplusTmp, &FcrossTmp,
         detectors[ifoNumber]->response,params->rightAscension,
         params->declination,0.,XLALGreenwichMeanSiderealTime(&segStartTime));
      Fplus[j*LAL_NUM_IFO + ifoNumber] = (REAL4) FplusTmp;
      Fcross[j*LAL_NUM_IFO + ifoNumber] = (REAL4) FcrossTmp;
    }
    LALFree(detectors[ifoNumber]);
  }
  

  numPoints = floor( params->segmentDuration * params->sampleRate + 0.5 );

  /* Initialize some of the structures */
  ifoNumber = LAL_NUM_IFO;
  channel[ifoNumber] = NULL;
  invspec[ifoNumber] = NULL;
  segments[ifoNumber] = NULL;
  PTFM[ifoNumber] = NULL;
  PTFN[ifoNumber] = NULL;
  PTFqVec[ifoNumber] = NULL;

  /* Create the relevant structures that will be needed */
  fcInitParams = LALCalloc( 1, sizeof( *fcInitParams ));
  fcTmplt = LALCalloc( 1, sizeof( *fcTmplt ) );
  fcTmpltParams = LALCalloc ( 1, sizeof( *fcTmpltParams ) );
  fcTmpltParams->approximant = FindChirpPTF;
  fcTmplt->PTFQtilde =
      XLALCreateCOMPLEX8VectorSequence( 5, numPoints / 2 + 1 );
/*  fcTmplt->PTFBinverse = XLALCreateArrayL( 2, 5, 5 );
  fcTmplt->PTFB = XLALCreateArrayL( 2, 5, 5 );*/
  fcTmplt->PTFQ = XLALCreateVectorSequence( 5, numPoints );
  fcTmpltParams->PTFphi = XLALCreateVector( numPoints );
  fcTmpltParams->PTFomega_2_3 = XLALCreateVector( numPoints );
  fcTmpltParams->PTFe1 = XLALCreateVectorSequence( 3, numPoints );
  fcTmpltParams->PTFe2 = XLALCreateVectorSequence( 3, numPoints );
  fcTmpltParams->fwdPlan =
        XLALCreateForwardREAL4FFTPlan( numPoints, 1 );
  fcTmpltParams->deltaT = 1.0/params->sampleRate;
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( params->haveTrig[ifoNumber] )
    {
      PTFM[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
      PTFN[ifoNumber] = XLALCreateREAL8ArrayL( 2, 5, 5 );
      PTFqVec[ifoNumber] = XLALCreateCOMPLEX8VectorSequence ( 5, numPoints );
    }
  }
  /* Create an inverser FFT plan */
  invPlan = XLALCreateReverseCOMPLEX8FFTPlan( numPoints, 1 );

  /* Read in the tmpltbank xml file */
  numTmplts = InspiralTmpltBankFromLIGOLw( &PTFtemplate,bankFileName,
      startTemplate, stopTemplate );
  PTFbankhead = PTFtemplate;
  /*fake_template (PTFtemplate);*/

  for (i = 0; (i < numTmplts); PTFtemplate = PTFtemplate->next, i++)
  {
    /* Determine if we can model this template as non-spinning */
    spinTemplate = 1;
    if (PTFtemplate->chi < 0.1) 
      spinTemplate = 0;
    PTFtemplate->approximant = FindChirpPTF;
    PTFtemplate->order = LAL_PNORDER_TWO;
    PTFtemplate->fLower = 38.;
    /* Generate the Q freq series of the template */
    coh_PTF_template(fcTmplt,PTFtemplate,fcTmpltParams);

    verbose("Generated template %d at %ld \n", i, time(NULL)-startTime);

    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
    {
      if ( params->haveTrig[ifoNumber] )
      {
        memset( PTFM[ifoNumber]->data, 0, 25 * sizeof(REAL8) );
        memset( PTFN[ifoNumber]->data, 0, 25 * sizeof(REAL8) );
        memset( PTFqVec[ifoNumber]->data, 0,
                  5 * numPoints * sizeof(COMPLEX8) );
        coh_PTF_normalize(params,fcTmplt,invspec[ifoNumber],PTFM[ifoNumber],
              PTFN[ifoNumber],PTFqVec[ifoNumber],
              &segments[ifoNumber]->sgmnt[0],invPlan,1);
      }
    }

    spinTemplate = coh_PTF_spin_checker(PTFM,PTFN,params,singleDetector,Fplus,Fcross,numSegments/2); 
    if (spinTemplate)
      verbose("Template %d treated as spin at %ld \n",i,time(NULL)-startTime);
    else
      verbose("Template %d treated as non spin at %ld \n",i,time(NULL)-startTime);
    if (spinTemplate)
    {
      if ( !PTFSpinTmplt )
      {
        PTFSpinTmpltHead = conv_insp_tmpl_to_sngl_table(PTFtemplate,i);
        PTFSpinTmplt = PTFSpinTmpltHead;
      }
      else
      {
        PTFLastTmplt = PTFSpinTmplt;
        PTFSpinTmplt = conv_insp_tmpl_to_sngl_table(PTFtemplate,i);
        PTFLastTmplt->next = PTFSpinTmplt;
      }
    }
    else
    {
      if ( !PTFNoSpinTmplt )
      {
        PTFNoSpinTmpltHead = conv_insp_tmpl_to_sngl_table(PTFtemplate,i);
        PTFNoSpinTmplt = PTFNoSpinTmpltHead;
      }
      else
      {
        PTFLastTmplt = PTFNoSpinTmplt;
        PTFNoSpinTmplt = conv_insp_tmpl_to_sngl_table(PTFtemplate,i);
        PTFLastTmplt->next = PTFNoSpinTmplt;
      }
    }

    if (! PTFtemplate->event_id)
    {
      PTFtemplate->event_id = (EventIDColumn *)
            LALCalloc(1, sizeof(EventIDColumn) );
      PTFtemplate->event_id->id = i;
    }
  }

  coh_PTF_output_tmpltbank(params->spinBankName,PTFSpinTmpltHead,procpar,params);
  coh_PTF_output_tmpltbank(params->noSpinBankName,PTFNoSpinTmpltHead,procpar,params);

  verbose("Generated output xml files, cleaning up and exiting at %ld \n",
      time(NULL)-startTime);

  LALFree(timeSlideVectors);
  coh_PTF_cleanup(params,procpar,fwdplan,psdplan,revplan,invPlan,channel,
      invspec,segments,eventList,NULL,PTFbankhead,fcTmplt,fcTmpltParams,
      fcInitParams,PTFM,PTFN,PTFqVec,timeOffsets,NULL,Fplus,Fcross,NULL,NULL,\
      NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  while ( PTFSpinTmpltHead )
  {
    PTFSpinTmplt = PTFSpinTmpltHead;
    PTFSpinTmpltHead = PTFSpinTmplt->next;
    if ( PTFSpinTmplt->event_id )
    {
      LALFree( PTFSpinTmplt->event_id );
    }
    LALFree( PTFSpinTmplt );
  }
  while ( PTFNoSpinTmpltHead )
  {
    PTFNoSpinTmplt = PTFNoSpinTmpltHead;
    PTFNoSpinTmpltHead = PTFNoSpinTmplt->next;
    if ( PTFNoSpinTmplt->event_id )
    {
      LALFree( PTFNoSpinTmplt->event_id );
    }
    LALFree( PTFNoSpinTmplt );
  }

  LALCheckMemoryLeaks();
  return 0;
}

int coh_PTF_spin_checker(
    REAL8Array              *PTFM[LAL_NUM_IFO+1],
    REAL8Array              *PTFN[LAL_NUM_IFO+1],
    struct coh_PTF_params   *params,
    UINT4                   singleDetector,
    REAL4                   *Fplus,
    REAL4                   *Fcross,
    INT4                    segmentNumber
)

{
  UINT4 i, j, k, l, vecLength, vecLengthTwo, vecLengthSquare, UNUSED vecLengthTwoSquare;
  UINT4 UNUSED nsVecLength, nsVecLengthTwo, UNUSED nsVecLengthSquare,UNUSED nsVecLengthTwoSquare;

  if (singleDetector == 1)
  {
    vecLength = 5;
    vecLengthTwo = 5;
    vecLengthSquare = 25;
    vecLengthTwoSquare = 25;
    nsVecLength = 2;
    nsVecLengthTwo = 2;
    nsVecLengthSquare = 4;
    nsVecLengthTwoSquare = 4;
  }
  else
  {
    vecLength = 5;
    vecLengthTwo = 10;
    vecLengthSquare = 25;
    vecLengthTwoSquare = 100;
    nsVecLength = 2;
    nsVecLengthTwo = 4;
    nsVecLengthSquare = 4;
    nsVecLengthTwoSquare = 16;
  }

/*  REAL8Array  *B, *Binv;*/
  REAL4 a[LAL_NUM_IFO], b[LAL_NUM_IFO];
  REAL4 zh[vecLengthSquare],sh[vecLengthSquare],yu[vecLengthSquare];

  gsl_matrix *B2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *N2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *nsB2 = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *nsN2 = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *tempM = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *nsTempM = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *b2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *b2i = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *n2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *nsb2 = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *nsb2i = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *nsn2 = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *Binv2 = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_eigen_symmv_workspace *matTemp = gsl_eigen_symmv_alloc (vecLengthTwo);
  gsl_eigen_symmv_workspace *nsMatTemp = gsl_eigen_symmv_alloc (nsVecLengthTwo);
  gsl_matrix *eigenvecs = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_vector *eigenvals = gsl_vector_alloc(vecLengthTwo);
  gsl_matrix *nsEigenvecs = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_vector *nsEigenvals = gsl_vector_alloc(nsVecLengthTwo);
  gsl_matrix *MMM = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *MMN = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *NMN = gsl_matrix_alloc(vecLengthTwo,vecLengthTwo);
  gsl_matrix *nsMMM = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *nsMMN = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);
  gsl_matrix *nsNMN = gsl_matrix_alloc(nsVecLengthTwo,nsVecLengthTwo);

  for (i = 0; i < LAL_NUM_IFO; i++)
  {
    a[i] = Fplus[segmentNumber*LAL_NUM_IFO+i];
    b[i] = Fcross[segmentNumber*LAL_NUM_IFO+i];
  }

  /* Create the M matrix (Diego calls this B)*/
  for (i = 0; i < vecLength; i++ )
  {
    for (j = 0; j < vecLength; j++ )
    {
      zh[i*vecLength+j] = 0;
      sh[i*vecLength+j] = 0;
      yu[i*vecLength+j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          /* Note that PTFM is always 5x5 even for non-spin */
          zh[i*vecLength+j] += a[k]*a[k] * PTFM[k]->data[i*5+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFM[k]->data[i*5+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFM[k]->data[i*5+j];
        }
      }
    }
  }

  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      if ( i < vecLength && j < vecLength )
      {
        gsl_matrix_set(B2,i,j,zh[i*vecLength+j]);
      }
      else if ( i > (vecLength-1) && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j,sh[(i-vecLength)*vecLength + (j-vecLength)]);
      }
      else if ( i < vecLength && j > (vecLength-1))
      {
        gsl_matrix_set(B2,i,j, yu[i*vecLength + (j-vecLength)]);
      }
      else if ( i > (vecLength-1) && j < vecLength)
      {
        gsl_matrix_set(B2,i,j,yu[j*vecLength + (i-vecLength)]);
      }
      else
        fprintf(stderr,"BUGGER! Something went wrong.");
      gsl_matrix_set(Binv2,i,j,gsl_matrix_get(B2,i,j));
      /*fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));*/
    }
    /*fprintf(stdout,"\n");*/
  }

  /* Create the Nmatrix */
  for (i = 0; i < vecLength; i++ )
  {
    for (j = 0; j < vecLength; j++ )
    {
      zh[i*vecLength+j] = 0;
      sh[i*vecLength+j] = 0;
      yu[i*vecLength+j] = 0;
      for( k = 0; k < LAL_NUM_IFO; k++)
      {
        if ( params->haveTrig[k] )
        {
          /* Note that PTFM is always 5x5 even for non-spin */
          zh[i*vecLength+j] += a[k]*a[k] * PTFN[k]->data[i*5+j];
          sh[i*vecLength+j] += b[k]*b[k] * PTFN[k]->data[i*5+j];
          yu[i*vecLength+j] += a[k]*b[k] * PTFN[k]->data[i*5+j];
        }
      }
    }
  }

  for (i = 0; i < vecLengthTwo; i++ )
  {
    for (j = 0; j < vecLengthTwo; j++ )
    {
      if ( i < vecLength && j < vecLength )
      {
        gsl_matrix_set(N2,i,j,zh[i*vecLength+j]);
      }
      else if ( i > (vecLength-1) && j > (vecLength-1))
      {
        gsl_matrix_set(N2,i,j,sh[(i-vecLength)*vecLength + (j-vecLength)]);
      }
      else if ( i < vecLength && j > (vecLength-1))
      {
        gsl_matrix_set(N2,i,j, yu[i*vecLength + (j-vecLength)]);
      }
      else if ( i > (vecLength-1) && j < vecLength)
      {
        gsl_matrix_set(N2,i,j,yu[j*vecLength + (i-vecLength)]);
      }
      else
        fprintf(stderr,"BUGGER! Something went wrong.");
      /*fprintf(stdout,"%f ",gsl_matrix_get(B2,i,j));*/
    }
    /*fprintf(stdout,"\n");*/
  }

  /*fprintf(stdout,"\n \n");*/

  /* Here we compute the eigenvalues and eigenvectors of B2 */
  gsl_eigen_symmv (B2,eigenvals,eigenvecs,matTemp);

  /* We also compute the eigenvalues and eigenvectors of non spin B2 */
  for (i = 0; i < nsVecLengthTwo; i++)
  {
    for (j = 0; j < nsVecLengthTwo; j++)
    {
      k = i;
      l = j;
      if (i > 1)
        k += 3;
      if ( j > 1)
        l += 3;
      gsl_matrix_set(nsB2,i,j,gsl_matrix_get(Binv2,k,l));
      gsl_matrix_set(nsN2,i,j,gsl_matrix_get(N2,k,l));
    }
  }

  gsl_eigen_symmv (nsB2,nsEigenvals,nsEigenvecs,nsMatTemp);

  for (i = 0; i < nsVecLengthTwo; i++)
  {
    for (j = 0; j < nsVecLengthTwo; j++)
    {
      k = i;
      l = j;
      if (i > 1)
        k += 3;
      if ( j > 1)
        l += 3;
      gsl_matrix_set(nsB2,i,j,gsl_matrix_get(Binv2,k,l));
    }
  }

  /* Determine the maximum eigenvalue */
  REAL8 max_eigenvalue = 0;
  for (i = 0; i < vecLengthTwo ; i++)
  {
    if (gsl_vector_get(eigenvals,i) > max_eigenvalue)
      max_eigenvalue = gsl_vector_get(eigenvals,i);
  }

  /* Rotate the M(B) and N matrices */

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,Binv2,eigenvecs,0.,tempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,eigenvecs,tempM,0.,b2);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,N2,eigenvecs,0.,tempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,eigenvecs,tempM,0.,n2);

  gsl_matrix_scale(b2,1./max_eigenvalue);
  gsl_matrix_scale(n2,1./max_eigenvalue);

  /* And repeat for the non spin matrices */

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsB2,nsEigenvecs,0.,nsTempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,nsEigenvecs,nsTempM,0.,nsb2);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsN2,nsEigenvecs,0.,nsTempM);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1,nsEigenvecs,nsTempM,0.,nsn2);

  gsl_matrix_scale(nsb2,1./max_eigenvalue);
  gsl_matrix_scale(nsn2,1./max_eigenvalue);

  for (i = 0; i < vecLengthTwo; i++)
  {
    for (j = 0; j < vecLengthTwo; j++)
    {
      if (i == j && gsl_matrix_get(b2,i,j) > 1E-5)
        gsl_matrix_set(b2i,i,j,1./gsl_matrix_get(b2,i,j));
      else
        gsl_matrix_set(b2i,i,j,0.);
    }
  }

  for (i = 0; i < nsVecLengthTwo; i++)
  {
    for (j = 0; j < nsVecLengthTwo; j++)
    {
      if (i == j)
      {
        gsl_matrix_set(nsb2i,i,j,1./gsl_matrix_get(nsb2,i,j));
      }
      else
      {
        gsl_matrix_set(nsb2i,i,j,0.);
      }
    }
  }


  /* Calculate all the MMM type matrices */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,b2i,b2,0.,tempM);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,b2,tempM,0.,MMM);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,b2i,n2,0.,tempM);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,n2,tempM,0.,NMN);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,b2,tempM,0.,MMN);

  /* And again for the non spin */

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsb2i,nsb2,0.,nsTempM);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsb2,nsTempM,0.,nsMMM);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsb2i,nsn2,0.,nsTempM);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsn2,nsTempM,0.,nsNMN);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1,nsb2,nsTempM,0.,nsMMN);

  /* Now we begin to loop over the simulations. Some definitions first */
  unsigned int iseed = (unsigned int)time(NULL);
  srand (iseed);
  UINT4 Ntrials = 500000;
  REAL4 distance,psi,beta,lamda,theta,varphi,phi;
  REAL4 cphi,sphi,aresp,bresp,Qcomp[9],Scomp[5],Tcomp[5],Pcomp[vecLengthTwo];
  REAL4 Pscaled[vecLengthTwo];
  REAL4 nsPcomp[nsVecLengthTwo],nsPscaled[nsVecLengthTwo];
  REAL4 c2varphi,s2varphi,svarphi,cvarphi;
  REAL4 c2theta,ctheta,s2theta,stheta,calpha,salpha;
  REAL4 numPi = 3.14159265;
  REAL4 pMMMp,pMMNp,pNMNp;
  REAL4 AdotA,AdotB,BdotB;
  REAL4 normScale;
  REAL4 ptfSNR,normSNR;
  UINT4 ptfCount=0;
  UINT4 normCount=0;
  UINT4 passCheck = 0;

  normScale = gsl_matrix_get(B2,0,0)/max_eigenvalue;
/*  FILE *tempFP;
  tempFP = fopen("ccode_output.dat","w");*/

  for (i = 0; i < Ntrials; i++)
  {
    distance = pow(rand()/(float)RAND_MAX,(1./3.))*0.2 ;
    if (! singleDetector)
      distance *= 2.;
    psi = rand()/(float)RAND_MAX * 2 * numPi;
    beta = acos(rand()/(float)RAND_MAX*2 - 1);
    lamda = rand()/(float)RAND_MAX * 2 * numPi;
    theta = acos(rand()/(float)RAND_MAX*2 - 1);
    varphi = rand()/(float)RAND_MAX * 2 * numPi;
    phi = rand()/(float)RAND_MAX * 2 * numPi;
    cphi = cos(phi);
    sphi = sin(phi);
    cvarphi = cos(varphi);
    c2varphi = cos(2*varphi);
    svarphi = sin(varphi);
    s2varphi = sin(2*varphi);
    ctheta = cos(theta);
    stheta = sin(theta);
    c2theta = cos(2*theta);
    s2theta = sin(2*theta);
    calpha = cos(psi);
    salpha = sin(psi);
    aresp = 0.5 * (1. + pow(cos(beta),2.))*cos(2.*lamda);
    bresp = cos(beta) * sin(2.*lamda);
    Qcomp[0] = -0.25*c2varphi*(3+c2theta);
    Qcomp[1] = ctheta*s2varphi;
    Qcomp[2] = -0.25 * s2varphi * (3+c2theta);
    Qcomp[3] = - ctheta*c2varphi;
    Qcomp[4] = 0.5 * cvarphi * s2theta;
    Qcomp[5] = - stheta*svarphi;
    Qcomp[6] = 0.5*svarphi*s2theta;
    Qcomp[7] = stheta*cvarphi;
    Qcomp[8] = pow(3,0.5) * 0.25 * (1- c2theta);
    Scomp[0] = Qcomp[0]*calpha + Qcomp[1]*salpha;
    Scomp[1] = Qcomp[2]*calpha + Qcomp[3]*salpha;
    Scomp[2] = Qcomp[4]*calpha + Qcomp[5]*salpha;
    Scomp[3] = Qcomp[6]*calpha + Qcomp[7]*salpha;
    Scomp[4] = Qcomp[8]*calpha;
    Tcomp[0] = -Qcomp[0]*salpha + Qcomp[1]*calpha;
    Tcomp[1] = -Qcomp[2]*salpha + Qcomp[3]*calpha;
    Tcomp[2] = -Qcomp[4]*salpha + Qcomp[5]*calpha;
    Tcomp[3] = -Qcomp[6]*salpha + Qcomp[7]*calpha;
    Tcomp[4] = -Qcomp[8]*salpha;
    if (! singleDetector)
    {
      for (j = 0; j < 10; j++ )
      {
        Pscaled[j] = 0;
        if (j < 5)
          Pcomp[j] = Scomp[j]/distance;
        else
          Pcomp[j] = Tcomp[j-5]/distance;
      }
      for (j = 0; j < 4; j++ )
      {
        nsPscaled[j] = 0;
        if (j < 2)
          nsPcomp[j] = Scomp[j]/distance;
        else
          nsPcomp[j] = Tcomp[j-2]/distance;
      }
    }
    else
    {
      for (j = 0; j < 5; j++ )
      {
        Pcomp[j] = (aresp * Scomp[j] + bresp * Tcomp[j])/distance;
        Pscaled[j] = 0;
      }
    }
    for (j = 0; j < vecLengthTwo; j++ )
    {
      for (k = 0; k < vecLengthTwo; k++ )
      {
        Pscaled[j] += gsl_matrix_get(eigenvecs,k,j)*Pcomp[k];
      }
    }
    for (j = 0; j < nsVecLengthTwo; j++ )
    {
      for (k = 0; k < nsVecLengthTwo; k++ )
      {
        nsPscaled[j] += gsl_matrix_get(nsEigenvecs,k,j)*nsPcomp[k];
      }
    }
    pMMMp = 0;
    pMMNp = 0;
    pNMNp = 0;
    for (j = 0; j < vecLengthTwo; j++ )
    {
      for (k = 0; k < vecLengthTwo; k++ )
      {
        pMMMp += gsl_matrix_get(MMM,j,k)*Pscaled[j]*Pscaled[k];
        pMMNp += gsl_matrix_get(MMN,j,k)*Pscaled[j]*Pscaled[k];
        pNMNp += gsl_matrix_get(NMN,j,k)*Pscaled[j]*Pscaled[k];
      }
    }

    /* Calculate dot products */
    AdotA = cphi*cphi * pMMMp - 2.*sphi*cphi*pMMNp;
    AdotA -= sphi*sphi*pNMNp;
    AdotB = cphi*cphi*pMMNp + cphi*sphi*pMMMp;
    AdotB += sphi*cphi*pNMNp - sphi*sphi*pMMNp;
    BdotB = - cphi*cphi*pNMNp + 2* cphi*sphi*pMMNp;
    BdotB += sphi*sphi*pMMMp;

    ptfSNR = AdotA + BdotB + pow(pow(AdotA-BdotB,2) + 4*AdotB*AdotB,0.5);
    ptfSNR *= 0.5;


    if (singleDetector)
      normSNR = (Pcomp[0]*Pcomp[0]+Pcomp[1]*Pcomp[1])*normScale;
    else
    {
      pMMMp = 0;
      pMMNp = 0;
      pNMNp = 0;
      for (j = 0; j < nsVecLengthTwo; j++ )
      {
        for (k = 0; k < nsVecLengthTwo; k++ )
        {
          pMMMp += gsl_matrix_get(nsMMM,j,k)*nsPscaled[j]*nsPscaled[k];
          pMMNp += gsl_matrix_get(nsMMN,j,k)*nsPscaled[j]*nsPscaled[k];
          pNMNp += gsl_matrix_get(nsNMN,j,k)*nsPscaled[j]*nsPscaled[k];
        }
      }

      /* Calculate dot products */
      AdotA = cphi*cphi * pMMMp - 2.*sphi*cphi*pMMNp;
      AdotA -= sphi*sphi*pNMNp;
      AdotB = cphi*cphi*pMMNp + cphi*sphi*pMMMp;
      AdotB += sphi*cphi*pNMNp - sphi*sphi*pMMNp;
      BdotB = - cphi*cphi*pNMNp + 2* cphi*sphi*pMMNp;
      BdotB += sphi*sphi*pMMMp;

      normSNR = AdotA + BdotB + pow(pow(AdotA-BdotB,2) + 4*AdotB*AdotB,0.5);
      normSNR *= 0.5;
    }

    /* Decide if bigger than thresholds */
    if (ptfSNR > params->spinSNR2threshold)
      ptfCount++;
    if (normSNR > params->nonspinSNR2threshold)
      normCount++;

  }

  if (ptfCount > normCount)
    passCheck = 1;

  gsl_matrix_free(B2);
  gsl_matrix_free(N2);
  gsl_matrix_free(tempM);
  gsl_matrix_free(b2);
  gsl_matrix_free(n2);
  gsl_matrix_free(Binv2);
  gsl_eigen_symmv_free(matTemp);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
  gsl_matrix_free(MMM);
  gsl_matrix_free(MMN);
  gsl_matrix_free(NMN);
  gsl_matrix_free(nsB2);
  gsl_matrix_free(nsN2);
  gsl_matrix_free(nsTempM);
  gsl_matrix_free(nsb2);
  gsl_matrix_free(nsn2);
  gsl_eigen_symmv_free(nsMatTemp);
  gsl_matrix_free(nsEigenvecs);
  gsl_vector_free(nsEigenvals);
  gsl_matrix_free(nsMMM);
  gsl_matrix_free(nsMMN);
  gsl_matrix_free(nsNMN);


  return passCheck;
}
