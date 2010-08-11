/* functions to create the likelihood for a pulsar search to be used with the
LALInference tools */

#define NUM_FACT 7
static const REAL8 inv_fact[NUM_FACT] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0),
(1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* a function to sum over the data */
REAL8Vector * sum_data(COMPLEX16Vect *data, UINT4Vector *chunkLength){
  INT4 chunkLength=0, length=0, i=0, j=0, count=0;
  COMPLEX16 B;
  REAL8Vector *sumData=NULL; 
  
  length = data->length + 1 - chunkLengths->data[chunkLengths->length-1];

  sumData = XLALCreateVector( length );
  
  for( i = 0 ; i < length ; i+= chunkLength ){
    chunkLength = chunkLengths->data[count];
    sumData->data[count] = 0.;

    for( j = i ; j < i + chunkLength ; j++){
      B.re = data->data[j].re;
      B.im = data->data[j].im;

      /* sum up the data */
      sumData->data[count] += (B.re*B.re + B.im*B.im);
    }

    count++;
  }
  
  sumData = XLALResizeREAL8Vector( sumData, count );
  
  return sumData;
}

/* detector response lookup table function  - this function will output a lookup
table of points in time and psi, covering a sidereal day from the start time
(t0) and from -pi/4 to pi/4 in psi */
void response_lookup_table(REAL8 t0, LALDetAndSource detAndSource,
  INT4 timeSteps, INT4 psiSteps, gsl_matrix *LUfplus, gsl_matrix *LUfcross){ 
  LIGOTimeGPS gps;
  REAL8 T=0;

  REAL8 fplus=0., fcross=0.;
  REAL8 psteps = (REAL8)lookupTable->psiSteps;
  REAL8 tsteps = (REAL8)lookupTable->timeSteps;

  INT4 i=0, j=0;
  
  for( i = 0 ; i < psiSteps ; i++ ){
    detAndSource.pSource->orientation = -(LAL_PI/4.) +
        (REAL8)i*(LAL_PI/2.) / ( psteps - 1. );

    for( j = 0 ; j < timeSteps ; j++ ){
      T = t0 + (REAL8)j*LAL_DAYSID_SI / tsteps;

      gps = XLALGPSSetREAL8(&gps, T);

      XLALComputeDetAMResponse(&fplus, &fcross,
        detAndSource.pDetector->response,
        detAndSource.pSource->equatorialCoords.longitude,
        detAndSource.pSource->equatorialCoords.latitude,
        detAndSource.pSource->orientation,
        XLALGreenwichMeanSiderealTime(&gps));
        
      gsl_matrix_set(LUfplus, i, j, fplus);
      gsl_matrix_set(LUfcross, i, j, fcross);
    }
  }
}


/* function to calculate the log(likelihood) given some data and a set of
   particular pulsar parameters */
LALLikelihoodFunction pulsar_log_likelihood( LALVariables *vars, 
  LALIFOData *data, LALTemplateFunction *get_pulsar_model ){
  INT4 i=0, j=0, count=0, k=0, cl=0;
  INT4 length=0;
  REAL8 chunkLength=0.;

  REAL8 tstart=0., T=0.;

  COMPLEX16 model;
  INT4 psibin=0, timebin=0;

  REAL8 plus=0., cross=0.;
  REAL8 sumModel=0., sumDataModel=0.;
  REAL8 chiSquare=0.;
  COMPLEX16 B, M;

  REAL8 exclamation[data.chunkMax+1]; /* all factorials up to chunkMax */
  REAL8 logOf2=log(2.);

  REAL8 loglike=0.; /* the log likelihood */

  INT4 first=0, through=0;

  REAL8 phiIni=0., phi=0.;

  /* copy model parameters to data parameters */
  data.modelParams = vars;
  
  /* get pulsar model */
  get_pulsar_model( data );
  
  /* to save time get all log factorials up to chunkMax */
  for( i = 0 ; i < data.chunkMax+1 ; i++ )
    exclamation[i] = log_factorial(i);
  
  for( i = 0 ; i < length ; i += chunkLength ){
    chunkLength = (REAL8)data->chunkLengths->data[count];

    if( chunkLength < data.chunkMin ){
      count++;

      if( through == 0 ) first = 0;

      continue;
    }

    through = 1;

    sumModel = 0.;
    sumDataModel = 0.;

    cl = i + (INT4)chunkLength;
    
    for( j = i ; j < cl ; j++){
      B.re = data->compTimeData->data->data[j].re;
      B.im = data->compTimeData->data->data[j].im;

      M.re = data->compModelData->data->data[j].re;
      M.im = data->compModelData->data->data[j].im;
      
      /* sum over the model */
      sumModel += M.re*M.re + M.im*M.im;

      /* sum over that data and model */
      sumDataModel += B.re*M.re + B.im*M.im;
    }

    chiSquare = data->sumData->data[count];
    chiSquare -= 2.*sumDataModel;
    chiSquare += sumModel;

    if( first == 0 ){
      loglike = (chunkLength - 1.)*logOf2;
      first++;
    }
    else loglike += (chunkLength - 1.)*logOf2;

    loglike += exclamation[(INT4)chunkLength];
    loglike -= chunkLength*log(chiSquare);

    count++;
  }

  return loglike;
}

/* function to creates a vector of the phase differences between that used for
an initial set of phase parameters (i.e. those used to heterodyne the data),
and a new set of phase parameters */
void get_pulsar_model( LALIFOData *data ){
  INT4 i=0, length=0;
  
  BinaryPulsarParams pars;
  
  REAL8Vector *dphi=NULL;
  REAL8Vector *amp=NULL;
  
  /* set model parameters */
  pars.h0 = *(REAL8*)getVariable( data->modelParams, "h0");
  pars.cosiota = *(REAL8*)getVariable( data->modelParams, "cosiota");
  pars.psi = *(REAL8*)getVariable( data->modelParams, "psi");
  pars.phi0 = *(REAL8*)getVariable( data->modelParams, "phi0");
  
  /* set the potentially variable parameters */
  pars.pepoch = *(REAL8*)getVariable( data->modelParams, "pepoch");
  pars.posepoch = *(REAL8*)getVariable( data->modelParams, "posepoch");
  
  pars.ra = *(REAL8*)getVariable( data->modelParams, "ra");
  pars.pmra = *(REAL8*)getVariable( data->modelParams, "pmra");
  pars.dec = *(REAL8*)getVariable( data->modelParams, "dec");
  pars.pmdec = *(REAL8*)getVariable( data->modelParams, "pmdec");
  
  pars.f0 = *(REAL8*)getVariable( data->modelParams, "f0");
  pars.f1 = *(REAL8*)getVariable( data->modelParams, "f1");
  pars.f2 = *(REAL8*)getVariable( data->modelParams, "f2");
  pars.f3 = *(REAL8*)getVariable( data->modelParams, "f3");
  pars.f4 = *(REAL8*)getVariable( data->modelParams, "f4");
  pars.f5 = *(REAL8*)getVariable( data->modelParams, "f5");
  
  pars.model = *(CHAR*)getVariable( data->modelParams, "model");
  
  /* binary parameters */
  if( pars.model != NULL ){
    pars.e = *(REAL8*)getVariable( data->modelParams, "e");
    pars.w0 = *(REAL8*)getVariable( data->modelParams, "w0");
    pars.Pb = *(REAL8*)getVariable( data->modelParams, "Pb");
    pars.x = *(REAL8*)getVariable( data->modelParams, "x");
    pars.T0 = *(REAL8*)getVariable( data->modelParams, "T0");
    
    pars.e2 = *(REAL8*)getVariable( data->modelParams, "e2");
    pars.w02 = *(REAL8*)getVariable( data->modelParams, "w02");
    pars.Pb2 = *(REAL8*)getVariable( data->modelParams, "Pb2");
    pars.x2 = *(REAL8*)getVariable( data->modelParams, "x2");
    pars.T02 = *(REAL8*)getVariable( data->modelParams, "T02");
    
    pars.e3 = *(REAL8*)getVariable( data->modelParams, "e3");
    pars.w03 = *(REAL8*)getVariable( data->modelParams, "w03");
    pars.Pb3 = *(REAL8*)getVariable( data->modelParams, "Pb3");
    pars.x3 = *(REAL8*)getVariable( data->modelParams, "x3");
    pars.T03 = *(REAL8*)getVariable( data->modelParams, "T03");
    
    pars.xpbdot = *(REAL8*)getVariable( data->modelParams, "xpbdot");
    
    pars.eps1 = *(REAL8*)getVariable( data->modelParams, "eps1");    
    pars.eps2 = *(REAL8*)getVariable( data->modelParams, "eps2");       
    pars.eps1dot = *(REAL8*)getVariable( data->modelParams, "eps1dot");
    pars.eps2dot = *(REAL8*)getVariable( data->modelParams, "eps2dot");
    pars.Tasc = *(REAL8*)getVariable( data->modelParams, "Tasc");
    
    pars.wdot = *(REAL8*)getVariable( data->modelParams, "wdot"); 
    pars.gamma = *(REAL8*)getVariable( data->modelParams, "gamma");
    pars.Pbdot = *(REAL8*)getVariable( data->modelParams, "Pbdot");  
    pars.xdot = *(REAL8*)getVariable( data->modelParams, "xdot");   
    pars.edot = *(REAL8*)getVariable( data->modelParams, "edot");
    
    pars.s = *(REAL8*)getVariable( data->modelParams, "s"); 
    pars.dr = *(REAL8*)getVariable( data->modelParams, "dr");
    pars.dth = *(REAL8*)getVariable( data->modelParams, "dth");   
    pars.a0 = *(REAL8*)getVariable( data->modelParams, "a0");
    pars.b0 = *(REAL8*)getVariable( data->modelParams, "b0"); 

    pars.M = *(REAL8*)getVariable( data->modelParams, "M"); 
    pars.m2 = *(REAL8*)getVariable( data->modelParams, "m2");
  }

  amp = get_amplitude_model( pars, data );
  length = amp->length;
  
  /* assume that timeData vector within the LALIFOData structure contains the
     phase calculated using the initial (heterodyne) values of the phase
     parameters */
   
  /* get difference in phase and perform extra heterodyne with it */ 
  if ( (dphi = get_phase_model( pars, data )) != NULL){
    for( i=0; i<length; i++ ){
      COMPLEX16 M;
      REAL8 dphit;
      REAL4 sp, cp;
    
      dphit = -fmod(dphi->data[i] - data->timeData->data[i], 1.);
    
      sin_cos_2PI_LUT( &sp, &cp, dphit );
    
      M.re = data->compModelData->data->data[i].re;
      M.im = data->compModelData->data->data[i].im;
    
      /* heterodyne */
      data->compModelData->data->data[i].re = M.re*cp - M.im*sp;
      data->compModelData->data->data[i].im = M.im*cp + M.re*sp;
    }
  }
}

REAL8Vector *get_phase_model( BinaryPulsarParams params, LALIFOData *data ){
  static LALStatus status;

  INT4 i=0, length;

  REAL8 T0=0., DT=0., DTplus=0., deltat=0., deltat2=0.;
  REAL8 interptime = 1800.; /* calulate every 30 mins (1800 secs) */

  EarthState earth, earth2;
  EmissionTime emit, emit2;
  REAL8 emitdt=0.;

  BinaryPulsarInput binput;
  BinaryPulsarOutput boutput;

  BarycenterInput bary;
  EphemerisData edat;
  
  REAL8Vector *phis=NULL;

  /* if edat is NULL then return a NULL poniter */
  if( data->ephem == NULL )
    return NULL;
  
  /* copy barycenter and ephemeris data */
  memcpy(bary, data->bary, sizeof(data->bary));
  memcpy(edat, data->ephem, sizeof(data->ephem));

   /* set the position and frequency epochs if not already set */
  if(params.pepoch == 0. && params.posepoch != 0.)
    params.pepoch = params.posepoch;
  else if(params.posepoch == 0. && params.pepoch != 0.)
    params.posepoch = params.pepoch;

  length = data->dataTimes.length;
  
  /* allocate memory for phases */
  phis = XLALCreateREAL8Vector( length );

  /* set 1/distance if parallax or distance value is given (1/sec) */
  if( params.px != 0. )
    bary.dInv = params.px*1e-3*LAL_C_SI/LAL_PC_SI;
  else if( params.dist != 0. )
    bary.dInv = LAL_C_SI/(params.dist*1e3*LAL_PC_SI);
  else
    bary.dInv = 0.;

  for( i=0; i<length; i++){
    REAL8 realT = XLALGPSGetREAL8(data->dataTimes->data[i]);
    
    T0 = params.pepoch;

    DT = realT - T0;

    /* only do call the barycentring routines every 30 minutes, otherwise just
       linearly interpolate between them */
    if( i==0 || DT > DTplus ){
      bary.tgps.gpsSeconds = data->dataTimes->data[i].gpsSeconds;
      bary.tgps.gpsNanoSeconds = data->dataTimes->data[i].gpsNanoSeconds;

      bary.delta = params.dec + (realT-params.posepoch) * params.pmdec;
      bary.alpha = params.ra + (realT-params.posepoch) *
         params.pmra/cos(bary.delta);

      /* call barycentring routines */
      LAL_CALL( LALBarycenterEarth(&status, &earth, &bary.tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit, &bary, &earth), &status );

      /* add interptime to the time */
      DTplus = DT + interptime;
      bary.tgps = XLALGPSAdd(bary.tgps, interptime);

      /* No point in updating the positions as difference will be tiny */
      LAL_CALL( LALBarycenterEarth(&status, &earth2, &bary.tgps, edat),
        &status );
      LAL_CALL( LALBarycenter(&status, &emit2, &bary, &earth2), &status );
    }

    /* linearly interpolate to get emitdt */
    emitdt = emit.deltaT + (DT - (DTplus - interptime)) *
      (emit2.deltaT - emit.deltaT)/interptime;

    /* check if need to perform binary barycentring */
    if( params.model != NULL ){
      binput.tb = realT + emitdt;

      XLALBinaryPulsarDeltaT( &boutput, &binput, &params );

      deltat = DT + emitdt + boutput.deltaT;
    }
    else
      deltat = DT + emitdt;

    /* work out phase */
    deltat2 = deltat*deltat;
    phis->data[i] = 2.*deltat*(params.f0 + 
      inv_fact[2]*params.f1*deltat +
      inv_fact[3]*params.f2*deltat2 +
      inv_fact[4]*params.f3*deltat*deltat2 +
      inv_fact[5]*params.f4*deltat2*deltat2 +
      inv_fact[6]*params.f5*deltat2*deltat2*deltat);
  }

  return phis;
}

REAL8Vector *get_amplitude_model( BinaryPulsarParams pars, LALIFOData *data ){
  INT4 i=0, length;
  
  REAL8Vector *amp=NULL;
  
  REAL8 psteps, tsteps;
  INT4 psibin, timebin;
  REAL8 tstart;
  REAL8 plus, cross;

  REAL8 Xplus, Xcross;
  REAL8 Xpcosphi_2, Xccosphi_2, Xpsinphi_2, Xcsinphi_2;
  REAL8 sinphi, cosphi;
  
  length = data->dataTimes.length;
  
  /* set lookup table parameters */
  psteps = *(INT4*)getVariable( data->modelParams, "psiSteps" );
  tsteps = *(INT4*)getVariable( data->modelParams, "timeSteps" );
  
  /* allocate memory for amplitudes */
  amp = XLALCreateREAL8Vector( length );
  
  sin_cos_LUT( &sinphi, &cosphi, pars.phi0 );
  
  /************************* CREATE MODEL *************************************/
  /* This model is a complex heterodyned time series for a triaxial neutron
     star emitting at twice its rotation frequency (as defined in Dupuis and
     Woan, PRD, 2005):
       real = (h0/2) * ((1/2)*F+*(1+cos(iota)^2)*cos(phi0) 
         + Fx*cos(iota)*sin(phi0))
       imag = (h0/2) * ((1/2)*F+*(1+cos(iota)^2)*sin(phi0)
         - Fx*cos(iota)*cos(phi0))
   ****************************************************************************/
  
  Xplus = 0.25*(1.+pars.cosiota*pars.cosiota)*pars.h0;
  Xcross = 0.5*pars.cosiota*pars.h0;
  Xpsinphi = Xplus*sinphi;
  Xcsinphi = Xcross*sinphi;
  Xpcosphi = Xplus*cosphi;
  Xccosphi = Xcross*cosphi;
 
  for( i=0; i<length; i++ ){
    /* set the psi bin for the lookup table */
    psibin = (INT4)ROUND( ( pars.psi + LAL_PI/4. ) * ( psteps-1. )/LAL_PI_2 );

    tstart = XLALGPSGetREAL8( data->dataTimes->data[0] ); /*time of first B_k*/
  
    /* set the time bin for the lookup table */
    /* sidereal day in secs*/
    T = fmod(XLALGPSGetREAL8(data->dataTimes->data[i]) - tstart,
      LAL_DAYSID_SI);
    timebin = (INT4)fmod( ROUND(T*tsteps/LAL_DAYSID_SI), tsteps );

    plus = gsl_matrix_get(data->lookupTablePlus, psibin, timebin);
    cross = gsl_matrix_get(data->lookupTableCross, psibin, timebin);
    
    /* create the complex signal amplitude model */
    data->compModelData->data->data[i].re = plus*Xpcosphi + cross*Xcsinphi;
    data->compModelData->data->data[i].im = plus*Xpsinphi - cross*Xccosphi;
  }
}

void add_initial_variables( LALVariables *ini, BinaryPulsarParams pars ){
  /* amplitude model parameters */
  addVariable(&ini, "h0", &pars.h0, REAL8_t);
  addVariable(&ini, "phi0", pars.phi0, REAL8_t);
  addVariable(&ini, "cosiota", pars.cosiota, REAL8_t);
  addVariable(&ini, "psi", pars.psi, REAL8_t);
  
  /* phase model parameters */
  
  /* frequency */
  addVariable(&ini, "f0", &pars.f0, REAL8_t);
  addVariable(&ini, "f1", &pars.f1, REAL8_t);
  addVariable(&ini, "f2", &pars.f2, REAL8_t);
  addVariable(&ini, "f3", &pars.f3, REAL8_t);
  addVariable(&ini, "f4", &pars.f4, REAL8_t);
  addVariable(&ini, "f5", &pars.f5, REAL8_t);
  addVariable(&ini, "pepoch", &pars.pepoch, REAL8_t);
  
  /* sky position */
  addVariable(&ini, "ra", &pars.ra, REAL8_t);
  addVariable(&ini, "pmra", &pars.pmra, REAL8_t);
  addVariable(&ini, "dec", &pars.dec, REAL8_t);
  addVariable(&ini, "pmdec", &pars.pmdec, REAL8_t);
  addVariable(&ini, "posepoch", &pars.posepoch, REAL8_t);
  
  /* binary system parameters */
  addVariable(&ini, "model", &pars.model, string);
  
  addVariable(&ini, "Pb", &pars.Pb, REAL8_t);
  addVariable(&ini, "e", &pars.e, REAL8_t);
  addVariable(&ini, "eps1", &pars.eps1, REAL8_t);
  addVariable(&ini, "eps2", &pars.eps2, REAL8_t);
  addVariable(&ini, "T0", &pars.T0, REAL8_t);
  addVariable(&ini, "Tasc", &pars.Tasc, REAL8_t);
  addVariable(&ini, "x", &pars.x, REAL8_t);
  addVariable(&ini, "w0", &pars.w0, REAL8_t);

  addVariable(&ini, "Pb2", &pars.Pb2, REAL8_t);
  addVariable(&ini, "e2", &pars.e2, REAL8_t);
  addVariable(&ini, "T02", &pars.T02, REAL8_t);
  addVariable(&ini, "x2", &pars.x2, REAL8_t);
  addVariable(&ini, "w02", &pars.w02, REAL8_t);
  
  addVariable(&ini, "Pb3", &pars.Pb3, REAL8_t);
  addVariable(&ini, "e3", &pars.e3, REAL8_t);
  addVariable(&ini, "T03", &pars.T03, REAL8_t);
  addVariable(&ini, "x3", &pars.x3, REAL8_t);
  addVariable(&ini, "w03", &pars.w03, REAL8_t);
  
  addVariable(&ini, "xpbdot", &pars.xpbdot, REAL8_t);
  addVariable(&ini, "eps1dot", &pars.eps1dot, REAL8_t);
  addVariable(&ini, "eps2dot", &pars.eps2dot, REAL8_t);
  addVariable(&ini, "wdot", &pars.wdot, REAL8_t);
  addVariable(&ini, "gamma", &pars.gamma, REAL8_t);
  addVariable(&ini, "Pbdot", &pars.Pbdot, REAL8_t);
  addVariable(&ini, "xdot", &pars.xdot, REAL8_t);
  addVariable(&ini, "edot", &pars.edot, REAL8_t);
  
  addVariable(&ini, "s", &pars.s, REAL8_t);
  addVariable(&ini, "dr", &pars.dr, REAL8_t);
  addVariable(&ini, "dth", &pars.dth, REAL8_t);
  addVariable(&ini, "a0", &pars.a0, REAL8_t);
  addVariable(&ini, "b0", &pars.b0, REAL8_t);
  addVariable(&ini, "M", &pars.M, REAL8_t);
  addVariable(&ini, "m2", &pars.m2, REAL8_t);
}

/* things I need to pass to the likelihood function via the LALInferenceRunState
  - 
  - a REAL8 value giving the initial time stamp for the above detector response 


/* FOR REFERENCE - using LIGOTimeGPSVector requires the
XLALCreateTimestampVector function in SFTutils */

/* functions to read in, or set, specific entries from an array */
/* REAL8 get_array_value( REAL8Array *array, UINT4Vector *entry ){
  UINT4 idx=0;
  UINT4 i=0; */
  
  /* check that the length of the vector containin the required array entry
and the array dimension length are the same */
  /* if ( array->dimLength->length != entry->length ){
    fprintf(stderr, "Error... entries not the same length!\n");
    exit(3);
  } */
  
  /* check that required entries are within the dimesion lengths */
  /* for( i=0; i<entry->length; i++ ){
    if (entry->data[i] > array->dimLength->data[i]){
      fprintf(stderr, "Error... entries greater than maximum value!\n");
      exit(4);
    }
  } */
  
  /* get the index of the entry in the array */
/*  for ( i=0; i<entry->length - 1; i++ )
    idx += entry->data[i] * array->dimLength->data[0];
  
  idx += entry->data[i];
  
  return array->data[idx];
} */

/* void set_array_value( REAL8Array *array, UINT4Vector *entry, REAL8 val ){
  UINT4 idx=0;
  UINT4 i=0; */
  
  /* check that the length of the vector containin the required array entry
and the array dimension length are the same */
/*  if ( array->dimLength->length != entry->length ){
    fprintf(stderr, "Error... entries not the same length!\n");
    exit(3);
  } */
  
  /* check that required entries are within the dimesion lengths */
  /*for( i=0; i<entry->length; i++ ){
    if (entry->data[i] > array->dimLength->data[i]){
      fprintf(stderr, "Error... entries greater than maximum value!\n");
      exit(4);
    }
  }*/
  
  /* get the index of the entry in the array */
  /*for ( i=0; i<entry->length - 1; i++ )
    idx += entry->data[i] * array->dimLength->data[0];
  
  idx += entry->data[i];
  
  array->data[idx] = val;
}*/