#include <lal/CoherentEstimation.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/AVFactories.h>
#include <strings.h>
#include <math.h>

#define EPS 1.0e-7

#define Cosh(x) cosh(x)
#define ACosh(x) acosh(x)

NRCSID( COHERENTESTIMATIONC, "$Id$" );

double cosh(double);
double acosh(double);

static INT4 jacobi(float **a, int n, float d[], float **v, int *nrot);

void
LALCoherentEstimation ( LALStatus          *stat,
			REAL4TimeSeries *output,
			CoherentEstimation *params,
			DetectorsData      *in) {
  /* 
     NOTES:
     o destroys input (in)
     o order of in must be same as order of params->filters
     o output time wrt center of Earth
  */

  INT4 /* Sret, */
    i, j, k, /* counters */
    iPad, ePad, /* used for padding */
    del;  /* integer time delay */
  REAL4 y; /* dummy */
  REAL8 p1,p2; /* weights in interpolation */
  REAL8 *tDelays; /* time delays wrt Earth center */
  REAL8 mDelay; /* max of time delays */
  REAL8 tmid; /* mid time point */
  LALDetAMResponse *F; /* response functions */
  DetTimeAndASource dtS; /* for time delays */
  LALPlaceAndGPS pGPS; /* for time delays */
  LALDetAndSource dAs; /* for responses */
  LALSource source; /* for responses */
  LIGOTimeGPS t0; /* start time of data */

  REAL8 *alpha; /* scale factor */
  INT4 nrot;
  REAL8 maxLambda, tmpLambda;
  REAL4 *lambda;
  REAL4 **Hmat;
  REAL4 **v;
  REAL8 **Cmat;

  /*
    {REAL8 tmp = Cosh(ACosh((double)3.02993)/3.0);
    printf("%g\n",tmp);}
  */

  /***********************************************************************/
  /* initialize status & validate input                                  */
  /***********************************************************************/
  INITSTATUS( stat, "LALCoherentEstimation", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );


  ASSERT ( in, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( in->data, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( in->Ndetectors > 0 && in->Ndetectors < 10, stat, COHERENTESTIMATIONH_E0DEC, COHERENTESTIMATIONH_MSGE0DEC );

  ASSERT ( in->Ndetectors == 3, stat, COHERENTESTIMATIONH_EUIMP, COHERENTESTIMATIONH_MSGEUIMP );

  for(i=0;i<(INT4)in->Ndetectors;i++) {
    ASSERT ( in->data[i].data, stat, COHERENTESTIMATIONH_E0DEC, COHERENTESTIMATIONH_MSGE0DEC );
  }

  t0 = in->data[0].epoch;

  for(i=1;i<(INT4)in->Ndetectors;i++) {
    ASSERT ( in->data[i].epoch.gpsSeconds == t0.gpsSeconds &&
	     in->data[i].epoch.gpsNanoSeconds == t0.gpsNanoSeconds, stat, COHERENTESTIMATIONH_EDST, COHERENTESTIMATIONH_MSGEDST );
  }

  ASSERT ( params, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->detectors, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->filters, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->position, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->Ndetectors == in->Ndetectors, stat, COHERENTESTIMATIONH_EICE, COHERENTESTIMATIONH_MSGEICE );

  /***********************************************************************/
  /* if params->preProcessed is 0, need to pre-processed by applying the */
  /* IIR filters twice                                                   */
  /***********************************************************************/
  if(!(params->preProcessed) && params->nPreProcessed) {
  
    REAL8Vector *tmpR8 = NULL;
    UINT4 jind;

    for(i=0;i<(INT4)in->Ndetectors;i++) {
      ASSERT ( params->filters[i], stat, COHERENTESTIMATIONH_EICE, COHERENTESTIMATIONH_MSGEICE );
    }

    /* loop over all detectors, and filter twice */
    for(i=0; i<(INT4)params->Ndetectors; i++) {

      TRY ( LALDCreateVector( stat->statusPtr,
			      &tmpR8,
			      in->data[i].data->length ), stat );

      for(jind = 0; jind<in->data[i].data->length; jind++) {
	tmpR8->data[jind] = (REAL8)(in->data[i].data->data[jind]);
      }
     
      for(j=0; j<params->nPreProcessed; j++) { 
        TRY ( LALIIRFilterREAL8Vector( stat->statusPtr,
	       			       tmpR8,
				       params->filters[i]), stat );
      }

      for(jind = 0; jind<in->data[i].data->length; jind++) {
	in->data[i].data->data[jind] = (REAL4)(tmpR8->data[jind]);
      }

      TRY ( LALDDestroyVector (  stat->statusPtr,
				 &tmpR8 ) , stat );

      /* set first 1/16 s to zero to avoid transient */
      bzero(in->data[i].data->data, sizeof(REAL4) * (UINT4)ceil(0.0635 / in->data[i].deltaT));

    }

    params->preProcessed = 1;

  }


  /* update output parameters */
  output->epoch = in->data[0].epoch;
  output->deltaT = in->data[0].deltaT;
  output->f0 = in->data[0].f0;
  memcpy(&(output->sampleUnits), &(in->data[0].sampleUnits), sizeof(LALUnit));

  /* make sure output is zero */
  bzero(output->data->data, output->data->length * sizeof(REAL4));



  /***********************************************************************/
  /* Compute time delays and response functions                          */
  /***********************************************************************/
  if(!(tDelays = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8)))) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }
  dtS.p_source = params->position;

  /* delays are computed wrt to center of data stretch */
  dtS.p_det_and_time = &pGPS;
  pGPS.p_gps = &t0;
  tmid = 0.5 * (REAL8)(in->data[0].data->length) * in->data[0].deltaT;
  pGPS.p_gps->gpsSeconds += (INT4)floor(tmid);
  pGPS.p_gps->gpsNanoSeconds += (INT4)floor(1E9*(tmid-floor(tmid)));
  if(pGPS.p_gps->gpsNanoSeconds >= 1000000000) {
    pGPS.p_gps->gpsSeconds += (INT4)floor((REAL8)(pGPS.p_gps->gpsNanoSeconds) / 1E9);
    pGPS.p_gps->gpsNanoSeconds -= 1000000000 * (INT4)floor((REAL8)(pGPS.p_gps->gpsNanoSeconds) / 1E9);
  }


  if(!(F = (LALDetAMResponse *)LALMalloc(params->Ndetectors * sizeof(LALDetAMResponse)))) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  dAs.pSource = &source;
  source.equatorialCoords = *(params->position);
  source.orientation = params->polAngle;

  for(i=0; i<(INT4)params->Ndetectors; i++) {

    LALGPSandAcc gpsAndAcc;

    pGPS.p_detector = dAs.pDetector = params->detectors + i;

    /* tDelays = arrival time at detector - arrival time a center of Earth */
    TRY ( LALTimeDelayFromEarthCenter( stat->statusPtr, 
				       tDelays + i, 
				       &dtS ), stat );


    /* JC: isnan is not allowed 
    if(isnan(tDelays[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }
    */

    gpsAndAcc.gps = *pGPS.p_gps;
    gpsAndAcc.accuracy = LALLEAPSEC_LOOSE; /* FIXME ??? */
    TRY ( LALComputeDetAMResponse ( stat->statusPtr, 
				    F + i,
				    &dAs,
				    &gpsAndAcc ), stat );

    /* JC: isnan is not allowed
    if(isnan(F[i].cross) || isnan(F[i].plus)) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }
    */
  }

  /***********************************************************************/
  /* Compute and store estimated data                                    */
  /***********************************************************************/
  /* set time origine on detector with largest delay */
  /*
  mDelay = tDelays[0];
  for(i=1; i<params->Ndetectors; i++) {
	if(tDelays[i] > mDelay) {
		mDelay = tDelays[i];
	}
  }
  */
  mDelay = tDelays[0];

  for(i=0; i<(INT4)params->Ndetectors; i++) {
    tDelays[i] -= mDelay;
  }


  alpha = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8));
  if(!alpha) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  lambda = (REAL4 *)LALMalloc(params->Ndetectors * sizeof(REAL4));
  if(!lambda) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  Hmat = (REAL4 **)LALMalloc(params->Ndetectors * sizeof(REAL4 *));
  if(!Hmat) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }
  for(i=0; i<(INT4)params->Ndetectors; i++) {
    Hmat[i] = (REAL4 *)LALMalloc(params->Ndetectors * sizeof(REAL4));
    if(!(Hmat[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
    }
  }

  v = (REAL4 **)LALMalloc(params->Ndetectors * sizeof(REAL4 *));
  if(!v) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }
  for(i=0; i<(INT4)params->Ndetectors; i++) {
    v[i] = (REAL4 *)LALMalloc(params->Ndetectors * sizeof(REAL4));
    if(!(v[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
    }
  }

  Cmat = params->CMat;

  for(i=0; i<(INT4)params->Ndetectors; i++) {
    for(j=0; j<(INT4)params->Ndetectors; j++) {

	Hmat[i][j] = (F[i].plus*F[j].plus*params->plus2cross + (F[i].plus*F[j].cross + F[j].plus*F[i].cross)*params->plusDotcross + F[i].cross*F[j].cross/params->plus2cross) / Cmat[i][i];

    }
  }

  if(jacobi(Hmat, params->Ndetectors, lambda, v, &nrot)) {
    ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
  }
  /*
 printf("%g\t%g\t%g\n",lambda[0],lambda[1],lambda[2]);
 */
  maxLambda = -1e30;
  for(k=0;k<(INT4)params->Ndetectors;k++) {

    tmpLambda = 0.0;

    for(i=0; i<(INT4)params->Ndetectors; i++) {
      for(j=0; j<(INT4)params->Ndetectors; j++) {
	tmpLambda += v[i][k]*v[j][k]*(F[i].plus*F[j].plus*params->plus2cross + (F[i].plus*F[j].cross + F[j].plus*F[i].cross)*params->plusDotcross + F[i].cross*F[j].cross/params->plus2cross);
      }
    }

    if(tmpLambda > maxLambda) {
      for(i=0; i<(INT4)params->Ndetectors; i++) {
	alpha[i] = v[i][k];
      }
      maxLambda = tmpLambda;
    }
/*
    printf("%u\t%g\t%g\t%g\t%g\n",k,lambda[k],alpha[0],alpha[1],alpha[2]);
    */
  }

  /* loop */
  for(i=0; i<(INT4)params->Ndetectors; i++) {

    /* setup padding and weights */
    if(tDelays[i] < 0.0) { 
      /* need padding at beginning */
      iPad = (INT4)floor(-tDelays[i]/output->deltaT);
      ePad = 0;

      /* set integer delay (for p1 weight) */
      del = -iPad;

      /* set weights */
      p1 = ceil(tDelays[i] / output->deltaT) - tDelays[i] / output->deltaT;
      p2 = 1.0 - p1;
    } else { 
      /* need padding at end */ 
      iPad = 0;
      ePad = (INT4)ceil(tDelays[i]/output->deltaT);

      /* integer delay */
      del = ePad;

      /* weights */
      p1 = ceil(tDelays[i] / output->deltaT) - tDelays[i] / output->deltaT;
      p2 = 1.0 - p1;
    }


    /* interpolate using time delays */
    for(j=iPad+1; j<(INT4)output->data->length - (INT4)ePad; j++) {

      y = p1 * in->data[i].data->data[del+j-1] + p2 * in->data[i].data->data[del+j];

      output->data->data[j] += y * alpha[i];
    }
  }

  /*
  {
    FILE *out;
    out = fopen("test.dat","w");
    for(j=0;j<output->data->length;j++) {
      fprintf(out,"%g\t%g\n",output->data->data[j].re,output->data->data[j].im);
    }
    fclose(out);

    exit(0);
  }
  */

  /***********************************************************************/
  /* clean up and normal return                                          */
  /***********************************************************************/
  LALFree(tDelays);
  LALFree(F);

  LALFree(alpha);
  LALFree(lambda);

  for(i=0; i<(INT4)params->Ndetectors; i++) {
    LALFree(Hmat[i]);
    LALFree(v[i]);
  }
  LALFree(Hmat);
  LALFree(v);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}


void 
LALClearCoherentData (
		      LALStatus     *stat,
		      DetectorsData *dat
		      ) {

  UINT4 i;
  
  INITSTATUS( stat, "LALClearCoherentData", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );

  if(!dat) {
    ABORT ( stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  }

  for(i=0; i<dat->Ndetectors; i++) {
    if(dat->data[i].data) {
      TRY ( LALDestroyVector(stat->statusPtr, &(dat->data[i].data)), stat );
    }
  }

  LALFree(dat->data);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}



void 
LALClearCoherentInfo (
		      LALStatus     *stat,
		      CoherentEstimation *dat
		      ) {

  UINT4 i;
  
  INITSTATUS( stat, "LALClearCoherentInfo", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );

  if(!dat) {
    ABORT ( stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  }

  LALFree(dat->detectors);

  for(i=0; i<dat->Ndetectors; i++) {
    if(dat->filters[i]) {
      TRY ( LALDestroyREAL8IIRFilter(stat->statusPtr, dat->filters + i), stat );
    }
    LALFree(dat->CMat[i]);
  }

  LALFree(dat->filters);

  LALFree(dat->CMat);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

int jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j,iq,ip,i/*,k*/;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b = (REAL4 *)LALMalloc(n * sizeof(REAL4));
	if(!b) {
	  return 1;
	}
	b -= 1;

	z = (REAL4 *)LALMalloc(n * sizeof(REAL4));
	if(!z) {
	  return 1;
	}
	z -= 1;

	d -= 1;

	for(i=0;i<n;i++) {
	  v[i] -= 1;
	  a[i] -= 1;
	}
	v -= 1;
	a -= 1;

	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
		  LALFree(z + 1);
		  LALFree(b + 1);
		  d += 1;
		  a += 1;
		  v += 1;
		  for(i=0;i<n;i++) {
		    v[i] += 1;
		    a[i] += 1;
		  }
		  return 0;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}

	fprintf(stderr,"Too many iterations in routine jacobi");
	return 1;
}
#undef ROTATE

