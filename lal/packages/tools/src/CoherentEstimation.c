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

  static INT4
    SolveHomogeneous(
		     REAL8 **Hmat,
		     REAL8 *tmpAlpha,
		     UINT4 N
		     );


  INT4 Sret,
    i, j, /* counters */
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

  REAL8 *alpha, *tmpAlpha; /* scale factor */
  UINT4 Nlambda;
  REAL8 a0, a1, a2, A0, A1, A2, A3, an, bn, cn, ab, ac, bc, p, q, C, maxLambda, tmpLambda;
  REAL8 *lambda;
  REAL8 **Hmat, **Cmat;

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

  for(i=0;i<in->Ndetectors;i++) {
    ASSERT ( in->data[i].data, stat, COHERENTESTIMATIONH_E0DEC, COHERENTESTIMATIONH_MSGE0DEC );
  }

  t0 = in->data[0].epoch;

  for(i=1;i<in->Ndetectors;i++) {
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

    for(i=0;i<in->Ndetectors;i++) {
      ASSERT ( params->filters[i], stat, COHERENTESTIMATIONH_EICE, COHERENTESTIMATIONH_MSGEICE );
    }

    /* loop over all detectors, and filter twice */
    for(i=0; i<params->Ndetectors; i++) {

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

  for(i=0; i<params->Ndetectors; i++) {

    pGPS.p_detector = dAs.pDetector = params->detectors + i;

    /* tDelays = arrival time at detector - arrival time a center of Earth */
    TRY ( LALTimeDelayFromEarthCenter( stat->statusPtr, 
				       tDelays + i, 
				       &dtS ), stat );


    if(isnan(tDelays[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }

    TRY ( LALComputeDetAMResponse ( stat->statusPtr, 
				    F + i,
				    &dAs,
				    pGPS.p_gps ), stat );

    if(isnan(F[i].cross) || isnan(F[i].plus)) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }
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

  for(i=0; i<params->Ndetectors; i++) {
    tDelays[i] -= mDelay;
  }


  alpha = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8));
  if(!alpha) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  tmpAlpha = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8));
  if(!tmpAlpha) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  lambda = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8));
  if(!lambda) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  Hmat = (REAL8 **)LALMalloc(params->Ndetectors * sizeof(REAL8 *));
  if(!Hmat) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }
  for(i=0; i<params->Ndetectors; i++) {
    Hmat[i] = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8));
    if(!(Hmat[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
    }
  }

  Cmat = params->CMat;

  an = F[0].plus*F[0].plus*params->plus2cross + F[0].cross*F[0].cross/params->plus2cross + 2.0*F[0].plus*F[0].cross*params->plusDotcross;
  bn = F[1].plus*F[1].plus*params->plus2cross + F[1].cross*F[1].cross/params->plus2cross + 2.0*F[1].plus*F[1].cross*params->plusDotcross;
  cn = F[2].plus*F[2].plus*params->plus2cross + F[2].cross*F[2].cross/params->plus2cross + 2.0*F[2].plus*F[2].cross*params->plusDotcross;

  ab = F[0].plus*F[1].plus*params->plus2cross + (F[0].plus*F[1].cross + F[0].cross*F[1].plus)*params->plusDotcross + F[0].cross*F[1].cross/params->plus2cross;
  ac = F[0].plus*F[2].plus*params->plus2cross + (F[0].plus*F[2].cross + F[0].cross*F[2].plus)*params->plusDotcross + F[0].cross*F[2].cross/params->plus2cross;
  bc = F[1].plus*F[2].plus*params->plus2cross + (F[1].plus*F[2].cross + F[1].cross*F[2].plus)*params->plusDotcross + F[1].cross*F[2].cross/params->plus2cross;

  A0 = -ac*ac*bn-bc*bc*an-ab*ab*cn+2.0*ab*ac*bc+an*bn*cn;
  A1 = -bc*bc*Cmat[0][0]+bn*cn*Cmat[0][0]+2.0*ac*bc*Cmat[0][1]-2.0*ab*cn*Cmat[0][1]-2.0*ac*bn*Cmat[0][2]+2.0*ab*bc*Cmat[0][2]-ac*ac*Cmat[1][1]+an*cn*Cmat[1][1]+2.0*ab*ac*Cmat[1][2]-2.0*an*bc*Cmat[1][2]-ab*ab*Cmat[2][2]+an*bn*Cmat[2][2];
  A2 = -cn*Cmat[0][1]*Cmat[0][1]+2.0*bc*Cmat[0][1]*Cmat[0][2]-bn*Cmat[0][2]*Cmat[0][2]+cn*Cmat[0][0]*Cmat[1][1]-2.0*ac*Cmat[0][2]*Cmat[1][1]-2.0*bc*Cmat[0][0]*Cmat[1][2]+2.0*ac*Cmat[0][1]*Cmat[1][2]+2.0*ab*Cmat[0][2]*Cmat[1][2]-an*Cmat[1][2]*Cmat[1][2]+bn*Cmat[0][0]*Cmat[2][2]-2.0*ab*Cmat[0][1]*Cmat[2][2]+an*Cmat[1][1]*Cmat[2][2];
  A3 = -Cmat[0][2]*Cmat[0][2]*Cmat[1][1]+2.0*Cmat[0][1]*Cmat[0][2]*Cmat[1][2]-Cmat[0][0]*Cmat[1][2]*Cmat[1][2]-Cmat[0][1]*Cmat[0][1]*Cmat[2][2]+Cmat[0][0]*Cmat[1][1]*Cmat[2][2];

  a0 = A0/A3;
  a1 = A1/A3;
  a2 = A2/A3;

  p = (3.0*a1-a2*a2)/3.0;
  q = (9.0*a1*a2-27.0*a0-2.0*a2*a2*a2)/27.0;
  C = 0.5*q*pow(3.0/fabs(p),1.5);

  if(C>=1.0) {
    Nlambda = 1;
    lambda[0] = 2.0*sqrt(fabs(p)/3.0)*Cosh(ACosh(C)/3.0)-a2/3.0;
  }
  if(C<=-1.0) {
    Nlambda = 1;
    lambda[0] = -2.0*sqrt(fabs(p)/3.0)*Cosh(ACosh(fabs(C))/3.0)-a2/3.0;

    /*
    {REAL8 tmp = Cosh(ACosh(3.02993)/3.0);
    printf("%g\n",tmp);}
    */

    /*
    printf("%g\t%g\t%g\t%g\t%g\n",p,C,a2,lambda[0],-2.0*sqrt(fabs(p)/3.0)*Cosh(ACosh(fabs(C))/3.0)-a2/3.0);
    */

  }
  if(fabs(C)<1.0) {
    Nlambda = 3;
    lambda[0] = 2.0*sqrt(fabs(p)/3.0)*cos(acos(C)/3.0)-a2/3.0;
    lambda[1] = 2.0*sqrt(fabs(p)/3.0)*cos((acos(C)+2.0*LAL_PI)/3.0)-a2/3.0;
    lambda[2] = 2.0*sqrt(fabs(p)/3.0)*cos((acos(C)+4.0*LAL_PI)/3.0)-a2/3.0;
  }

  maxLambda = -1e30;
  for(i=0;i<Nlambda;i++) {
    Hmat[0][0] = an + lambda[i]*Cmat[0][0];
    Hmat[0][1] = Hmat[1][0] = ab + lambda[i]*Cmat[0][1];
    Hmat[0][2] = Hmat[2][0] = ac + lambda[i]*Cmat[0][2];

    Hmat[1][1] = bn + lambda[i]*Cmat[1][1];
    Hmat[1][2] = Hmat[2][1] = bc + lambda[i]*Cmat[1][2];

    Hmat[2][2] = cn + lambda[i]*Cmat[2][2];

    /* note: following instruction modifies Hmat */
    if((Sret = SolveHomogeneous(Hmat,tmpAlpha,params->Ndetectors))) {
      if(Sret == 1) {
	ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
      } else {
	bzero(alpha,params->Ndetectors * sizeof(REAL8));
	break;
      }
    }

    tmpLambda = (tmpAlpha[0]*tmpAlpha[0]*an+tmpAlpha[1]*tmpAlpha[1]*bn+tmpAlpha[2]*tmpAlpha[2]*cn+2.0*tmpAlpha[0]*tmpAlpha[1]*ab+2.0*tmpAlpha[0]*tmpAlpha[2]*ac+2.0*tmpAlpha[1]*tmpAlpha[2]*bc)/(tmpAlpha[0]*tmpAlpha[0]*Cmat[0][0]+tmpAlpha[1]*tmpAlpha[1]*Cmat[1][1]+tmpAlpha[2]*tmpAlpha[2]*Cmat[2][2]+2.0*tmpAlpha[0]*tmpAlpha[1]*Cmat[0][1]+2.0*tmpAlpha[0]*tmpAlpha[2]*Cmat[0][2]+2.0*tmpAlpha[1]*tmpAlpha[2]*Cmat[1][2]);

    if(tmpLambda > maxLambda) {
      memcpy(alpha,tmpAlpha,params->Ndetectors * sizeof(REAL8));
      maxLambda = tmpLambda;
    }

    /*
    printf("%u\t%g\t%g\t%g\t%g\n",i,lambda[i],alpha[0],alpha[1],alpha[2]);
    */
  }

  /* loop */
  for(i=0; i<params->Ndetectors; i++) {

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
    for(j=iPad+1; j<output->data->length - ePad; j++) {

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
  LALFree(tmpAlpha);
  LALFree(lambda);

  for(i=0; i<params->Ndetectors; i++) {
    LALFree(Hmat[i]);
  }
  LALFree(Hmat);

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


static REAL4 sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static REAL4 maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static INT4 iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static INT4 SolveHomogeneous(
		 REAL8 **ar,
		 REAL8 *Alpha,
		 UINT4 m
		 ) {
  UINT4 n = m;
  REAL4 *w;
  REAL4 **v;
  REAL4 **a;

  static REAL4 pythag(REAL4 a, REAL4 b);
  INT4 flag,i,its,j,jj,k,l,nm;
  REAL4 anorm,c,f,g,h,s,scale,x,y,z,*rv1, wmin;

  /*
  printf("det = %g\n",a[0][0]*(a[1][1]*a[2][2]-a[2][1]*a[1][2])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]));
  */

  a = (REAL4 **)LALMalloc(n * sizeof(REAL4 *));
  if(!a) {
    return 1;
  }
  for(i=0;i<n;i++) {
    a[i] = (REAL4 *)LALMalloc(n * sizeof(REAL4));
    if(!(a[i])) {
      return 1;
    }
    for(j=0;j<n;j++) {
      a[i][j] = (REAL4)ar[i][j];
    }
    a[i] -= 1;
  }

  a -= 1;

  w = (REAL4 *)LALMalloc(n * sizeof(REAL4));
  if(!w) {
    return 1;
  }
  w -= 1;

  v = (REAL4 **)LALMalloc(n * sizeof(REAL4 *));
  if(!v) {
    return 1;
  }
  for(i=0;i<n;i++) {
    v[i] = (REAL4 *)LALMalloc(n * sizeof(REAL4));
    if(!(v[i])) {
      return 1;
    }
    v[i] -= 1;
  }
  v -= 1;

  rv1 = (REAL4 *)LALMalloc(n * sizeof(REAL4));
  if(!rv1) {
    return 1;
  }
  rv1 -= 1;

  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((float)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((float)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((float)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30) {
	fprintf(stderr,"no convergence in 30 svdcmp iterations. Setting all weights to zero.");
	return 2;
      }
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  

  /******************************************/
            
  i = 1;
  wmin = fabs(w[i]);

  for(j=2;j<=n;j++) {
    if(fabs(w[j]) < wmin) {
      i = j;
      wmin = fabs(w[j]);
    }
  }

  for(j=0;j<n;j++) {
    Alpha[j] = (REAL8)(v[j+1][i]);
  }

  LALFree(rv1 + 1);
  LALFree(w + 1);
  for(i=1;i<=n;i++) {
    LALFree(v[i] + 1);
    LALFree(a[i] + 1);
  }
  LALFree(v + 1);
  LALFree(a + 1);

  return 0;
}


float pythag(float a, float b)
{
        float absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


#undef SIGN
#undef FMAX
#undef IMIN
#undef SQR
