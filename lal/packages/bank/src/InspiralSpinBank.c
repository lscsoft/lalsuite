/**** <lalVerbatim file="InspiralSpinBankCV">
 * Authors: Hanna, C. R. and Owen, B. J.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{InspiralSpinBank.c}}
 *
 * This module creates a bank of templates to search for precessing
 * binaries.
 *
 * \subsubsection*{Prototypes}
 * \input{InspiralSpinBankCP}
 * %% \idx{LALInspiralSpinBank()}
 *
 * \subsubsection*{Description}
 *
 * This function creates a bank of templates to search for precessing
 * binaries.
 *
 * \subsubsection*{Algorithm}
 *
 * The target region of parameter space is a distorted box in the
 * coordinates $(x=\psi_0, y=\psi_3, z=\beta)$. The metric at high values of
 * $\beta$ is flat. It is convenient to rotate to coordinates $(x',y',z')$
 * which lie along eigenvectors of the metric.
 *
 * The algorithm first draws a rectilinear box in the primed coordinates
 * which includes the distorted box, then steps through along the directions
 * of the primed coordinates.  At each point it tests if the point lies
 * within the distorted box. If the point is inside the distorted box, the
 * algorithm adds a template to the linked list. If not, it continues.
 *
 * The tiling is done with a body-centered cubic lattice. People usually
 * solve the non-overlapping bcc problem rather than the overlapping one
 * here, so it's worth mentioning how we do it. I don't have time to stick
 * in the 3D figures you need to show it properly, but you can figure out
 * the spacing by finding the smallest sphere that contains the Wigner-Seitz
 * cell. When you do that you find that the lattice constant (spacing
 * between templates in the plane, in proper distance) is
 * $(4/3)\sqrt{2\mu}$. So the coordinate spacing is that divided by the
 * square root of the corresponding eigenvalue of the metric. (The vertical
 * spacing in the bcc lattice is multiplied by a further 1/2.)
 *
 * \subsubsection*{Uses}
 *
 * \begin{verbatim}
 * LALCalloc()
 * LALFree()
 * LALGetInspiralMoments()
 * LALSSymmetricEigenVectors()
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * Currently we use a static function for the metric based on an
 * approximation that is good only for large $\beta$. We should update it
 * and put it out in the LAL namespace.
 *
 * The metric relies on approximations that make it valid only for a binary
 * system with a total mass $<15M\odot$ where the larger body's minimum mass
 * is at least twice the smaller body's maximum mass.  Using mass values
 * that violate these conditions will result in an error message.   
 *
 * \vfill{\footnotesize\input{InspiralSpinBankCV}}
 *
 **** </lalLaTeX> */


#include <math.h>
#include <lal/AVFactories.h>
#include <lal/FlatMesh.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALMalloc.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/MatrixUtils.h>
#include <lal/SeqFactories.h>


#define INSPIRALSPINBANKC_ENOTILES 5 
#define INSPIRALSPINBANKC_MSGENOTILES "No templates were generated"




NRCSID(INSPIRALSPINBANKC, "$Id$");

/* Internal structures and functions --------------------------------------- */

static void cleanup(LALStatus *,
    REAL4Array 	**, 
    UINT4Vector **, 
    REAL4Vector **, 
    SnglInspiralTable 	*, 
    SnglInspiralTable 	*,
    INT4 	*
    ); /* cleanup() prototype */

static REAL4 calculateX(
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    INT4,
    REAL4
    ); /* calculateX() prototype */

static REAL4 calculateY(
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    INT4,
    REAL4
    ); /* calculateY() prototype */

static REAL4 calculateZ(
    REAL4,
    REAL4,
    REAL4
    ); /* calculateZ() prototype */

static INT4 test(
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4,
    REAL4 
    ); /* test() prototype */

/* END - Internal structures and functions --------------------------------- */


/* LALINSPIRALSPINBANKMETRIC() --------------------------------------------- */
static void
LALInspiralSpinBankMetric(
   LALStatus		*status,
   REAL4Array		*metric,
   InspiralMomentsEtc	*moments,
   InspiralTemplate	*inspiralTemplate,
   REAL4		*f0
   )
{
  INT4 loop = 0;
  REAL8 J1  = 0.0;
  REAL8 J4  = 0.0;
  REAL8 J6  = 0.0;
  REAL8 J9  = 0.0;
  REAL8 J11 = 0.0;
  REAL8 J12 = 0.0;
  REAL8 J14 = 0.0;
  REAL8 J17 = 0.0;

  INITSTATUS( status, "LALInspiralSpinBank", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( status );

  if (!metric){
    ABORT(status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
    }
  if (!moments){
    ABORT(status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
    }
  
  
  /* Rescale the moments to F0 = Noise Curve Minimum */
  for(loop = 1; loop <=17; loop++){
    moments->j[loop] *= pow((inspiralTemplate->fLower/(*f0)), ((7.0-(REAL4) loop)/3.0));
    }

/* This just copies the noise moment data from *moments */
  J1  = moments->j[1];
  J4  = moments->j[4];
  J6  = moments->j[6];
  J9  = moments->j[9];
  J11 = moments->j[11];
  J12 = moments->j[12];
  J14 = moments->j[14];
  J17 = moments->j[17];
                                                                                                                                                
  /* Set metric components as functions of moments. */

  metric->data[0] = (REAL4) 0.5*(J17-J12*J12-(J9-J4*J12)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[1] = (REAL4) 0.5*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[2] = (REAL4) 0;
  metric->data[3] = (REAL4) 0.5*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[4] = (REAL4) 0.5*(J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4));

  metric->data[5] = (REAL4) 0.0;
  metric->data[6] = (REAL4) 0.0;
  metric->data[7] = (REAL4) 0.0;
  metric->data[8] = (REAL4) 0.5*(J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4));
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALInspiralSpinBankMetric */


/* Allocate a template at these coordinate values. */
/* Does no error checking, so check immediately after calling. */
static void
allocate(
    REAL4               x,
    REAL4               y,
    REAL4               z,
    REAL4               f0,
    SnglInspiralTable **tmplt,
    INT4               *ntiles )
{
  REAL4 mass, eta, m1, m2;

  *tmplt = (*tmplt)->next = (SnglInspiralTable *) LALCalloc( 1,
           sizeof(SnglInspiralTable) );

  mass = -y/x / (16.0*LAL_PI*LAL_PI*f0);
  eta = 16.0457 * pow( -x*x/y/y/y/y/y, 0.3333333 );
  m1 = 0.5*mass* (1 + sqrt(1 - 4*eta));
  m2 = 0.5*mass* (1 - sqrt(1 - 4*eta));

  (*tmplt)->mass1 = m1;
  (*tmplt)->mass2 = m2;
  (*tmplt)->eta = eta;
  (*tmplt)->mchirp = pow(m1*m2,0.6)/pow(m1+m2,0.2);
  (*tmplt)->psi0 = x*pow(f0,5./3);
  (*tmplt)->psi3 = y*pow(f0,2./3);
  (*tmplt)->beta = z*pow(f0,2./3);
  ++(*ntiles);
} /* allocate() */


/* <lalVerbatim file="InspiralSpinBankCP"> */
void
LALInspiralSpinBank(
    LALStatus         	 *status,
    SnglInspiralTable   **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn *coarseIn
    )
/* </lalVerbatim> */
{
  SnglInspiralTable *tmplt = 	  NULL; /* loop counter */
  REAL4Array *metric = 		  NULL; /* parameter-space metric */
  UINT4Vector *metricDimensions = NULL;	/* contains the dimension of metric */
  REAL4Vector *eigenval =  	  NULL; /* eigenvalues of metric */
  InspiralMomentsEtc moments; 		/* Added for LALGetInspiralMoments() */
  InspiralTemplate inspiralTemplate; 	/* Added for LALGetInspiralMoments() */
  REAL4 x, y, z;             		/* psi0, psi3, beta coordinates */
  REAL4 x0, myy0, z0;          		/* minimum values of x, y, z */
  REAL4 x1, myy1, z1;          		/* maximum values of x, y, z */
  REAL4 xp, yp, zp;          		/* metric eigenvector coordinates */
  REAL4 xp0, yp0, zp0;       		/* minimum values of xp, yp, zp */
  REAL4 xp1, yp1, zp1;       		/* maximum values of xp, yp, zp */
  REAL4 dxp, dyp, dzp;       		/* step sizes in xp, yp, zp */
  REAL4 theta;               		/* angle of rotation for xp and yp */
  REAL4 m1Min, m1Max;       		/* range of m1 to search */
  REAL4 m2Min, m2Max;        		/* range of m2 to search */
  REAL4 f0 = -1;  			/* frequency of minimum of noise curve */
  INT2 bccFlag = 0;      		/* determines offset for bcc tiling */
  INT4 cnt = 0;				/* loop counter set to value of ntiles */
  REAL4 shf0 = 1;			/* used to find minimum of shf */
  
  /* Set up status pointer. */
  INITSTATUS( status, "LALInspiralSpinBank", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( status );
  
  
  ASSERT( coarseIn, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );
  /* Check that minimal match is OK. */
  ASSERT( coarseIn->mmCoarse > 0, status, LALINSPIRALBANKH_ECHOICE,
          LALINSPIRALBANKH_MSGECHOICE );
  ASSERT( coarseIn->mmCoarse < 1, status, LALINSPIRALBANKH_ECHOICE,
          LALINSPIRALBANKH_MSGECHOICE );
  /* Another mass bound needed, or go to psi bounds instead? */
  ASSERT( coarseIn->mMin > 0, status, LALINSPIRALBANKH_ECHOICE,
          LALINSPIRALBANKH_MSGECHOICE );
  ASSERT( coarseIn->MMax > 0, status, LALINSPIRALBANKH_ECHOICE,
          LALINSPIRALBANKH_MSGECHOICE );
  ASSERT( coarseIn->MMax > 2*coarseIn->mMin, status,
          LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE );
  /* Check that noise curve exists. */
  ASSERT( coarseIn->shf.data, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );
  ASSERT( coarseIn->shf.data->data, status, LALINSPIRALBANKH_ENULL,
          LALINSPIRALBANKH_MSGENULL );

  /*These parameters have not been added to InspiralCoarseBankIn yet, but when they are the will need to be checked */
  /*
    if (coarseIn->betaMax < 0) 
      ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  */

  /* Allocate memory for 3x3 parameter-space metric. */
  /* BEN: mess creating all these structures & adding TRYs etc */
  /* BEN: do it by hand, since it's so simple? */
  LALU4CreateVector( status->statusPtr, &metricDimensions, (UINT4) 2 );
  BEGINFAIL(status)
    cleanup(status->statusPtr, &metric, &metricDimensions, &eigenval, *tiles, tmplt, ntiles);
  ENDFAIL(status);
  metricDimensions->data[0] = 3;
  metricDimensions->data[1] = 3;
  LALSCreateArray( status->statusPtr, &metric, metricDimensions );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt, ntiles);
  ENDFAIL(status);

  /* Set f0 to frequency of minimum of noise curve. */
  for( cnt = 0; cnt < (INT4) coarseIn->shf.data->length; cnt++ )
  {
    if( coarseIn->shf.data->data[cnt] > 0 && coarseIn->shf.data->data[cnt] <
        shf0 )
    {
      f0 = (REAL4) coarseIn->shf.deltaF * cnt;
      shf0 = coarseIn->shf.data->data[cnt];
    }
  }
  /* Compute noise moments. */
  inspiralTemplate.fLower = coarseIn->fLower;
  inspiralTemplate.fCutoff = coarseIn->fUpper;
  LALGetInspiralMoments( status->statusPtr, &moments, &(coarseIn->shf),
                         &inspiralTemplate );
  BEGINFAIL(status)                                                           
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
  ENDFAIL(status);

  /* Call the metric */
  LALInspiralSpinBankMetric(status->statusPtr, metric, &moments, &inspiralTemplate, &f0);	
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt, ntiles);
  ENDFAIL(status);

  /* Find eigenvalues and eigenvectors of metric. */
  eigenval = NULL;
  LALSCreateVector( status->statusPtr, &eigenval, 3 );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric, &metricDimensions,&eigenval,*tiles,tmplt, ntiles);
  ENDFAIL(status);
  LALSSymmetricEigenVectors( status->statusPtr, eigenval, metric );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric, &metricDimensions,&eigenval,*tiles,tmplt, ntiles);
  ENDFAIL(status);
    
  /* Set stepsizes and xp-yp rotation angle from metric. */
  if( eigenval->data[0] <= 0 || eigenval->data[1] <=0 || eigenval->data[2]
      <= 0 )
  {
    ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  }
  dxp = 1.333333*sqrt(2*(1-coarseIn->mmCoarse)/eigenval->data[0]);
  dyp = 1.333333*sqrt(2*(1-coarseIn->mmCoarse)/eigenval->data[1]);
  dzp = 0.6666667*sqrt(2*(1-coarseIn->mmCoarse)/eigenval->data[2]);
  theta = atan2( -metric->data[3], -metric->data[0] );

  /* Hardcode mass range on higher mass for the moment. */
  m2Min = coarseIn->mMin*LAL_MTSUN_SI;
  m2Max = coarseIn->MMax*LAL_MTSUN_SI/2;
  m1Min = 2.0*m2Max;
  m1Max = 15.0*LAL_MTSUN_SI - m2Max;

  /* Set box on unprimed coordinates including region. */
  x0 = 0.9*(3.0/128) / (pow(LAL_PI*f0*(m1Max+m2Max),1.666667)*(m1Max*m2Max/pow(m1Max+m2Max,2)));
  myy0 = 1.1*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Max+m2Min),0.6666667)*(m1Max*m2Min/pow(m1Max+m2Min,2)));
  z0 = 0;
  x1 = 1.1*(3.0/128) / (pow(LAL_PI*f0*(m1Min+m2Min),1.666667)*(m1Min*m2Min/pow(m1Min+m2Min,2)));
  myy1 = .9*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Min+m2Max),0.6666667)*(m1Min*m2Max/pow(m1Min+m2Max,2)));
  z1 = 3.8* LAL_PI/29.961432 * (1+0.75*m2Max/m1Min) * (m1Max/m2Min) * pow(LAL_MTSUN_SI*100.0/(m1Min+m2Min), 0.6666667);

  /* Set boundaries of box in primed coordinates. */
  xp0 = x0 + sin(theta)*sin(theta) * (x1 - x0);
  yp0 = myy0 - cos(theta)*sin(theta) * (x1 - x0);
  yp1 = sin(theta) * (x1 - x0) + cos(theta) * (myy1 - myy0);
  xp1 = sin(theta) * (myy1 - myy0) + cos(theta) * (x1 - x0);
  zp0 = z0;
  zp1 = z1;
    
  /* Allocate first template, which will remain blank. */
  *tiles = tmplt = (SnglInspiralTable *) LALCalloc( 1,
                   sizeof(SnglInspiralTable) );
  if (*tiles == NULL) {
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
    ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }
  *ntiles = 0;

  /* This loop generates the template bank. */
  for (zp = 0; zp <= zp1; zp += dzp)
  {
    bccFlag++;
    for (yp = 0; yp<= yp1; yp += dyp)
    {
      for (xp = 0; xp <= xp1; xp += dxp)
      {
        /* Calculate Coordinate values at this point */
        x = calculateX(0, xp0, xp, dxp, yp, dyp, bccFlag, theta);
        y = calculateY(0, yp0, xp, dxp, yp, dyp, bccFlag, theta);
        z = calculateZ(0, zp, dzp);
        /* If the point is in the target region, allocate a template. */
        if( test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0) )
        {
          allocate( x, y, z, f0, &tmplt, ntiles );
          if( tmplt == NULL )
          {
            cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
            ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
          }
        }
        
        /* Test a spot dx behind */
        x = calculateX(-1, xp0, xp, dxp, yp, dyp, bccFlag, theta);
        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          x = calculateX(-0.5, xp0, xp, dxp, yp, dyp, bccFlag, theta);
          /* If its not in the range check 1/2 dx behind. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if( tmplt == NULL )
            {
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
            }
          }
        }

        /* Test a spot dy behind */
        x = calculateX(0, xp0, xp, dxp, yp, dyp, bccFlag, theta);
        y = calculateY(-1, yp0, xp, dxp, yp, dyp, bccFlag, theta);
        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          y = calculateY(-0.5, yp0, xp, dxp, yp, dyp, bccFlag, theta);
          /* If its not in the range check 1/2 dy behind. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if (!tmplt){
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
              }
          }
        }

        /* Test a spot dz behind */
        y = calculateY(0, yp0, xp, dxp, yp, dyp, (bccFlag+1), theta);
        z = calculateZ(-1, zp, dzp);
        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          z = calculateZ(-0.5, zp, dzp);
          /*  if its not in the range check 1/2 dz behind. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if (!tmplt){
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
              }
          }
        }

        /* Test a spot dx ahead */
        x = calculateX(1, xp0, xp, dxp, yp, dyp, bccFlag, theta);
        y = calculateY(0, yp0, xp, dxp, yp, dyp, (bccFlag+1), theta);
        z = calculateZ(0, zp, dzp);

        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          x = calculateX(0.5, xp0, xp, dxp, yp, dyp, bccFlag, theta);
          /* If its not in the range check 1/2 dx. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if (!tmplt){
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
              }
          }
        }

        /* Test a spot dy ahead */
        x = calculateX(0, xp0, xp, dxp, yp, dyp, bccFlag, theta);
        y = calculateY(1, yp0, xp, dxp, yp, dyp, bccFlag, theta);
        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          y = calculateY(0.5, yp0, xp, dxp, yp, dyp, bccFlag, theta);
          /* If its not in the range check 1/2 dy. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if (!tmplt){
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
              }
          }
        }

        /* Test a spot dz ahead */
        y = calculateY(0, yp0, xp, dxp, yp, dyp, (bccFlag+1), theta);
        z = calculateZ(1, zp, dzp);
        if(!test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
          z = calculateZ(0.5, zp, dzp);
          /*  if its not in the range check 1/2 dz. */
          if (test(x,y,z,m1Min,m1Max,m2Min,m2Max,f0)){
            allocate( x, y, z, f0, &tmplt, ntiles );
            if (!tmplt){
              cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
              }
          }
        }
      } /* for (zp...) */
    } /* for (yp...) */
  } /* for (zp...) */

  /* Trim the first template, which was left blank. */
  tmplt = (*tiles)->next;
  LALFree( *tiles );
  *tiles = tmplt;

  /* What if no templates were allocated? ABORT */
  if (*tiles == NULL) 
  {  
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,*tiles,tmplt,ntiles);
    ABORT(status, INSPIRALSPINBANKC_ENOTILES, INSPIRALSPINBANKC_MSGENOTILES);
  }
  
  /* free memory allocated for the vectors and arrays */
  cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,NULL,NULL,&cnt);
  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALInspiralSpinBank() */


static void cleanup(
    LALStatus *s, 
    REAL4Array **m, 
    UINT4Vector **md, 
    REAL4Vector **e, 
    SnglInspiralTable *f,
    SnglInspiralTable *t,
    INT4 *nt)
{
  INITSTATUS( s, "LALInspiralSpinBank-cleanup", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( s );

  if (m){
    TRY(LALU4DestroyVector(s->statusPtr,md),s);
    }
  if (md){
    TRY(LALSDestroyVector(s->statusPtr, e),s);
    }
  if (e){
    TRY(LALSDestroyArray(s->statusPtr, m),s); 
    }
  if (t && f)
  {
    t = f;
    while ((t->next) && (*nt > 0))
    {
      f = t;
      t = t->next;
      LALFree(f);
      --(*nt);
    }/* while(tmplt) */
  LALFree(t);
  --(*nt);
  }
  DETATCHSTATUSPTR( s );
  RETURN( s );
}

static REAL4 calculateX(REAL4 n,
               REAL4 xp0,
               REAL4 xp,
               REAL4 dxp,
               REAL4 yp,
               REAL4 dyp,
               INT4 bccFlag, 
               REAL4 theta)
  {
  return (xp0 + (n*dxp + xp+dxp/2.0*((bccFlag)%2))*cos(theta) - 
                (yp+dyp/2.0*((bccFlag)%2))*sin(theta));
  } /* REAL4 x(); */

static REAL4 calculateY(REAL4 n,
               REAL4 yp0,
               REAL4 xp,
               REAL4 dxp,
               REAL4 yp,
               REAL4 dyp,
               INT4 bccFlag,
               REAL4 theta)
  {
  return (yp0 + (xp+dxp/2.0*((bccFlag)%2))*sin(theta) + 
                (n*dyp + yp+dyp/2.0*((bccFlag)%2))*cos(theta));
  } /* REAL4 y(); */

static REAL4 calculateZ(REAL4 n,
               REAL4 zp,
               REAL4 dzp)
  {
  return (zp + n*dzp);
  } /* REAL4 z(); */

static INT4 test(REAL4 x, 
                 REAL4 y, 
                 REAL4 z,
                 REAL4 m1Min,
                 REAL4 m1Max,
                 REAL4 m2Min,
                 REAL4 m2Max,
                 REAL4 f0){

  REAL4 mass, eta, m1, m2, betaMax;
  mass = -y/x / (16.0*LAL_PI*LAL_PI*f0);
  eta = 16.0457 * pow( -x*x/y/y/y/y/y, 0.3333333 );
  if (eta > 0.25 || eta < 0)
    return 0;
  m1 = 0.5*mass* (1 + sqrt(1 - 4*eta));
  m2 = 0.5*mass* (1 - sqrt(1 - 4*eta));
  if (m1 > m1Max || m1 < m1Min || m2 > m2Max || m2 < m2Min)
    return 0;
  betaMax = 3.8*LAL_PI/29.961432 * (1+0.75*m2/m1)*(m1/m2) * pow((LAL_MTSUN_SI*100.0/mass),0.6666667);
  if (z > betaMax)
    return 0;
  return 1;
}
