/**** <lalVerbatim file="InspiralSpinBankCV">
 * Authors: Hanna, C. and Owen, B. J.
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
 * The tiling is done with a body-centered cubic lattice.
 *
 * At the end it copies the linked list into the inspiral package's array
 * format.
 *
 * \subsubsection*{Uses}
 *
 * \begin{verbatim}
 * LALCalloc()
 * LALFree()
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * Currently we use dummy functions for the metric and noise moments. These
 * should be updated, especially to account for real noise spectra.
 *
 * The issue of linked lists vs.\ arrays needs to be seriously addressed. As
 * our experience with this code shows, multidimensional tiling of
 * complicated parameter regions demands the flexibility of linked lists.
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

#define INSPIRALSPINBANKC_ENONPOSITIVEMM 1
#define INSPIRALSPINBANKC_MSGENONPOSITIVEMM "Minimum match value is not positive."



NRCSID(INSPIRALSPINBANKC, "$Id$");


static void tmpmetric( REAL4Array *metric );


/* <lalVerbatim file="InspiralSpinBankCP"> */
void
LALInspiralSpinBank(
    LALStatus         	 *status,
    InspiralTemplateList **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn  coarseIn
    )
/* </lalVerbatim> */
{
  SnglInspiralTable *tmplt;  /* loop counter */
  SnglInspiralTable *output; /* head of output linked list */
  REAL4Array *metric;        /* parameter-space metric */
  UINT4Vector *metricDimensions;
  REAL4Vector *eigenval;     /* eigenvalues of metric */
  REAL4 x, y, z;             /* psi0, psi3, beta coordinates */
  REAL4 x0, y0, z0;          /* minimum values of x, y, z */
  REAL4 x1, y1, z1;          /* maximum values of x, y, z */
  REAL4 xp, yp, zp;          /* metric eigenvector coordinates */
  REAL4 xp0, yp0, zp0;       /* minimum values of xp, yp, zp */
  REAL4 xp1, yp1, zp1;       /* maximum values of xp, yp, zp */
  REAL4 dxp, dyp, dzp;       /* step sizes in xp, yp, zp */
  REAL4 theta;               /* angle of rotation for xp and yp */
  REAL4 m1;                  /* greater binary component mass */
  REAL4 m1Min, m1Max;        /* range of m1 to search */
  REAL4 m2;                  /* lesser binary component mass */
  REAL4 m2Min, m2Max;        /* range of m2 to search */
  REAL4 mass;                /* total mass of binary */
  REAL4 eta;                 /* symmetric mass ratio of binary */
  REAL4 betaMax;             /* maximum spin parameter of binary */
  REAL4 f0;                  /* frequency of minimum of noise curve */
  INT2 bccFlag = 0;          /* determines offset for bcc tiling */

  /* Set up status pointer. */
  INITSTATUS( status, "LALInspiralSpinBank", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( status );

  if (coarseIn.mmCoarse <= 0) ABORT(status, INSPIRALSPINBANKC_ENONPOSITIVEMM, INSPIRALSPINBANKC_MSGENONPOSITIVEMM);

  /* Get noise power moments and trig moments. */
  /* Hardcode this for the moment. */
  f0 = 164;

  /* Get 3x3 parameter-space metric. */
  /* BEN: mess creating all these structures & adding TRYs etc */
  /* BEN: do it by hand, since it's so simple? */
  metricDimensions = NULL;

  LALU4CreateVector( status->statusPtr, &metricDimensions, 2 );
  BEGINFAIL(status)
    TRY(LALU4DestroyVector(status->statusPtr,&metricDimensions),status);
  ENDFAIL(status);
  
  metricDimensions->data[0] = 3;
  metricDimensions->data[1] = 3;
  metric = NULL;

  LALSCreateArray( status->statusPtr, &metric, metricDimensions );
  BEGINFAIL(status)
  {
    TRY(LALSDestroyArray(status->statusPtr, &metric),status);
    TRY(LALU4DestroyVector(status->statusPtr,&metricDimensions),status);
  }
  ENDFAIL(status);

  tmpmetric( metric );

  /* Find eigenvalues and eigenvectors of metric. */
  eigenval = NULL;

  LALSCreateVector( status->statusPtr, &eigenval, 3 );
  BEGINFAIL(status)
  {
    TRY(LALU4DestroyVector(status->statusPtr,&metricDimensions),status);
    TRY(LALSDestroyVector(status->statusPtr,&eigenval),status);
    TRY(LALSDestroyArray(status->statusPtr, &metric),status);  
  }
  ENDFAIL(status);

  LALSSymmetricEigenVectors( status->statusPtr, eigenval, metric );
  BEGINFAIL(status)
  {
    TRY(LALU4DestroyVector(status->statusPtr,&metricDimensions),status);
    TRY(LALSDestroyVector(status->statusPtr,&eigenval),status);
    TRY(LALSDestroyArray(status->statusPtr, &metric),status);  
  }
  ENDFAIL(status);

  /* Set stepsizes and xp-yp rotation angle from metric. */
  dxp = 1.333333*sqrt(2*coarseIn.mmCoarse/eigenval->data[0]);
  dyp = 1.333333*sqrt(2*coarseIn.mmCoarse/eigenval->data[1]);
  dzp = 0.6666667*sqrt(2*coarseIn.mmCoarse/eigenval->data[2]);
  theta = atan2( -metric->data[3], -metric->data[0] );

  /* Hardcode mass range etc for the moment. */
  m1Min = 5*LAL_MTSUN_SI;
  m1Max = 10*LAL_MTSUN_SI;
  m2Min = 1*LAL_MTSUN_SI;
  m2Max = 2*LAL_MTSUN_SI;

  /* Set box on unprimed coordinates including region. */
  x0 = 0.9*(3.0/128) / (pow(LAL_PI*f0*(m1Max+m2Max),1.666667)*(m1Max*m2Max/pow(m1Max+m2Max,2)));
  y0 = 1.1*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Max+m2Min),0.6666667)*(m1Max*m2Min/pow(m1Max+m2Min,2)));
  z0 = 0;
  x1 = 1.1*(3.0/128) / (pow(LAL_PI*f0*(m1Min+m2Min),1.666667)*(m1Min*m2Min/pow(m1Min+m2Min,2)));
  y1 = .9*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Min+m2Max),0.6666667)*(m1Min*m2Max/pow(m1Min+m2Max,2)));
  z1 = 3.8* LAL_PI/29.961432 * (1+0.75*m2Max/m1Min) * (m1Max/m2Min) * pow(LAL_MTSUN_SI*100.0/(m1Min+m2Min), 0.6666667);

  /* Set boundaries of box in primed coordinates. */
  xp0 = x0 + sin(theta)*sin(theta) * (x1 - x0);
  yp0 = y0 - cos(theta)*sin(theta) * (x1 - x0);
  yp1 = sin(theta) * (x1 - x0) + cos(theta) * (y1 - y0);
  xp1 = sin(theta) * (y1 - y0) + cos(theta) * (x1 - x0);
  zp0 = z0;
  zp1 = z1;
    
  /* Allocate first template, which will remain blank. */
  output = tmplt = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );
  /* BEN: error check here */

  /* This loop generates the template bank. */
  *ntiles = 0;
  for (zp = 0; zp <= zp1; zp += dzp)
  {
    bccFlag++;
    for (yp = 0; yp<= yp1; yp += dyp)
    { 
      for (xp = 0; xp <= xp1; xp += dxp)
      {
        x = xp0 + (xp+dxp/2.0*(bccFlag%2))*cos(theta) - (yp+dyp/2.0*(bccFlag%2))*sin(theta);
	y = yp0 + (xp+dxp/2.0*(bccFlag%2))*sin(theta) + (yp+dyp/2.0*(bccFlag%2))*cos(theta);
	z = zp;
        mass = -y/x / (16.0*LAL_PI*LAL_PI*f0);
        eta = 16.0457 * pow( -x*x/y/y/y/y/y, 0.3333333 );
        if (eta > 0.25 || eta < 0)
          continue;
        m1 = 0.5*mass* (1 + sqrt(1 - 4*eta));
        m2 = 0.5*mass* (1 - sqrt(1 - 4*eta));
        if (m1 > m1Max || m1 < m1Min || m2 > m2Max || m2 < m2Min)
          continue;
        betaMax = 3.8*LAL_PI/29.961432 * (1+0.75*m2/m1)*(m1/m2) * pow((LAL_MTSUN_SI*100.0/mass),0.6666667);
        if (z > betaMax)
          continue;
        tmplt = tmplt->next = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );
        tmplt->mass1 = m1;
        tmplt->mass2 = m2;
        tmplt->eta = eta;
        tmplt->mchirp = pow(m1*m2,0.6)/pow(m1+m2,0.2);
        tmplt->psi0 = x;            
        tmplt->psi3 = y; 
/*        tmplt->beta = z;*/
        ++(*ntiles);
      } /* for (zp...) */
    } /* for (yp...) */
  } /* for (zp...) */

  /* Trim the first template, which was left blank. */
  tmplt = output->next;
  LALFree( output );
  /* BEN: error check here */
  output = tmplt;

  /* What if no templates were allocated? ABORT or what? */
  if (!output) { }
  
  
  /* Convert output to communicate with LALInspiralCreateCoarseBank(). */
  *tiles = (InspiralTemplateList *) LALCalloc( *ntiles, sizeof(InspiralTemplateList) );
  for (tmplt = output; tmplt; tmplt = tmplt->next)
  {
    (*tiles)->params.mass1 = tmplt->mass1;
    (*tiles)->params.mass2 = tmplt->mass2;
    (*tiles)->params.psi0 = tmplt->psi0;
    (*tiles)->params.psi3 = tmplt->psi3;
/*    (*tiles)->params.beta = tmplt->beta;*/
  } /* for(tmplt...) */
  
  /* Free the memory allocated for the linked list. */
  tmplt = output;
  while (tmplt->next)
  {
    output = tmplt;
    tmplt = tmplt->next;
    LALFree(output);
  }/* while(tmptl) */
  
/*  ASSERT(!tmplt->next) to make sure the deallocation worked
    LALFree(tmplt)       free the final node */

  /* Free the memory for the vectors and arrays allocated */
  TRY(LALU4DestroyVector(status->statusPtr,&metricDimensions),status);
  TRY(LALSDestroyVector(status->statusPtr,&eigenval),status);
  TRY(LALSDestroyArray(status->statusPtr, &metric),status);  
  
  /* Clean up and leave. */
  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALInspiralSpinBank() */


/* Temporary metric function. OK for high beta, LIGO-I SRD noise. */

static void tmpmetric( REAL4Array *metric )
{
  /* These noise moment values are from Ben's Mathematica notebook. */
  REAL4 J1  = 1.20978;
/*  REAL4 J2  = 1.05754;*/
/*  REAL4 J3  = 0.971883;*/
  REAL4 J4  = 0.931058;
/*  REAL4 J5  = 0.924365;*/
  REAL4 J6  = 0.947410;	
/*  REAL4 J7  = 1;*/
/*  REAL4 J8  = 1.08541;*/
  REAL4 J9  = 1.21048;	
/*  REAL4 J10 = 1.38642;*/
  REAL4 J11 = 1.63039;	
  REAL4 J12 = 1.96795;
/*  REAL4 J13 = 2.43713;*/
  REAL4 J14 = 3.09454;
/*  REAL4 J15 = 4.02482;*/
/*  REAL4 J16 = 5.35540;*/
  REAL4 J17 = 7.27921;

  /* Set metric components as functions of moments. */
  metric->data[0] = (REAL4) (1.5)*(J17-J12*J12-(J9-J4*J12)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[1] = (REAL4) (1.5)*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));	
  metric->data[2] = (REAL4) 0;
  metric->data[3] = (REAL4) (1.5)*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));		
  metric->data[4] = (REAL4) (1.5)*(J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4));
  metric->data[5] = (REAL4) 0.0;
  metric->data[6] = (REAL4) 0.0;
  metric->data[7] = (REAL4) 0.0;
  metric->data[8] = (REAL4) J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4);
} /* tmpmetric() */
