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


#define INSPIRALSPINBANKC_ENOTILES 5 
#define INSPIRALSPINBANKC_MSGENOTILES "No templates were generated"




NRCSID(INSPIRALSPINBANKC, "$Id$");

typedef struct Node
{
  REAL4 psi0;
  REAL4 psi3;
  REAL4 beta;
  REAL4 eta;
  REAL4 chirpMass;
  REAL4 mass1;
  REAL4 mass2;
  struct Node *next;
} Node;

  
static void tmpmetric( REAL4Array *metric );
static void cleanup(LALStatus *s,
    REAL4Array **m, 
    UINT4Vector **md, 
    REAL4Vector **e, 
    Node *f, 
    Node *t,
    INT4 *nt);


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
  Node *tmplt = NULL;  	/* loop counter */
  Node *output = NULL; 	/* head of output linked list */
  REAL4Array *metric = NULL;        	/* parameter-space metric */
  UINT4Vector *metricDimensions = NULL;
  REAL4Vector *eigenval =  NULL;     	/* eigenvalues of metric */
  REAL4 x, y, z;             		/* psi0, psi3, beta coordinates */
  REAL4 x0, y0, z0;          		/* minimum values of x, y, z */
  REAL4 x1, y1, z1;          		/* maximum values of x, y, z */
  REAL4 xp, yp, zp;          		/* metric eigenvector coordinates */
  REAL4 xp0, yp0, zp0;       		/* minimum values of xp, yp, zp */
  REAL4 xp1, yp1, zp1;       		/* maximum values of xp, yp, zp */
  REAL4 dxp, dyp, dzp;       		/* step sizes in xp, yp, zp */
  REAL4 theta;               		/* angle of rotation for xp and yp */
  REAL4 m1;                  		/* greater binary component mass */
  REAL4 m1Min, m1Max;       		/* range of m1 to search */
  REAL4 m2;                  		/* lesser binary component mass */
  REAL4 m2Min, m2Max;        		/* range of m2 to search */
  REAL4 mass;                		/* total mass of binary */
  REAL4 eta;                 		/* symmetric mass ratio of binary */
  REAL4 betaMax;             		/* maximum spin parameter of binary */
  REAL4 f0;                  		/* frequency of minimum of noise curve */
  INT2 bccFlag = 0;          		/* determines offset for bcc tiling */
  INT4 cnt = 0;		     		/* loop counter set to value of ntiles */

  /* Set up status pointer. */
  INITSTATUS( status, "LALInspiralSpinBank", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( status );
  
  
  /* Check to make sure that all the parameters are okay */
  if (coarseIn.mmCoarse <= 0) 
    ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  
  /* These parameters have not been added to InspiralCoarseBankIn yet, but when they are the will need to be checked */
  /*  if ((coarseIn.m1Min <= 0) || (coarseIn.m2Min <= 0) || (coarseIn.m1Max <= 0) || (coarseIn.m2Max <= 0))
      ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);

    if (coarseIn.betaMax < 0) 
      ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);

  */
    
  /* Get noise power moments and trig moments. */
  /* Hardcode this for the moment. */
  f0 = 164;

  /* Get 3x3 parameter-space metric. */
  /* BEN: mess creating all these structures & adding TRYs etc */
  /* BEN: do it by hand, since it's so simple? */
  metricDimensions = NULL;

  LALU4CreateVector( status->statusPtr, &metricDimensions, 2 );
  BEGINFAIL(status)
    cleanup(status->statusPtr, &metric, &metricDimensions, &eigenval, output, tmplt, ntiles);
  ENDFAIL(status);
  
  metricDimensions->data[0] = 3;
  metricDimensions->data[1] = 3;
  metric = NULL;

  LALSCreateArray( status->statusPtr, &metric, metricDimensions );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,output,tmplt, ntiles);
  ENDFAIL(status);

  tmpmetric( metric );

  /* Find eigenvalues and eigenvectors of metric. */
  eigenval = NULL;

  LALSCreateVector( status->statusPtr, &eigenval, 3 );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric, &metricDimensions,&eigenval,output,tmplt, ntiles);
  ENDFAIL(status);
  
  LALSSymmetricEigenVectors( status->statusPtr, eigenval, metric );
  BEGINFAIL(status)
    cleanup(status->statusPtr,&metric, &metricDimensions,&eigenval,output,tmplt, ntiles);
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
  output = tmplt = (Node *) LALCalloc( 1, sizeof(Node) );
  if (!output) 
  {
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,output,tmplt,ntiles);
    ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }



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
        tmplt = tmplt->next = (Node *) LALCalloc( 1, sizeof(Node) );
        /* check to see if calloc worked */
        if (!tmplt) 
        {
          cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,output,tmplt,ntiles);
          ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
        }
        tmplt->mass1 = m1;
        tmplt->mass2 = m2;
        tmplt->eta = eta;
        tmplt->chirpMass = pow(m1*m2,0.6)/pow(m1+m2,0.2);
        tmplt->psi0 = x;            
        tmplt->psi3 = y; 
        tmplt->beta = z;
        ++(*ntiles);
      } /* for (zp...) */
    } /* for (yp...) */
  } /* for (zp...) */

  /* Trim the first template, which was left blank. */
  tmplt = output->next;
  LALFree( output );
  /* BEN: error check here */
  output = tmplt;

  /* What if no templates were allocated? ABORT */
  if (!output) 
  {  
    cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,output,tmplt,ntiles);
    ABORT(status, INSPIRALSPINBANKC_ENOTILES, INSPIRALSPINBANKC_MSGENOTILES);
  }
  
  /* Convert output to communicate with LALInspiralCreateCoarseBank(). */
  *tiles = (InspiralTemplateList *) LALCalloc( *ntiles, sizeof(InspiralTemplateList));
  cnt = 0;
  for (tmplt = output; tmplt; tmplt = tmplt->next)
  {
    (*tiles)[cnt].params.mass1 = tmplt->mass1;
    (*tiles)[cnt].params.mass2 = tmplt->mass2;
    (*tiles)[cnt].params.psi0 = tmplt->psi0;
    (*tiles)[cnt].params.psi3 = tmplt->psi3;
    (*tiles)[cnt].params.eta = tmplt->eta;
    (*tiles)[cnt].params.chirpMass = tmplt->chirpMass;
    (*tiles)[cnt].params.beta = tmplt->beta;
    ++cnt;
  } /* for(tmplt...) */

  
  
  /* prepare the link list to be freed by copying the number of tiles to cnt */
  tmplt = output;
  cnt = *ntiles;
  
  /* free memory allocated for the linked list, vectors and arrays */
  cleanup(status->statusPtr,&metric,&metricDimensions,&eigenval,output,tmplt,&cnt);
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

static void cleanup(
    LALStatus *s, 
    REAL4Array **m, 
    UINT4Vector **md, 
    REAL4Vector **e, 
    Node *f,
    Node *t,
    INT4 *nt)
{
  INITSTATUS( s, "LALInspiralSpinBank-cleanup", INSPIRALSPINBANKC );
  ATTATCHSTATUSPTR( s );

  if (m)
    TRY(LALU4DestroyVector(s->statusPtr,md),s);
  if (md)
    TRY(LALSDestroyVector(s->statusPtr, e),s);
  if (e)
    TRY(LALSDestroyArray(s->statusPtr, m),s); 
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
