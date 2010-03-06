/*
*  Copyright (C) 2007 Chad Hanna, Benjamin Owen
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

/**** <lalVerbatim file="InspiralSpinBankwNDTemplateBankCV">
 * Authors: Hanna, C. R. and Owen, B. J.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{InspiralSpinBankwNDTemplateBank.c}}
 *
 * This module creates a bank of templates to search for precessing
 * binaries.
 *
 * \subsubsection*{Prototypes}
 * \input{InspiralSpinBankwNDTemplateBankCP}
 * %% \idx{LALInspiralSpinBankwNDTemplateBank()}
 *
 * \subsubsection*{Description}
 *
 * This function creates a bank of templates to search for precessing
 * binaries. It uses the general tiling algorithm LALNDTemplateBank()
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
 *
 * The metric relies on approximations that make it valid only for a binary
 * system with a total mass $<15M\odot$ where the larger body's minimum mass
 * is at least twice the smaller body's maximum mass.  Using
 * that violate these conditions will result in an error message.
 *
 * The issue of linked lists vs.\ arrays needs to be seriously addressed. As
 * our experience with this code shows, multidimensional tiling of
 * complicated parameter regions demands the flexibility of linked lists.
 *
 * \vfill{\footnotesize\input{InspiralSpinBankwNDTemplateBankCV}}
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
#include <lal/TemplateBankGeneration.h>


NRCSID(INSPIRALSPINBANKWNDTEMPLATEBANKC, "$Id$");


/* LALINSPIRALSPINBANKMETRIC() --------------------------------------------- */
void
LALInspiralSpinBankMetric(
   LALStatus		*status,
   NDTemplateBankInput  *input,
   REAL4Array		*metric
   )
{
  InspiralMomentsEtc moments;           /* Added for LALGetInspiralMoments() */
  InspiralTemplate inspiralTemplate;    /* Added for LALGetInspiralMoments() */

  INT2 loop = 0;
  REAL4 f0 = input->f0;
  REAL8 J1  = 0.0;
  REAL8 J4  = 0.0;
  REAL8 J6  = 0.0;
  REAL8 J9  = 0.0;
  REAL8 J11 = 0.0;
  REAL8 J12 = 0.0;
  REAL8 J14 = 0.0;
  REAL8 J17 = 0.0;

  INITSTATUS( status, "LALInspiralSpinBankMetric", INSPIRALSPINBANKWNDTEMPLATEBANKC );
  ATTATCHSTATUSPTR( status );


  inspiralTemplate.fLower  = 30;        /* These are arbitrarily chosen for now */
  inspiralTemplate.fCutoff = 2000;      /* They are necessary for LALInspiralGetMoments() */

  LALGetInspiralMoments( status->statusPtr, &moments, input->PSD, &inspiralTemplate );

  /* Rescale the moments to F0 = Noise Curve Minimum */
  for(loop = 1; loop <=17; loop++){
    moments.j[loop] *= pow((inspiralTemplate.fLower/(f0)), ((7.0-(REAL4) loop)/3.0));
    }

/* This just copies the noise moment data from *moments */
  J1  = moments.j[1];
  J4  = moments.j[4];
  J6  = moments.j[6];
  J9  = moments.j[9];
  J11 = moments.j[11];
  J12 = moments.j[12];
  J14 = moments.j[14];
  J17 = moments.j[17];

  /* Set metric components as functions of moments. */
  metric->data[0] = (REAL4) (1.5)*(J17-J12*J12-(J9-J4*J12)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[1] = (REAL4) (1.5)*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[2] = (REAL4) 0.0;
  metric->data[3] = (REAL4) (1.5)*(J14-J9*J12-(J6-J4*J9)*(J9-J4*J12)/(J1-J4*J4));
  metric->data[4] = (REAL4) (1.5)*(J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4));
  metric->data[5] = (REAL4) 0.0;
  metric->data[6] = (REAL4) 0.0;
  metric->data[7] = (REAL4) 0.0;
  metric->data[8] = (REAL4) J11-J9*J9-(J6-J4*J9)*(J6-J4*J9)/(J1-J4*J4);

  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALInspiralSpinBankMetric */


/* LALInspiralSpinBankwNDTemplateBank() --------------------------------------------------- */
/* <lalVerbatim file="LALInspiralSpinBankwNDTemplateBankCP"> */
void
LALInspiralSpinBankwNDTemplateBank(
    LALStatus         	 *status,
    InspiralTemplateList **tiles,
    INT4      		 *ntiles,
    InspiralCoarseBankIn  coarseIn
    )
/* </lalVerbatim> */
{

  /* Initialize variables */
  REAL4Array *metric = 		  NULL; /* parameter-space metric */
  UINT4Vector *metricDimensions = NULL;	/* contains the dimension of metric */
  NDTemplateBankInput NDinput;
  NDTemplateBankOutput *NDoutput = NULL;
  NDTemplateBankOutput *NDFirst  = NULL;
  NDTemplateBankFunctionPtrs NDFunctionPtrs;
  INT4 cnt = 0;
  REAL4 f0, m1Max, m1Min, m2Max, m2Min;

  /* Set up status pointer. */
  INITSTATUS( status, "LALInspiralSpinBankwNDTemplateBank", INSPIRALSPINBANKWNDTEMPLATEBANKC );
  ATTATCHSTATUSPTR( status );

  /* Check to make sure that all the parameters are okay */
  if (coarseIn.mmCoarse <= 0){
    ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
    }

  if ((coarseIn.mMin <= 0) || (coarseIn.MMax <= 0) ||
      (coarseIn.mMin >= coarseIn.MMax) || (3.0*coarseIn.MMax >= 15.0)){
    ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
    }

  /*These parameters have not been added to InspiralCoarseBankIn yet, but when they are the will need to be checked */
  /*
    if (coarseIn.betaMax < 0)
      ABORT(status, LALINSPIRALBANKH_ECHOICE, LALINSPIRALBANKH_MSGECHOICE);
  */

  /* Get 3x3 parameter-space metric. */
  /* BEN: mess creating all these structures & adding TRYs etc */
  /* BEN: do it by hand, since it's so simple? */

  /* allocate memory for the metric */
  LALU4CreateVector( status->statusPtr, &metricDimensions, (UINT4) 2 );
  metricDimensions->data[0] = 3;
  metricDimensions->data[1] = 3;
  LALSCreateArray( status->statusPtr, &metric, metricDimensions );


  /* Hardcode mass range etc for the moment. */
  m2Min = NDinput.minParameters[1] = coarseIn.mMin*LAL_MTSUN_SI;
  m2Max = NDinput.maxParameters[1] = coarseIn.MMax*LAL_MTSUN_SI;
  m1Min = NDinput.minParameters[0] = 2.0*m2Max;
  m1Max = NDinput.maxParameters[0] = 15.0*LAL_MTSUN_SI - m2Max;
  f0 = NDinput.f0 = 153.0; /*FIX THIS FIX THIS FIX THIS FIX THIS FIX THIS */

  /* Set NDinput parameters */
  NDinput.mm = coarseIn.mmCoarse;
  NDinput.type = PrecessingType;
  NDinput.dimension = 3;
  NDinput.PSD = &coarseIn.shf;

  /* Set box on unprimed coordinates including region. */
  /*psi0*/ NDinput.minCoordinates[0] = 0.9*(3.0/128) / (pow(LAL_PI*f0*(m1Max+m2Max),1.666667)*(m1Max*m2Max/pow(m1Max+m2Max,2)));
  /*psi3*/ NDinput.minCoordinates[1] = 1.1*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Max+m2Min),0.6666667)*(m1Max*m2Min/pow(m1Max+m2Min,2)));
  /*beta*/ NDinput.minCoordinates[2] = 0;
  /*psi0*/ NDinput.maxCoordinates[0] = 1.1*(3.0/128) / (pow(LAL_PI*f0*(m1Min+m2Min),1.666667)*(m1Min*m2Min/pow(m1Min+m2Min,2)));
  /*psi3*/ NDinput.maxCoordinates[1] = .9*(-.375*LAL_PI) / (pow(LAL_PI*f0*(m1Min+m2Max),0.6666667)*(m1Min*m2Max/pow(m1Min+m2Max,2)));
  /*beta*/ NDinput.maxCoordinates[2] = 3.8* LAL_PI/29.961432 * (1+0.75*m2Max/m1Min) * (m1Max/m2Min) * pow(LAL_MTSUN_SI*100.0/(m1Min+m2Min), 0.6666667);

  /* Set the function pointers for NDTemplate */
  NDFunctionPtrs.metric = LALInspiralSpinBankMetric;
  NDFunctionPtrs.test = LALInspiralSpinBankBoundary;

  printf("\ncalling LALNDTemplateBank()...\n");
  LALNDTemplateBank(status->statusPtr, &NDinput, &NDFunctionPtrs, &NDFirst);

  printf("\nconverting output...\n");
  NDoutput = NDFirst;
  while(NDoutput->next){
    (*ntiles)++;
    NDoutput = NDoutput->next;
    }

  printf("\n...counter is %i\n", *ntiles);

  *tiles = (InspiralTemplateList *) LALCalloc( *ntiles, sizeof(InspiralTemplateList));
  NDoutput = NDFirst;
  cnt = 0;

  for (cnt = 0; cnt < *ntiles; cnt++)
  {
    (*tiles)[cnt].params.mass1 = NDoutput->parameterVals[0];
    (*tiles)[cnt].params.mass2 = NDoutput->parameterVals[1];
    (*tiles)[cnt].params.psi0 = NDoutput->coordinateVals[0];
    (*tiles)[cnt].params.psi3 = NDoutput->coordinateVals[1];
    (*tiles)[cnt].params.beta = NDoutput->coordinateVals[2];
    /*(*tiles)[cnt].params.eta = tmplt->eta;
    (*tiles)[cnt].params.chirpMass = tmplt->chirpMass;*/
    NDoutput = NDoutput->next;
  } /* for(tmplt...) */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralSpinBankBoundary(
   LALStatus            *status,
   NDTemplateBankInput  *input,
   NDTemplateBankOutput *output,
   INT2 		*testFlag
   )
  {
  REAL4 f0, mass, eta, m1, m2, m1Min, m1Max, m2Min, m2Max, betaMax, x, y, z;
  x = output->coordinateVals[0];
  y = output->coordinateVals[1];
  z = output->coordinateVals[2];
  f0 = input->f0;
  m1Min = input->minParameters[0];
  m1Max = input->maxParameters[0];
  m2Min = input->minParameters[1];
  m2Max = input->maxParameters[1];
  betaMax = input->maxParameters[2];
  *testFlag =  1;

  mass = -y/x / (16.0*LAL_PI*LAL_PI*f0);
  eta = 16.0457 * pow( -x*x/y/y/y/y/y, 0.3333333 );
  if (eta > 0.25 || eta < 0)
    *testFlag = 0;
  output->parameterVals[0] = m1 = 0.5*mass* (1 + sqrt(1 - 4*eta));
  output->parameterVals[1] = m2 = 0.5*mass* (1 - sqrt(1 - 4*eta));
  if (m1 > m1Max || m1 < m1Min || m2 > m2Max || m2 < m2Min)
    *testFlag = 0;
  betaMax = 3.8*LAL_PI/29.961432 * (1+0.75*m2/m1)*(m1/m2) * pow((LAL_MTSUN_SI*100.0/mass),0.6666667);
  if (z > betaMax)
    *testFlag = 0;
  }






