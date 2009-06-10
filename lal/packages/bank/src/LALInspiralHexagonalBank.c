/*
*  Copyright (C) 2007 Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralHexagonalBankCV">
Author: Cokelaer Thomas
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralHexagonalBank.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralHexagonalBankCP}
\idx{LALInspiralHexagonalBank()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALInspiralParameterCalc()
LALInspiralComputeMetric()
\end{verbatim}

\subsubsection*{Notes}

\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralHexagonalBankCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>


/* Thomas:: definition for hexagonal grid. */


NRCSID(LALINSPIRALHEXAGONALBANKC, "$Id$");


void
LALInspiralCreatePNCoarseBankHexa(
    LALStatus            *status,
    InspiralTemplateList **list,
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    )
{
  INT4                  i;
  INT4 			firstId = 0;
  REAL4                 piFl;
  REAL4 		A0, A3;
  InspiralBankParams    bankPars;
  InspiralTemplate      *tempPars;
  InspiralMomentsEtc    moments;
  InspiralCell          *cells;
  HexaGridParam         gridParam;
  CellEvolution         cellEvolution;
  CellList 		*cellList = NULL;

  INITSTATUS( status, "LALInspiralHexagonalBank",
      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.mMin > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.mMax > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.MMax >= 2.*coarseIn.mMin, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  /* Set the elements of the metric and tempPars structures in  */
  /* conformity with the coarseIn structure                     */
  if ( !(tempPars = (InspiralTemplate *)
  			LALCalloc( 1, sizeof(InspiralTemplate))))
  {
    LALFree(tempPars);
    LALFree(cells);
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }


  LALInspiralSetParams( status->statusPtr, tempPars, coarseIn );
  CHECKSTATUSPTR( status );

  /* Identify the boundary of search and parameters for the     */
  /* first lattice point                                        */
  LALInspiralSetSearchLimits( status->statusPtr, &bankPars, coarseIn );
  CHECKSTATUSPTR( status );

  tempPars->totalMass   = coarseIn.MMax;
  tempPars->eta         = 0.25;
  tempPars->ieta        = 1.L;
  tempPars->fLower      = coarseIn.fLower;
  tempPars->massChoice  = m1Andm2;
  tempPars->mass1       = coarseIn.mMin;
  tempPars->mass2       = coarseIn.mMax;

  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );

  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric  once for all.   */
  LALGetInspiralMoments(
  		status->statusPtr,
  		&moments,
   		&coarseIn.shf,
   	 	tempPars );
  CHECKSTATUSPTR( status );

  /* Allocate memory for one cell */
  cells = (InspiralCell*)
      LALCalloc(1,   sizeof(InspiralCell) );

  /*define gridParam*/
  gridParam.mm 			= coarseIn.mmCoarse;
  gridParam.x0Min     	= bankPars.x0Min;
  gridParam.x0Max     	= bankPars.x0Max;
  gridParam.x1Min     	= bankPars.x1Min;
  gridParam.x1Max     	= bankPars.x1Max;
  gridParam.mMin      	= coarseIn.mMin;
  gridParam.mMax      	= coarseIn.mMax;
  gridParam.MMin      	= coarseIn.MMin;
  gridParam.MMax      	= coarseIn.MMax;
  gridParam.etaMin    	= coarseIn.etamin;
  gridParam.space     	= coarseIn.space;
  gridParam.massRange 	= coarseIn.massRange;
  gridParam.gridSpacing = coarseIn.gridSpacing;


  cellEvolution.nTemplate 		= 1;
  cellEvolution.nTemplateMax 	= 1;
  cellEvolution.fertile 		= 0;

  /* initialise that first cell */
  tempPars->massChoice  = t03;
  cells[0].t0           = tempPars->t0;
  cells[0].t3           = tempPars->t3;

  /* some aliases */
  piFl  = LAL_PI * tempPars->fLower;
  A0    = 5. / pow(piFl, 8./3.) / 256.;
  A3    = LAL_PI / pow(piFl, 5./3.)/8.;


  /* Initialise the first template */
  LALInitHexagonalBank(
  			status->statusPtr,
		       	&cells, firstId,
		       	&moments, tempPars,
		       	&gridParam, &cellEvolution,
		       	&cellList);
  CHECKSTATUSPTR( status );

  {
    INT4 k, kk; /*some indexes*/
    INT4 *list 		= NULL;
    CellList *ptr 	= NULL;
    INT4 length 	= 1; /* default size of the bank when we
    						start the bank generation. */

    /* we re-allocate an array which size equals the
     * template bank size. */
    if (! (list =  LALMalloc(length*sizeof(INT4))))
    {
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }

    /* while there are cells/template which can propagate, we carry on the loop.*/
    while (cellEvolution.fertile)
    {
      length = LALListLength(cellList);
      /*realloc some memory for the next template*/
      if (! (list =  LALRealloc(list, length*sizeof(INT4))))
      {
		ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
		/* freeing memory here ? */
      }
      ptr = cellList;
      /* we extract the ids which might change within the LALPopulateCell
       * function. Indeed the bank might grow and then we will lost track
       * of ids/bank size and so on. */
      for ( k = 0; k < length; k++)
      {
		list[k] = ptr->id;
		ptr = ptr->next;
      }
      /* look at all the template/ids in the current bank to search for fertile cells */
      for (kk = 0; kk < length; kk++)
	  {
		k = list[kk];
		if ( cells[k].status == Fertile)
		{
          LALPopulateCell(status->statusPtr, &moments, &cells,
              k,  tempPars, &gridParam, &cellEvolution, &cellList);
		  CHECKSTATUSPTR( status );
	  	  /* now the bank might have grown, but we only look at the
	  	   * template created before this for loop, when we entered
	  	   * in the while loop
	  	   * */
        }
      }
    }
    LALFree(list);
  }

  if (cellList != NULL)
    ABORT(status, LALINSPIRALBANKH_EHEXAINIT,LALINSPIRALBANKH_MSGEHEXAINIT);
	/* Here is the current number of template generated. Now, we need
	 * to clean some of them which might be redundant.
	 * */
  *nlist = cellEvolution.nTemplate;

  {
    INT4 k ;
    INT4 length;
    length = cellEvolution.nTemplate;

    for ( k = 0; k < length; k++)
    {
      REAL4 a;
      REAL4 b;
      REAL4 x0;
      REAL4 tempA3;
      SFindRootIn input;
      INT4 valid;

      PRIN  prin;

      tempA3              = pow(A3, -5./2.)/pow(0.25,-1.5);
      tempPars->t0        = cells[k].t0;
      tempPars->t3        = cells[k].t3;

      /* if non physical parameter i.e below eta=0.25*/
      if(cells[k].RectPosition[0] == Below )
      {
        INT4 above=0, below=0, in=0, out=0;

		/*first, we define the line which is along the long semi-axis of the
		 * ambiguity function, defined by the angle theta and the position of
		 * the template.
		 * */
		a = tan(cells[k].metric.theta);
		b = cells[k].t3 - a * cells[k].t0;
		/* and call a function to search for a solution along eta=1/4 */
		input.function 	= LALSPAF;
		input.xmin 		= cells[k].t3-1e-3;
		input.xmax 		= 1000;
		input.xacc 		= 1e-6;

		prin.ct = a * A0 * tempA3;
		prin.b = b;

		LALSBisectionFindRoot(status->statusPtr,
			&x0, &input, (void *)&prin);
		CHECKSTATUSPTR( status );

		tempPars->t3 = x0 + 1e-3; /* to be sure it is physical */
		tempPars->t0 = (tempPars->t3 - b)/a;
		if (tempPars->t0 > 0)
		{
	  	  LALInspiralParameterCalc(status->statusPtr, tempPars);
	  	  CHECKSTATUSPTR( status );
        }
		cells[k].t0  = tempPars->t0;
		cells[k].t3  = tempPars->t3;

		/* update its position values */
		valid = 1;
		GetPositionRectangle(status->statusPtr, &cells, k,  tempPars ,
			     &gridParam,
			     &cellEvolution,
			     &cellList,
			     &valid);

		{
	  		switch (cells[k].RectPosition[1]){
			  case In:    in    +=1; break;
			  case Below: below +=1; break;
			  case Above: above +=1; break;
			  case Out:   out   +=1; break;
			  }
			  switch (cells[k].RectPosition[2]){
			  case In:    in    +=1; break;
			  case Below: below +=1; break;
			  case Above: above +=1; break;
			  case Out:   out   +=1; break;
			  }
			  switch (cells[k].RectPosition[3]){
			  case In:    in    +=1; break;
			  case Below: below +=1; break;
			  case Above: above +=1; break;
			  case Out:   out   +=1; break;
			  }
			  switch (cells[k].RectPosition[4]){
			  case In:    in    +=1; break;
			  case Below: below +=1; break;
			  case Above: above +=1; break;
			  case Out:   out   +=1; break;
			  }
			}

		  if (above == 2 && cells[k].position == In)
		  {
		    cells[cells[k].child[0]].position = Out;
		  }
      }
    }
  }

  for (i=0; i<cellEvolution.nTemplate; i++) {
    if (cells[i].position == In ) {
      *nlist = *nlist +1;
    }
  }
;


  /* allocate appropriate memory and fill the output bank */
  *list = (InspiralTemplateList*)
    LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist+1) );
  if ( ! *list )
  {
    LALFree( tempPars );
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }
  memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );
  {
    *nlist = 0 ;
    for (i=0; i<cellEvolution.nTemplate; i++)
    {
      if (cells[i].position == In)
      {
        tempPars->t0  = cells[i].t0;
        tempPars->t3  = cells[i].t3;
        tempPars->massChoice = t03;
        tempPars->fLower = coarseIn.fLower;

        LALInspiralParameterCalc( status->statusPtr, tempPars );
        CHECKSTATUSPTR( status );

        (*list)[*nlist].ID            = *nlist;
        (*list)[*nlist].params        = *tempPars;
        (*list)[*nlist].metric        = cells[i].metric;
        ++(*nlist);
      }
    }
  }

  LALFree( cells );
  LALFree( tempPars );

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}





void
LALPopulateCell(
		LALStatus               *status,
		InspiralMomentsEtc      *moments,
		InspiralCell            **cell,
		INT4                     headId,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam,
		CellEvolution           *cellEvolution,
		CellList		**cellList
		)
{
  REAL4 dx0, dx1;
  REAL4 newt0, newt3;
  INT4 i, id1, id2;
  REAL4 theta, ctheta, stheta;
  INT4 offSpring;
  INT4 it;
  INT4 add = 0;

  INITSTATUS( status, "LALPopulateCell",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  /* aliases to get the characteristics of the parent template,
   * that we refer to its ID (headId) */
  dx0           = (*cell)[headId].dx0;
  dx1           = (*cell)[headId].dx1;
  theta         = (*cell)[headId].metric.theta;
  ctheta        = cos(theta);
  stheta        = sin(theta);
  offSpring     = cellEvolution->nTemplate;

   /* Around the parent, the offspring can be at most 6 (hexagonal grid).
   * By default the child are unset. If so it is created and have the
   * properties of its parents. However, a child migh have been created
   * earlier. In that case, we do not do anything.  */
  it = 0 ;

  for (i = 0; i < 6; i++)
  {
    if ((*cell)[headId].child[i] == -1)
    {
      add++;
      /* reallocate memory by set of 1000 cells if needed*/
      if ( (offSpring+add)>cellEvolution->nTemplateMax)
      {
        *cell = (InspiralCell*)
          LALRealloc( *cell,
           sizeof(InspiralCell) * (cellEvolution->nTemplateMax + 1000) );
        if ( !cell ) {
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        cellEvolution->nTemplateMax +=  1000;
      }

      /* creates the child connection if needed. A child heritates the
       * properties of its parent */
      switch ( i ){
      case 0:
		newt0   = dx0 ;
		newt3   = 0 ;
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 *ctheta + stheta* newt3;
		(*cell)[offSpring + it].t3   += newt0 *stheta - ctheta* newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
	  	break;
      case 1:
		newt0   =   dx0/2. ;
		newt3   =   -dx1 *sqrt(3./2) ;
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
		(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
		break;
      case 2:
		newt0   =  -dx0/2 ;
		newt3   =  -dx1 *sqrt(3./2);
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
		(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
		break;
      case 3:
		newt0   = -dx0 ;
		newt3   = 0;
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
		(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
		break;
      case 4:
		newt0   =  -dx0/2. ;
		newt3   =  dx1 *sqrt(3./2);
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
		(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
		break;
      case 5:
		newt0   = dx0/2. ;
		newt3   = dx1 *sqrt(3./2);
		(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
		(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
		(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
		(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
		LALInitHexagonalBank(status->statusPtr,  cell,  offSpring+it,
		    moments, paramsIn, gridParam, cellEvolution, cellList);
		break;
      }

      /* Now, tricky part, if a child has been creating, he must have a
       * connection with its parents and vice-versa.  */
      if ((*cell)[offSpring + it].child[(i+3)%6] == -1)
      {
		(*cell)[offSpring + it].child[(i+3)%6] = (*cell)[headId].ID;
		(*cell)[headId].child[i] = offSpring+it;
      }
      /* a new cell index */
      it += 1;
    }
  }

  cellEvolution->nTemplate +=it;


  /* Here, the parent has its 6 children set; he become sterile. */
  (*cell)[headId].status 	= Sterile;
  (cellEvolution->fertile) 	= cellEvolution->fertile-1;
  LALListDelete(cellList, headId);

  /* what shall we do with that parent. Is he valid ? inside the space,
   * outside since eta > 0.25 but close to the boundary .... */
  {
    if ((*cell)[headId].RectPosition[0] == Above && (*cell)[headId].in == 1)
    {
	  (*cell)[headId].RectPosition[0]=Out;
    }
  }

  /* propagate  connections to the brothers to avoid redundancies */
  for (i=0; i<6; i++)
  {
  	/* for each child*/
    id1 = (*cell)[headId].child[i%6];
    id2 = (*cell)[headId].child[(i+1)%6];
    (*cell)[id1].child[(i+2)%6] = (*cell)[id2].ID;
    (*cell)[id2].child[(i+4+1)%6] = (*cell)[id1].ID;
  }

  /* enfin trouver position[0] (In/out)? of the children. */
  for (i=0; i<6; i++)
  {/* for each child find position[0]*/
    id1 = (*cell)[headId].child[i%6];

    if ((*cell)[id1].status == Fertile)
    {
      LALSPAValidPosition(status->statusPtr, cell, id1,
			  moments, cellEvolution, cellList);
      CHECKSTATUSPTR( status );

      if ((*cell)[id1].position != In )
      {
        if ((*cell)[id1].status == Fertile)
        {
          (*cell)[id1].status= Sterile;
          cellEvolution->fertile=cellEvolution->fertile-1;
	  	  LALListDelete(cellList, id1);
        }
      }
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}



void
LALInitHexagonalBank(
	LALStatus               *status,
	InspiralCell            **cell,
	INT4                    id,
	InspiralMomentsEtc      *moments,
	InspiralTemplate        *paramsIn,
	HexaGridParam           *gridParam,
	CellEvolution           *cellEvolution,
	CellList **cellList)
{
  INT4          i;
  INT4 		valid;

  INITSTATUS( status, "LALInitHexagonalBank",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  /* a new cell is created; by default it can create new children,
     therefore it is fertile */
  cellEvolution->fertile = cellEvolution->fertile + 1;;
  (*cell)[id].status = Fertile;
  LALListAppend(cellList, id);


  /* all of whom are unset and do not have any id set yet*/
  for (i = 0; i < 6; i++)
  {
    (*cell)[id].child[i] = -1;
  }

  /* filled some values related to the space */
  (*cell)[id].ID        = id;
  (*cell)[id].position  = In;
  (*cell)[id].metric.space = gridParam->space;


  /* before any further computation, check that t0, t3 are positive.*/
  if ((*cell)[id].t0 > 0 && (*cell)[id].t3 > 0)
  {
    /* Get the metric at the position of the cell */
    paramsIn->t0 = (*cell)[id].t0;
    paramsIn->t3 = (*cell)[id].t3;

    LALInspiralComputeMetric( status->statusPtr,
			      &((*cell)[id].metric),
			      paramsIn,
			      moments);
    CHECKSTATUSPTR( status );

    /* let us store the dx0 and dx3 at that point. */
    (*cell)[id].dx0 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g00 );
    (*cell)[id].dx1 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g11 );

    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[0]), paramsIn, gridParam);
    CHECKSTATUSPTR( status );

    /* if outside, this is a sterile cell which can not propagate */
    if ((*cell)[id].RectPosition[0] == Out)
    {
      (*cell)[id].position      = Out;
      for (i = 0; i < 5; i++)
      {
        (*cell)[id].RectPosition[i] = Out;
      }
      (*cell)[id].status = Sterile;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      LALListDelete(cellList, id);

      DETATCHSTATUSPTR(status);
      RETURN(status);
    }
    else
    {
      valid = 1;
      GetPositionRectangle(status->statusPtr, &(*cell), id,  paramsIn ,
			   gridParam, cellEvolution, &(*cellList), &valid);
    }
  }
  else
  {/* if t0 or t3 < 0 , this is not a valid cell*/
    valid = 0;
  }

  /* If this is not a valid template, we remove it from the bank*/
  if (valid == 0)
  {
    for (i=0; i<5; i++)
    {
      (*cell)[id].RectPosition[i] = Out;
    }
    (*cell)[id].position 		= Out;
    (*cell)[id].status 			= Sterile;
    (cellEvolution->fertile)	=cellEvolution->fertile-1;
    LALListDelete(cellList, id);
  }




#if 1
  if (gridParam->gridSpacing == HybridHexagonal)
  {
    INT4 below=0, above=0;
    for (i=1; i<=4; i++){
      if ( (*cell)[id].RectPosition[i] == Below) below++;
      if ( (*cell)[id].RectPosition[i] == Above) above++;
    }
    if (below==2 && above == 2){
      (*cell)[id].status = Edge;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      LALListDelete(cellList, id);

    }


  }
#endif



  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/* Get the position of the rectangle corners which are inscribe within the ambiguity
 * function. Are they within the parameter space or not ?*/
void
GetPositionRectangle(
		LALStatus 		*status,
		InspiralCell 		**cell,
		INT4 			id,
		InspiralTemplate 	*params,
		    HexaGridParam 		*gridParam,
		    CellEvolution 		*cellEvolution,
		    CellList 			**cellList,
		    INT4 				*valid)
{
  RectangleIn   RectIn;
  RectangleOut  RectOut;
  InspiralTemplate paramsIn;

  INITSTATUS( status, "GetPosition",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  /* let us investigate this particular template : */
  RectIn.x0    = params->t0;
  RectIn.y0    = params->t3;
  RectIn.dx    = (*cell)[id].dx0 ;
  RectIn.dy    = (*cell)[id].dx1 ;
  RectIn.theta = (*cell)[id].metric.theta;

  /* what is the rectangle ? */
  LALRectangleVertices(status->statusPtr, &RectOut, &RectIn);
  CHECKSTATUSPTR( status );

  /* for each corner, let us decide where it lies in the parameter space */
  paramsIn = *params;
  paramsIn.t0 = RectOut.x1;
  paramsIn.t3 = RectOut.y1;

  if (RectOut.x1<0 || RectOut.y1<0
   || RectOut.x2<0 || RectOut.y2<0
   || RectOut.x3<0 || RectOut.y3<0
   || RectOut.x4<0 || RectOut.y4<0)
  {
   	*valid = 0;
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  if (RectOut.x1>0 && RectOut.y1>0)
  {
    LALFindPosition(status->statusPtr,(*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[1]),
		    &paramsIn,
		    gridParam);
    CHECKSTATUSPTR( status );
  }

  paramsIn.t0 = RectOut.x2;
  paramsIn.t3 = RectOut.y2;
  if (RectOut.x2>0 && RectOut.y2>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[2]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );

  }

  paramsIn.t0 = RectOut.x3;
  paramsIn.t3 = RectOut.y3;
  if (RectOut.x3>0 && RectOut.y3>0)
  {
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[3]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
  }

  paramsIn.t0 = RectOut.x4;
  paramsIn.t3 = RectOut.y4;
  if (RectOut.x4>0 && RectOut.y4>0)
  {
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[4]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
  }

  DETATCHSTATUSPTR( status );
  RETURN ( status );
}








void
LALSPAValidPosition(LALStatus *status,
		    InspiralCell **cell,
		    INT4 id1,
		    InspiralMomentsEtc *moments,
		    CellEvolution *cellEvolution,
		    CellList **cellList
		    )
{
  INT4 below = 0, in = 0, out = 0, above = 0;

  INITSTATUS( status, "LALSPAFindPosition",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );


  switch ((*cell)[id1].RectPosition[1]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[2]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[3]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[4]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }
  switch ((*cell)[id1].RectPosition[0]){
  case In:    in    +=1; break;
  case Below: below +=1; break;
  case Above: above +=1; break;
  case Out:   out   +=1; break;
  }

  (*cell)[id1].in = in;

  if ((*cell)[id1].RectPosition[0]==In)
    {
      (*cell)[id1].position = In;
      if ((*cell)[id1].status == Sterile)
	{
	  (*cell)[id1].status = Fertile;
	  (cellEvolution->fertile)=cellEvolution->fertile+1;;
	  LALListAppend(cellList, id1);
	}
      DETATCHSTATUSPTR(status);
      RETURN(status);
  }

  if ( above == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);

      }
  }
  else if ( below == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);

      }
  }
  else if ( out == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);

      }
  }
  else if (in >= 1){
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile+1;
	LALListAppend(cellList, id1);
      }

  }
  else if (above+below >= 5){
    if(out==1){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);

      }
    }
    else
    {
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);


      }
    }
  }
  else{
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	LALListDelete(cellList, id1);

      }
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALFindPosition(LALStatus       *status,
		REAL4                   dx0,
		REAL4                   dx1,
		Position                *position,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam
)
{
  REAL8 	mint3;
  REAL4   	eta;
  REAL4 	totalMass,ieta, oneby4, tiny, piFl, A0, A3;

  INITSTATUS( status, "LALFindPosition",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  ieta 	        = 1.;
  oneby4 	= 1./4.;
  tiny 	        = 1.e-10;
  piFl 	        = LAL_PI * paramsIn->fLower;
  A0    = 5. / pow(piFl, 8./3.) / 256.;
  A3    = LAL_PI / pow(piFl, 5./3.)/8.;

  ASSERT(paramsIn->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(paramsIn->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  /* given t0, t3 we get the totalMass and eta.
     We do not need to call ParameterCalc again and again here. */
  totalMass     = A0 * paramsIn->t3/(A3 * paramsIn->t0);
  eta           = A0/(paramsIn->t0 * pow(totalMass, fiveby3));

  /* be sure eta is inside the space if it is suppose to be */
  if (eta > oneby4) {
    eta-=tiny;
  }

  /* let us fill the param strucutre now : eta, Mass, mass1, mass2*/
  paramsIn->eta = eta;
  totalMass     = paramsIn->totalMass = totalMass/LAL_MTSUN_SI;
  if (eta <= oneby4) {
    paramsIn->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
    paramsIn->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
  }

  /* does t3 positive*/
  if ((paramsIn->t3-dx1)<0)
  {
    mint3 = 0;
  }
  else
  {
    mint3 = paramsIn->t3-dx1;
  }

  if ( (paramsIn->t0 <gridParam->x0Min - dx0)
      ||(paramsIn->t0 >gridParam->x0Max + dx0)
          || (paramsIn->t3 <= mint3))
  {
    *position = Out;
    DETATCHSTATUSPTR(status);
    RETURN(status);
  }

  switch ( gridParam->massRange )
  {
    case MinMaxComponentMass:
      if (
          paramsIn->mass1 >= gridParam->mMin &&
          paramsIn->mass2 >= gridParam->mMin &&
          paramsIn->mass1 <= gridParam->mMax &&
          paramsIn->mass2 <= gridParam->mMax &&
          paramsIn->eta <= 0.25 &&
          paramsIn->eta >= gridParam->etaMin
          )
        {
          *position = In;
        }
      else
        if (paramsIn->eta > .25){
          *position = Below;
        }
        else{
          *position = Above;
        }
      break;

    case MinComponentMassMaxTotalMass:
      if (
          paramsIn->mass1 >= gridParam->mMin &&
          paramsIn->mass2 >= gridParam->mMin &&
          paramsIn->totalMass <= gridParam->MMax &&
          paramsIn->eta <= 0.25 &&
          paramsIn->eta >= gridParam->etaMin
          )
        {
          *position = In;
        }
      else
        if (paramsIn->eta > .25){
          *position = Below;
        }
        else{
          *position = Above;
        }
      break;

    case MinMaxComponentTotalMass:
      if (
          paramsIn->mass1 >= gridParam->mMin &&
          paramsIn->mass2 >= gridParam->mMin &&
          paramsIn->totalMass <= gridParam->MMax &&
          paramsIn->totalMass >= gridParam->MMin &&
          paramsIn->eta <= 0.25 &&
          paramsIn->eta >= gridParam->etaMin
          )
        {
          *position = In;
        }
      else if (paramsIn->eta > .25 ){
          *position = Below;
        }
      else{
        *position = Above;
        }

      /* Now cut out unnecessary templates */
      if ( paramsIn->totalMass < gridParam->MMin )
      {
        REAL4 totalMass2 = A0 * (paramsIn->t3 - dx1)/(A3 * paramsIn->t0);
        totalMass2 = totalMass2 / LAL_MTSUN_SI;
        totalMass     = A0 * paramsIn->t3/(A3 * (paramsIn->t0 - dx0));
        totalMass = totalMass / LAL_MTSUN_SI;

        if ( totalMass < gridParam->MMin && totalMass2 < gridParam->MMin )
        {
          *position = Out;
        }
      }
      break;

    default:
      ABORT(status, 999, "Invalid choice for enum InspiralBankMassRange");
      break;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* This function corresponds to the eta=1/4 line? */
void LALSPAF(
	LALStatus 	*status,
	REAL4 		*result,
	REAL4 		t3,
	void 		*param)
{
  REAL4 ct, b;
  PRIN *prin;

  INITSTATUS( status, "LALSPAF",
	      LALINSPIRALHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  prin = (PRIN *)param;
  ct = prin->ct;
  b  = prin->b;

  *result = ct*pow(t3,5./2.) - t3 + b;

  DETATCHSTATUSPTR( status );
  RETURN(status);
}
