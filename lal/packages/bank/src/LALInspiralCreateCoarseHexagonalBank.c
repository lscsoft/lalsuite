/*  <lalVerbatim file="LALInspiralCreateCoarseBankCV">
Author: T.Cokelaer.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralCreateCoarseHexagonalBank.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralCreateCoarseBankCP}
\idx{LALInspiralCreateCoarseHexagonalBank()}

\begin{itemize}
   \item \texttt{list,} Output, an array containing the template bank parameters
   \item \texttt{nlist,} Output, the number of templates found by the function
   \item \texttt{coarseIn,} Input, specifies the search space, range of masses, etc.
\end{itemize}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
LALInspiralNextTemplate()
\end{verbatim}

\subsubsection*{Notes}
\clearpage


\vspace{0.1in}
\vfill{\footnotesize\input{LALInspiralCreateCoarseBankCV}}
</lalLaTeX>  */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>
#define Power(a,b) (pow((a),(b)))
#define Sqrt(a)    (sqrt((a)))
#define MT_SUN     (4.92549095e-6)

double solveForM (double x, double p, double q, double a);
double bisectionLine (double x, double fa, double mMin, double mMax);
double getSqRootArgument (double x, double p, double q, double a);


/* Thomas:: temporary definition for SPA hexagonal grid. */
static REAL4  A0;
static REAL4  A3;

typedef struct{
REAL4 ct;
REAL4 b;
}
PRIN;

void
LALPopulateNarrowEdge(LALStatus               *status,
		      InspiralMomentsEtc      *moments,
		      InspiralCell            **cell, 
		      INT4                     headId,
		      InspiralTemplate        *paramsIn,
		      HexaGridParam           *gridParam,
		      CellEvolution           *cellEvolution, 
		      CellList **cellList,
		      INT4 flag
		      );

static void LALSPAF(LALStatus *status,  REAL4 *result, REAL4 x, void *t3);
/* end comments ::Thomas */

NRCSID(LALINSPIRALCREATECOARSEHEXAGONALBANKC, "$Id$");

void 
LALInspiralCreatePNCoarseHexagonalBank(
    LALStatus            *status, 
    InspiralTemplateList **list, 
    INT4                 *nlist,
    InspiralCoarseBankIn coarseIn
    ) 
{  

  InspiralBankParams    bankPars;
  InspiralTemplate      *tempPars;
  InspiralMomentsEtc    moments;
  INT4                  i;
  InspiralCell          *cells;
  REAL4                 piFl;
  HexaGridParam         gridParam;
  CellEvolution         cellEvolution;
  INT4                  firstId=0;
  CellList              *cellList=NULL;

  INITSTATUS( status, "LALInspiralCreateCoarseHexagonalBank", 
      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
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
                LALCalloc( 1, sizeof(InspiralTemplate)))) {
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
  tempPars->mass1       = coarseIn.mMax;
  tempPars->mass2       = coarseIn.mMin;
  tempPars->approximant  = coarseIn.approximant;

  LALInspiralParameterCalc( status->statusPtr, tempPars );
  CHECKSTATUSPTR( status );
  
  /* Get the moments of the PSD integrand and other parameters */
  /* required in the computation of the metric  once for all.   */
  LALGetInspiralMoments( status->statusPtr, &moments, &coarseIn.shf, tempPars );
  CHECKSTATUSPTR( status );
  
  /* Allocate memory for one cell */
  cells = (InspiralCell*)
    LALCalloc(1,   sizeof(InspiralCell) );

  /*define gridParam*/
  gridParam.mm = coarseIn.mmCoarse;
  gridParam.x0Min     = bankPars.x0Min;
  gridParam.x0Max     = bankPars.x0Max;
  gridParam.x1Min     = bankPars.x1Min;
  gridParam.x1Max     = bankPars.x1Max;
  gridParam.mMin      = coarseIn.mMin;
  gridParam.mMax      = coarseIn.mMax;
  gridParam.etaMin    = coarseIn.etamin;
  gridParam.space     = coarseIn.space;
  gridParam.fLower    = coarseIn.fLower;

  cellEvolution.nTemplate = 0;
  cellEvolution.nTemplateMax = 1;
  cellEvolution.fertile = 0;

  /* initialise that first cell */

  tempPars->massChoice  = t03;
  cells[0].t0           = tempPars->t0;
  cells[0].t3           = tempPars->t3;

  /* some aliases set once for all (uses pow function...)*/
  piFl  = LAL_PI * tempPars->fLower;
  A0    = 5. / pow(piFl, 8./3.) / 256.;
  A3    = LAL_PI / pow(piFl, 5./3.)/8.;


  LALCellInit(status->statusPtr, 
	      &cells, firstId, 
	      &moments, tempPars,
	      &gridParam, &cellEvolution, 
	      &cellList, 0);
  CHECKSTATUSPTR( status );

/*  fprintf(stderr, "fertility : %d\n", cellEvolution.fertile);*/
  /* the main cell generator*/
  {
    INT4        k, kk;


    INT4       *sublist=NULL;
    CellList *ptr=NULL;
    INT4 length=1;

    if (! (sublist =  LALMalloc(length*sizeof(INT4))))
      {
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
	
      }

    while (cellEvolution.fertile) {
      length = Length(cellList);
      
      
      if (! (sublist =  LALRealloc(sublist, length*sizeof(INT4))))
      {
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
	/* freeing memory here ? */
      }
      ptr = cellList;

      for (k=0; k< length; k++){
	sublist[k]=ptr->id;	
	ptr = ptr->next;
      }      

      for (kk = 0; kk < length; kk++) 
	{	
	  k = sublist[kk];
	  
	  if ( cells[k].status == Fertile) {
	    
	    LALPopulateCell(status->statusPtr, &moments, &cells,
			    k,  tempPars, &gridParam, &cellEvolution, &cellList);          
	    CHECKSTATUSPTR( status );         	 	  
	  }
	}
    }
    LALFree(sublist);
  }
  


#if 1 

  if (cellList != NULL)
    printf("wierd behaviour here\n");

  {
    INT4 edge1=0, edge2=0; 
    i=0;
    while (i<cellEvolution.nTemplate){    
      if (cells[i].status == Edge){
	edge1 = i;
	cells[i].status = In;   
	i=cellEvolution.nTemplate;
      }
      i++;
    }
    i=0;
    while (i<cellEvolution.nTemplate){    
      if (cells[i].status == Edge){
	edge2=i;
	cells[i].status = In;   
	i=cellEvolution.nTemplate;
      }
      i++;
    }
  
    
    if (cells[edge1].t0 > cells[edge2].t0){
      LALPopulateNarrowEdge(status->statusPtr, &moments, &cells,
			    edge1,  tempPars, &gridParam, &cellEvolution, &cellList, 0);          
      CHECKSTATUSPTR( status );          
      LALPopulateNarrowEdge(status->statusPtr, &moments, &cells,
			    edge2,  tempPars, &gridParam, &cellEvolution, &cellList, 1);          
      CHECKSTATUSPTR( status );          
    }
    else{
      LALPopulateNarrowEdge(status->statusPtr, &moments, &cells,
			    edge1,  tempPars, &gridParam, &cellEvolution, &cellList, 1);          
      CHECKSTATUSPTR( status );          
      LALPopulateNarrowEdge(status->statusPtr, &moments, &cells,
			    edge2,  tempPars, &gridParam, &cellEvolution, &cellList, 0);          
      CHECKSTATUSPTR( status );          
      
    }
  }


#endif

  *nlist = cellEvolution.nTemplate;
  /*move back the template below the eta=1/4 boundary. i put it here as a single 
    block so it can be removed if one accept non physical template position */
  {
    INT4 k ;
    INT4 length;
    length = cellEvolution.nTemplate;

    for (k=0; k<length; k++)
      {  
	REAL4  a,b, x0, tempA3;
	SFindRootIn input;

	PRIN  prin;
	
	tempA3              = pow(A3, -5./2.)/pow(0.25,-1.5);
	tempPars->t0        = cells[k].t0;
	tempPars->t3        = cells[k].t3; 
               /* a safety factor to avoid numerical errors. Sometimes the 
                * template could be slighlty below the et=1/4 line and therefore
                * might be considerer as above. We do not want that so we arti-
                * ficially push the template slightlty below its actual position
                * before computing its position. */
	
	if(cells[k].RectPosition[0] == Below ) {
	  
	  a = tan(cells[k].metric.theta);
	  b = cells[k].t3 - a * cells[k].t0;
	  
	  input.function = LALSPAF;
	  input.xmin = cells[k].t3-1e-3;
	  input.xmax = 1000;
	  input.xacc = 1e-6;
	  
	  prin.ct = a * A0 * tempA3;
	  prin.b = b;


	  /* could fail. if so something is wrong before anyway.*/
	  LALSBisectionFindRoot(status->statusPtr,&x0, &input, (void *)&prin);
	  CHECKSTATUSPTR( status );         
	  
	  tempPars->t3 = x0 + 1e-3; /*I do not like that arbitrary number but 
				      it works right now. it should be related 
				      to the metric component though.*/ 
	  tempPars->t0 = (tempPars->t3 - b)/a;
	  if (tempPars->t0 > 0) {
	    LALInspiralParameterCalc(status->statusPtr, tempPars);
	    CHECKSTATUSPTR( status );         		  
	   
	    cells[k].t0  = tempPars->t0;
	    cells[k].t3  = tempPars->t3;    
	  }
	  else{
	    /*todo*/
	  }
	  
	} 
      }
  }
  
  for (i=0; i<cellEvolution.nTemplate; i++) {
    if (cells[i].position == In ) {
      *nlist = *nlist +1; 
    }
  }
  
  /* allocate appropriate memory */
  
  
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
    for (i=0; i<cellEvolution.nTemplate; i++) {
      if (cells[i].position == In) {
	tempPars->t0  = cells[i].t0;
	tempPars->t3  = cells[i].t3;
	tempPars->massChoice = t03;
	tempPars->fLower = coarseIn.fLower;
	
	LALInspiralParameterCalc( status->statusPtr, tempPars );	
	CHECKSTATUSPTR( status );

	while (tempPars->eta>.25){
	  tempPars->t3+=1e-4;
	  /* let's put a warning or info here*/
	  /*printf("%f %f %f",	tempPars->t0,	tempPars->t3 , tempPars->eta);*/
	  LALInspiralParameterCalc( status->statusPtr, tempPars );	
	  CHECKSTATUSPTR( status );
	}


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
LALPopulateCell(LALStatus               *status,
		InspiralMomentsEtc      *moments,
		InspiralCell            **cell, 
		INT4                     headId,
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam,
		CellEvolution           *cellEvolution, 
		CellList **cellList
		)
{
  REAL4 dx0, dx1,  newt0, newt3;  
  INT4 i, id1, id2;
  REAL4 theta, ctheta,stheta;
  INT4 offSpring;
  INT4 it;
  INT4 add=0;

  INITSTATUS( status, "LALPopulateCell", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  /* here are the distance to be used and orientation of the ellipses*/
  dx0 = (*cell)[headId].dx0/sqrt(2.);
  dx1 = (*cell)[headId].dx1/sqrt(2.);
  theta = (*cell)[headId].metric.theta ;

  /* aliases */
  ctheta        = cos(theta);
  stheta        = sin(theta);
  offSpring     = cellEvolution->nTemplate;

  
  /* 
     The input arguments provided the properties of the genitor. 
     Around the genitor, there are 6 connections allowed which 
     corresponds to an offspring to be populated. Since the 
     offspring populated an hexagonal space. Some of the children 
     might already exists (flag set) and therefore the offspring 
     can be less than 6.

  */
     
  it = 0 ; 
      
  for (i = 0; i < 6; i++) {
    if ((*cell)[headId].child[i] == -1) {
      add++;
      /* reallocate memory by set of 1000 cells if needed*/
      if ( (offSpring+add)>cellEvolution->nTemplateMax){
        *cell = (InspiralCell*) 
          LALRealloc( *cell, sizeof(InspiralCell) * (cellEvolution->nTemplateMax + 1000) );
        if ( ! cell ) {
          ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
        }
        cellEvolution->nTemplateMax += 1000;
      }
      
      /* since we are dealing with ellipses and not circles, the 
	 hexagonal packing is not unique. We have done tests starting 
	 from t0=dx0, t3=0  and then t0 = dx0*1.5 and t3=dx1*sqrt(3)/2. 
	 Results are quite similar as expected but slightlty better in 
	 the second case that we are using here. 
      */

      switch ( i ){
      case 0:
	newt0   = dx0 ;
	newt3   = 0 ;
	newt0   = dx0 *1.5;
	newt3   = dx1*sqrt(3.)/2. ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 *ctheta + stheta* newt3;
	(*cell)[offSpring + it].t3   += newt0 *stheta - ctheta* newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      case 1:
	newt0   =   dx0/2. ;
	newt3   =   -dx1 *sqrt(3./2) ;
	newt0   =   dx0*1.5 ;
	newt3   =   -dx1 *sqrt(3.)/2. ;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      case 2:
	newt0   =  -dx0/2 ;
	newt3   =  -dx1 *sqrt(3./2);
	newt0   =   0;
	newt3   =  -dx1 *sqrt(3.);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      case 3:
	newt0   = -dx0 ;
	newt3   = 0;
	newt0   = -dx0 *1.5;
	newt3   = -dx1*sqrt(3.)/2.;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      case 4:
	newt0   =  -dx0/2. ;
	newt3   =  dx1 *sqrt(3./2);
	newt0   = -dx0 *1.5;
	newt3   = dx1*sqrt(3.)/2.;
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      case 5:
	newt0   = dx0/2. ;
	newt3   = dx1 *sqrt(3./2);
	newt0   =   0;
	newt3   =  dx1 *sqrt(3.);
	(*cell)[offSpring + it].t0   = (*cell)[headId].t0;
	(*cell)[offSpring + it].t3   = (*cell)[headId].t3;
	(*cell)[offSpring + it].t0   += newt0 * ctheta + stheta * newt3;
	(*cell)[offSpring + it].t3   += newt0 * stheta - ctheta * newt3;
	LALCellInit(status->statusPtr,  cell,  offSpring+it, 
		    moments, paramsIn, gridParam, cellEvolution, cellList,0);
	break;
      }      
      
      /* Now, tricky part, if a child has been creating, he must have a
       * connection with its parents and vice-versa.  */
      if ((*cell)[offSpring + it].child[(i+3)%6] == -1){
	(*cell)[offSpring + it].child[(i+3)%6] = (*cell)[headId].ID;
	(*cell)[headId].child[i] = offSpring+it;
      }
      /* a new cell index */
      it += 1;
    }
  }
  

  /* Here, the parent has its 6 children set; he become sterile. */
  (*cell)[headId].status = Sterile;
  (cellEvolution->fertile)=cellEvolution->fertile-1;
  Delete(cellList, headId);


  /* propagate  connection annexe aux freres pour eviter redondance */  
  for (i=0; i<6; i++){/* for each child*/
    id1 = (*cell)[headId].child[i%6];
    id2 = (*cell)[headId].child[(i+1)%6];
    (*cell)[id1].child[(i+2)%6] = (*cell)[id2].ID;
    (*cell)[id2].child[(i+4+1)%6] = (*cell)[id1].ID;   
  }
  
  /* enfin trouver position[0] (In/out)? of the children. */
  /* not for the first case */
  for (i=0; i<6; i++){/* for each child find position[0]*/
    id1 = (*cell)[headId].child[i%6];

    if ((*cell)[id1].status == Fertile) {
      LALSPAValidPosition(status->statusPtr, cell, id1, 
			  moments, cellEvolution, cellList);
      CHECKSTATUSPTR( status );


      if ((*cell)[id1].position != In ) {
        if ((*cell)[id1].status == Fertile) {
          (*cell)[id1].status= Sterile;
          cellEvolution->fertile=cellEvolution->fertile-1;
	  Delete(cellList, id1);
        }
      }
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN ( status );

}


void
LALPopulateNarrowEdge(LALStatus               *status,
		      InspiralMomentsEtc      *moments,
		      InspiralCell            **cell, 
		      INT4                     headId,
		      InspiralTemplate        *paramsIn,
		      HexaGridParam           *gridParam,
		      CellEvolution           *cellEvolution, 
		      CellList **cellList, 
		      INT4                    flag
		      )
{
  REAL4 dx0, dx1;  
  REAL4 theta, ctheta,stheta;
  INT4 offSpring;
  INT4 next, iteration;
  REAL4 x_int, y_int,xr_int, yr_int, c,s, dy,theta_min, theta_max, theta_int, a, b, t0, t3;

  INITSTATUS( status, "LALPopulateNarrowEdge",
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  /* aliases to get the characteristics of the parent template, that we refer
   * to its ID (headId) */  
  


  while ( (*cell)[headId].t0 < gridParam->x0Max && (*cell)[headId].t0 > gridParam->x0Min) {
    dx0           = (*cell)[headId].dx0/sqrt(2.);
    dx1           = (*cell)[headId].dx1/sqrt(2.);
    theta         = (*cell)[headId].metric.theta;
    ctheta        = cos(theta); /*aliases*/
    stheta        = sin(theta); /*aliases*/
    offSpring     = cellEvolution->nTemplate;
    

    /* reallocate memory by set of 1000 cells if needed*/
    if ( cellEvolution->nTemplate  >= cellEvolution->nTemplateMax){
      *cell = (InspiralCell*) 
	LALRealloc( *cell, sizeof(InspiralCell) * (cellEvolution->nTemplateMax + 1000) );
      if ( ! cell ) {
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      cellEvolution->nTemplateMax +=  1000;
    }
    

    next = cellEvolution->nTemplate;
    

    theta_min = 0.+.1 +LAL_PI/2.;
    theta_max = 2*LAL_PI-.1+LAL_PI/2.;

    theta_min = 0.+.1 ;
    theta_max = 2*LAL_PI-.1;

    t0 = (*cell)[headId].t0;
    t3 = (*cell)[headId].t3;
    a  = dx0 * sqrt(3.);
    b  = dx1 * sqrt(3.);
    c = cos(theta);
    s = sin(theta);    


    iteration = 1;
    while (fabs(theta_max-theta_min)>(.1/180.*LAL_PI) && iteration<20){

      theta_int = (theta_max + theta_min)/2.;

      xr_int = a*cos(theta_int);
      yr_int = b*sin(theta_int);
      
      x_int = xr_int *c - yr_int * s +t0;
      y_int = xr_int *s + yr_int * c +t3;
      
      dy = y_int -  bisectionLine(x_int, gridParam->fLower, gridParam->mMin,gridParam->mMax);
      /*printf("# %f %f\n",x_int, y_int);*/

      if (flag==0){
	if (dy>0 )
	  theta_max = theta_int;
	else
	  theta_min = theta_int;
      }
      else{
	if (dy>0 )
	  theta_min = theta_int;
	else
	  theta_max = theta_int;
      }
      iteration++;
    }

    
    (*cell)[next].t0   = x_int;
    (*cell)[next].t3   = y_int;

    if ( (*cell)[next].t0  > gridParam->x0Max ){
      (*cell)[next].t0 = gridParam->x0Max;
      (*cell)[next].t3 = bisectionLine(gridParam->x0Max, gridParam->fLower, gridParam->mMin,gridParam->mMax);
    }

    if ( (*cell)[next].t3  > gridParam->x1Max ){
      (*cell)[next].t0 = gridParam->x0Max;
      (*cell)[next].t3 = bisectionLine(gridParam->x0Max, gridParam->fLower, gridParam->mMin,gridParam->mMax);
    }

    if ( (*cell)[next].t0  < gridParam->x0Min ){
	  (*cell)[next].t0 = gridParam->x0Min;
	  (*cell)[next].t3 = bisectionLine(gridParam->x0Min, gridParam->fLower, gridParam->mMin,gridParam->mMax);
    }

    LALCellInit(status->statusPtr,  cell,  next, 
		moments, paramsIn, gridParam, cellEvolution, cellList,1);



    (*cell)[next].status = Sterile;
    (cellEvolution->fertile)=cellEvolution->fertile-1;
    Delete(cellList, next);
    headId=next;
    
    
  }
 (*cell)[headId].status = Sterile;
    (cellEvolution->fertile)=cellEvolution->fertile-1;
    Delete(cellList, headId);	
  DETATCHSTATUSPTR( status );
  RETURN ( status );

}


void 
LALCellInit(    LALStatus               *status,
                InspiralCell            **cell, 
                INT4                    id,
                InspiralMomentsEtc      *moments, 
                InspiralTemplate        *paramsIn, 
		HexaGridParam           *gridParam, 
		CellEvolution           *cellEvolution,
		CellList                **cellList, 
                INT4                    type)
{
  
  INT4          i;
  INT4 valid;  
  REAL4 totalMass;
  INITSTATUS( status, "LALCellInit", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );
  
  /* a new cell is created; by default it can create new children, 
     therefore it is fertile */
  cellEvolution->fertile = cellEvolution->fertile + 1;
  (*cell)[id].status = Fertile;  
  append(cellList, id);


  /* all of whom are unset and do not have any id set for the time being */
  for (i = 0; i < 6; i++) {
    (*cell)[id].child[i] = -1;
  } 
 
  /* filled some values related to the space */
  (*cell)[id].ID        = id;  /* here is a unique id number which*/
  (*cell)[id].position  = In;  /* we assume its position inside the space */
  (*cell)[id].metric.space = gridParam->space; /* heritates the gridParam.space*/


  /* before any further computation, check that t0, t3 are positive.*/
  if ((*cell)[id].t0 > 0 && (*cell)[id].t3 > 0){
    /* Get the metric at the position of the cell */ 
    paramsIn->t0 = (*cell)[id].t0;
    paramsIn->t3 = (*cell)[id].t3;

    LALInspiralComputeMetric( status->statusPtr, 
			      &((*cell)[id].metric),
			      paramsIn,
			      moments);
    CHECKSTATUSPTR( status );
  

    
    totalMass = A0 * paramsIn->t3/(A3 * paramsIn->t0);
    totalMass = totalMass/LAL_MTSUN_SI;

    if (type==0){ /*we are in the hexa grid placement*/
            if (gridParam->mMax >3 ) /*bbh and bhns case*/
                    {        
                      (*cell)[id].metric.g00 /= 0.80;
                      (*cell)[id].metric.g11 /= 0.80;
                      (*cell)[id].metric.theta += 15./180.*LAL_PI*0;
                    }
            else /*bns case*/
                    {
                      (*cell)[id].metric.g00 /= 0.93;
                      (*cell)[id].metric.g11 /= 0.93;
                      (*cell)[id].metric.theta += 15./180.*LAL_PI*0;    
                    }
    }
    else
    {
        /* we are in the narrow edges*/
          /* and look at BBH*/
            if ((gridParam->mMax >3 )&&(gridParam->mMin>=3)) /*bbh and bhns case*/
            {
              (*cell)[id].metric.g00 /= 0.85;
              (*cell)[id].metric.g11 /= 0.85;
              
            }
    }

      
    (*cell)[id].dx0 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g00 );
    (*cell)[id].dx1 = sqrt(2.L * (1.L - gridParam->mm)/(*cell)[id].metric.g11 );

    
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[0]), paramsIn, gridParam);
    CHECKSTATUSPTR( status );

    valid = 1;
      GetPositionRectangle(status->statusPtr, &(*cell), id,  paramsIn , 
			   gridParam,  &valid);
   
  }
  else{/* if t0 or t3 < 0 , this is not a valid cell*/
    valid = 0;   
  }


  
/* fprintf(stdout, "valid of id %d = %d", id, valid);
 fflush(stdout);*/
  if (valid == 0){
    for (i=0; i<5; i++){(*cell)[id].RectPosition[i] = Out;}
    (*cell)[id].position = Out;
    (*cell)[id].status = Sterile;
    (cellEvolution->fertile)=cellEvolution->fertile-1;
    Delete(cellList, id);
  }

  /* figure out if the template already covers the whole parameter space.*/
#if 1 
  {
    INT4 below=0, above=0;
    
    for (i=1; i<=4; i++){
      if ( (*cell)[id].RectPosition[i] == Below) below++;
      if ( (*cell)[id].RectPosition[i] == Above) above++;
    }
    
    if (below==2 && above == 2){
      (*cell)[id].status = Edge;
      (cellEvolution->fertile)=cellEvolution->fertile-1;
      Delete(cellList, id);
    }
  }
#endif

  cellEvolution->nTemplate++;
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}








void
GetPositionRectangle(LALStatus *status, 
		     InspiralCell **cell,
		     INT4 id,
		     InspiralTemplate *params, 
		     HexaGridParam *gridParam, 
		     INT4 *valid)
{
  RectangleIn   RectIn;
  RectangleOut  RectOut;
  InspiralTemplate paramsIn;

  INITSTATUS( status, "GetPositionRectangle", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  RectIn.x0    = params->t0;
  RectIn.y0    = params->t3;
  RectIn.dx    = (*cell)[id].dx0 *1.1; /*safety factor to be kept */
  RectIn.dy    = (*cell)[id].dx1 *1.1; /*safety factor to be kept */
  RectIn.theta = (*cell)[id].metric.theta;
  
  LALRectangleVertices(status->statusPtr, &RectOut, &RectIn);
  CHECKSTATUSPTR( status );

  paramsIn = *params;

  paramsIn.t0 = RectOut.x1;
  paramsIn.t3 = RectOut.y1;

  /*default values in case of any problem below.*/
  (*cell)[id].RectPosition[1]=Below;
  (*cell)[id].RectPosition[2]=Below;
  (*cell)[id].RectPosition[3]=Below;
  (*cell)[id].RectPosition[4]=Below;
  
  if (RectOut.x1>0 && RectOut.y1>0){
    LALFindPosition(status->statusPtr,(*cell)[id].dx0, (*cell)[id].dx1, 
		    &((*cell)[id].RectPosition[1]), 
		    &paramsIn, 
		    gridParam);    
    CHECKSTATUSPTR( status );
  }
  else {
    *valid = 0; 
    DETATCHSTATUSPTR(status);
    RETURN(status);     
  }
  
  paramsIn.t0 = RectOut.x2;
  paramsIn.t3 = RectOut.y2;
  if (RectOut.x2>0 && RectOut.y2>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[2]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
    
  }
  else
    {
      *valid = 0;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
  
  paramsIn.t0 = RectOut.x3;
  paramsIn.t3 = RectOut.y3; 
  if (RectOut.x3>0 && RectOut.y3>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[3]), &paramsIn, gridParam);
    CHECKSTATUSPTR( status );
  }
  else
    {
      *valid = 0 ;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
  
  paramsIn.t0 = RectOut.x4;
  paramsIn.t3 = RectOut.y4;
  if (RectOut.x4>0 && RectOut.y4>0){
    LALFindPosition(status->statusPtr, (*cell)[id].dx0, (*cell)[id].dx1,
		    &((*cell)[id].RectPosition[4]), &paramsIn, gridParam); 
    CHECKSTATUSPTR( status );
  }
  else
    {
      *valid = 0;
      DETATCHSTATUSPTR(status);
      RETURN(status);     
    }
  DETATCHSTATUSPTR( status );
  RETURN ( status );
}

/* that function returns the position of the template inside the parameter
 * space and update the status of the cell (fertile or sterile)*/
void
LALSPAValidPosition(LALStatus *status, 
		    InspiralCell **cell,
		    INT4 id1,
		    InspiralMomentsEtc *moments, 
		    CellEvolution *cellEvolution, 
		    CellList **cellList
		    )
{
  INT4 below=0, in=0, out=0, above=0;
  
  INITSTATUS( status, "LALSPAFindPosition", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
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
	  append(cellList, id1);
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
	Delete(cellList, id1);

      }
  }
  else if ( below == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }  
  else if ( out == 5){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }
  else if (in >= 1){
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile+1;
	append(cellList, id1);

  
      }

  }
  else if (above+below >= 5){
    if(out==1){
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
    }
    else
    {
    (*cell)[id1].position = In;
    if ((*cell)[id1].status == Sterile)
      {
	(*cell)[id1].status = Fertile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);


      }
    }  
  }
  else{
    (*cell)[id1].position = Out;
    if ((*cell)[id1].status == Fertile)
      {
	(*cell)[id1].status = Sterile;
	(cellEvolution->fertile)=cellEvolution->fertile-1;
	Delete(cellList, id1);

      }
  }

	    

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*return position of the template position with respect to the 
  parameter space boundaries. Inside, above or below by using 
  its position t0 and t3 as well as the lower cutoff frequency
  and parameters of min and max mass components.*/
void
LALFindPosition(LALStatus               *status, 
		REAL4                   dx0, 
		REAL4                   dx1,
		Position                *position, 
		InspiralTemplate        *paramsIn,
		HexaGridParam           *gridParam
)

{
  REAL8 mint3;  
  REAL4   eta, totalMass,ieta, oneby4, tiny, piFl;

  INITSTATUS( status, "LALFindPosition", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );


  ieta 	        = 1.;
  oneby4 	= 1./4.;
  tiny 	        = 1.e-4;
  piFl 	        = LAL_PI * paramsIn->fLower;
  
  ASSERT(paramsIn->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(paramsIn->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  /* given t0, t3 we get the totalMass and eta. 
     We do not need to call ParameterCalc again and again here.
     it gains 50% of computation time !!! */
  paramsIn->t3-=tiny;  /* be sure that if close to eta=1/4 it wont be 
                        * considered as inside. Just push it outside 
                        * to be sure to keep it. 
                        */
  totalMass = A0 * paramsIn->t3/(A3 * paramsIn->t0);
  eta       = A0/(paramsIn->t0 * pow(totalMass, fiveby3));
  paramsIn->t3 +=tiny; /* get back original value.*/

  
  /* let us fill the param structure now : eta, Mass, mass1, mass2*/
  paramsIn->eta = eta;
  totalMass     = paramsIn->totalMass = totalMass/LAL_MTSUN_SI;
  if (eta <= oneby4) {
    paramsIn->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
    paramsIn->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
  }
  
  /* does t3 positive*/  
  if ((paramsIn->t3-dx1)<0){ 
    mint3 = 0;
  }
  else{
    mint3 = paramsIn->t3-dx1;
  }
  
  if ( 
      (paramsIn->t0 <gridParam->x0Min - dx0)
      ||(paramsIn->t0 >gridParam->x0Max + dx0) 
          || (paramsIn->t3 <= mint3))
    {
      *position = Out;
      DETATCHSTATUSPTR(status);
      RETURN(status);
    }   

  if (
      paramsIn->mass1 >= gridParam->mMin-tiny &&
      paramsIn->mass2 >= gridParam->mMin-tiny &&
      paramsIn->mass1 <= gridParam->mMax+tiny &&
      paramsIn->mass2 <= gridParam->mMax+tiny &&
      paramsIn->eta <= 0.25+tiny && 
      paramsIn->eta >= gridParam->etaMin-tiny
      ) 
    {
      *position = In;
    }
  else
    if (paramsIn->eta > .25+tiny){
      *position = Below; 
    }
    else{
      *position = Above;
    }    
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}




/* function which computes the eta=1/4 boundary line.*/
static void LALSPAF(LALStatus *status,  
		    REAL4 *result, 
		    REAL4 t3, 
		    void *param)
{
  REAL4 ct, b;
  PRIN *prin;

  INITSTATUS( status, "LALSPAF", 
	      LALINSPIRALCREATECOARSEHEXAGONALBANKC );
  ATTATCHSTATUSPTR( status );

  prin = (PRIN *)param;
  ct = prin->ct;
  b  = prin->b;

  *result = ct*pow(t3,5./2.) - t3 + b;

  DETATCHSTATUSPTR( status );  
  RETURN(status);
}



/* The following function are  not "lalified" if someone wants to do it ... */

/* The following small functions are an implementation and nested list usefule
 * to generate the template bank interactively. */

void
print_list(CellList *head)
{
  if (head==NULL){
    printf("\n");
  }
  else { 
    printf(" %d", head->id);
    print_list(head->next);
  }
}

/* return the length of the list.*/
int Length(CellList *list)
{

  int count = 0; 

  while (list!=NULL){
    count++;
    list = list->next;
  }
  return count;
}

/* Append a cell to the list*/
void append(CellList **headRef, INT4 id)
{
  CellList *current;

  if ((current = malloc(sizeof(*current))) == NULL) {
    {
      printf("Error with malloc\n");
      exit(0);
    }
  }

  current->id = id;
  current->next = *headRef;
  *headRef = current;

}

/* delete the entire list*/
void DeleteList(CellList **headRef)
{
  CellList *tmp;
  
  while (headRef !=NULL){
    tmp = (*headRef)->next;
    free(headRef);
    (*headRef) = tmp;
  }
  
}


/*delete one particulat list*/
void Delete(CellList **headRef, INT4 id)
{
  CellList *ptr  = NULL;
  CellList *prev = NULL;


  for (ptr = *headRef; ptr != NULL; ptr= ptr->next){
    if (id == ptr->id)
      break;

    prev = ptr;
  }

  if (ptr==NULL)
    return ;

  if (prev!=NULL){
    prev->next = ptr->next;
  }
  else{
    *headRef = ptr->next;
  }

  /* free the data here if needed such as free(ptr->id); */
  free(ptr);

  return ;


 
}

/* code from Anand Sengupta which  estimate the bissectrice line to the binary parameter space
 * boundaries. Useful to place template optimally when only one template can
 * cover the eta=1/4 line as well as the (m1=min, m2)or (m1=max, m2) lines.  */
/* it has to be written properly at some points, but works pretty well at the
 * moment. */
/******************************************************************/
double bisectionLine (double x, double fa, double mMin, double mMax)
{

    double pi, alpha, beta, piFa;
    double y1, y2, p, q, M, S, eta, xbndry;

    double A, B;

    A=5./256/LAL_PI/fa;
    B=1./8./fa;

    /*    return pow(.25, -3./5.)*B*pow(x/A, 2./5.);*/

    pi   = 4.0 * atan(1.0);
    piFa = pi * fa;

    alpha = (5.0 / 256.0) * pow(piFa, (-8.0/3.0));
    beta  = 1.0 / (8.0 * fa * pow(piFa, (2.0/3.0)));

    /* First we solve for the lower (equal mass) limit */
    y1 = 4.0*beta;
    y1 *= (pow(x/(4.0*alpha),2.0/5.0));

    /* Figure out the boundary between m1 = mMin and m1 = mMax */
    M   = mMin + mMax;
    eta = (mMin*mMax)/pow(M, 2.0);
    xbndry = alpha * pow(M*MT_SUN, -5.0/3.0) / eta;

    /* Next we solve for the upper part */
    p = pow(x*mMin/alpha, -3.0/5.0);
    if (x >= xbndry )
          q = mMin;
    else
          q = mMax;
    M = solveForM (x, MT_SUN, q, alpha);
    S = getSqRootArgument (x, MT_SUN, q, alpha);
    if (x >= xbndry ) 
          eta = mMin*(M-mMin)/ pow(M,2.0);
    else
          eta = mMax*(M-mMax)/ pow(M,2.0);

    y2 = beta * (pow(M*MT_SUN, -2.0/3.0)) / eta;

    /*    fprintf (stderr, "M = %e SqRootArg = %e eta = %e y1 = %e y2 = %e "
	  "xbndry = %e y = %e\n", M, S, eta, y1, y2, xbndry, 0.5*(y1+y2));*/

    return 0.5*(y1+y2);
}


double solveForM (double x, double p, double q, double a)
    /* x = tau_0
     * p = solarmass in sec
     * q = mMin in solar masses
     * a = alpha
     */
{
    double ans;
    double temp;


    temp = (-4.*Power(a,9.)*Power(p,15.)*Power(q,9.)*
	    Power(x,9.) + 
	    27.*Power(a,6.)*Power(p,20.)*Power(q,14.)*
	    Power(x,12.));
    
    if (temp <0)
      temp = 0 ;
    else
      temp = sqrt(temp);
    
    
    
    ans = q + (Power(0.6666666666666666,0.3333333333333333)*
	       Power(a,3.0))/
      Power(9.0*Power(a,3.0)*Power(p,10.0)*Power(q,7.0)*
	    Power(x,6.0) + Sqrt(3.)* temp ,0.3333333333333333) + 
      Power(9.*Power(a,3.)*Power(p,10.)*Power(q,7.)*
	    Power(x,6.) + Sqrt(3.)* temp ,0.3333333333333333)/
(Power(2.,0.3333333333333333)*
 Power(3.,0.6666666666666666)*Power(p,5.)*Power(q,3.)*
 Power(x,3.)); 

    return ans;
}

double getSqRootArgument (double x, double p, double q, double a)
{
    double ans;

    ans = -4*Power(a,9)*Power(p,15)*Power(q,9)*Power(x,9) + 
   27*Power(a,6)*Power(p,20)*Power(q,14)*Power(x,12);

    return ans;

}
