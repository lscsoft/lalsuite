/* <lalVerbatim file="GenerateInspiralCV">
 * Author: Thomas Cokelaer
 * $Id$
 * </lalVerbatim>  */

/*
<lalLaTeX>
	\subsection{Module \texttt{GenerateInspiral.c}}
	\label{ss:GenerateInspiral.c}
	\noindent Generates a CoherentGW inspiral waveform for injection.

	\subsubsection*{Prototypes}
      	\vspace{0.1in}
      	\input{LALGenerateInspiralCP}
        \input{LALGenerateInspiralGetApproxFromStringCP}
	\input{LALGenerateInspiralGetOrderFromStringCP}
	\input{LALGenerateInspiralGetModelFromStringCP}
	\input{LALGenerateInspiralPopulatePPNCP}
	\input{LALGenerateInspiralPopulateInspiralCP}


     	\idx{LALGenerateInspiral}
        \idx{LALGenerateInspiralGetApproxFromString}
	\idx{LALGenerateInspiralGetOrderFromString}
	\idx{LALGenerateInspiralGetModelFromString}
	\idx{LALGenerateInspiralPopulatePPN}
	\idx{LALGenerateInspiralPopulateInspiral}


	\begin{description}
 	\item[\texttt{LALGenerateInspiral()}] create an inspiral binary 
	waveform generated either by the \texttt{inspiral} package (EOB,
	PadeT1, TaylorT1, TaylorT2, TaylorT3, SpinTaylor) or the 
	\texttt{inject} package	(GeneratePPN).	It is used in the module 
	\texttt{FindChirpSimulation} in \texttt{findchirp} package.
 
	There are three  parsed arguments
 	\begin{itemize}
		 \item a \texttt{CoherentGW}  structure which stores amplitude, 
		 frequency and phase of the  waveform (output)
		 \item a \texttt{thisEvent}  structure which provides some 
		 waveform parameters (input)
		 \item a \texttt{PPNParamStruc} which gives some input 
		 parameters needed by the PPN waveform  generation. That 
		 arguments is also used as an output by all the different
		 approximant  (output/input).
	 \end{itemize}

	The input must be composed of a valid thisEvent structure as well as 
	the  variable deltaT of the  PPNparamsStruct. All others variables 
	of the  PPNParamStruc are populated within that function.

   	\item[\texttt{LALGenerateInspiralGetModelFromString()}] convert a string  
   	provided by the \texttt{CoherentGW} structure in order to retrieve the 
	order and approximant of the waveform to generate. 

   	\item[\texttt{LALGenerateInspiralPopulatePPN()}] Populate the 
	PPNParamsStruc	with the input argument \texttt{thisEvent}. That 
	structure is used by both inspiral waveforms inject waveforms.

  	\item[\texttt{LALGenerateInspiralPopulateInspiral()}]  Populate the 
	InspiralTemplate strucuture if the model chosen belongs to the 
	inspiral package.

	\end{description}
 
   	\subsubsection*{Algorithm}
   	\noindent None.
   
   	\subsubsection*{Notes}
   	Inject only time-domain waveforms for the time being such as PPN, 
   	TaylorT1, TaylorT2, TaylorT3, PadeT1 and EOB , Spintaylor..
   	\subsubsection*{Uses}
   	\begin{verbatim}
   	
 LALInspiralWaveLength()            
 LALInspiralParameterCalc()            
 LALInspiralWave1ForInjection()
 LALInspiralWave2ForInjection()
 LALInspiralWave3ForInjection()
 LALEOBWaveformForInjection()
 LALSTPNWaveformForInjection()
 LALCalloc()
 LALFree()
      
 *  	\end{verbatim}
 *  
 *  	\vfill{\footnotesize\input{GenerateInspiralCV}}
 *  	</lalLaTeX> 
 * */


#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/Units.h>


void 
LALInspiralDummyWaveformForInjection (
				      LALStatus        *status,
				      CoherentGW       *waveform,
				      InspiralTemplate *params,
				      PPNParamStruc    *ppnParams
				      ) ;

NRCSID( GENERATEINSPIRALC, 
"$Id$" );



/* <lalVerbatim file="LALGenerateInspiralCP"> */
void LALGenerateInspiral(
			LALStatus		*status,
			CoherentGW            *waveform,
			SimInspiralTable      *thisEvent,
			PPNParamStruc         *ppnParams
						)
/* </lalVerbatim> */
{
  UINT4			order		= -1; 	/* Order of the model        */
  UINT4			approximant	= -1;	/* And its approximant value */
  /* If we need to use the inspiral package, then we need to fill 
   * the following structure */
  InspiralTemplate      inspiralParams;        										
  CHAR msg[1024]; 
  
  INITSTATUS(status, "LALGenerateInspiral",GENERATEINSPIRALC);
  ATTATCHSTATUSPTR(status);

  /* NO ASSERT FOR THE TIME BEING. WE ASSUME THAT ASSERTS ARE PRESENT 
     IN THE FINDCHIRP FILE OR INSPIRAL PACKAGE. */
  
  /* First, we need to know which approximant and which order to inject.
   * They are provided by the xml injection file and thisEvent structure 
   * as a string such as "EOBtwoPN". However, we need the order  as an 
   * integer to be used in the switch. So,  I wrote two functions to do 
   * the conversion.
   * */
  LALGenerateInspiralGetModelFromString(status->statusPtr,
		  			thisEvent->waveform,
					&order,
					&approximant);
  CHECKSTATUSPTR(status);

  /* Well a priori, we do not need that anymore. Indeed it seems that 
   * it has been decided to put that switch within FindChirpSimulation. 
   * It means that we do enter in that file only for non PPN cases....
   * Any way I keep it for the time being. It was working fine like
   * that though.*/
  /*   --- switch for the model approximant  --- */
  if  (approximant == PPN )
    {
      LALGenerateInspiralPopulatePPN(	status->statusPtr, 
		      			ppnParams, 
					thisEvent);   
      CHECKSTATUSPTR(status);    
      
      LALGeneratePPNInspiral( 	status->statusPtr, 
		      		waveform,
				ppnParams );  
      CHECKSTATUSPTR(status);    
    }
  else
    {
      inspiralParams.approximant = approximant;
      inspiralParams.order       = order;
      /* We fill ppnParams. Redundant with FindChirpSimulation.c but at 
       * least we are sure to fill the structure. Later we can delete
       * the equivalent part within FindChirpSimulation ...  */      
      LALGenerateInspiralPopulatePPN(	status->statusPtr, 
		      			ppnParams, 
					thisEvent);
      CHECKSTATUSPTR(status);      
      /* we fill inspiralParams structure as well.*/
      LALGenerateInspiralPopulateInspiral(	status->statusPtr, 
		      				&inspiralParams, 
						thisEvent, 
						ppnParams);
      CHECKSTATUSPTR(status);      

      /* the waveform generation itself*/
      LALInspiralWaveForInjection(	status->statusPtr,
		      			waveform, 
					&inspiralParams, 
					ppnParams);
      CHECKSTATUSPTR(status);      
      /* If nothing has been computed, then we create a dummy waveform. 
       * Do we need to do that ? 
       * maybe there is a better way to send a null waveform.
       * or just abort ? but it is a pity to abort a job which has
       * run for ages just becuase one waveform can not be injected ? */
      if (waveform->a == NULL) {	
	LALWarning(status, "no waveform generated");
	LALInspiralDummyWaveformForInjection (	status->statusPtr, 
						waveform, 
						&inspiralParams, 
						ppnParams);				 
	CHECKSTATUSPTR(status);      	      
      }      
    } 
/*  fprintf(stderr,"%e %e\n", waveform->phi->data->data[0],
	  waveform->phi->data->data[waveform->phi->data->length-1]);*/

 LALSnprintf( msg, sizeof(msg)/sizeof(*msg),
        "Injected waveform parameters:\n"
        "ppnParams->mTot\t= %e\n"
        "ppnParams->eta\t= %e\n"
        "ppnParams->d\t= %e\n"
        "ppnParams->inc\t= %e\n"
        "ppnParams->phi\t= %e\n"
	"ppnParams->psi\t= %e\n"
        "ppnParams->fStartIn\t= %e\n"
	"ppnParams->fStopIn\t= %e\n"
        "ppnParams->position.longitude\t= %e\n"
        "ppnParams->position.latitude\t= %e\n"
	"ppnParams->position.system\t= %e\n"
	"ppnParams->epoch.gpsSeconds\t= %e\n"       
	"ppnParams->epoch.gpsNanoSeconds\t= %e\n"
        " ppnParams->tC\t= %e\n"
        " ppnParams->dfdt\t =%e\n", 
        ppnParams->mTot, 
        ppnParams->eta, 
        ppnParams->d,
        ppnParams->inc,
        ppnParams->phi,
        ppnParams->psi, 
	ppnParams->fStartIn,
	ppnParams->fStopIn,	
        ppnParams->position.longitude,
        ppnParams->position.latitude,
	ppnParams->position.system,
	ppnParams->epoch.gpsSeconds,	
  	ppnParams->epoch.gpsNanoSeconds,
	ppnParams->tc,
       	ppnParams->dfdt );     
  LALInfo( status, msg );
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}




/* <lalVerbatim file="LALGenerateInspiralGetModelFromStringCP"> */
void LALGenerateInspiralGetModelFromString(LALStatus 	*status,
					   CHAR 	*thisEvent,
					   UINT4 	*order,
					   UINT4 	*model) 
/* </lalVerbatim> */
{
  INITSTATUS(	status, 
		"LALGenerateInspiralGetOrderFromString", 
		GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( thisEvent, status, 
		  GENERATEINSPIRALH_ENULL, 
		  GENERATEINSPIRALH_MSGENULL);
  
  /* Get the order */  
  LALGenerateInspiralGetOrderFromString( status->statusPtr,
					 thisEvent,
					 order);
  CHECKSTATUSPTR(status);      
  /* Get the approximant name */
  LALGenerateInspiralGetApproxFromString( status->statusPtr,
					  thisEvent,
					  model);
  CHECKSTATUSPTR(status);      
  
  DETATCHSTATUSPTR( status );
  RETURN( status );  
}



/* <lalVerbatim file="LALGenerateInspiralGetOrderFromStringCP"> */
void LALGenerateInspiralGetOrderFromString(LALStatus *status,
					   CHAR *thisEvent,
					   UINT4 *order)
/* </lalVerbatim> */
{
  
  /* Function to transform the string of the PN order into an 
   * integer value. */
  CHAR  *ptr 	= NULL;
  
  INITSTATUS( 	status, 
		"LALGenerateInspiralGetOrderFromString", 
		GENERATEINSPIRALC );
  
  ATTATCHSTATUSPTR( status );
  
  ASSERT( thisEvent, status, 
		  GENERATEINSPIRALH_ENULL, 
		  GENERATEINSPIRALH_MSGENULL);
  
  
  if (  (ptr = strstr(thisEvent, 	"newtonian") ) ){
    {*order = newtonian;}}
  else if (  (ptr = strstr(thisEvent, 	"oneHalfPN") ) ){
    {*order = oneHalfPN;}}
  else if (  (ptr = strstr(thisEvent,  	"onePN") ) ){
    {*order = onePN;}}
  else if (  (ptr =  strstr(thisEvent,  "onePointFivePN") ) ){
    {*order = onePointFivePN;}}
  else if (  (ptr = strstr(thisEvent, 	"twoPN") ) ){
    {*order = twoPN;}}
  else if (  (ptr = strstr(thisEvent,  	"twoPointFivePN") ) ){
    {*order = twoPointFivePN;}}
  else if (  (ptr = strstr(thisEvent, 	"threePN") ) ){
    {*order = threePN;}}
  else if (  (ptr = strstr(thisEvent, 	"threePointFivePN") ) ){
    {*order = threePointFivePN;}}
  else {
    ABORT(status, LALINSPIRALH_EORDER, LALINSPIRALH_MSGEORDER);
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="LALGenerateInspiralGetApproxFromStringCP"> */
void LALGenerateInspiralGetApproxFromString(LALStatus *status,
					    CHAR *thisEvent,
					    UINT4 *approximant)  
/* </lalVerbatim> */
{
  /* Function to search for the approximant into a string */
  CHAR  *ptr 	= NULL;
  
  INITSTATUS( 	status, 
		" LALGenerateInspiralGetApproxFromString", 
		GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( 	thisEvent, 
		status, 
		GENERATEINSPIRALH_ENULL, 
		GENERATEINSPIRALH_MSGENULL);
 
  if (  (ptr = strstr(thisEvent, 	"TaylorT1" ) ) ){
    {*approximant = TaylorT1;}}
  else if (  (ptr = strstr(thisEvent, 	"TaylorT2" ) ) ){
    {*approximant = TaylorT2;}}
  else if (  (ptr = strstr(thisEvent, 	"TaylorT3" ) ) ){
    {*approximant = TaylorT3;}}
  else if (  (ptr = strstr(thisEvent, 	"EOB" ) ) ){
    {*approximant = EOB;}}
  else if (  (ptr = strstr(thisEvent, 	"SpinTaylor" ) ) ){
    {*approximant = SpinTaylor;}}
  else if (  (ptr = strstr(thisEvent, 	"PadeT1" ) ) ){
    {*approximant = PadeT1;}}
 else if (  (ptr = strstr(thisEvent, 	"GeneratePPN" ) ) ){
    {*approximant = PPN;}}
  else {
    ABORT(status, LALINSPIRALH_EAPPROXIMANT, LALINSPIRALH_MSGEAPPROXIMANT);
  }
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}





/* That function is just a copy and paste of FindChirpSimulation.c 
 * code related to the injection of PPN event.
 */
/* <lalVerbatim file="LALGenerateInspiralPopulatePPNCP"> */
void  LALGenerateInspiralPopulatePPN(LALStatus             *status,
				     PPNParamStruc         *ppnParams,
				     SimInspiralTable      *thisEvent)
/* </lalVerbatim> */
{

  INITSTATUS( 	status, 
		" LALGenerateInspiralPopulatePPN", 
		GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );

  /* input fields */
  ppnParams->mTot = thisEvent->mass1 + thisEvent->mass2;
  ppnParams->eta  = thisEvent->eta;
  ppnParams->d    = thisEvent->distance* 1.0e6 * LAL_PC_SI; /*in Mpc*/
  ppnParams->inc  = thisEvent->inclination;
  ppnParams->phi  = thisEvent->coa_phase;
  
  /* frequency cutoffs */
  if (thisEvent->fLower >0 ){
    ppnParams->fStartIn = thisEvent->fLower;   
  }
  else {
    ppnParams->fStartIn = 40.;
  }
  ppnParams->fStopIn  = -1.0 /  /* fCutoff of the inspiral package ?? why negative ?*/
    (6.0 * sqrt(6.0) * LAL_PI * ppnParams->mTot * LAL_MTSUN_SI);
  
  /* passed fields */
  ppnParams->position.longitude   = thisEvent->longitude;
  ppnParams->position.latitude    = thisEvent->latitude;
  ppnParams->position.system      = COORDINATESYSTEM_EQUATORIAL;
  ppnParams->psi                  = thisEvent->polarization;
  ppnParams->epoch.gpsSeconds     = 0;  /* Is it correct ? */ 
  ppnParams->epoch.gpsNanoSeconds = 0;  /* Is it correct ? */
  

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}





/* <lalVerbatim file="LALGenerateInspiralPopulateInspiralCP"> */
void LALGenerateInspiralPopulateInspiral(LALStatus		*status,
					 InspiralTemplate      	*inspiralParams,
					 SimInspiralTable      	*thisEvent,
					 PPNParamStruc        	*ppnParams)
     
/* </lalVerbatim> */
{
  INITSTATUS( 	status, 
		" LALGenerateInspiralPopulateInspiral",
		GENERATEINSPIRALC );
  ATTATCHSTATUSPTR( status );

  /* --- Let's fill the inspiral structure now --- */
  inspiralParams->mass1	     	=  thisEvent->mass1;  	/* masses 1 */
  inspiralParams->mass2	     	=  thisEvent->mass2;  	/* masses 2 */
  inspiralParams->fLower	=  ppnParams->fStartIn; /* lower cutoff 
							   frequency */
  inspiralParams->fCutoff	= 1./ (ppnParams->deltaT)/2.-1;  
				/* 	-1 to be  in agreement with the 
					inspiral assert.*/
  inspiralParams->tSampling	  = 1./ (ppnParams->deltaT); /* sampling*/
  inspiralParams->signalAmplitude = 1.; 
  inspiralParams->distance	  =  thisEvent->distance * LAL_PC_SI * 1e6;
  					/*distance in Mpc*/
  inspiralParams->startTime	  =  0.0;
  inspiralParams->startPhase	  =  thisEvent->coa_phase;
  inspiralParams->startPhase      = 0.0;
    
  inspiralParams->OmegaS = GENERATEINSPIRAL_OMEGAS;/* EOB 3PN contribution */
  inspiralParams->Theta	 = GENERATEINSPIRAL_THETA; /* EOB 3PN contribution */
  inspiralParams->Zeta2	 = GENERATEINSPIRAL_ZETA2; /* EOB 3PN contribution */

  inspiralParams->alpha	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi0	 = -1.;      /* bcv useless for the time being */
  inspiralParams->psi3	 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha1 = -1.;      /* bcv useless for the time being */
  inspiralParams->alpha2 = -1.;      /* bcv useless for the time being */
  inspiralParams->beta	 = -1.;      /* bcv useless for the time being */

  inspiralParams->inclination =  thisEvent->inclination ; 
  					/*inclination of the binary*/

  inspiralParams->ieta	    =  1;
  inspiralParams->nStartPad =  0;
  inspiralParams->nEndPad   =  0;
  
  inspiralParams->massChoice  = m1Andm2;      
  

  inspiralParams->sourceTheta = GENERATEINSPIRAL_SOURCETHETA;
  inspiralParams->sourcePhi   = GENERATEINSPIRAL_SOURCEPHI;
  inspiralParams->spin1[0]    = thisEvent->spin1x;
  inspiralParams->spin2[0]    = thisEvent->spin2x;
  inspiralParams->spin1[1]    = thisEvent->spin1y;
  inspiralParams->spin2[1]    = thisEvent->spin2y;
  inspiralParams->spin1[2]    = thisEvent->spin1z;
  inspiralParams->spin2[2]    = thisEvent->spin2z;

  /*put an assert here, theta0 can not be equal to zero.*/
  /*  if (orbitTheta0 == 0){
    ABORT()
    }
  */
  inspiralParams->orbitTheta0 = thisEvent->theta0;
  inspiralParams->orbitPhi0   = thisEvent->phi0;

  DETATCHSTATUSPTR( status );
  RETURN( status );  
}

/*  This is just a dummy function which fill waveform structure even 
 *  though values provided by thisEvent are wrong. I've done that 
 *  function sothat the inspiral code do not stop. Sould add warning 
 *   or something equivalent though.
 */
void 
LALInspiralDummyWaveformForInjection (LALStatus        *status,
				      CoherentGW       *waveform,
				      InspiralTemplate *params,
				      PPNParamStruc    *ppnParams) 
{
  INT4 length = 2;
  CreateVectorSequenceIn in;
  
  INITSTATUS(status, "LALDummyWaveformForInjection", GENERATEINSPIRALC);
  ATTATCHSTATUSPTR(status);

  /** -- Now we can fill the coherent GW strucuture for injection -- */  
  if ( ( waveform->a = (REAL4TimeVectorSeries *)
	 LALCalloc(length,  sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( waveform->f = (REAL4TimeSeries *)
	 LALCalloc(length, sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( waveform->a ); waveform->a = NULL;
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( waveform->phi = (REAL8TimeSeries *)
	 LALCalloc(length, sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( waveform->a ); waveform->a = NULL;
    LALFree( waveform->f ); waveform->f = NULL;
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );

  in.length = length;
  in.vectorLength = 2;

  LALSCreateVectorSequence( status->statusPtr, &( waveform->a->data ), &in );
  CHECKSTATUSPTR(status);      
  LALSCreateVector( status->statusPtr,  &( waveform->f->data ), length);
  CHECKSTATUSPTR(status);      
  LALDCreateVector( status->statusPtr,  &( waveform->phi->data ), length );
  CHECKSTATUSPTR(status);   

  waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
    = 1./params->tSampling;
  
  waveform->a->sampleUnits = lalStrainUnit;
  waveform->f->sampleUnits = lalHertzUnit;
  waveform->phi->sampleUnits = lalDimensionlessUnit;

  LALSnprintf( waveform->a->name,   LALNameLength, "Dummy waveform" );
  LALSnprintf( waveform->f->name,   LALNameLength, "Dummy waveform" ); 
  LALSnprintf( waveform->phi->name, LALNameLength, "Dummy waveform" );

  /* --- fill some output --- */
  ppnParams->tc     	= 0.;
  ppnParams->length 	= 2;
  ppnParams->dfdt   	= 1.; /*waht should we put here ?*/
  ppnParams->fStop  	= params->fLower;
  ppnParams->termCode   = GENERATEPPNINSPIRALH_EFBAD;
  ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFBAD; 
  ppnParams->fStart   	= ppnParams->fStartIn;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
