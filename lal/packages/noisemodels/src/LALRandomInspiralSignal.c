/*  <lalVerbatim file="LALRandomInspiralSignalCV">
Author: Sathyaprakash, B. S., Thomas Cokelaer
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALRandomInspiralSignal.c}}
Module to generate 

(a) inspiral signals with random masses or chirp times 
that have values within the parameter space specified by an input struct,

(b) simulated Gaussian noise of PSD expected in a given interferometer 

(c) inspiral signal as in (a) but of a specified amplitude added 
to simulated Gaussian noise as in (b).

In all cases the returned vector is the Fourier transform of the relevant signal.  

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALRandomInspiralSignalCP}
\idx{LALRandomInspiralSignal()}

\subsubsection*{Description}
The function receives input struct of type RandomInspiralSignalIn 
whose members are
\begin{verbatim}
typedef struct 
tagRandomInspiralSignalIn
{
   INT4 useed;
   INT4 type;

   REAL8 mMin;
   REAL8 mMax;
   REAL8 MMax;
   REAL8 SignalAmp;
   REAL8 NoiseAmp;
   REAL8 etaMin;
   REAL8 t0Min;
   REAL8 t0Max;
   REAL8 tnMin;
   REAL8 tnMax;

   InspiralTemplate param;
   REAL8Vector psd;                 
   RealFFTPlan *fwdp;
} RandomInspiralSignalIn;
\end{verbatim}

Depending on the value of the parameter (\texttt{randIn.type}) this
code returns the Fourier transform of 

(a) a pure inspiral signal of a given type
(\texttt{randIn.type=0}),

(b) simulated noise expected
in a chosen interferometer \texttt{randIn.type=1} or 

(c) 
$\mathtt{SignalAmp}\times s+\mathtt{NoiseAmp}\times n$ (\texttt{randIn.type=2}), 
where $s$ is normalised signal and $n$ random Gaussian noise whose PSD is
that expected in a given interferometer with zero mean and unit rms.

User must specify the following quantities in the input structure
\begin{table}
\begin{tabular}{lcl}
\hline
\hline
     Parameter         &   i/o  &   Comment \\
\hline
\texttt {INT4 useed}        &  input &   Seed for the random number generator    \\
\texttt {INT4 type}         &  input &   Type of signal required to be generated    \\\\
\texttt {InspiralTemplate p}  &  i/o   &   user must input certain params; others will be output\\
\texttt {p.startTime}         &        &   usually 0.\\
\texttt {p.startPhase}        &        &   $[0,\pi/2]$\\
\texttt {p.nStartPad}         &        &   number of zeros in the vector before the signal begins\\
\texttt {p.nEndPad}           &        &   number of zeros in the vector after the signal ends\\
\texttt {p.signalAmplitude}   &        &   usually 1\\
\texttt {p.ieta}              &        &   1 for comparable mass systems 0 for test mass model\\
\texttt {p.fLower}            &        &   lower frequency cutoff in Hz\\
\texttt {p.fCutoff}           &        &   upper frequency cutoff in Hz\\
\texttt {p.tSampling}         &        &   sampling rate in Hz\\
\texttt {p.order}             &        &   order of the PN approximant of the signal \\
\texttt {p.approximant}       &        &   PN approximation to be used for inspiral signal generation\\
\texttt {InputMasses massChoice}        & input &   space in which parameters are chosen; \texttt{m1Andm2, totalMassAndEta, totalMassUAndEta, t02, t03, bhns}\\\\
\texttt {REAL8Vector psd}   &  input &   pre-computed power spectral density used for coloring the noise \\
\texttt {RealFFTPlan *fwdp} &  input &   pre-computed fftw plan to compute forward Fourier transform   \\\\
\texttt {REAL8 mMin}        &  input &   smallest component mass allowed   \\
\texttt {REAL8 mMax}        &  input &   largest component mass allowed   {\bf OR} \\
\texttt {REAL8 MMax}        &  input &   largest total mass allowed   \\
\texttt {REAL8 SignalAmp}   &  input &   amplitude of the signal (relevant only when \texttt{type=2})   \\
\texttt {REAL8 NoiseAmp}    &  input &   amplitude of noise (relevant only when \texttt{type=2})   \\
\texttt {REAL8 etaMin}      &  input &   smallest value of the symmetric mass ratio    \\
\hline
\multicolumn{3}{c}{Following chirp times are needed 
only if \texttt{param.massChoice {\rm is} t02 {\rm or} t03}} \\
\hline
\texttt {REAL8 t0Min}       &  input &   smallest Newtonian chirp time   \\
\texttt {REAL8 t0Max}       &  input &   largest Newtonian chirp time   \\
\texttt {REAL8 tnMin}       &  input &   smallest 1 chirp time if \texttt{param.massChoice=t02}\\
                            &        &   smallest 1.5 chirp time if \texttt{param.massChoice=t03}\\
\texttt {REAL8 tnMax}       &  input &   largest 1 chirp time  if \texttt{param.massChoice=t02}\\
                            &        &   largest 1.5 chirp time  if \texttt{param.massChoice=t03}\\
\hline
\end{tabular}
\caption{Input structure needed for the function \texttt{LALRandomInspiralSignal}}.
\end{table}
When repeatedly called, the parameters of the signal will be 
uniformly distributed in the space of 

(a) component masses in the range \texttt{[randIn.mMin, randIn.mMax]} if  
\texttt{param.massChoice=m1Andm2},
  
(b) component masses greater than \texttt{randIn.mMin} and total mass 
less than \texttt{randIn.MMax} if  \texttt{param.massChoice=totalMassAndEta},
  
(c) component masses greater than \texttt{randIn.mMin} and \texttt{uniform} total mass 
less than \texttt{randIn.MMax} if  \texttt{param.massChoice=totalMassUAndEta},
  
(d) Newtonian and first post-Newtonian chirp times if 
\texttt{param.massChoice=t02},

(e) Newtonian and 1.5 post-Newtonian chirp times if 
\texttt{param.massChoice=t03} and.

(f) component masses in the range \texttt{[randIn.mMin, randIn.mMax]} one of them being a neutron start and the other a black hole (one above 3 solar mass and one below) if  \texttt{param.massChoice=bhns}. The function therefore checks the mass range validity i.e. randIn.mMin must be less than 3 and randIn.mMax greater than 3.

\subsubsection*{Algorithm}
No special algorithm, only a series of calls to pre-existing functions.
\subsubsection*{Uses}
\begin{verbatim}
random
LALInspiralParameterCalc
LALInspiralWave
LALREAL4VectorFFT
LALInspiralWaveNormaliseLSO
LALCreateRandomParams
LALNormalDeviates
LALDestroyRandomParams
LALREAL4VectorFFT
LALColoredNoise
LALAddVectors
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRandomInspiralSignalCV}}
%Laldoc Closed at: Wed Jan 16 08:39:35 2002

</lalLaTeX>  */
#include <stdlib.h>
#include <lal/LALNoiseModels.h>
#include <lal/Random.h>

#define random() rand()
#define srandom( seed ) srand( seed )

NRCSID (LALRANDOMINSPIRALSIGNALC, "$Id$");

/*  <lalVerbatim file="LALRandomInspiralSignalCP"> */

void
LALRandomInspiralSignal
   (
   LALStatus              *status, 
   REAL4Vector            *signal,
   RandomInspiralSignalIn *randIn
   )
{  /*  </lalVerbatim>  */

  REAL8 maxTemp; /* temporary variable*/
  INT4 iMax;     /* temporary indice */
  UINT4 indice;
   REAL8 epsilon1, epsilon2, norm;
   REAL4Vector noisy, buff;
   AddVectorsIn addIn;
   INT4 valid;
   static RandomParams *randomparams;
   InspiralWaveNormaliseIn normin;
   
   INITSTATUS (status, "LALRandomInspiralSignal", LALRANDOMINSPIRALSIGNALC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (randIn->psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (randIn->mMin > 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->MMax > 2*randIn->mMin, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->type >= 0, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (randIn->type <= 2, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   buff.length = signal->length;
   if (!(buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length))) {
      ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
   }
   /* Use the seed to initialize random(). */
   srandom(randIn->useed);
   /* use the random number so generated as the next seed */
   randIn->useed = random();
   /* we need random parameters only if we need to generate a signal (i.e. type 0/2) */
   if (randIn->type==0 || randIn->type==2)
     {
       valid = 0;
       /* Keep generating random parameters until they 
	* are located within the specified region */
       while (!valid) 
	 {
	   epsilon1 = (float) random()/(float)RAND_MAX;
	   epsilon2 = (float) random()/(float)RAND_MAX;
	   switch (randIn->param.massChoice) 
	     {
	       
	     case bhns:
	       ASSERT(randIn->mMin<=3 && randIn->mMax>=3, status, 10,
		      "if bhns requested, mass1 must be <= 3 solar mass and mass2 >= 3 solar mass\n");
	       randIn->param.mass1 = randIn->mMin + 
		 (3 - randIn->mMin) * epsilon1;
	       randIn->param.mass2 = 3  + 
		 (randIn->mMax - 3) * epsilon2;
	       randIn->param.massChoice=m1Andm2;
	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
	       break;

	     case m1Andm2: 
	       /* 
		* restriction is on the minimum and maximum individual 
		* masses of the two component stars. 
		*/
	       randIn->param.mass1 = randIn->mMin + 
		 (randIn->mMax - randIn->mMin) * epsilon1;
	       randIn->param.mass2 = randIn->mMin  + 
		 (randIn->mMax - randIn->mMin) * epsilon2;
	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));	       
               CHECKSTATUSPTR(status);
               break;
               
	     case totalMassAndEta: 
	       /*
		* restriction is on the total mass of the binary 
		* and the minimum mass of the component stars
		*/
	       randIn->param.mass1 = randIn->mMin 
		 + (randIn->MMax - 2.*randIn->mMin) * epsilon1;
	       randIn->param.mass2 = randIn->mMin 
		 + (randIn->MMax - randIn->param.mass1 - randIn->mMin) * epsilon2;
	       randIn->param.totalMass = randIn->param.mass1 + randIn->param.mass2 ;
	       randIn->param.eta = (randIn->param.mass1*randIn->param.mass2)/pow(randIn->param.totalMass,2.L); 
	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
	       break;
	     case totalMassUAndEta: 
	       /*
		* restriction is on the total mass of the binary 
		* and the etaMin which depends on max total mass.
		*/
               {
                 REAL4 etaMin;
                 
                 randIn->param.totalMass =  2*randIn->mMin  +  epsilon1 * (randIn->MMax - 2 * randIn->mMin) ;
                
		 if (randIn->param.totalMass < (randIn->mMin+ randIn->mMax)) {
                   etaMin = (randIn->mMin / randIn->param.totalMass);		 
                   etaMin = etaMin - etaMin *etaMin;
                 }
		 else {
		   etaMin = (randIn->mMax / randIn->param.totalMass);		 
		   etaMin = etaMin - etaMin *etaMin;
		 }
                 randIn->param.eta = etaMin + epsilon2 * (.25 - etaMin);      
                 
                 LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
                 CHECKSTATUSPTR(status);
	       }               
	       break;  

	     case fixedMasses: /* the user has already given individual masses*/
 	       randIn->param.massChoice = m1Andm2;
 	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
 	       break;
 	     case fixedPsi: /* the user has already given psi0/psi3*/
 	       randIn->param.massChoice = psi0Andpsi3;
 	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
 	       break;
 	     case fixedTau: /* the user has already given tau0/tau3*/
 	       randIn->param.massChoice = t03;
 	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
 	       break;
      	     case t02: 
	       /* chirptimes t0 and t2 are required in a specified range */
	       randIn->param.t0 = randIn->t0Min + (randIn->t0Max - randIn->t0Min)*epsilon1;
	       randIn->param.t2 = randIn->tnMin + (randIn->tnMax - randIn->tnMin)*epsilon2;
	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
	       break;
	     case t03: 
	       /* chirptimes t0 and t3 are required in a specified range */
	       randIn->param.t0 = randIn->t0Min + (randIn->t0Max - randIn->t0Min)*epsilon1;
	       randIn->param.t3 = randIn->tnMin + (randIn->tnMax - randIn->tnMin)*epsilon2;
	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
	       break;
	     case psi0Andpsi3: 
	       /* BCV parameters are required in a specified range */
	       randIn->param.psi0 = randIn->psi0Min + (randIn->psi0Max - randIn->psi0Min)*epsilon1;
	       randIn->param.psi3 = randIn->psi3Min + (randIn->psi3Max - randIn->psi3Min)*epsilon2;
	       break;
	     case t04: 
	     default:
	       /* if the choice of parameters is wrong abort the run */
	       ABORT (status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
	       break;
	     }
	   /*
	     printf("%e %e %e %e\n", randIn->param.t0, randIn->param.t3, randIn->param.mass1, randIn->param.mass2);
	   */
	   switch (randIn->param.massChoice) 
	     {
	     case bhns:
	       valid=1;
	       break;
	     case m1Andm2: 
	     case t03: 
	     case t02:
	       /*
		* The following imposes a range in which min and
		* max of component masses are restricted. 
		*/
	       if (
		   randIn->param.mass1 >= randIn->mMin &&
		   randIn->param.mass2 >= randIn->mMin &&
		   randIn->param.mass1 <= randIn->mMax &&
		   randIn->param.mass2 <= randIn->mMax &&
		   randIn->param.eta <= 0.25 &&
		   randIn->param.eta >= randIn->etaMin
		   )
		 {
		   valid = 1;
		 }
	       break;
	     case totalMassAndEta: 
	     case totalMassUAndEta: 
	       
	       /*
		* The following imposes a range in which min of 
		* component masses and max total mass are restricted. 
		*/
	       if (
		   randIn->param.mass1 >= randIn->mMin &&
		   randIn->param.mass2 >= randIn->mMin &&
		   randIn->param.totalMass <= randIn->MMax &&
		   randIn->param.eta <= 0.25 &&
		   randIn->param.eta >= randIn->etaMin &&
		   randIn->param.mass1 <= randIn->mMax && 
		   randIn->param.mass2 <= randIn->mMax 
		   )
		 
		 {
		   valid = 1;
		 }
	       break;
	       
	     case psi0Andpsi3: 
	       /* 
		* the following makes sure that the BCV has
		* a well defined end-frequency
		*/

	       LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
               CHECKSTATUSPTR(status);
	       valid = 1; 
	       /*	       if (randIn->param.totalMass > 0.)
		 {
		   REAL8 fLR, fLSO, fend;
		   epsilon1 = (float) random()/(float)RAND_MAX;
		   fLR = 1.L/(LAL_PI * pow (3.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI);
		   fLSO = 1.L/(LAL_PI * pow (6.L,1.5) * randIn->param.totalMass * LAL_MTSUN_SI);
		   fend = fLSO + (fLR - fLSO) * epsilon1;
		   if (fend > randIn->param.tSampling/2. || fend < randIn->param.fLower) break;
		   randIn->param.fFinal = fend;
		   valid = 1;
		   }*/
	       break;
	     case t04: 
	     default:
	       ABORT (status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
	       break; 
	     }
	 }
     }

   
   /* set up the structure for normalising the signal */
   normin.psd = &(randIn->psd);
   normin.df = randIn->param.tSampling / (REAL8) signal->length;
   normin.fCutoff = randIn->param.fCutoff;
   normin.samplingRate = randIn->param.tSampling;

   switch (randIn->type) 
     {
     case 0:
       
	      /* First deal with the signal only case:
	       * if the signal is generated in the Fourier domain no
	       * need for Fourier transform
	       */
       if (randIn->param.approximant == BCV ||
	   randIn->param.approximant == TaylorF1 ||
	   randIn->param.approximant == TaylorF2 ||
	   randIn->param.approximant == PadeF1)
	 { 
	   LALInspiralWave(status->statusPtr, signal, &randIn->param);
	   CHECKSTATUSPTR(status);
	 }
       else
	 {
	   /* force to comput fFinal is it really necessary  ? */
	   randIn->param.fFinal=0;
	   LALInspiralWave(status->statusPtr, &buff, &randIn->param);
	   CHECKSTATUSPTR(status);
	   LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
	   CHECKSTATUSPTR(status);
	 }
	

       /* we might want to know where is the signal injected*/
       maxTemp  = 0;
       iMax     = 0;
       for ( indice = 0 ; indice< signal->length; indice++)
	 {
	   if (fabs(signal->data[indice]) > maxTemp){
	     iMax = indice;
	     maxTemp = fabs(signal->data[indice]);
	   }
	 }
       randIn->coalescenceTime = iMax;

       normin.fCutoff = randIn->param.fFinal;

       LALInspiralWaveNormaliseLSO(status->statusPtr, signal, &norm, &normin);
       CHECKSTATUSPTR(status);
       break;
     case 1:
       /* 
	* next deal with the noise only case:
	*/
       /*
	 Old method of generating Gaussian noise
	 LALGaussianNoise(status->statusPtr, &buff, &randIn->useed);
        */
       LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
       CHECKSTATUSPTR(status);
       LALNormalDeviates(status->statusPtr, &buff, randomparams);
       CHECKSTATUSPTR(status);
       LALDestroyRandomParams(status->statusPtr, &randomparams);
       CHECKSTATUSPTR(status);
       LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
       CHECKSTATUSPTR(status);
       LALColoredNoise(status->statusPtr, signal, randIn->psd);
       CHECKSTATUSPTR(status);
	      
       /* multiply the noise vector by the correct normalisation factor */
       {
	       double a2 = randIn->NoiseAmp * sqrt (randIn->param.tSampling)/2.L; 
	       UINT4 i;
	       for (i=0; i<signal->length; i++) signal->data[i] *= a2;
       }
       break;
     default:
       /* 
	* finally deal with the noise+signal only case:
	*/
         noisy.length = signal->length;
         if (!(noisy.data = (REAL4*) LALMalloc(sizeof(REAL4)*noisy.length))) 
	   {
	     if (buff.data != NULL) LALFree(buff.data);
	     buff.data = NULL;
	     ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
	   }
         LALCreateRandomParams(status->statusPtr, &randomparams, randIn->useed);
         CHECKSTATUSPTR(status);
         LALNormalDeviates(status->statusPtr, &buff, randomparams);
         CHECKSTATUSPTR(status);
         LALDestroyRandomParams(status->statusPtr, &randomparams);
         CHECKSTATUSPTR(status);
         LALREAL4VectorFFT(status->statusPtr, &noisy, &buff, randIn->fwdp);
         CHECKSTATUSPTR(status);
         LALColoredNoise(status->statusPtr, &noisy, randIn->psd);
         CHECKSTATUSPTR(status);

	 if (randIn->param.approximant == BCV ||
	     randIn->param.approximant == TaylorF1 ||
	     randIn->param.approximant == TaylorF2 ||
	     randIn->param.approximant == PadeF1)
	 {
		 LALInspiralWave(status->statusPtr, &buff, &randIn->param);
		 CHECKSTATUSPTR(status);
	 }
	 else
	 {
		 LALInspiralWave(status->statusPtr, signal, &randIn->param);
		 CHECKSTATUSPTR(status);
		 LALREAL4VectorFFT(status->statusPtr, &buff, signal, randIn->fwdp);
		 CHECKSTATUSPTR(status);
	 }
	 maxTemp  = 0; 
	 iMax     = 0;
	 for ( indice = 0 ; indice< signal->length; indice++)
	   {
	     if (fabs(signal->data[indice]) > maxTemp){
	       iMax = indice;
	       maxTemp = fabs(signal->data[indice]);
	     }
	   }
	 randIn->coalescenceTime = iMax;

         normin.fCutoff = randIn->param.fFinal;
	 LALInspiralWaveNormaliseLSO(status->statusPtr, &buff, &norm, &normin);
         CHECKSTATUSPTR(status);

         addIn.v1 = &buff;
         addIn.a1 = randIn->SignalAmp;
         addIn.v2 = &noisy;
	 /* After experimenting we found that the following factor SamplingRate*sqrt(2)
	  * is needed for noise amplitude to get the output of the signalless correlation
	  * equal to unity; the proof of this is still needed 
	  */
         addIn.a2 = randIn->NoiseAmp * sqrt (randIn->param.tSampling)/2.L; 
         LALAddVectors(status->statusPtr, signal, addIn);
         CHECKSTATUSPTR(status);
         if (noisy.data != NULL) LALFree(noisy.data);
         break;
   }
   if (buff.data != NULL) LALFree(buff.data);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

