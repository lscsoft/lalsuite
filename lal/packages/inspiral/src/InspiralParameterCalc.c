#include <lal/LALStdlib.h>
#include <lal/Inspiral.h>



NRCSID (INSPIRALPARAMETERCALCC, "$Id$");


void LALInspiralParameterCalc (LALStatus *status,
			     InspiralParamsOutput *output,
			     InspiralParamsInput *params)
{

  REAL8 m1, m2, totalMass, eta, mu, fLower, piFl;
  REAL8 dumm_const1, dumm_const2, dumm_const3;


  INITSTATUS (status, "LALInspiralParameterCalc", INSPIRALPARAMETERCALCC );

  ASSERT(output, status, INSPIRALPARAMETERCALCALC_ENULL,INSPIRALPARAMETERCALCALC_MSGENULL);
  ASSERT(params, status, INSPIRALPARAMETERCALCALC_ENULL,INSPIRALPARAMETERCALCALC_MSGENULL);


   switch(params->massChoice) {

 	case m1Andm2:

  		output->m1 = m1 = params->m1;
  		output->m2 = m2 = params->m2;
  		output->totalMass = totalMass = m1+m2;
  		output->eta = eta = m1*m2/pow(totalMass,2);
  		output->mu = mu = m1*m2/totalMass;
  		output->chirpMass = pow(mu,0.6)*pow(totalMass,0.4);

  		break;

  	case totalMassAndEta:

  		output->totalMass = totalMass = params->totalMass;
  		output->eta = eta = params->eta;
  		output->m1 = 0.5*totalMass + (totalMass*sqrt(0.25-eta));
  		output->m2 = 0.5*totalMass - (totalMass*sqrt(0.25-eta));
  		output->mu = eta*totalMass;
  		output->chirpMass = pow(eta,0.6)*totalMass;

  		break;

  	case totalMassAndMu:

  		output->totalMass = totalMass = params->totalMass;
  		output->mu = params->mu;
  		output->eta = eta =  (output->mu)/totalMass;
  		output->m1 = 0.5*totalMass + (totalMass*sqrt(0.25-eta));
  		output->m2 = 0.5*totalMass - (totalMass*sqrt(0.25-eta));
  		output->chirpMass = pow(eta,0.6)*totalMass;

  	default:
                fprintf(stderr, "No choice of masses in LALInspiralParameterCalc ... exiting\n");
                exit(0);
   }
   totalMass = totalMass*LAL_MTSUN_SI;
   fLower = params->fLower;
   piFl = LAL_PI*fLower;
   output->tau0 = (5.0*pow(eta,-1.0)*pow(totalMass,(-5.0/3.0))
                  *pow(piFl,(-8.0/3.0)))/256.0;
   output->tau2 = (3715.0 + (4620.0*eta))/(64512.0*eta*totalMass*pow(piFl,2.0));
   output->tau3 = LAL_PI/(8.0*eta*pow(totalMass,(2.0/3.0))*pow(piFl,(5.0/3.0)));
   dumm_const1 = 3058673.0/1016064.0;
   dumm_const2 = (5429.0*eta)/1008.0;
   dumm_const3 = (617.0*pow(eta,2.0))/144.0;
   output->tau4 = (5.0/(128.0*eta*pow(totalMass,(1.0/3.0))*pow(piFl,(4.0/3.0))))
   *(dumm_const1 + dumm_const2 + dumm_const3);
   output->tauC = output->tau0 + output->tau2 - output->tau3 + output->tau4;

   switch(params->order) {
			
           case newtonian:
	   case oneHalfPN:
	        fprintf(stderr,"Choosing Newtonian waveforms \n");
                output->tau2=0.0;
         	output->tau3=0.0;
         	output->tau4=0.0;
         	output->tauC = output->tau0;

         	break;

	   case onePN:
	        fprintf(stderr,"Choosing 1 PN waveforms \n");
		output->tau3=0.0;
         	output->tau4=0.0;
         	output->tauC = output->tau0 + output->tau2;

         	break;

           case onePointFivePN:
	        fprintf(stderr,"Choosing 1.5 PN waveforms \n");
         	output->tau4=0.0;
         	output->tauC = output->tau0 + output->tau2 - output->tau3;

		break;

	   case twoPN:
	        fprintf(stderr,"Choosing 2 PN waveforms \n");

		break;

	   default:
		fprintf(stderr,"You haven't chosen a PN order in LALInspiralParameterCalc\n");
   }
   RETURN (status);

}
