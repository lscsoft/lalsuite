/* choosemodel.c
   February 21, 00.
   To set pointers to functions dEnergy and flux to the
   desired energy and flux functions.
*/

#include <math.h>
#include "Inspiral.h"
#include "LALStdlib.h"

NRCSID (CHOOSEMODELC, "$Id$");





static REAL8 et0(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, x;
   x = v*v;
   energy = ak->eTaN * x;
   return (energy);
}

static REAL8 et2(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, x;
   x = v*v;
   energy = ak->eTaN * x * (1. + ak->eTa1*x);
   return (energy);
}

static REAL8 et4(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, x;
   x = v*v;
   energy = ak->eTaN * x * (1. + ak->eTa1*x + ak->eTa2*x*x);
   return (energy);
}

static REAL8 Et0(REAL8 v, expnCoeffs *ak) {
   REAL8 Energy, x;
   x = v*v;
   Energy = ak->ETaN * x;
   return (Energy);
}

static REAL8 Et2(REAL8 v, expnCoeffs *ak) {
   REAL8 Energy, x;
   x = v*v;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x);
   return (Energy);
}

static REAL8 Et4(REAL8 v, expnCoeffs *ak) {
   REAL8 Energy, x;
   x = v*v;
   Energy = ak->ETaN * x * (1. + ak->ETa1*x + ak->ETa2*x*x);
   return (Energy);
}

static REAL8 dEt0(REAL8 v, expnCoeffs *ak) {
   REAL8 dEnergy;
   dEnergy = ak->dETaN * v;
   return (dEnergy);
}

static REAL8 dEt2(REAL8 v, expnCoeffs *ak) {
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x);
   return (dEnergy);
}

static REAL8 dEt4(REAL8 v, expnCoeffs *ak) {
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x);
   return (dEnergy);
}





static REAL8 Ft0(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->FTaN * pow(v,10.);
   return (flux);
}

static REAL8 Ft2(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->FTaN * pow(v,10.) * (1. + ak->FTa2*pow(v,2.));
   return (flux);
}

static REAL8 Ft3(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->FTaN * pow(v,10.) * (1. + ak->FTa2*pow(v,2.) + ak->FTa3*pow(v,3.));
   return (flux);
}

static REAL8 Ft4(REAL8 v, expnCoeffs *ak) {
   REAL8 flux; 
   flux = ak->FTaN * pow(v,10.) * (1. + ak->FTa2*pow(v,2.) + ak->FTa3*pow(v,3.) + ak->FTa4*pow(v,4.));
   return (flux);
}

static REAL8 Ft5(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->FTaN * pow(v,10.) * (1.+ ak->FTa2*pow(v,2.) + ak->FTa3*pow(v,3.) + ak->FTa4*pow(v,4.) + ak->FTa5*pow(v,5.));
   return (flux);
}














static REAL8 ep0(REAL8 v, expnCoeffs *ak) {
   REAL8 x, energy;
   x = v*v;
   energy = -x;
   return (energy);
}

static REAL8 ep2(REAL8 v, expnCoeffs *ak) {
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1 * x);
   return (energy);
}

static REAL8 ep4(REAL8 v, expnCoeffs *ak) {
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x));
   return (energy);
}

static REAL8 dEp0(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, denergy, Energy, dEnergy, y;
   energy = ep0(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   denergy = -1;
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}

static REAL8 dEp2(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep2(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   denergy = -1 / pow(1.+ak->ePa1 * x, 2.);
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}

static REAL8 dEp4(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep4(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   denergy = (1. + 2.*ak->ePa2*x + (ak->ePa1 + ak->ePa2) * ak->ePa2 * pow(x,2.))/pow(1. + (ak->ePa1 + ak->ePa2) * x ,2.);
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


/*
static REAL8 Ep0(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, Energy, y;
   energy = ep0(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   return(Energy);
}

static REAL8 Ep2(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, Energy, x, y;
   x = v*v;
   energy = ep2(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   return(Energy);
}

static REAL8 Ep4(REAL8 v, expnCoeffs *ak) {
   REAL8 energy, Energy, x, y;
   x = v*v;
   energy = ep4(v, ak);
   y = pow(1.+energy, 0.5);
   Energy = pow (1. + 2.* ak->eta * (y - 1.), 0.5) - 1.;
   return(Energy);
}

*/













static REAL8 Fp0(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.);
   return (flux);
}

static REAL8 Fp1(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.)/ ((1.+ak->fPa1*v) * (1.-v/ak->vpoleP4));
   return (flux);
}

static REAL8 Fp2(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.)/ ((1.+ak->fPa1*v / (1.+ak->fPa2*v)) * (1.-v/ak->vpoleP4));
   return (flux);
}

static REAL8 Fp3(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.)/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v))) * (1.-v/ak->vpoleP4));
   return (flux);
}

static REAL8 Fp4(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.)/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v / (1.+ak->fPa4*v)))) * (1.-v/ak->vpoleP4));
   return (flux);
}

static REAL8 Fp5(REAL8 v, expnCoeffs *ak) {
   REAL8 flux;
   flux = ak->fPaN * pow(v,10.)/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v / (1.+ak->fPa4*v / (1.+ak->fPa5*v))))) * (1.-v/ak->vpoleP4));
   return (flux);
}






void ChooseModel(Status *status,
		 expnFunc *f,
		 expnCoeffs *ak,
		 InspiralTemplate *params)
{
   REAL8 vn, vlso;
   TofVIn in1;
   REAL8 tofv;
   void *in2;

   INITSTATUS (status, "ChooseModel", CHOOSEMODELC);
   ATTATCHSTATUSPTR(status);

   ASSERT (f,  status, CHOOSEMODEL_ENULL, CHOOSEMODEL_MSGENULL);
   ASSERT (ak,  status, CHOOSEMODEL_ENULL, CHOOSEMODEL_MSGENULL);
   ASSERT (params,  status, CHOOSEMODEL_ENULL, CHOOSEMODEL_MSGENULL);


   switch (params->approximant) {

	case taylor:

		switch (params->order) {
			case newtonian:
	         		ak->vn = ak->vlso = vlso = ak->vlsoT0;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
				ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEt0;
	   			f->flux = Ft0;
	   			fprintf(stderr, "Choosing Newtonian Taylor waveforms\n");
	   			break;
			case oneHalfPN:
	   			fprintf(stderr, "You can't choose 1/2 PN T-approximants\n");
	   			break;
			case onePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoT2;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEt2;
	   			f->flux = Ft2;
	   			fprintf(stderr, "Choosing 1 PN T-approximants\n");
	   			break;
			case onePointFivePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoT2;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEt2;
	   			f->flux = Ft3;
	   			fprintf(stderr, "Choosing 1.5 PN T-approximants \n");
	   			break;
			case twoPN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoT2;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEt4;
	   			f->flux = Ft4;
	   			fprintf(stderr, "Choosing 2 PN T-approximants \n");
	   			break;
			case twoPointFivePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoT2;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEt4;
	   			f->flux = Ft5;
	   			fprintf(stderr, "Choosing 2.5 PN T-approximants \n");
  				}

	break;

	case pade:

		switch (params->order) {
			case newtonian:			
	   			ak->vn = ak->vlso = vlso = ak->vlsoP0;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp0;
	   			f->flux = Fp0;
	   			fprintf(stderr, "Choosing Newtonian pade waveforms\n");
	   			break;
			case oneHalfPN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoP0;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp0;
	   			f->flux = Fp1;
	   			fprintf(stderr, "Choosing 1/2 PN P-approximants\n");
	   			break;
			case onePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoP0;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp2;
	   			f->flux = Fp2;
	   			fprintf(stderr, "Choosing 1 PN P-approximants\n");
	   			break;
			case onePointFivePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoP0;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp2;
	   			f->flux = Fp3;
	   			fprintf(stderr, "Choosing 1.5 PN P-approximants \n");
	   			break;
			case twoPN:
	  			ak->vn = ak->vlso = vlso = ak->vlsoP4;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp4;
	   			f->flux = Fp4;
	   			fprintf(stderr, "Choosing 2 PN P-approximants \n");
	   			break;
			case twoPointFivePN:
	   			ak->vn = ak->vlso = vlso = ak->vlsoP4;
	   			if (ak->fn) {
	      			vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
	      			ak->vn = (vn < vlso) ? vn :  vlso;
	   			} 
	   			f->dEnergy = dEp4;
	   			f->flux = Fp5;
	   			fprintf(stderr, "Choosing 2.5 PN P-approximants \n");
      				}
			}

      in1.t=0.0;
      in1.v0=ak->v0;
      in1.t0=ak->t0;
      in1.totalmass = ak->totalmass;
      in1.dEnergy = f->dEnergy;
      in1.flux = f->flux;
      in1.coeffs = ak;

      in2 = (void *) &in1;      

      TofV(status->statusPtr, &tofv, ak->vn, in2);
      CHECKSTATUSPTR(status);


      ak->tn = -tofv - ak->samplinginterval;
      ak->fn = pow(ak->vn, 3.)/(LAL_PI * ak->totalmass);

/*      fprintf(stderr, "%e %e %e %e %e %e\n", ak->tn, ak->fn, ak->vn, ak->t0, ak->f0, ak->v0);
*/


   DETATCHSTATUSPTR(status);
   RETURN (status); 

}

