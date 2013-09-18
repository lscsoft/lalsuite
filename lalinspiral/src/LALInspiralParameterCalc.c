/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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

/**
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Given a pair of masses (or other equivalent parameters) compute
 * related chirp parameters.
 *
 * ### Prototypes ###
 *
 * <tt>XLALInspiralParameterCalc()</tt>
 * <ul>
 * <li>\c params: Input/Output, given a pair of binary parameters and a lower
 * frequency cutoff, other equivalent parameters are computed by this function.</li>
 * </ul>
 *
 * ### Description ###
 *
 * The code takes as its input <tt>params->fLower</tt> in Hz and
 * a pair of masses (in units of \f$M_\odot\f$) or chirptimes (in seconds measured from <tt>params->fLower</tt>)
 * and computes all the other {\em mass} parameters in the \c params structure.
 * Users choice of input pair of {\em masses} should be specified by appropriately setting
 * the variable <tt>params->massChoice</tt> as described in the Table below:
 *
 * <table align="center" class="doxtable">
 * <caption align="top" style="text-align: left; font-weight: normal;">
 * Table I. For a given <tt>params->massChoice</tt> in column 1 the user should specify the
 * parameters as in column 2, in units as in column 3. Column 4 gives the conventional meaning
 * of the parameters. Chirp times are measured from a lower frequency cutoff given
 * in <tt>params->fLower.</tt></caption>
 * <tr><th><tt>params->massChoice</tt></th><th>User should set</th><th>in units</th><th>which means</th></tr>
 * <tr><td>m1Andm2</td><td>(<tt>mass1, mass2</tt>)</td><td>\f$(M_\odot, M_\odot)\f$</td><td>\f$(m_1,m_2)\f$</td></tr>
 * <tr><td>totalMassAndEta</td><td>(<tt>totalmass, eta</tt>)</td><td>\f$(M_\odot, 0 < \eta \le 1/4)\f$</td><td>\f$(m, \eta)\f$</td></tr>
 * <tr><td>totalMassAndMu</td><td>(<tt>totalmass, mu</tt>)</td><td>\f$(M_\odot, M_\odot)\f$</td><td>\f$(m, \mu)\f$</td></tr>
 * <tr><td>t02</td><td>(<tt>t0, t2</tt>)</td><td>(sec, sec)</td><td>\f$(\tau_0, \tau_2)\f$</td></tr>
 * <tr><td>t03</td><td>(<tt>t0, t3</tt>)</td><td>(sec, sec)</td><td>\f$(\tau_0, \tau_3)\f$</td></tr>
 * <tr><td>t04</td><td>(<tt>t0, t4</tt>)</td><td>(sec, sec)</td><td>\f$(\tau_0, \tau_4)\f$</td></tr>
 * </table>
 *
 * If \c massChoice is not set properly an error condition will occur and
 * the function is aborted with status code defined by
 * #LALINSPIRALH_EMASSCHOICE in \ref LALInspiral.h.
 * In the above list \f$m_{1}\f$ and \f$m_{2}\f$ are the masses of
 * the two compact objects, \f$m=m_{1}+m_{2}\f$ is the total
 * mass, \f$\eta=m_{1}m_{2}/(m_{1}+m_{2})^{2}\f$ is the
 * symmetric mass ratio, \f$\mu=m_{1}m_{2}/(m_{1}+m_{2})\f$ is
 * the reduced mass and \f$\tau\f$'s are the chirptimes
 * defined in terms of \f$f_{a}\f$=\c fLower by:
 * \f{eqnarray}{
 * \tau_{0} = \frac{5}{256 \eta m^{5/3} (\pi f_{a})^{8/3}}, \ \ \
 * \tau_{2} = \frac{(3715 + 4620 \eta)}{64512 \eta m (\pi f_{a})^{2}}, \ \ \
 * \tau_{3} = \frac{\pi}{8 \eta m^{2/3} (\pi f_{a})^{5/3}} \\
 * \tau_{4} = \frac{5}{128 \eta m^{1/3} (\pi f_{a})^{4/3}} \left[ \frac{3058673}{1016064} +
 * \frac{5429}{1008} \eta + \frac{617}{144} \eta^{2} \right],\ \ \
 * \tau_5 = \frac {5}{256\eta f_a}  \left (\frac {7729}{252} + \eta \right ).
 * \f}
 * %% Beyond 2.5 PN order, chirp times do not have an
 * %% explicit expression in terms of the masses and \f$f_a.\f$
 * Whichever pair of parameters is given to the function as an input, the function
 * calculates the rest.  Apart from the various masses and chirptimes the function
 * also calculates the chirp mass \f$\mathcal{M}=(\mu^{3} m^{2})^{1/5}\f$ and
 * the total chirp time \f$\tau_C\f$ consistent with the approximation chosen:
 *
 * <table align="center" class="doxtable">
 * <caption align="top" style="text-align: left; font-weight: normal;">
 * Table II: \f$t_C\f$ will be set according to the PN order chosen in <tt>params->approximant.</tt>
 * </caption>
 * <tr><th></th><th>Newtonian</th><th>onePN</th><th>onePointFivePN</th><th>twoPN</th><th>twoPointFivePN</th></tr>
 * <tr><td>\f$\tau_C\f$</td><td>\f$\tau_0\f$</td><td>\f$\tau_0 + \tau_2\f$</td><td>\f$\tau_0 + \tau_2-\tau_3\f$
 * </td><td>\f$\tau_0 + \tau_2-\tau_3 + \tau_4\f$</td><td>\f$\tau_0 + \tau_2-\tau_3 + \tau_4 - \tau_5\f$</td></tr>
 * </table>
 *
 * ### Algorithm ###
 *
 * Root finding by bisection method is used to solve for mass ratio \f$\eta\f$ when
 * chirptimes \f$(\tau_0,\, \tau_2)\f$ or \f$(\tau_0,\, \tau_4)\f$ is input.
 *
 * ### Uses ###
 *
 * When appropriate this function calls:
 * <code>
 * XLALDBisectionFindRoot()
 * XLALEtaTau02()
 * XLALEtaTau04()
 * </code>
 *
 * ### Notes ###
 *
 */



#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

void
LALInspiralParameterCalc (
   LALStatus        *status,
   InspiralTemplate *params
   )
{
   XLALPrintDeprecationWarning("LALInspiralParameterCalc", "XLALInspiralParameterCalc");

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   XLALInspiralParameterCalc(params);
   if (xlalErrno)
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralParameterCalc (
   InspiralTemplate *params
   )
{
   REAL8 m1, m2, totalMass, eta, mu, piFl, etamin, tiny, ieta;
   REAL8 x1, x2, A0, A2, A3, A4, B2, B4, C4,v,tN;
   REAL8 theta = -11831.L/9240.L;
   REAL8 lambda = -1987.L/3080.L;
   static REAL8 oneby4;
   void *pars;
   REAL8 (*rootfunction)(REAL8, void *);
   REAL8 xmin, xmax, xacc;
   EtaTau02In Tau2In;
   EtaTau04In Tau4In;

   if (params == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if ((INT4)params->massChoice < 0)
      XLAL_ERROR(XLAL_EDOM);
   if ((INT4)params->massChoice > 15)
      XLAL_ERROR(XLAL_EDOM);

   totalMass 	= 0.0;
   ieta 	= params->ieta;
   ieta 	= 1.;
   oneby4 	= 1./4.;
   etamin 	= 1.e-10;
   tiny 	= 1.e-10;
   piFl 	= LAL_PI * params->fLower;

   switch(params->massChoice)
   {
      case massesAndSpin:
      /*case spinOnly:*/
      case minmaxTotalMass:
      case m1Andm2:
      case fixedMasses:

         if (params->mass1 <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->mass2 <= 0)
            XLAL_ERROR(XLAL_EDOM);

         m1 = params->mass1;
         m2 = params->mass2;
         params->totalMass = totalMass = m1+m2;
         params->eta = eta = m1*m2/(totalMass*totalMass);
         if (params->eta > oneby4) {
      		 params->eta -= tiny;
         }
         params->mu = mu = m1*m2/totalMass;
         params->chirpMass = pow(mu,0.6)*pow(totalMass,0.4);
         params->psi0 = 3./128./params->eta
	                * 1. * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-5./3.) ;
         params->psi3 = -3./128./params->eta
	                * (16 * LAL_PI) * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-2./3.);

      break;

      case totalMassAndEta:
      case totalMassUAndEta:

         if (params->totalMass <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->eta <= 0)
            XLAL_ERROR(XLAL_EDOM);

         if (params->eta > oneby4) {
		params->eta -= tiny;
         }
         if (params->eta > oneby4)
            XLAL_ERROR(XLAL_EDOM);

         totalMass = params->totalMass;
         eta = params->eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->mu = eta*totalMass;
         params->chirpMass = pow(eta,0.6)*totalMass;

      break;

      case totalMassAndMu:

         if (params->totalMass <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->mu <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->mu >= params->totalMass)
            XLAL_ERROR(XLAL_EDOM);

         totalMass = params->totalMass;
         mu = params->mu;
         eta =  (params->mu)/totalMass;
         if (eta > oneby4) {
		 eta -= tiny;
	 }
            params->eta = eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t02:

         if (params->t0 <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->t2 <= 0)
            XLAL_ERROR(XLAL_EDOM);

         A0 = 5./ pow(piFl, (8./3.))/256.;
         A2 = 3715.0/(64512.0*pow(piFl,2.0));
         B2 = 4620.0/3715 * ieta;
         Tau2In.t2 = params->t2;
         Tau2In.A2 = A2 * pow(params->t0/A0, 0.6);
         Tau2In.B2 = B2;

	 pars = (void *) &Tau2In;
         rootfunction = &XLALEtaTau02;
         xmax = oneby4+tiny;
         xmin = etamin;
         xacc = 1.e-8;
         x1 = XLALEtaTau02(xmax, pars);
         if (XLAL_IS_REAL8_FAIL_NAN(x1))
            XLAL_ERROR(XLAL_EFUNC);
         x2 = XLALEtaTau02(xmin, pars);
         if (XLAL_IS_REAL8_FAIL_NAN(x2))
            XLAL_ERROR(XLAL_EFUNC);

         if (x1*x2 > 0) {
            params->eta = 0.;
            return XLAL_SUCCESS;
         } else {
            eta = XLALDBisectionFindRoot(rootfunction, xmin, xmax, xacc, pars);
            if (XLAL_IS_REAL8_FAIL_NAN(eta))
               XLAL_ERROR(XLAL_EFUNC);
         }
         if (eta > oneby4) {
		 eta-=tiny;
   	}
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t03:

         if (params->t0 <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->t3 <= 0)
            XLAL_ERROR(XLAL_EDOM);

         A0 = 5./ pow(piFl, (8./3.))/256.;
         A3 = LAL_PI / pow(piFl, (5./3.))/8.;
         totalMass = A0 * params->t3/(A3 * params->t0);
         eta = A0/(params->t0 * pow(totalMass, (5./3.)));

	 if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t04:

         if (params->t0 <= 0)
            XLAL_ERROR(XLAL_EDOM);
         if (params->t4 <= 0)
            XLAL_ERROR(XLAL_EDOM);

	 A0 = 5./(256. * pow(piFl, (8./3.)));
         A4 = 5./(128.0 * pow(piFl,(4./3.))) * 3058673./1016064.;
         B4 = 5429./1008 * 1016064./3058673. * ieta;
         C4 = 617./144. * 1016064./3058673. * ieta;
         Tau4In.t4 = params->t4;
         Tau4In.A4 = A4 * pow(params->t0/A0, 0.2);
         Tau4In.B4 = B4;
         Tau4In.C4 = C4;

	 pars = (void *) &Tau4In;
         rootfunction = &XLALEtaTau04;
         xmax = oneby4+tiny;
         xmin = etamin;
         xacc = 1.e-8;
         x1 = XLALEtaTau04(xmax, pars);
         if (XLAL_IS_REAL8_FAIL_NAN(x1))
            XLAL_ERROR(XLAL_EFUNC);
         x2 = XLALEtaTau04(xmin, pars);
         if (XLAL_IS_REAL8_FAIL_NAN(x2))
            XLAL_ERROR(XLAL_EFUNC);

	 if (x1*x2 > 0) {
            params->eta = 0.;
            return XLAL_SUCCESS;
         } else {
            eta = XLALDBisectionFindRoot(rootfunction, xmin, xmax, xacc, pars);
            if (XLAL_IS_REAL8_FAIL_NAN(eta))
               XLAL_ERROR(XLAL_EFUNC);
         }
         if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;


      case psi0Andpsi3:
      if (params->psi0 > 0 && params->psi3 < 0)
      {
	      params->totalMass = totalMass = -params->psi3/(16.L * LAL_PI * LAL_PI * params->psi0)/LAL_MTSUN_SI;
	      params->eta = eta = 3.L/(128.L * params->psi0 * pow (LAL_PI * totalMass*LAL_MTSUN_SI, (5./3.)));

	      /* if eta < 1/4 amd M > 0 then physical values*/
	      if (eta <= oneby4)
	      {
		      params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
		      params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
		      params->mu = eta*totalMass;
		      params->chirpMass = pow(eta,0.6)*totalMass;
	      }
      }
      else
      {
	      params->eta = 0.;
	      return XLAL_SUCCESS;
      }
      break;

     default:
      XLALPrintError("XLAL Error - %s: Improper choice for massChoice\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
      break;
   }

   if (params->eta > oneby4) {
	   params->eta-=tiny;
	}
   totalMass 	= totalMass*LAL_MTSUN_SI;

   /* Should use the coefficients from LALInspiraSetup.c to avoid errors.
    * */
   v = pow(piFl * totalMass, 1.L/3.L);
   tN = 5.L/256.L / eta * totalMass / pow(v,8.L);

   params->t0 	= 5.0L/(256.0L*eta*pow(totalMass,(5./3.))*pow(piFl,(8./3.)));
   params->t2 	= (3715.0L + (4620.0L*ieta*eta))/(64512.0*eta*totalMass*pow(piFl,2.0));
   params->t3 	= LAL_PI/(8.0*eta*pow(totalMass,(2./3.))*pow(piFl,(5./3.)));
   params->t4 	= (5.0/(128.0*eta*pow(totalMass,(1./3.))*pow(piFl,(4./3.))))
              	* (3058673./1016064. + 5429.*ieta*eta/1008. +617.*ieta*eta*eta/144.);
   params->t5 	= -5.*(7729./252. - 13./3.*ieta*eta)/(256.*eta*params->fLower);
   /* This is a ddraft. t6 and t7 need to be checked propely*/
   params->t6 =  -10052469856691./23471078400. + 128./3.*LAL_PI*LAL_PI
     +(15335597827.L/15240960.L-451.L/12.L*LAL_PI*LAL_PI+352./3.*theta-2464.L/9.L*lambda)*ieta*eta
     +6848.L/105.L* LAL_GAMMA
     -15211.L/1728.L*ieta*eta*eta+25565.L/1296.L*eta*eta*eta*ieta;
   params->t6 = tN * (params->t6  + 6848.L/105.L*log(4.*v)) * pow(v,6);
   params->t7 = (-15419335.L/127008.L-75703.L/756.L*ieta*eta+14809.L/378.L*ieta*eta*eta) * LAL_PI * tN * pow(v,7);

   params->psi0 = 3.L/(128.L * eta * pow(LAL_PI * totalMass, (5./3.)));
   params->psi3 = -3.L * LAL_PI/(8.L * eta * pow(LAL_PI * totalMass, (2./3.)));

   switch (params->order) {

      case LAL_PNORDER_NEWTONIAN:
      case LAL_PNORDER_HALF:
         params->t2=0.0;
/*       params->t3=0.0;*/
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0;
      break;

      case LAL_PNORDER_ONE:
         params->t3=0.0;
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2;
      break;

      case LAL_PNORDER_ONE_POINT_FIVE:
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3;
      break;

      case LAL_PNORDER_TWO:
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4;
      break;

      case LAL_PNORDER_TWO_POINT_FIVE:
         params->t6 = 0.0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5;

      case LAL_PNORDER_THREE:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6;

      case LAL_PNORDER_THREE_POINT_FIVE:
      default:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0 and t7=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6 - params->t7;
      break;
   }

   return XLAL_SUCCESS;

}
