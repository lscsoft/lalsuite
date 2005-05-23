#include <lal/LALStdlib.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALConstants.h>


NRCSID (LALESTIMATEEFFECTIVEDISTANCEC, "$Id$");

void LALEstimateEffectiveDistance (
        LALStatus               *status,
        InspiralTemplate        param,
        REAL8                   df,
        REAL8Vector             *psd,
        REAL8                   snr,
        REAL8                   *effDistance
        )
{
    REAL8   msevenby3, powerNorm, totalMass, eta, flso, distanceNorm;
    REAL8   f, ins_amp;
    INT4    i;
    REAL8   dynRange = 1.0;

    INITSTATUS (status, "LALEstimateEffectiveDistance", LALESTIMATEEFFECTIVEDISTANCEC);
    ATTATCHSTATUSPTR (status);

    msevenby3 = -7.L/3.L;
    powerNorm = 0.;
    totalMass = param.totalMass*LAL_MTSUN_SI;
    eta       = param.eta;

    flso      = 1.L/(pow(6.L,1.5L)*totalMass*LAL_PI);

    /* Integrate over frequency - we can't start from i=0 (f=0) because f^-7/3 will become infinity */
    for (i=1; i<(INT4)psd->length; i++)  {
        f = i*df;
        if (f > flso) break;
        if (psd->data[i]) powerNorm += pow(f,msevenby3)/psd->data[i];
    } 

    /* I am not sure if there should be a factor of 2.0 here inside the sqrt ()
       i.e distanceNorm = sqrt(2.0*powerNorm * df);
       Multiplying by 2.0 makes dist agree with HW injections */
    distanceNorm = 2.*sqrt(2.*powerNorm * df);

    ins_amp = (LAL_MTSUN_SI * LAL_C_SI / (1.0e6 *  LAL_PC_SI)) 
            * sqrt( 5.0*param.mu / 96.0 ) 
            * ( pow( param.totalMass/(LAL_PI*LAL_PI) , 0.33333 ) / pow(LAL_MTSUN_SI, 1.0 / 6.0) ) ;

    distanceNorm *= (ins_amp * sqrt(dynRange));

    /* We need to calculate randIn.SignalAmp = distanceNorm / deff (in
     * Mpc)*/
    (*effDistance) = (distanceNorm / snr); 

    /* Normal exit */
    DETATCHSTATUSPTR (status);
    RETURN (status);

}
