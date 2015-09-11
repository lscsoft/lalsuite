#include <stdlib.h>
#include <gsl/gsl_vector.h>



/*#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lal/LALSimIMR.h>
#include <lal/Date.h>

#include <lal/SeqFactories.h>*/

#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimFindAttachTime.h"
#include "LALSimIMRSpinEOBHamiltonian.c"





double  XLALSimLocateOmegaTime(
    REAL8Array *dynamicsHi,
    unsigned int numdynvars,
    unsigned int retLenHi,
    SpinEOBParams   seobParams,
    SpinEOBHCoeffs  seobCoeffs,
    REAL8 m1,
    REAL8 m2,
    REAL8 *radiusVec,
    int *found,
    REAL8* tMaxOmega
        )
{        
    /* 
    * Locate merger point (max omega), 
    * WaveStep 1.1: locate merger point */
    int debugPK = 0;
    int debugRD = 0;
    FILE *out = NULL; 
    gsl_spline    *spline = NULL;
    gsl_interp_accel *acc = NULL;
    
    if (debugPK) {debugRD = 0;}
    
    unsigned int peakIdx, i, j;
    REAL8Vector *values = NULL;
    REAL8Vector *dvalues = NULL;
    REAL8Vector *omegaHi = NULL;

    
    if ( !(values = XLALCreateREAL8Vector( numdynvars )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    if ( !(dvalues = XLALCreateREAL8Vector( numdynvars )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    if ( !(omegaHi = XLALCreateREAL8Vector( retLenHi )) )
    {
        XLAL_ERROR(  XLAL_ENOMEM );
    }
    REAL8 rdotvec[3] = {0,0,0};
    REAL8 rvec[3] = {0,0,0};
    REAL8 rcrossrdot[3] = {0,0,0};
    REAL8Vector timeHi;
  
    timeHi.length = retLenHi;
    timeHi.data = dynamicsHi->data;
    
    double dt = timeHi.data[1] - timeHi.data[0];
    double ddradiusVec[timeHi.length - 2];
    unsigned int k;
    for (k = 1; k < timeHi.length-1; k++) {
        ddradiusVec[k] = (radiusVec[k+1] - 2.*radiusVec[k] + radiusVec[k-1])/dt/dt;
        //        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
    }
    //    for (k = timeHi->length-3; k>=1; k--) {
    //        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
    //        if (ddradiusVec[k] < 0) {
    //            break;
    //        }
    //    }
    for (k = 0; k < timeHi.length-2; k++) {
        //        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
        if (dt*k > dt*( timeHi.length-2)-20 && ddradiusVec[k] > 0) {
            break;
        }
    }
    double minoff = dt*( timeHi.length-2 - k) > 0.2 ? dt*( timeHi.length-2 - k) : 0.2;
    if (debugPK) {
        printf("Change of sign in ddr %3.10f M before the end\n", minoff );
    }
    // First we search for the maximum (extremum) of amplitude
    double maxoff = 20.0;
    unsigned int Nps = timeHi.length;
    // this definesthe search interval for maximum (we might use min0ff= 0.051 instead)
    double tMin  = timeHi.data[Nps-1] - maxoff;
    // FIXME
    // minoff = 0.0;
    double tMax = timeHi.data[Nps-1] - minoff;
    *tMaxOmega = tMax;
    tMin = tMax - 20.;
    if ( debugPK ) {
        printf("tMin, tMax = %3.10f %3.10f\n", tMin, tMax);
    }
    
    double omega  = 0.0;

    double magR;
    double time1, time2, omegaDeriv, omegaDerivMid, tPeakOmega;
  
    if(debugPK) {
        out = fopen( "omegaHi.dat", "w" );
        printf("length of values = %d, retLenHi = %d\n", values->length, retLenHi);
        fflush(NULL);
    }
    if(debugRD) {
        out = fopen( "omegaHi.dat", "w" );
    }
  
    for ( i = 0; i < retLenHi; i++ )
    {
        for ( j = 0; j < values->length; j++ )
            { values->data[j] = *(dynamicsHi->data+(j+1)*retLenHi+i); }
    
        /* Calculate dr/dt */
        memset( dvalues->data, 0, numdynvars*sizeof(dvalues->data[0]));
        if( XLALSpinHcapRvecDerivative( 0, values->data, dvalues->data, 
            &seobParams) != XLAL_SUCCESS )
        {
                printf(
                    " Calculation of dr/dt failed while computing omegaHi time series\n");
                XLAL_ERROR( XLAL_EFUNC );
        }
    
        /* Calculare r x dr/dt */
        for (j=0; j<3; j++){
            rvec[j] = values->data[j];
            rdotvec[j] = dvalues->data[j];
        }
    
        //memcpy(rdotvec, dvalues->data, 3*sizeof(REAL8));
        //rvec[0] = posVecxHi.data[i]; rvec[1] = posVecyHi.data[i]; 
        //rvec[2] = posVeczHi.data[i];
        cross_product( rvec, rdotvec, rcrossrdot );        
   
        /* Calculate omega = |r x dr/dt| / r*r */
        magR = sqrt(inner_product(rvec, rvec));
        omegaHi->data[i] = sqrt(inner_product(rcrossrdot, rcrossrdot)) / (magR*magR); 
        if(debugPK || debugRD){
            fprintf( out, "%.16e\t%.16e\n", timeHi.data[i], omegaHi->data[i]);
        }
    }


    // Searching for crude omega_max (extremum)
    peakIdx = 0;
    *found = 0;
    for ( i = 1, peakIdx = 0; i < retLenHi-1; i++ ){
        omega = omegaHi->data[i];

        if (omega >= omegaHi->data[i-1] && omega > omegaHi->data[i+1] && tMax>=timeHi.data[i] && timeHi.data[i]>=tMin){
            peakIdx = i;
            *found = 1;
            if (debugPK){
                printf("PK: Crude peak of Omega is at idx = %d. t = %f,  OmegaPeak = %.16e\n", 
                    peakIdx, timeHi.data[peakIdx], omega);
                fflush(NULL);

            }
        }  
    }
    
    if(debugPK) {
        fclose(out);
        if (peakIdx ==0){
            printf("Stas: peak of orbital frequency was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, retLenHi, i);
            fflush(NULL);
        }
    }
    if(debugRD) {
        fclose(out);
    }
  
    // refining the omega_max search (if it is found)
    tPeakOmega = 0.0;
    if(peakIdx != 0){
        spline = gsl_spline_alloc( gsl_interp_cspline, retLenHi );
        acc    = gsl_interp_accel_alloc();
        time1 = timeHi.data[peakIdx-2];
        gsl_spline_init( spline, timeHi.data, omegaHi->data, retLenHi );
        omegaDeriv = gsl_spline_eval_deriv( spline, time1, acc );
   
        if ( omegaDeriv > 0. ) { time2 = timeHi.data[peakIdx+2]; }
        else{
            time2 = time1;
            peakIdx = peakIdx-2;
	        time1 = timeHi.data[peakIdx-2];	      
	        omegaDeriv = gsl_spline_eval_deriv( spline, time1, acc );
        }
   
        do
        {
            tPeakOmega = ( time1 + time2 ) / 2.;
	        omegaDerivMid = gsl_spline_eval_deriv( spline, tPeakOmega, acc );
	   
	        if ( omegaDerivMid * omegaDeriv < 0.0 ) { time2 = tPeakOmega; }
	        else
	        {
		        omegaDeriv = omegaDerivMid;
		        time1 = tPeakOmega;
		    }
            if (debugPK){
                printf("Stas: searching for orbital max: %f, %f, %f, %f \n", time1, time2, omegaDeriv, omegaDerivMid);
            }
        } while ( time2 - time1 > 1.0e-5 );
        if(debugPK) {
          printf( "Estimation of the orbital peak is now at time %.16e \n", tPeakOmega);
          fflush(NULL);
        }
    }
  
    if(*found == 0 || debugRD || debugPK){
        if(debugPK){
            printf("Stas: We couldn't find the maximum of orbital frequency, search for maximum of A(r)/r^2 \n");
        }   
        REAL8 rad, rad2, m1PlusetaKK, bulk, logTerms, deltaU, u, u2, u3, u4, u5;
        REAL8 listAOverr2[retLenHi];
        REAL8 Aoverr2;
        REAL8Vector *sigmaStar = NULL;
        REAL8Vector *sigmaKerr = NULL;
        if ( !(sigmaStar = XLALCreateREAL8Vector( 3 )) )
        {
          XLALDestroyREAL8Vector( sigmaStar );
          XLAL_ERROR( XLAL_ENOMEM );
        }
        if ( !(sigmaKerr = XLALCreateREAL8Vector( 3 )) )
        {
          XLALDestroyREAL8Vector( sigmaStar );
          XLAL_ERROR( XLAL_ENOMEM );
        }
        REAL8Vector s1Vec, s2Vec;
        s1Vec.length = s2Vec.length = 3;
        REAL8 s1Data[3], s2Data[3];
        REAL8 mTotal = m1 + m2;
        REAL8 a;
        REAL8 eta = m1*m2/(mTotal*mTotal);
        
        if(debugPK || debugRD){ 
            out = fopen( "OutAofR.dat", "w" );
        }
        for ( i = 0; i < retLenHi; i++ )
        {
            for ( j = 0; j < values->length; j++ )
            {
                values->data[j] = *(dynamicsHi->data+(j+1)*retLenHi+i);
            }
            for( j = 0; j < 3; j++ )
            {
                //s1DataNorm[k] = values->data[k+6];
                //s2DataNorm[k] = values->data[k+9];
                s1Data[j] = values->data[j+6] * mTotal * mTotal;
                s2Data[j] = values->data[j+9] * mTotal * mTotal;
            }
            s1Vec.data = s1Data;
            s2Vec.data = s2Data;
            XLALSimIMRSpinEOBCalculateSigmaStar( sigmaStar, m1, m2, &s1Vec, &s2Vec );
            XLALSimIMRSpinEOBCalculateSigmaKerr( sigmaKerr, m1, m2, &s1Vec, &s2Vec );
            
            seobParams.a = a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
            m1PlusetaKK = -1. + eta * seobCoeffs.KK;
            rad2 =  values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2];
            rad = sqrt(rad2);
            u = 1./rad;
            u2 = u*u;
            u3 = u2*u;
            u4 = u2*u2;
            u5 = u4*u;
            bulk = 1./(m1PlusetaKK*m1PlusetaKK) + (2.*u)/m1PlusetaKK + a*a*u2;
            logTerms = 1. + eta*seobCoeffs.k0 + eta*log(1. + seobCoeffs.k1*u + seobCoeffs.k2*u2 + seobCoeffs.k3*u3 + seobCoeffs.k4*u4 + seobCoeffs.k5*u5 + seobCoeffs.k5l*u5*log(u));
            deltaU = bulk*logTerms;
            listAOverr2[i] = deltaU / rad2;
            if(debugPK || debugRD){
                fprintf(out, "%3.10f %3.10f\n", timeHi.data[i], listAOverr2[i]);
            }
            
        }
        if(debugPK || debugRD ) fclose(out);
        if (*found == 0){
            // searching formaximum of A(r)/r^2
            peakIdx = 0;
            //*found = 0;
            for ( i = 1, peakIdx = 0; i < retLenHi-1; i++ ){
                Aoverr2 = listAOverr2[i];
                if (Aoverr2 >= listAOverr2[i-1] && Aoverr2 > listAOverr2[i+1]){
                    if (timeHi.data[i] > tMin){
                        peakIdx = i;
                        tPeakOmega = timeHi.data[i];
                        *found = 1;
                        if (debugPK){
                            printf("PK: Peak of A(r)/r^2 is at idx = %d. t = %f, Peak ampl. = %.16e\n", 
                                peakIdx, timeHi.data[peakIdx], Aoverr2);
                            fflush(NULL);
                        }
                        break;
                    }
                }  
            }
        }
    
        if(debugPK) {
            if (peakIdx ==0){
                printf("Stas: peak of A(r)/r^2 was not found. \
                    peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, retLenHi, i);
                fflush(NULL);
            }
        }
        XLALDestroyREAL8Vector(sigmaStar);
        XLALDestroyREAL8Vector(sigmaKerr);
    }
    if (spline != NULL)
        gsl_spline_free(spline);
    if (acc != NULL)
        gsl_interp_accel_free(acc);
    XLALDestroyREAL8Vector( values );
    XLALDestroyREAL8Vector( dvalues );
    XLALDestroyREAL8Vector( omegaHi );
    if (*found == 0){
        return(timeHi.data[retLenHi-1]);
    }
    else{
        return(tPeakOmega);
    }
}
  
double XLALSimLocateAmplTime(
    REAL8Vector *timeHi, 
    COMPLEX16Vector *hP22,
    REAL8 *radiusVec,
    int *found,
    REAL8* tMaxAmp)
{
    int debugPK = 0;
    int debugRD = 0;
    FILE *out = NULL; 
    gsl_spline    *spline = NULL;
    gsl_interp_accel *acc = NULL;
    if (debugPK) {debugRD = 0;}
    
    double dt = timeHi->data[1] - timeHi->data[0];
    double ddradiusVec[timeHi->length - 2];
    unsigned int k;
    for (k = 1; k < timeHi->length-1; k++) {
        ddradiusVec[k] = (radiusVec[k+1] - 2.*radiusVec[k] + radiusVec[k-1])/dt/dt;
//        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
    }
//    for (k = timeHi->length-3; k>=1; k--) {
//        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
//        if (ddradiusVec[k] < 0) {
//            break;
//        }
//    }
    for (k = 0; k < timeHi->length-2; k++) {
//        printf("%3.10f %3.10f\n", timeHi->data[k], ddradiusVec[k]);
        if (dt*k > dt*( timeHi->length-2)-20 && ddradiusVec[k] > 0) {
            break;
        }
    }
    double minoff = dt*( timeHi->length-2 - k) > 0.2 ? dt*( timeHi->length-2 - k) : 0.2;
    if (debugPK) {
        printf("Change of sign in ddr %3.10f M before the end\n", minoff );
    }
    // First we search for the maximum (extremum) of amplitude
    unsigned int i, peakIdx; 
    double maxoff = 20.0;
    unsigned int Nps = timeHi->length; 
    // this definesthe search interval for maximum (we might use min0ff= 0.051 instead)
    
    //FIXME
    //minoff = 0.0; 
    double tMin  = timeHi->data[Nps-1] - maxoff;
    double tMax = timeHi->data[Nps-1] - minoff;
    *tMaxAmp = tMax;
    tMin = tMax - maxoff;
    if ( debugPK ) {
        printf("tMin, tMax = %3.10f %3.10f \n", tMin, tMax);
    }
    unsigned int iMin = ceil(tMin/dt);
    unsigned int iMax = floor(tMax/dt);
    unsigned int NpsSmall = iMax - iMin + 1;


    double AmplN, AmplO;
    double tAmpMax, AmpMax, tAmp;
    tAmpMax = 0.;
    REAL8 tSeries[NpsSmall], Ampl[NpsSmall];
    
    if(debugPK || debugRD) {
            out = fopen( "AmpPHi.dat", "w" );
    }
    AmplO = sqrt(creal(hP22->data[iMin + 0])*creal(hP22->data[iMin + 0]) + cimag(hP22->data[iMin + 0])*cimag(hP22->data[iMin + 0]));
    Ampl[0] = AmplO;
    peakIdx = 0;
    for (i=0; i<NpsSmall-1; i++){
        tSeries[i] = timeHi->data[iMin + i];
        AmplN = sqrt(creal(hP22->data[iMin + i+1])*creal(hP22->data[iMin + i+1]) + cimag(hP22->data[iMin + i+1])*cimag(hP22->data[iMin + i+1]));
        //Ampl = sqrt(hreP22->data[i]*hreP22->data[i] + himP22->data[i]*himP22->data[i]);
        if(debugPK || debugRD){
            fprintf(out, "%3.10f %3.10f\n", tSeries[i], Ampl[i]);
        }
        if (Ampl[i] >= AmplO && Ampl[i] >AmplN){
            if (*found !=1){
                tAmp = timeHi->data[iMin + i];
                if (tAmp >=tMin && tAmp <= tMax ){
                    *found = 1;
                    tAmpMax = tAmp;
                    AmpMax = Ampl[i];
                    peakIdx = iMin + i;
                }else{
                    if (debugPK){
                        printf("Stas dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                            tAmp, tMin, tMax);
                    }
                }
            }        
        }
        AmplO = Ampl[i];
        Ampl[i+1] = AmplN;                
    }
    
    if (debugPK) 
    {
        fclose(out);
        if (*found ==0){
            printf("Stas: peak of 2,2 mode in P-frame was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, Nps, i);
            fflush(NULL);
        }else{
            printf("Stas: we have found maximum of amplitude %3.10f at t = %3.10f \n", AmpMax, tAmpMax);
        }
    }
    if (debugRD) 
    {
        fclose(out);
    }

    if (*found ==0 || debugRD || debugPK){
        // we haven't found the maximum of amplitude -> search for minimum of derivative (extremum)
//        spline = gsl_spline_alloc( gsl_interp_cspline, NpsSmall );
//        acc    = gsl_interp_accel_alloc();
//        gsl_spline_init( spline, tSeries, Ampl, NpsSmall );
        
        REAL8 AmpDot[NpsSmall];
        REAL8 AmpDDot[NpsSmall];
        
        for (i=1; i<NpsSmall-2; i++){
//            AmpDot[i] = gsl_spline_eval_deriv(spline, tSeries[i] , acc);
            AmpDot[i] = (Ampl[i+1] - Ampl[i-1])/2./dt;
            AmpDDot[i] = (Ampl[i+1] - 2.0*Ampl[i] + Ampl[i-1])/(dt*dt);
        }
        AmpDot[0] = AmpDot[1]; 
        AmpDot[NpsSmall -2] = AmpDot[NpsSmall-3];
        AmpDot[NpsSmall -1] = AmpDot[NpsSmall-2];
        AmpDDot[0] = AmpDDot[1]; 
        AmpDDot[NpsSmall -2] = AmpDDot[NpsSmall-3];
        AmpDDot[NpsSmall -1] = AmpDDot[NpsSmall-2];
        //printf("Stas, check AmDot %f, %f, %f \n", AmpDot[NpsSmall -3],  AmpDot[NpsSmall -2],  AmpDot[NpsSmall -1]);
        //printf("Stas, check AmDDot %f, %f, %f \n", AmpDDot[NpsSmall -3],  AmpDDot[NpsSmall -2],  AmpDDot[NpsSmall -1]);



        REAL8 AmpDotSmooth[NpsSmall];
        // Compute moving average over 7 points
        unsigned int win = 3;
       // unsigned int win = 5;
        //int j;
        double norm = 1.0/(2.0*win+1.0);

        AmpDotSmooth[win] = 0;
        for (i=0; i<(2*win +1); i++){
            AmpDotSmooth[win] += AmpDot[i];
        }
        AmpDotSmooth[win] *= norm;
        for (i=0; i<win; i++){
            AmpDotSmooth[i] = AmpDotSmooth[win];
        }

        for (i=win+1; i<NpsSmall-1 -win; i++){
            AmpDotSmooth[i] = AmpDotSmooth[i-1] + norm*(AmpDot[i+win] - AmpDot[i-win-1]);             
        }
        for (i=0; i<win; i++){
            AmpDotSmooth[NpsSmall-win-1+i] = AmpDotSmooth[NpsSmall-win-2];
        }
        
        // second deriv (in case)
        REAL8 AmpDDotSmooth[NpsSmall];
        unsigned int win2 = 100;
        // unsigned int win = 5;
        //int j;
        norm = 1.0/(2.0*win2+1.0);

        AmpDDotSmooth[win2] = 0;
        for (i=0; i<(2*win2 +1); i++){
            AmpDDotSmooth[win2] += AmpDDot[i];
        }
        AmpDDotSmooth[win2] *= norm;
        for (i=0; i<win2; i++){
            AmpDDotSmooth[i] = AmpDDotSmooth[win2];
        }

        for (i=win2+1; i<NpsSmall-1 -win2; i++){
            AmpDDotSmooth[i] = AmpDDotSmooth[i-1] + norm*(AmpDDot[i+win2] - AmpDDot[i-win2-1]);             
        }
        for (i=0; i<win2; i++){
            AmpDDotSmooth[NpsSmall-win2-1+i] = AmpDDotSmooth[NpsSmall-win2-2];
        }

        
        if(debugPK || debugRD) {
            out = fopen( "DotAmpPHi.dat", "w" );
            for (i=0; i<NpsSmall - 1; i++){
                fprintf(out, "%3.10f %3.10f %3.10f %3.10f %3.10f\n", tSeries[i], AmpDot[i], AmpDotSmooth[i], AmpDDot[i], AmpDDotSmooth[i]);
            } 

        }
        if (*found ==0){
            if (debugPK || debugRD){
                printf("Max of Amplitude is not found, looking for min of dot{Ampl} %d \n", iMin);
            }
            for (i=1; i<NpsSmall-1-win; i++){
                   if (AmpDotSmooth[i] < AmpDotSmooth[i-1] && AmpDotSmooth[i] < AmpDotSmooth[i+1]){
                        tAmp = tSeries[i];
                        //tAmp = tSeries[i-iMin];
                        //printf("Stas check i = %d, tSeries = %3.10f, tAmp = %3.10f \n", i, tSeries[i], tAmp);
                        if (tAmp >=tMin && tAmp <= tMax  && *found==0){
                            *found = 1;
                            tAmpMax = tAmp;
                            AmpMax = AmpDotSmooth[i];
                            //AmpMax = AmpDotSmooth[i-iMin];
                            peakIdx = i;
                            if (debugPK || debugRD){
                                printf("we have found min of Adot at t= %f\n", tAmpMax);
                            }
                            //break;
                        }else{
                            if (debugPK){
                                printf("Stas, AmplDot - dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                                    tAmp, tMin, tMax);
                            }
                        }                              
                   
                   }  
            }
        }
        if (*found ==0){
            if (debugPK || debugRD)
                printf("Min of Adot is not found, looking for min of ddot{Ampl} \n");
            for (i=win2*2; i<NpsSmall-1-win2; i++){
                   if (AmpDDotSmooth[i] < AmpDDotSmooth[i-1] && AmpDDotSmooth[i] < AmpDDotSmooth[i+1]){
                        tAmp = tSeries[i];
                        //printf("Stas check i = %d, tSeries = %3.10f, tAmp = %3.10f,  %3.10f, %3.10f, %3.10f \n", i, tSeries[i], tAmp,  AmpDDotSmooth[i-1], AmpDDotSmooth[i], AmpDDotSmooth[i+1]);
                        //if (tAmp >=tMin && tAmp <= tMax  && *found==0)
                        if (tAmp >=tMin && tAmp <= tMax){
                            *found = 1;
                            tAmpMax = tAmp;
                            AmpMax = AmpDDotSmooth[i];
                            peakIdx = i;
                            //break;
                        }else{
                            if (debugPK){
                                printf("Stas, AmplDDot - dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                                    tAmp, tMin, tMax);
                            }
                        }                              
                   }  
            }
        }
        
//        if(*found==0){
//            REAL8 hRe[NpsSmall], hIm[NpsSmall];
//            REAL8 dhRe[NpsSmall], dhIm[NpsSmall];
//            for (i=0; i<NpsSmall - 1; i++) {
//                hRe[i] = creal(hP22->data[iMin + i]);
//                hIm[i] = cimag(hP22->data[iMin + i]);
//            }
//            for (i=1; i<NpsSmall - 2; i++) {
//                dhRe[i] = (hRe[i+1] - hRe[i-1])/2./dt;
//                dhIm[i] = (hIm[i+1] - hIm[i-1])/2./dt;
//            }
//            dhRe[0]=dhRe[1];
//            dhIm[0]=dhIm[1];
//            dhRe[NpsSmall-1]=dhRe[NpsSmall-2];
//            dhIm[NpsSmall-1]=dhIm[NpsSmall-2];
//            
//            REAL8 OmegaWave[NpsSmall], dOmegaWave[NpsSmall];
//            double hNorm2;
//            for (i=0; i<NpsSmall - 1; i++) {
//                hNorm2 = hRe[i]*hRe[i] + hIm[i]*hIm[i];
//                OmegaWave[i] = (-hRe[i]*dhIm[i] + hIm[i]*dhRe[i])/hNorm2;
//            }
//            for (i=1; i<NpsSmall - 2; i++) {
//                hNorm2 = hRe[i]*hRe[i] + hIm[i]*hIm[i];
//                dOmegaWave[i] = (OmegaWave[i+1] - OmegaWave[i-1])/2./dt;
//            }
//            dOmegaWave[0]=dOmegaWave[1];
//            dOmegaWave[NpsSmall-1]=dOmegaWave[NpsSmall-2];
//            
//            if(debugPK || debugRD) {
//                out = fopen( "OmegaW.dat", "w" );
//                for (i=0; i<NpsSmall - 1; i++){
//                    fprintf(out, "%3.10f %3.10f %3.10f\n", tSeries[i], OmegaWave[i], dOmegaWave[i]);
//                }
//            }
//            fclose(out);
//
//        }
        
        /*
        for (i=1; i<Nps-3; i++){
            AmplDerivN2 = gsl_spline_eval_deriv(spline, timeHi->data[i+2], acc);
            if(debugPK || debugRD){
                fprintf(out, "%3.10f %3.10f\n", timeHi->data[i], AmplDeriv);
            }

            if (*found == 0){
                if ((AmplDeriv  <= AmplDerivO1 && AmplDeriv < AmplDerivN1) && (AmplDerivO2 > AmplDerivO1 && AmplDerivN2>=AmplDerivN1) ){
                    printf("check %.16e, %.16e, %.16e, %.16e, %.16e \n", AmplDerivO2, AmplDerivO1, AmplDeriv, AmplDerivN1, AmplDerivN2);
                    tAmp = timeHi->data[i];
                    if (tAmp >=tMin && tAmp <= tMax  && *found==0){
                        *found = 1;
                        tAmpMax = tAmp;
                        AmpMax = AmplDeriv;
                        peakIdx = i;
                        //break;
                    }else{
                        if (debugPK){
                            printf("Stas dismissing time %3.10f outside limits %3.10f, %3.10f \n", 
                                tAmp, tMin, tMax);
                        }
                    }                              
                }
            }
            AmplDerivO2 = AmplDerivO1;
            AmplDerivO1 = AmplDeriv;
            AmplDeriv = AmplDerivN1;
            AmplDerivN1 = AmplDerivN2;
        }
        
        if (debugPK) 
        {
            fclose(out);
            if (*found ==0){
                printf("Stas: peak of 2,2 mode in P-frame was not found. peakIdx = %d, retLenHi = %d, i at exit = %d\n", peakIdx, Nps, i);
                fflush(NULL);
            }else{
                printf("Stas: we have found maximum of amplitude %3.10f at t = %3.10f \n", AmpMax, tAmpMax);
            }
        }
         if (debugRD) 
        {
            fclose(out);
        }
        */
        
    }
    

    
    if (spline != NULL)
        gsl_spline_free(spline);
    if (acc != NULL)
        gsl_interp_accel_free(acc);
    if (*found == 0){
        return(timeHi->data[Nps-1]);
    }
    else{
        return(tAmpMax);
    }
     
}   


INT4 XLALSimCheckRDattachment(
    REAL8Vector * signal1,	/**<< Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2,	/**<< Imag of inspiral waveform to which we attach ringdown */
    REAL8* ratio,           /**<< output ratio  */
    const REAL8 tAtt,       /**<< time of RD attachment */
    const INT4 l,	/**<< Current mode l */
    const INT4 m,	/**<< Current mode m */
    const REAL8 dt,	/**<< Sample time step (in seconds) */
    const REAL8 mass1,	/**<< First component mass (in Solar masses) */
    const REAL8 mass2,	/**<< Second component mass (in Solar masses) */
    const REAL8 spin1x,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1y,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1z,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin2x,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2y,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2z,	/**<<The spin of the second object; only needed for spin waveforms */
    REAL8Vector * timeVec,	/**<< Vector containing the time values */
    REAL8Vector * matchrange,	/**<< Time values chosen as points for performing comb matching */
    Approximant approximant,	/**<<The waveform approximant being used */
    const REAL8 JLN           /**<< cosine of the angle between J and LN at the light ring */
    )
{
    int debugPK = 0;
    unsigned int i;
    unsigned int i_att = 0;
    REAL8 Amp[signal1->length];
    // sanity check 
    int ind_att = (int) matchrange->data[1]*(((mass1 + mass2) * LAL_MTSUN_SI / dt)) + 1;
    if (debugPK){
        printf("attach_ind = %d, t =%f, %f \n", ind_att, matchrange->data[1], timeVec->data[ind_att]); 
    } 
    if (signal1->data[ind_att] == 0.0 && signal2->data[ind_att] == 0.0){
        printf("Opyat' signal = 0 \n");
        //FILE *out1 = fopen( "Andrea1.dat","w");
        //for (i = 0; i < timeVec->length; i++) {
        //    fprintf(out1, "%.16e   %.16e   %.16e\n", timeVec->data[i], signal1->data[i], signal2->data[i]);
            //printf("%.16e %.16e\n", timeVec->data[j], y[j]);
        //}
        //fclose(out1);
        //exit(0);
        XLAL_ERROR(XLAL_EFAILED);


    }
 
    
    if ( XLALSimIMREOBHybridAttachRingdown( signal1, signal2, l, m,
                dt, mass1, mass2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                timeVec, matchrange, approximant, JLN )
            == XLAL_FAILURE )
    {
        XLAL_ERROR( XLAL_EFUNC );
    }
   
    Amp[0] = sqrt(signal1->data[0]*signal1->data[0] + signal2->data[0]*signal2->data[0]);
    for (i=1; i<signal1->length; i++){
        Amp[i] = sqrt(signal1->data[i]*signal1->data[i] + signal2->data[i]*signal2->data[i]);
        if (timeVec->data[i-1] <= tAtt && timeVec->data[i] > tAtt){
            i_att = i;
            //if (debugPK){
            //    printf(" the attachment index = %d, time = %f, %f \n", i_att, timeVec->data[i_att], tAtt);
            //}
        } 
    }

    REAL8 maxL = Amp[0]; 
    for (i=0; i<i_att; i++){
        if (Amp[i] >= maxL){
           maxL = Amp[i];
        } 
    }  
    REAL8 maxR = Amp[i_att];
    for (i=i_att; i<signal1->length; i++){
        if (Amp[i] >= maxR){
           maxR = Amp[i];
        } 
    }  
     
    if (debugPK){
        printf(" the ratio of amplitudes = %f  , ampls = %f, %f \n", maxR/maxL, maxR, maxL);
    }

    *ratio = maxR/maxL;
    if (maxR/maxL != maxR/maxL){
        //this is nan
        *ratio = 1000.0;
    }


    return XLAL_SUCCESS;
   
}


int XLALSimAdjustRDattachmentTime( 
    REAL8Vector * signal1,	/**<< Output Real of inspiral waveform to which we attach ringdown */
    REAL8Vector * signal2,	/**<< Output Imag of inspiral waveform to which we attach ringdown */
    COMPLEX16TimeSeries* h22,   /**<< input time series (inspiral) */
    COMPLEX16TimeSeries* h2m2,  /**<< input time series (inspiral) */
    REAL8* ratio22,      /**<< output ratio for 2,2 mode */
    REAL8* ratio2m2,     /**<< output ratio  for 2,-2 mode*/
    REAL8* tAttach,       /**<< output/input time of RD attachment */
    const REAL8 thr,        /**<< threshold on the ratio */
    const REAL8 dt,	/**<< Sample time step (in seconds) */
    const REAL8 m1,	/**<< First component mass (in Solar masses) */
    const REAL8 m2,	/**<< Second component mass (in Solar masses) */
    const REAL8 spin1x,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1y,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin1z,	/**<<The spin of the first object; only needed for spin waveforms */
    const REAL8 spin2x,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2y,	/**<<The spin of the second object; only needed for spin waveforms */
    const REAL8 spin2z,	/**<<The spin of the second object; only needed for spin waveforms */
    REAL8Vector * timeVec,	/**<< Vector containing the time values */
    REAL8Vector * matchrange,	/**<< Time values chosen as points for performing comb matching */
    Approximant approximant,	/**<<The waveform approximant being used */
    const REAL8 JLN,            /**<< cosine of the angle between J and LN at the light ring */
    const REAL8 combSize,        /**<< combsize for RD attachment */
    const REAL8 tMaxOmega,
    const REAL8 tMaxAmp
    )
{
    int debugPK = 0;
    unsigned int retLenHi = h22->data->length;
    unsigned int i;
    int pass = 0;
    REAL8 tAtt;
    tAtt = *tAttach; // replace with the loop
    REAL8 maxDeltaT = 10.0;
    REAL8 thrStore22L = 0., thrStore2m2L = 0., thrStore22R = 0., thrStore2m2R = 0., tLBest = *tAttach, tRBest = *tAttach;

    REAL8 mTScaled = (retLenHi-1)*dt/matchrange->data[2];
    REAL8 tMax = timeVec->data[retLenHi - 2] - 0.5 ;
    double hNorm2, dsignal1, dsignal2;    double omegaVec[retLenHi - 1];

//    printf("tMaxOmega, tMaxAmp %f %f %f\n", tMaxOmega, tMaxAmp, tMax);

    if ( tMaxAmp < tMax) {
        tMax = tMaxAmp;
    }
    if ( tMaxOmega < tMax) {
        tMax = tMaxOmega;
    }
    tMax = tMax - 3.;
    if(tMax > tAtt + 5.0){
        tMax = tAtt + 5.0;
    }
    if (debugPK){
        printf("tmax = %f, tAtt = %f, tmaxAmp = %f, tmaxOm = %f\n", tMax, tAtt, tMaxAmp, tMaxOmega);
    }  
//    printf("tAtt, tMax = %f %f\n", tAtt, tMax);
    while(pass == 0 && (tAtt >= *tAttach - maxDeltaT)){
        tAtt = tAtt - 0.5;
        memset( signal1->data, 0, signal1->length * sizeof( signal1->data[0] ));
        memset( signal2->data, 0, signal2->length * sizeof( signal2->data[0] ));
        for ( i = 0; i < retLenHi; i++ )
        {
            signal1->data[i] = creal(h22->data->data[i]);
            signal2->data[i] = cimag(h22->data->data[i]);
        }
       
        matchrange->data[0] = combSize < tAtt ? tAtt - combSize : 0;
        matchrange->data[1] = tAtt;
        matchrange->data[0] -= fmod( matchrange->data[0], dt/mTScaled );
        matchrange->data[1] -= fmod( matchrange->data[1], dt/mTScaled );
        if (debugPK) printf("left 2,2 mode tAtt = %f     ", tAtt); 
        if( XLALSimCheckRDattachment(signal1, signal2, ratio22, tAtt, 2, 2,
                        dt, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                        timeVec, matchrange, approximant, JLN ) == XLAL_FAILURE )
        {
              XLAL_ERROR( XLAL_EFUNC );
        }
     
        memset( signal1->data, 0, signal1->length * sizeof( signal1->data[0] ));
        memset( signal2->data, 0, signal2->length * sizeof( signal2->data[0] ));
        for ( i = 0; i < retLenHi; i++ )
        {
            signal1->data[i] = creal(h2m2->data->data[i]);
            signal2->data[i] = cimag(h2m2->data->data[i]);
        }
       
        if (debugPK) printf("left 2,-2 mode tAtt = %f     ", tAtt); 
        if( XLALSimCheckRDattachment(signal1, signal2, ratio2m2, tAtt, 2, -2,
                        dt, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                        timeVec, matchrange, approximant, JLN ) == XLAL_FAILURE )
        {
              XLAL_ERROR( XLAL_EFUNC );
        }
        
        if ( thrStore22L != 0. && thrStore2m2L != 0. ) {
            if ( (*ratio22 - thr)*(*ratio22 - thr) + (*ratio2m2 - thr)*(*ratio2m2 - thr) < (thrStore22L - thr)*(thrStore22L - thr) + (thrStore2m2L - thr)*(thrStore2m2L - thr)  ) {
                thrStore22L = *ratio22;
                thrStore2m2L = *ratio2m2;
                tLBest = tAtt;
                if(debugPK)printf("tLBest is now %f %f %f %f\n", tLBest, *ratio22 ,*ratio2m2, (*ratio22 - thr)*(*ratio22 - thr) + (*ratio2m2 - thr)*(*ratio2m2 - thr));
            }
        }
        else {
            thrStore22L = *ratio22;
            thrStore2m2L = *ratio2m2;
            tLBest = tAtt;
        }

//        if (*ratio22 <= thr && *ratio2m2 <= thr){
//            pass = 1;
//        }


    }
    memset( signal1->data, 0, signal1->length * sizeof( signal1->data[0] ));
    memset( signal2->data, 0, signal2->length * sizeof( signal2->data[0] ));
    if(debugPK){
        if (pass == 1){
            printf("Going left, we have found better attachment point: new tAtt = %f, old = %f, ratios = %f, %f \n", tAtt, *tAttach, *ratio22, *ratio2m2);
        }else{
            printf("Going left did nto find the best attachment point\n");
        }
                
    }
    //REAL8 left_r22 = *ratio22;
    //REAL8 left_r2m2 = *ratio2m2;
    //REAL8 left_tAtt = tLBest;
    int pass_left = pass;
    int iBad = retLenHi - 1;

    pass = 0;
    tAtt = *tAttach;
    while(pass == 0 && (tAtt < tMax-0.5)){
//        printf("tAtt tMax %f %f\n", tAtt, tMax);
        tAtt = tAtt + 0.5;
        memset( signal1->data, 0, signal1->length * sizeof( signal1->data[0] ));
        memset( signal2->data, 0, signal2->length * sizeof( signal2->data[0] ));
        matchrange->data[0] = combSize < tAtt ? tAtt - combSize : 0;
        matchrange->data[1] = tAtt;
        matchrange->data[0] -= fmod( matchrange->data[0], dt/mTScaled );
        matchrange->data[1] -= fmod( matchrange->data[1], dt/mTScaled );
        for ( i = 0; i < retLenHi; i++ )
        {
            signal1->data[i] = creal(h22->data->data[i]);
            signal2->data[i] = cimag(h22->data->data[i]);
        }
        for ( i = 0; i < retLenHi-1; i++ )
        {
            dsignal1 = (signal1->data[i+1] - signal1->data[i]);
            dsignal2 = (signal2->data[i+1] - signal2->data[i]);
            hNorm2 = signal1->data[i]*signal1->data[i] + signal2->data[i]*signal2->data[i];
            omegaVec[i] = (-dsignal1*signal2->data[i] + dsignal2*signal1->data[i])/hNorm2;
        }
        for ( i = 0; i < retLenHi-2; i++ )
        {
            if (omegaVec[i]*omegaVec[i+1] < 0) {
                iBad = i;
                break;
            }
        }
        //for ( i = iBad; i < retLenHi; i++ )
        //{
        //    signal1->data[i] = 0.;
        //    signal2->data[i] = 0.;
        //}

        
        if (debugPK) printf("right 2,2 mode tAtt = %f     ", tAtt);
        if (matchrange->data[1] < iBad){
            if( XLALSimCheckRDattachment(signal1, signal2, ratio22, tAtt, 2, 2,
                            dt, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                            timeVec, matchrange, approximant, JLN ) == XLAL_FAILURE )
            {
                  XLAL_ERROR( XLAL_EFUNC );
            }

     
            memset( signal1->data, 0, signal1->length * sizeof( signal1->data[0] ));
            memset( signal2->data, 0, signal2->length * sizeof( signal2->data[0] ));
            for ( i = 0; i < retLenHi; i++ )
            {
                signal1->data[i] = creal(h2m2->data->data[i]);
                signal2->data[i] = cimag(h2m2->data->data[i]);
            }
       
            if (debugPK) printf("right 2,-2 mode tAtt = %f     ", tAtt); 
            if( XLALSimCheckRDattachment(signal1, signal2, ratio2m2, tAtt, 2, -2,
                            dt, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
                            timeVec, matchrange, approximant, JLN ) == XLAL_FAILURE )
            {
                  XLAL_ERROR( XLAL_EFUNC );
            }
        
            if ( thrStore22R != 0. && thrStore2m2R != 0. ) {
                if ( tRBest < timeVec->data[iBad] && (*ratio22 - thr)*(*ratio22 - thr) + (*ratio2m2 - thr)*(*ratio2m2 - thr) < (thrStore22R - thr)*(thrStore22R - thr) + (thrStore2m2R - thr)*(thrStore2m2R - thr)  ) {
                    thrStore22R = *ratio22;
                    thrStore2m2R = *ratio2m2;
                    tRBest = tAtt;
                    if(debugPK)printf("tRBest is now %f %f %f %f\n", tRBest, *ratio22 , *ratio2m2, (*ratio22 - thr)*(*ratio22 - thr) + (*ratio2m2 - thr)*(*ratio2m2 - thr));
                }
            }
            else {
                thrStore22R = *ratio22;
                thrStore2m2R = *ratio2m2;
                tRBest = tAtt;
            }
        }
        
        
//        if (*ratio22 <= thr && *ratio2m2 <= thr){
//            pass = 1;
//        }

    }
    if(debugPK){
        if (pass == 1){
            printf("Going right, we have found better attachment point: new tAtt = %f, old = %f, ratios = %f, %f \n", tAtt, *tAttach, *ratio22, *ratio2m2);
        }else{
            printf("Going right did nto find the best attachment point\n");
        }
                
    }
    int pass_right = pass;
    
    if (1==1 || (pass_right == 0 && pass_left == 0) ) {
        if ( debugPK ) {
            printf("Cannot go below required threshold on RD/insp amplitude\n");
        }
        if ( (thrStore22L - thr)*(thrStore22L - thr) + (thrStore2m2L - thr)*(thrStore2m2L - thr) < (thrStore22R - thr)*(thrStore22R - thr) + (thrStore2m2R - thr)*(thrStore2m2R - thr)) {
            *tAttach = tLBest;
            if ( debugPK ) {
                printf("tLBest %f\n", tLBest);
            }
        }
        else {
            *tAttach = tRBest;
            if ( debugPK ) {
                printf("tRBest %f\n", tRBest);
            }

        }
            return(2);
    }

    /*if( pass_right == 1 && pass_left == 0){
        *tAttach = tAtt;
        return(1);
    }
    if( pass_left == 1 && pass_right == 0){
       *tAttach =left_tAtt;
       *ratio22 = left_r22;
       *ratio2m2 = left_r2m2;
       return(1);
    }
    if (pass_left == 1 && pass_right == 1){
   
        // hard choice
        if (left_r22 <= 1.0 || left_r2m2 <= 1.0){
            if (*ratio22 <= 1.0 || *ratio2m2 <= 1.0){
                *tAttach =left_tAtt;
                *ratio22 = left_r22;
                *ratio2m2 = left_r2m2;
                return(1); // choose left
            }else{
                *tAttach = tAtt;
                return(1); // choose right
            }
        }else{
            if (*ratio22 <= 1.0 || *ratio2m2 <= 1.0){
                *tAttach =left_tAtt;
                *ratio22 = left_r22;
                *ratio2m2 = left_r2m2;
                return(1); // choose left
            }else{ // booth looks ok so far
                if (left_r22*left_r2m2 < (*ratio22)*(*ratio2m2)){
                    *tAttach =left_tAtt;
                    *ratio22 = left_r22;
                    *ratio2m2 = left_r2m2;
                    return(1); // choose left
                }else{
                    *tAttach = tAtt;
                    return(1); // choose right
                }

            }
        }
    }*/
    return(0);

} 


