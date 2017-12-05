//
//  LALInferenceFVectorMultiBanding.c:
//  
//
//  Created by Serena Vinciguerra on 17/02/2015.
//
//

#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/LALInferenceMultibanding.h>


/** F(t) and T(f) for newtonian waveform */
static double LALInferenceTimeFrequencyRelation(double mc, double inPar, UINT4 flag_f);


static double LALInferenceTimeFrequencyRelation(double mc, double inPar, UINT4 flag_f)
{
    double outPar = 0.0;
    if (!flag_f){
        //outPar is time
        inPar=inPar/1.1; /* dividing the frequency by a safety factor 1.1*/
        outPar=-5.0*pow((8.0*LAL_PI*inPar),-8./3.)*pow(mc, -5./3.);
    }
    else{
        //outPar is frequency
        outPar=pow(-inPar/5.,-3./8.)*pow(mc,-5./8.)/(8.*LAL_PI);
    }
    return (outPar);
}

REAL8Sequence *LALInferenceMultibandFrequencies(int NBands, double f_min, double f_max, double deltaF0, double mc)
{
    
    if (0) printf("f_max %f", f_max);
    if (NBands==0){
        XLALPrintError("ERROR: ask for Multi Banding procedure with 0 bands!\n");
        return NULL;
    }

    double mc_sec=mc*LAL_MTSUN_SI;

    f_min=deltaF0*ceil(f_min/deltaF0);
    f_max=deltaF0*floor(f_max/deltaF0);
    

    double fi = f_min;
    int n_max = floor(- log2(deltaF0*2.1));
    float t_nmax= 1./(pow(2.,n_max)*deltaF0) - 2.1;
    float f_nmax= LALInferenceTimeFrequencyRelation(mc_sec,-t_nmax,1);
    
    int n_min = floor(- log2(deltaF0*(2.1 - LALInferenceTimeFrequencyRelation(mc_sec,f_min,0))));
    if (n_min < 0 ) n_min = 0;
    int max_NBands = 0;
    max_NBands = n_max+1 - n_min;
    if (NBands==-1)
        NBands = max_NBands;
        if (NBands > max_NBands){
        XLALPrintError(" WORNING required number of bands bigger than the maximum allowed, new number of bands %i \n",max_NBands);
        printf("WORNING required number of bands bigger than the maximum allowed, new number of bands %i \n",max_NBands);
        NBands = max_NBands;
    }
    
    REAL8Vector *F_inf = XLALCreateREAL8Vector(NBands);
    REAL8Vector *F_sup = XLALCreateREAL8Vector(NBands);
    REAL8Vector *deltaF = XLALCreateREAL8Vector(NBands);
    REAL8Vector *desired_df = XLALCreateREAL8Vector(NBands);
    REAL8Sequence *Frequencies=NULL;
    
    if (F_sup== NULL || F_inf==NULL || deltaF ==NULL) {
        XLALPrintError(" ERROR F_inf/sup/deltaF MBfile == NULL.\n");
        printf(" ERROR new/sup/deltaF (MB file)== NULL.\n");
        exit(1);
    }
    
    F_inf->data[0]=f_min;
    F_sup->data[NBands-1]=f_max;
    
    int nC = 1;
    if (NBands == 1) n_max = n_min;
    deltaF->data[NBands - 1] = (pow(2.,n_max)*deltaF0);
    F_inf->data[NBands - 1] = f_nmax;
    F_sup->data[NBands - 1] = f_max;
    fi = f_max;
    int n_i = n_max;
    
    for (int i=NBands - 1; i>=0; i--) {
       
      if (i< NBands - 1){
        n_i = n_i - 1;
        if (n_i<0) n_i = 0 ;
        deltaF->data[i] = (pow(2.,n_i)*deltaF0);
        if (i< NBands - 1) F_sup->data[i] = F_inf->data[i+1];
        double t_ni= 1./(pow(2.,n_i)*deltaF0) - 2.1;
        F_inf->data[i] = LALInferenceTimeFrequencyRelation(mc_sec,-t_ni,1);
        if(i == 0 && F_inf->data[i]<= f_min ) F_inf->data[i]= f_min;
        if (i== 0 && F_inf->data[i] > f_min ) {
            n_i = n_min;
            deltaF->data[i] = (pow(2.,n_i)*deltaF0);
            F_inf->data[i] = f_min;
        }
      }
      while((fi-deltaF->data[i])>=F_inf->data[i]){
            fi=fi-deltaF->data[i];
            nC++;
      }
    }
    
    UINT4 NFreq = nC;
    Frequencies = XLALCreateREAL8Sequence(NFreq);
    Frequencies->data[nC - 1]= f_max;
    int k= nC - 1;
    printf("PHASE INTERPOLATION ACTIVATED: %i bands used to define the %u frequencies at which compute the waveforms \n",NBands,NFreq);


    for (int i=NBands - 1; i>=0; i--) {
         while((Frequencies->data[k]-deltaF->data[i])>=F_inf->data[i]){
            k--;
            Frequencies->data[k] = Frequencies->data[k+1]-deltaF->data[i];
        }
    }

    XLALDestroyREAL8Vector(F_inf);
    XLALDestroyREAL8Vector(F_sup);
    XLALDestroyREAL8Vector(deltaF);
    XLALDestroyREAL8Vector(desired_df);
    return(Frequencies);
    
}
