/************************************************************************************/
/* The functions below are all associated with reading in source iniformation from  */
/* a source file and filling in the source structure.                               */ 
/*                                                                                  */
/*			           C. Messenger                                     */
/*                                                                                  */
/*                         BIRMINGHAM UNIVERISTY -  2005                            */
/************************************************************************************/

#include "ReadSourceFile_v1.h"

LALStatus status;

int ReadSource(char *sourcefile, char *sourcename, LIGOTimeGPS *obsstart, binarysource *source) 
{

  FILE *fpsource;
  INT4 foundsource;
  char line[1024];
  char sourceline[1024];
  char source_name_one[256];
  char source_name_two[256];
  INT4 ra_hr;
  INT4 ra_min;
  REAL8 ra_sec;
  REAL8 ra_err;
  INT4 dec_deg;
  INT4 dec_arcmin;
  REAL8 dec_arcsec;
  REAL8 dec_err;
  REAL8 freq_one;
  REAL8 freq_one_err;
  REAL8 freq_two;
  REAL8 freq_two_err;
  REAL8 period;
  REAL8 period_err;
  INT4 period_epoch;
  LIGOTimeGPS period_epochGPS;
  REAL8 sma;
  REAL8 sma_err;
  INT4 tperi;
  REAL8 tperi_err;
  REAL8 argp;
  REAL8 argp_err;
  REAL8 ecc;
  REAL8 ecc_err;
  REAL8 ra_rad;
  REAL8 ra_rad_min;
  REAL8 ra_rad_max;
  REAL8 dec_rad;
  REAL8 dec_rad_min;
  REAL8 dec_rad_max;
  REAL8 freq_one_min;
  REAL8 freq_one_max;
  REAL8 freq_two_min;
  REAL8 freq_two_max;
  REAL8 period_min;
  REAL8 period_max;
  REAL8 sma_min;
  REAL8 sma_max;
  LIGOTimeGPS tperi_dummy;
  LALTimeInterval interval;
  REAL8 period_epochdiff;
  INT4 n_period;
  REAL8 extra_err;
  REAL8 argp_min;
  REAL8 argp_max;
  REAL8 ecc_min;
  REAL8 ecc_max;
  INT4 sign;
  LIGOTimeGPS tperi_min;
  LIGOTimeGPS tperi_max;
  INT4 nband;
  

  /* this function simply reads the source file and extracts the information relating to */
  /* a particular source */

  /* opening the source file */
  fpsource=fopen(sourcefile,"r");
  if (fpsource==NULL) {
    fprintf(stderr,"Unable to open source file %s\n",sourcefile);
    return 1;
  }

  foundsource=0;
  /* loop over the sources */
  while ((fgets(line,1023,fpsource)!=NULL)&&(foundsource==0)) {

    /* read in the name fields */
    sscanf(line,"%s%s",source_name_one,source_name_two); 
    
    /* check if the source is in the list */
    if ((strcmp(sourcename,source_name_one)==0)||(strcmp(sourcename,source_name_two)==0)) {
      /* if it is in the list then record the line */ 
      foundsource=1;
      strcpy(sourceline,line);
    }
    
  }

  /* if we have not found the source in the file */
  if (foundsource==0) {
    fprintf(stderr,"ERROR : source (%s) not found in file %s\n",sourcename,sourcefile);
    exit(1);
  }

  /* read in the source variables */
  sscanf(sourceline,"%s%s %d%d%lf%lf %d%d%lf%lf %lf%lf%lf%lf %lf%lf%d %lf%lf %d%lf %lf%lf %lf%lf",
	 source_name_one,source_name_two,
	 &ra_hr,&ra_min,&ra_sec,&ra_err,
	 &dec_deg,&dec_arcmin,&dec_arcsec,&dec_err,
	 &freq_one,&freq_one_err,&freq_two,&freq_two_err,
	 &period,&period_err,&period_epoch,
	 &sma,&sma_err,
	 &tperi,&tperi_err,
	 &argp,&argp_err,
	 &ecc,&ecc_err); 

  /* do some error checking here */

  /* calculate some intermediate variables */

  /* could add proper motions later */
  /* first convert sky positions to radians */
  ra_rad=(LAL_TWOPI/24.0)*((REAL8)ra_hr+((REAL8)ra_min/60.0)+(ra_sec/3600.0));
  sign=(INT4)(abs(dec_deg)/dec_deg);
  dec_rad=(sign*LAL_TWOPI/360.0)*((REAL8)(abs(dec_deg))+((REAL8)dec_arcmin/60.0)+(dec_arcsec/3600.0));
  ra_rad_min=ra_rad-((LAL_TWOPI*3600.0*ra_err)/24.0);
  ra_rad_max=ra_rad+((LAL_TWOPI*3600.0*ra_err)/24.0);
  dec_rad_min=dec_rad-((LAL_TWOPI*3600.0*dec_err)/360.0);
  dec_rad_max=dec_rad+((LAL_TWOPI*3600.0*dec_err)/360.0);

  /* sort out the frequency ranges */
  freq_one_min=freq_one-freq_one_err;
  freq_one_max=freq_one+freq_one_err;
  freq_two_min=freq_two-freq_two_err;
  freq_two_max=freq_two+freq_two_err;

  /* sort out period ranges */
  period_min=period-period_err;
  period_max=period+period_err;

  /* sort out the sma ranges */
  sma_min=sma-sma_err;
  sma_max=sma+sma_err;

  /* in the future we need to sort out correct errors */
  /* sort out the tperi ranges */
  tperi_dummy.gpsSeconds=tperi;
  tperi_dummy.gpsNanoSeconds=0;
  period_epochGPS.gpsSeconds=period_epoch;
  period_epochGPS.gpsNanoSeconds=0;
  
  /* lets calculate the extra errors due to accumulating time (its a bit simple at present) */
  if (obsstart!=NULL) {
    LALDeltaGPS(&status,&interval,obsstart,&period_epochGPS);
    LALIntervalToFloat(&status,&period_epochdiff,&interval);
    n_period=period_epochdiff/period;
    extra_err=n_period*period_err;
    tperi_err=sqrt((tperi_err*tperi_err)+(extra_err*extra_err));
  }

  /* now add the errors to find range */
  LALFloatToInterval(&status,&interval,&tperi_err);
  LALDecrementGPS(&status,&tperi_min,&tperi_dummy,&interval);
  LALIncrementGPS(&status,&tperi_max,&tperi_dummy,&interval);

  /* sort out the argp ranges */
  argp_min=argp-argp_err;
  argp_max=argp+argp_err;

  /* sort out the ecc ranges */
  ecc_min=ecc-ecc_err;
  ecc_max=ecc+ecc_err;

  /* allocate memory to the structure */
  nband=1;
  if (freq_two!=0.0) nband=2;
  source->freq.f_min=(REAL8 *)LALMalloc(nband*sizeof(REAL8));
  source->freq.f_max=(REAL8 *)LALMalloc(nband*sizeof(REAL8));
  source->freq.f_err=(REAL8 *)LALMalloc(nband*sizeof(REAL8));
  
  /* fill in the source structure */
  strcpy(source->name,sourcename);
  /*printf("source name = %s\n",sourcename);*/
  source->skypos.ra=ra_rad;
  /*printf("ra_rad = %f\n",ra_rad);*/
  source->skypos.dec=dec_rad;
  /*printf("dec_rad = %f\n",dec_rad);*/
  source->skypos.ra_min=ra_rad_min;
  /*printf("ra_rad_min = %f\n",ra_rad_min);*/
  source->skypos.ra_max=ra_rad_max;
  /*printf("ra_rad_max = %f\n",ra_rad_max);*/
  source->skypos.dec_min=dec_rad_min;
  /*printf("dec_rad_min = %f\n",dec_rad_min);*/
  source->skypos.dec_max=dec_rad_max;
  /*printf("dec_rad_max = %f\n",dec_rad_max);*/
  source->skypos.ra_err=ra_err;
  /*printf("ra_err = %f\n",ra_err);*/
  source->skypos.dec_err=dec_err;
  /*printf("dec_err = %f\n",dec_err);*/
  source->freq.f_min[0]=freq_one_min;
  /*printf("freq_one_min = %f\n",freq_one_min);*/
  source->freq.f_max[0]=freq_one_max;
  /*printf("freq_one_max = %f\n",freq_one_max);*/
  source->freq.f_err[0]=freq_one_err;
  /*printf("freq_one_err = %f\n",freq_one_err);*/
  if (nband==2) {
    source->freq.f_min[1]=freq_two_min;
    /*printf("freq_two_min = %f\n",freq_two_min);*/
    source->freq.f_max[1]=freq_two_max;
    /*printf("freq_two_max = %f\n",freq_two_max);*/
    source->freq.f_err[1]=freq_two_err;
    /*printf("freq_two_err = %f\n",freq_two_err);*/
  }
  source->freq.nband=nband;
  /*printf("nband = %d\n",nband);*/
  source->orbit.period=period;
  /*printf("period = %f\n",period);*/
  source->orbit.period_min=period_min;
  /*printf("period_min = %f\n",period_min);*/
  source->orbit.period_max=period_max;
  /*printf("period_max = %f\n",period_max);*/
  source->orbit.period_err=period_err;
  /*printf("period_err = %f\n",period_err);*/
  source->orbit.sma=sma;
  /*printf("sma = %f\n",sma);*/
  source->orbit.sma_min=sma_min;
  /*printf("sma_min = %f\n",sma_min);*/
  source->orbit.sma_max=sma_max;
  /*printf("sma_max = %f\n",sma_max);*/
  source->orbit.sma_err=sma_err;
  /*printf("sma_err = %f\n",sma_err);*/
  source->orbit.tperi.gpsSeconds=tperi_dummy.gpsSeconds;
  /*printf("tperi_sec = %d\n",tperi_dummy.gpsSeconds);*/
  source->orbit.tperi.gpsNanoSeconds=tperi_dummy.gpsNanoSeconds;
  /*printf("tperi_nan = %d\n",tperi_dummy.gpsNanoSeconds);*/
  source->orbit.tperi_min.gpsSeconds=tperi_min.gpsSeconds;
  /*printf("tperi_min_sec = %d\n",tperi_min.gpsSeconds);*/
  source->orbit.tperi_min.gpsNanoSeconds=tperi_min.gpsNanoSeconds;
  /*printf("tperi_min_nan = %d\n",tperi_min.gpsNanoSeconds);*/
  source->orbit.tperi_max.gpsSeconds=tperi_max.gpsSeconds;
  /*printf("tperi_max_sec = %d\n",tperi_max.gpsSeconds);*/
  source->orbit.tperi_max.gpsNanoSeconds=tperi_max.gpsNanoSeconds;
  /*printf("tperi_max_nan = %d\n",tperi_max.gpsNanoSeconds);*/
  source->orbit.tperi_err=tperi_err;
  /*printf("tperi_err = %f\n",tperi_err);*/
  source->orbit.argp=argp;
  /*printf("argp = %f\n",argp);*/
  source->orbit.argp_min=argp_min;
  /*printf("argp_min = %f\n",argp_min);*/
  source->orbit.argp_max=argp_max;
  /*printf("argp_max = %f\n",argp_max);*/
  source->orbit.argp_err=argp_err;
  /*printf("argp_err = %f\n",argp_err);*/
  source->orbit.ecc=ecc;
  /*printf("ecc = %f\n",ecc);*/
  source->orbit.ecc_min=ecc_min;
  /*printf("ecc_min = %f\n",ecc_min);*/
  source->orbit.ecc_max=ecc_max;
  /*printf("ecc_max = %f\n",ecc_max);*/
  source->orbit.ecc_err=ecc_err;
  /*printf("ecc_err = %f\n",ecc_err);*/


  return 0;

}

/*********************************************************************************************/

