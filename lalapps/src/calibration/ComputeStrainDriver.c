/*
*  Copyright (C) 2007 Xavier Siemens
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

/*********************************************************************************/
/*                         Driver for calibration function                       */
/*                                                                               */
/*         X. Siemens ideas & implementation, J. Burguet-Castell bugs            */
/*                                                                               */
/*                               UWM - August 2004                               */
/*********************************************************************************/

#include <config.h>
#if !defined HAVE_LIBGSL || !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no gsl or no lal frame library support.\n", stderr);return 1;}
#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <glob.h>
#include <errno.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include <pwd.h>
#include <time.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/LALCache.h>
#include <lal/FrameStream.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/ComputeDataQualityVector.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>
#include <lal/ReadFiltersFile.h>
#include <lal/TimeSeries.h>
#include <lal/LALVersion.h>
#include <lal/LALFrameIO.h>
#include <lal/LALDetectors.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALString.h>
#include <lal/FrequencySeries.h>
#include <lal/LALVCSInfo.h>
#include <lalapps.h>
#include <LALAppsVCSInfo.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)



#define CVS_ID "$Id$"
#define CVS_HEADER "$Header$"
#define CVS_AUTHOR "$Author$"
#define CVS_NAME "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define PROGRAM_NAME "ComputeStrainDriver"

/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 To;                /* factors integration time */
  INT4 GPSStart;           /* Start GPS time for the segment to be calibrated */
  INT4 GPSEnd;             /* End GPS time for the segment to be calibrated */
  INT4 testsensing;
  INT4 testactuation;
  char *FrCacheFile;       /* Frame cache file for corresponding time */
  char *ifo;               /* interferometer name (H1, H2, L1) */
  char *filterfile;        /* file with filter coefficients */
  char *frametype;
  char *strainchannel;
  char *datadirL1;         /* (deprecated) */
  char *datadirL2;
  INT4 checkfilename;
} CommandLineArgs;

/***************************************************************************/

/* GLOBAL VARIABLES */

StrainIn InputData;
StrainOut OutputData;
INT4 duration;
LIGOTimeGPS gpsStartepoch;

INT4TimeSeries OutputDQ;  /* data quality */

static LALStatus status;
LALCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;
char sv_cname[] = "Xn:IFO-SV_STATE_VECTOR",                      /* channel names */
    lax_cname[] = "Xn:LSC-LA_PTRX_NORM", lay_cname[] = "Xn:LSC-LA_PTRY_NORM",
    asq_cname[] = "Xn:LSC-DARM_ERR",  /* temporary hack: set name of (unused) asq to darm_err */
    dctrl_cname[] = "Xn:LSC-DARM_CTRL",
    derr_cname[] = "Xn:LSC-DARM_ERR", exc_cname[] = "Xn:LSC-DARM_CTRL_EXC_DAQ";


/***************************************************************************/
/* to avoid a warning */
int gethostname(char *name, size_t len);
/* int getdomainname(char *name, size_t len); */


/* FUNCTION PROTOTYPES (defined static so they cannot be used outside) */
/* Reads the command line */
static int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads T seconds of AS_Q, EXC, and DARM_CTRL data */
static int ReadData(struct CommandLineArgsTag CLA);

/* Writes frame file */
static int WriteFrame(int argc,char *argv[],struct CommandLineArgsTag CLA);

/* Frees the memory */
static int FreeMem(void);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

  /* Check whether files exist before proceeding */
  if (CommandLineArgs.checkfilename)
    {
      /* char fname[FILENAME_MAX]; */  /* *** DEPRECATED *** */
      char fname2[FILENAME_MAX];
      char site;
      INT4 t0, dt;

      site = CommandLineArgs.ifo[0];

      t0 = CommandLineArgs.GPSStart+InputData.wings;
      dt = CommandLineArgs.GPSEnd-CommandLineArgs.GPSStart-2*InputData.wings;

      /* snprintf( fname, sizeof( fname ), "%s/%c-%s%s-%d-%d.gwf", CommandLineArgs.datadirL1,site, CommandLineArgs.frametype,"_L1", t0, dt ); */  /* *** DEPRECATED *** */
      snprintf( fname2, sizeof( fname2 ), "%s/%c-%s%s-%d-%d.gwf", CommandLineArgs.datadirL2,site, CommandLineArgs.frametype,"_L2", t0, dt );

    /*   if (access(fname, F_OK) == 0) { */
	/* fprintf(stdout, "Frame file %s exists. Exiting.\n",fname); */
	/* return 0; */
    /*   } */  /* *** DEPRECATED *** */
      if (access(fname2, F_OK) == 0) {
	fprintf(stdout, "Frame file %s exists. Exiting.\n",fname2);
	return 0;
      }
    }

  if (XLALReadFiltersFile(CommandLineArgs.filterfile, &InputData)) return 3;

  if (ReadData(CommandLineArgs)) return 2;

  LALComputeStrain(&status, &OutputData, &InputData);
  TESTSTATUS( &status );

  XLALComputeDQ(InputData.StateVector.data->data, InputData.StateVector.data->length/OutputDQ.data->length,
                InputData.LAX.data->data, InputData.LAY.data->data, InputData.LAX.data->length/OutputDQ.data->length,
                OutputData.alphabeta.data->data, OutputData.alphabeta.data->length/OutputDQ.data->length,
                0, 0, InputData.wings,
                0,   /* how can I actually know if it is missing in the DMT or not?? */
                OutputDQ.data->data, OutputDQ.data->length);

  if (WriteFrame(argc,argv,CommandLineArgs)) return 4;

  if(FreeMem()) return 5;

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/

/*  FUNCTIONS */

/*******************************************************************************/

int WriteFrame(int argc,char *argv[],struct CommandLineArgsTag CLA)
{
  /* This is mostly ripped off some code Jolien sent me a while back,
     and calibration frame writing code */

  FrFile *frfile;
  FrameH *frame;
  char fname[FILENAME_MAX];
  char tmpfname[FILENAME_MAX];
  char site;
  INT4 t0;
  INT4 dt;
  INT4 FrDuration;
  int detectorFlags;
  char hostnameanduser[4096];
  char hostname[1024];
  char domainname[1024];
  char allargs[16384];
  char lalappsconfargs[16384];
  char lalconfargs[16384];
  char headerinfo[16384];
  int i;
  REAL4TimeSeries *alphare = NULL;
  REAL4TimeSeries *alphaim = NULL;
  REAL4TimeSeries *gammare = NULL;
  REAL4TimeSeries *gammaim = NULL;
  char alphareName[] = "Xn:CAL-CAV_FAC";
  char gammareName[] = "Xn:CAL-OLOOP_FAC";
  char alphaimName[] = "Xn:CAL-CAV_FAC_Im";
  char gammaimName[] = "Xn:CAL-OLOOP_FAC_Im";
  char dqName[] = "Xn:LSC-DATA_QUALITY_VECTOR";
  char freqInfo[] = "Frequency validity range: 40Hz-5kHz.";

  char *cnames[] = { alphareName, gammareName, alphaimName, gammaimName, dqName };

  for (i = 0; i < 5; i++)
      memcpy(cnames[i], CLA.ifo, 2);  /* set the proper name of the channels */

  /*re-size h(t) and data time series*/
  if(!XLALResizeREAL8TimeSeries(&(OutputData.h), (int)(InputData.wings/OutputData.h.deltaT),
			   OutputData.h.data->length-2*(UINT4)(InputData.wings/OutputData.h.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  strncpy( OutputData.h.name,  CLA.strainchannel, sizeof( OutputData.h.name  ) );

  /* resize alpha and gamma time series */
  if(!XLALResizeCOMPLEX16TimeSeries(&(OutputData.alpha), (int)(InputData.wings/OutputData.alpha.deltaT),
			   OutputData.alpha.data->length-2*(UINT4)(InputData.wings/OutputData.alpha.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  if(!XLALResizeCOMPLEX16TimeSeries(&(OutputData.alphabeta), (int)(InputData.wings/OutputData.alphabeta.deltaT),
			   OutputData.alphabeta.data->length-2*(UINT4)(InputData.wings/OutputData.alphabeta.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);

  /* Resize State Vector and Data Quality time series */
  if(!XLALResizeREAL4TimeSeries(&(InputData.StateVector), (int)(InputData.wings/InputData.StateVector.deltaT),
               InputData.StateVector.data->length-2*(UINT4)(InputData.wings/InputData.StateVector.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  if(!XLALResizeINT4TimeSeries(&(OutputDQ), (int)(InputData.wings/OutputDQ.deltaT),
			   OutputDQ.data->length-2*(UINT4)(InputData.wings/OutputDQ.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  strncpy(OutputDQ.name, dqName, sizeof(OutputDQ.name));   /* also set the name of the channel */

  /* Resize DARM_CTRL, DARM_ERR, EXC and AS_Q*/
  if(!XLALResizeREAL4TimeSeries(&(InputData.DARM), (int)(InputData.wings/InputData.DARM.deltaT),
			   InputData.DARM.data->length-2*(UINT4)(InputData.wings/InputData.DARM.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  if(!XLALResizeREAL4TimeSeries(&(InputData.DARM_ERR), (int)(InputData.wings/InputData.DARM_ERR.deltaT),
			   InputData.DARM_ERR.data->length-2*(UINT4)(InputData.wings/InputData.DARM_ERR.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  if(!XLALResizeREAL4TimeSeries(&(InputData.EXC), (int)(InputData.wings/InputData.EXC.deltaT),
			   InputData.EXC.data->length-2*(UINT4)(InputData.wings/InputData.EXC.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);
  if(!XLALResizeREAL4TimeSeries(&(InputData.AS_Q), (int)(InputData.wings/InputData.AS_Q.deltaT),
			   InputData.AS_Q.data->length-2*(UINT4)(InputData.wings/InputData.AS_Q.deltaT)))
    XLAL_ERROR(XLAL_EFUNC);

  /* based on IFO name, choose the correct detector */
  if ( 0 == strcmp(CLA.ifo, "H2") )
    detectorFlags = LAL_LHO_2K_DETECTOR_BIT;
  else if ( 0 == strcmp( CLA.ifo, "H1") )
    detectorFlags = LAL_LHO_4K_DETECTOR_BIT;
  else if ( 0 == strcmp( CLA.ifo, "L1") )
    detectorFlags = LAL_LLO_4K_DETECTOR_BIT;
  else
    return 1;  /* Error: not a recognized name */
  site = CLA.ifo[0];

  /* based on series metadata, generate standard filename */
  FrDuration = OutputData.h.deltaT * OutputData.h.data->length;
  t0 = OutputData.h.epoch.gpsSeconds;
  dt = ceil( XLALGPSGetREAL8( &OutputData.h.epoch ) + FrDuration ) - t0;
  if ( t0 < 0 || dt < 1 )
    return 1;  /* Error: invalid time or FrDuration */
  snprintf( fname, sizeof( fname ), "%s/%c-%s%s-%d-%d.gwf", CLA.datadirL1,site, CLA.frametype,"_L1", t0, dt );
  snprintf( tmpfname, sizeof( tmpfname ), "%s.tmp", fname );

  /* Harwired numbers in call to XLALFrameNew: */
  /* Run number is set to 0, this is an annoying thing to have to
     change all the time and if you can't tell the run number from
     the GPS time you have a problem */
  /* number of frames in frame to 1 */
  frame = XLALFrameNew( &OutputData.h.epoch , FrDuration, "LIGO", 0, 1, detectorFlags );

  /* Here's where I need to add a bunch of things */
  /* Add cvs header */
  snprintf( headerinfo, sizeof( headerinfo), "Code header info: %s",CVS_HEADER);
  FrHistoryAdd( frame, headerinfo);

  /* Add lalapps info */
  snprintf( lalappsconfargs, sizeof( lalappsconfargs), "LALApps Info:\n                          LALApps Version: %s\n                          Git Tag: %s\n                          Git ID: %s\n                          Configure Date: %s\n                          Configure Arguments: %s",
	       LALAPPS_VERSION , lalAppsVCSInfo.vcsTag, lalAppsVCSInfo.vcsId, LALAPPS_CONFIGURE_DATE , LALAPPS_CONFIGURE_ARGS );
  FrHistoryAdd( frame, lalappsconfargs);

  /* Add lal info */
  snprintf( lalconfargs, sizeof( lalconfargs), "LAL Info:\n                          LAL Version: %s\n                          Git Tag: %s\n                          Git ID: %s\n                          Configure Date: %s\n                          Configure Arguments: %s",
	       LAL_VERSION , lalHeaderVCSInfo.vcsTag, lalHeaderVCSInfo.vcsId, LAL_CONFIGURE_DATE , LAL_CONFIGURE_ARGS );
  FrHistoryAdd( frame, lalconfargs);

  /* Create string with all command line arguments and add it to history */
  strcat(allargs, "Command line run: ");

  for(i = 0; i < argc; i++)
    {
      strcat(allargs,argv[i]);
      strcat(allargs, " ");
    }
  FrHistoryAdd( frame, allargs);

  /* hostname and user */
  gethostname(hostname,sizeof(hostname));
  if ( getdomainname(domainname,sizeof(domainname)) == -1 )
    XLALPrintError ("\ngetdomainname() failed!\n");
  snprintf( hostnameanduser, sizeof( hostnameanduser), "Made by user: %s. Made on machine: %s.%s",getlogin(),hostname,domainname);
  FrHistoryAdd( frame, hostnameanduser);

  /* Frequency range of validity (FIXME: This should be updated regularly somehow) */
  FrHistoryAdd( frame, freqInfo);

  /* Filters file checksum and cvs info (first 2 lines in filters file) */
  {
    char buffer[1024];
    snprintf(buffer, sizeof buffer, "Filters file checksum and header: %s\n%s",
             InputData.filter_chksum, InputData.filter_vc_info);
    FrHistoryAdd(frame, buffer);
  }

  /* Add in the h(t) data */
  XLALFrameAddREAL8TimeSeriesProcData( frame, &OutputData.h);

  /* Add in the state vector data */
  XLALFrameAddREAL4TimeSeriesProcData( frame, &InputData.StateVector);

  /* Add in the data quality data */
  XLALFrameAddINT4TimeSeriesProcData( frame, &OutputDQ);

  /* Add in the factors data */
  alphare = XLALCreateREAL4TimeSeries( alphareName, &OutputData.alpha.epoch, 0.0, OutputData.alpha.deltaT,
				     &lalDimensionlessUnit,  OutputData.alpha.data->length);
  gammare = XLALCreateREAL4TimeSeries( gammareName, &OutputData.alphabeta.epoch, 0.0, OutputData.alphabeta.deltaT,
				     &lalDimensionlessUnit,  OutputData.alphabeta.data->length);
  alphaim = XLALCreateREAL4TimeSeries( alphaimName, &OutputData.alpha.epoch, 0.0, OutputData.alpha.deltaT,
				       &lalDimensionlessUnit,  OutputData.alpha.data->length);
  gammaim = XLALCreateREAL4TimeSeries( gammaimName, &OutputData.alphabeta.epoch, 0.0, OutputData.alphabeta.deltaT,
				       &lalDimensionlessUnit,  OutputData.alphabeta.data->length);

  for (i=0; i < (int)OutputData.alpha.data->length; i++)
    {
      alphare->data->data[i]=creal(OutputData.alpha.data->data[i]);
      alphaim->data->data[i]=cimag(OutputData.alpha.data->data[i]);
      gammare->data->data[i]=creal(OutputData.alphabeta.data->data[i]);
      gammaim->data->data[i]=cimag(OutputData.alphabeta.data->data[i]);
    }

  /* Deprecating Stat Data Calibration...
  XLALFrameAddCalFac( frame, alphare, atoi(&CLA.frametype[9]) );
  XLALFrameAddCalFac( frame, alphaim, atoi(&CLA.frametype[9]) );
  XLALFrameAddCalFac( frame, gammare, atoi(&CLA.frametype[9]));
  XLALFrameAddCalFac( frame, gammaim, atoi(&CLA.frametype[9]));
  ... add it as Proc Data instead: */
  XLALFrameAddREAL4TimeSeriesProcData( frame, alphare );
  XLALFrameAddREAL4TimeSeriesProcData( frame, alphaim );
  XLALFrameAddREAL4TimeSeriesProcData( frame, gammare );
  XLALFrameAddREAL4TimeSeriesProcData( frame, gammaim );


  XLALDestroyREAL4TimeSeries( gammare );
  XLALDestroyREAL4TimeSeries( alphare );
  XLALDestroyREAL4TimeSeries( gammaim );
  XLALDestroyREAL4TimeSeries( alphaim );


  /* write level 2 frame */
  {
    char fname2[FILENAME_MAX];
    char tmpfname2[FILENAME_MAX];

    snprintf( fname2, sizeof( fname2 ), "%s/%c-%s%s-%d-%d.gwf", CLA.datadirL2,site, CLA.frametype,"_L2", t0, dt );
    snprintf( tmpfname2, sizeof( tmpfname2 ), "%s.tmp", fname2 );

    /* write first to tmpfile then rename it */
    frfile = FrFileONew( tmpfname2, -1); /* 1 = GZIP */
    if ( ! frfile )
      return 1;  /* Error: could not open frame file */

    FrameWrite( frame, frfile );
    FrFileOEnd( frfile );
    /* now rename */
    if ( rename( tmpfname2, fname2 ) < 0 )
      return 1; /* Error: system error */
  }

  return 0;

}

/*******************************************************************************/


int ReadData(struct CommandLineArgsTag CLA)
{

FrPos pos1;

static FrChanIn chanin_darm;
static FrChanIn chanin_darmerr;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;
static FrChanIn chanin_sv;   /* state vector */
static FrChanIn chanin_lax;  /* light in x-arm */
static FrChanIn chanin_lay;  /* light in y-arm */

  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_darmerr.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;
  chanin_sv.type  = ADCDataChannel;
  chanin_lax.type  = ADCDataChannel;
  chanin_lay.type  = ADCDataChannel;

  chanin_asq.name  = asq_cname;
  chanin_darm.name = dctrl_cname;
  chanin_darmerr.name = derr_cname;
  chanin_exc.name  = exc_cname;
  chanin_sv.name  = sv_cname;
  chanin_lax.name  = lax_cname;
  chanin_lay.name  = lay_cname;

  /* create Frame cache, open frame stream and delete frame cache */
  framecache = XLALCacheImport(CommandLineArgs.FrCacheFile);
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  XLALDestroyCache(framecache);

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&gpsStartepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.AS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM,&chanin_darm,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.EXC,&chanin_exc,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM_ERR,&chanin_darmerr,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.StateVector,&chanin_sv,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.LAX,&chanin_lax,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.LAY,&chanin_lay,framestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALSCreateVector(&status,&InputData.AS_Q.data,(UINT4)(duration/InputData.AS_Q.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.DARM.data,(UINT4)(duration/InputData.DARM.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.DARM_ERR.data,(UINT4)(duration/InputData.DARM_ERR.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.EXC.data,(UINT4)(duration/InputData.EXC.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.StateVector.data,(UINT4)(duration/InputData.StateVector.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.LAX.data,(UINT4)(duration/InputData.LAX.deltaT +0.5));
  TESTSTATUS( &status );
  LALSCreateVector(&status,&InputData.LAY.data,(UINT4)(duration/InputData.LAX.deltaT +0.5));
  TESTSTATUS( &status );

  /* Read in the data */
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.AS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM,&chanin_darm,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.DARM_ERR,&chanin_darmerr,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.EXC,&chanin_exc,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.StateVector,&chanin_sv,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.LAX,&chanin_lax,framestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&InputData.LAY,&chanin_lay,framestream);
  TESTSTATUS( &status );

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  /* Set the rest of the input variables */
  InputData.To=CLA.To;

  /* check input data epoch agrees with command line arguments */
  if ( InputData.DARM_ERR.epoch.gpsSeconds != CLA.GPSStart )
    {
      fprintf(stderr,"GPS start time of data (%d) does not agree with requested start time (%d). Exiting.",
	      InputData.DARM_ERR.epoch.gpsSeconds, CLA.GPSStart);
      return 1;
    }

  /* Allocate output data */
  OutputData.h.epoch=InputData.DARM_ERR.epoch;
  OutputData.h.deltaT=InputData.DARM_ERR.deltaT;
  LALDCreateVector(&status,&OutputData.h.data,(UINT4)(duration/OutputData.h.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.hC.epoch=InputData.DARM_ERR.epoch;
  OutputData.hC.deltaT=InputData.DARM_ERR.deltaT;
  LALDCreateVector(&status,&OutputData.hC.data,(UINT4)(duration/OutputData.hC.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.hR.epoch=InputData.DARM_ERR.epoch;
  OutputData.hR.deltaT=InputData.DARM_ERR.deltaT;
  LALDCreateVector(&status,&OutputData.hR.data,(UINT4)(duration/OutputData.hR.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.alpha.epoch=InputData.DARM_ERR.epoch;
  OutputData.alpha.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.alpha.data,(UINT4)(duration/OutputData.alpha.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.alphabeta.epoch=InputData.DARM_ERR.epoch;
  OutputData.alphabeta.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.alphabeta.data,(UINT4)(duration/OutputData.alphabeta.deltaT +0.5));
  TESTSTATUS( &status );

  OutputData.beta.epoch=InputData.DARM_ERR.epoch;
  OutputData.beta.deltaT=CLA.To;
  LALZCreateVector(&status,&OutputData.beta.data,(UINT4)(duration/OutputData.beta.deltaT +0.5));
  TESTSTATUS( &status );

  OutputDQ.epoch=InputData.StateVector.epoch;
  OutputDQ.deltaT=1;  /* Data Quality channel written at 1 Hz */
  LALI4CreateVector(&status,&OutputDQ.data,(UINT4)(duration/OutputDQ.deltaT +0.5));
  TESTSTATUS( &status );

  return 0;
}



/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 errflg=0;
  struct option long_options[] = {
    {"factors-time",        required_argument, NULL,  't'},
    {"filters-file",        required_argument, NULL,  'F'},
    {"frame-cache",         required_argument, NULL,  'C'},
    {"ifo",                 required_argument, NULL,  'i'},
    {"gps-start-time",      required_argument, NULL,  's'},
    {"gps-end-time",        required_argument, NULL,  'e'},
    {"wings",               required_argument, NULL,  'o'},
    {"test-sensing",        no_argument,       NULL,  'r'},
    {"test-actuation",      no_argument,       NULL,  'c'},
    {"delta",               no_argument,       NULL,  'd'},
    {"no-factors",          no_argument,       NULL,  'u'},
    {"td-fir",              no_argument,       NULL,  'x'},
    {"require-histories",   required_argument, NULL,  'H'},
    {"frame-type",          required_argument, NULL,  'T'},
    {"strain-channel",      required_argument, NULL,  'S'},
    {"data-dirL1",          required_argument, NULL,  'z'},
    {"data-dirL2",          required_argument, NULL,  'p'},
    {"check-file-exists",   no_argument,       NULL,  'v'},
    {"darm-err-only",       no_argument,       NULL,  'w'},
    {"gamma-fudge-factor",  required_argument, NULL,  'y'},
    {"help",                no_argument,       NULL,  'h'},
    {0, 0, 0, 0}
  };
  char args[] = "hrcdux:C:F:s:e:i:t:o:H:T:S:z:v:wy:p:";

  /* Initialize default values */
  CLA->To=0.0;
  CLA->FrCacheFile=NULL;
  CLA->filterfile=NULL;
  CLA->ifo=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->testsensing=0;
  CLA->testactuation=0;
  CLA->frametype=NULL;
  CLA->strainchannel=NULL;
  CLA->datadirL1=NULL;
  CLA->datadirL2=NULL;
  CLA->checkfilename=0;

  InputData.delta=0;
  InputData.wings=0;
  InputData.usefactors=1;
  InputData.fftconv=1;
  InputData.outalphas=0;
  InputData.darmctrl=1;
  InputData.gamma_fudgefactor=1.0;

  /* Scan through list of command line arguments */
  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
    case 't':
      /* factors integration time */
      CLA->To=atof(optarg);
      break;
    case 'C':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'F':
      /* name of filter cache file */
      CLA->filterfile=optarg;
      break;
    case 'i':
      /* name of interferometer */
      if (strcmp(optarg, "H1") != 0 && strcmp(optarg, "H2") != 0 &&
          strcmp(optarg, "L1") != 0) {
        fprintf(stderr, "Bad ifo: %s   (must be H1, H2 or L1)\n", optarg);
        exit(1);
      }
      CLA->ifo=optarg;
      {
        int i;
        char *cnames[] = { sv_cname, lax_cname, lay_cname, asq_cname,
                           dctrl_cname, derr_cname, exc_cname };
        for (i = 0; i < 7; i++)
          memcpy(cnames[i], CLA->ifo, 2);  /* set channel names appropiately */
      }
      break;
    case 's':
      /* GPS start */
      CLA->GPSStart=atof(optarg);
      break;
    case 'e':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'd':
      /*  use unit impulse*/
      InputData.delta=1;
      break;
    case 'r':
      /*  output residual signal */
      CLA->testsensing=1;
      break;
    case 'c':
      /* output control signal */
      CLA->testactuation=1;
      break;
    case 'o':
      /*  wing size (at each end have this many extra seconds) */
      InputData.wings=atoi(optarg);
      break;
    case 'u':
      /* don't use calibration factors in teh strain computation */
      InputData.usefactors=0;
      break;
    case 'x':
      /* use time domain FIR filtering instead of FFT convolution */
      InputData.fftconv=0;
      break;
    case 'T':
      CLA->frametype=optarg;
      break;
    case 'S':
      CLA->strainchannel=optarg;
      break;
    case 'z':
      CLA->datadirL1=optarg;
      break;
    case 'p':
      CLA->datadirL2=optarg;
      break;
    case 'v':
      CLA->checkfilename=1;
      break;
    case 'w':
      /* don't use calibration factors in the strain computation */
      InputData.darmctrl=0;
      break;
    case 'y':
      InputData.gamma_fudgefactor=atof(optarg);
      break;
    case 'h':
      /* print usage/help message */
      printf("Arguments are:\n"
"\t-t, --factors-time    FLOAT     Factors integration time in seconds.\n"
"\t-s, --gps-start-time  INT       GPS start time.\n"
"\t-e, --gps-end-time    INT       GPS end time.\n"
"\t-F, --filters-file    STRING    Name of file containing filters and histories.\n"
"\t-C, --frame-cache     STRING    Name of frame cache file.\n"
"\t-i, --ifo             STRING    Name of the interferometer (H1, H2, L1).\n"
"\t-o, --wings           INTEGER   Size of wings in seconds.\n"
"\t-r, --test-sensing    FLAG      Output residual strain only.\n"
"\t-c, --test-actuation  FLAG      Output control strain only.\n"
"\t-d, --delta           FLAG      Use unit impulse.\n"
"\t-u, --no-factors      FLAG      Do not use factors in strain computation.\n"
"\t-x, --td-fir          FLAG      Use time-domain FIR filtering (default is FFT convolution).\n"
"\t-y, --gamma-fudge-factor FLOAT  Fudge factor used to adjust factor values. Gamma is divided by that value.\n"
"\t-T, --frame-type      STRING    Frame type to be written (eg, H1_RDS_C01_LX).\n"
"\t-S, --strain-channel  STRING    Strain channel name in frame (eg, H1:LSC-STRAIN).\n"
"\t-z, --data-dirL1      STRING    Ouput L1 frame to this directory (eg, /tmp/S4/H1/) [DEPRECATED].\n"
"\t-p, --data-dirL2      STRING    Ouput L2 frame to this directory (eg, /tmp/S4/H1/).\n"
"\t-v, --check-file-exists FLAG    Checks frame files exist and if they do it exits gracefully.\n"
"\t-w, --darm-err-only   FLAG      Do darm_err only calibration. Default is to use darm_err and darm_ctrl. For first epoch of S5.\n"
"\t-h, --help            FLAG      This message.\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  }

  if(CLA->To == 0.0)
    {
      fprintf(stderr,"No integration time for the factors specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->filterfile == NULL)
    {
      fprintf(stderr,"No filter file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
   if(CLA->ifo == NULL)
    {
      fprintf(stderr,"No ifo specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
   if(CLA->frametype == NULL)
    {
      fprintf(stderr,"No frame type specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
   if(CLA->strainchannel == NULL)
    {
      fprintf(stderr,"No strain channel specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
   if(CLA->datadirL1 != NULL)
    {
      fprintf(stdout,"Warning: L1 frame directory specified, but Level 1 "
              "frames are not going to be produced anymore.\n");
    }
   if(CLA->datadirL2 == NULL)
    {
      fprintf(stderr,"No L2 frame directory specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
   if ( (InputData.wings < 2) || (InputData.wings%2 != 0) )
     {
       fprintf(stderr,"Overlap %d is either not >=2, or not exactly divisible by two. Exiting.",InputData.wings);
       return 1;
     }

   /* Set some global variables */
   duration = CLA->GPSEnd - CLA->GPSStart;
   gpsStartepoch.gpsSeconds = CLA->GPSStart;
   gpsStartepoch.gpsNanoSeconds = 0;

  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{

  /* Free input data */
  LALSDestroyVector(&status,&InputData.AS_Q.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.DARM.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.DARM_ERR.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.EXC.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.StateVector.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.LAX.data);
  TESTSTATUS( &status );
  LALSDestroyVector(&status,&InputData.LAY.data);
  TESTSTATUS( &status );

  /* Free filters */
  LALDDestroyVector(&status,&InputData.Cinv->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.Cinv->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.Cinv->history);
  TESTSTATUS( &status );
  LALFree(InputData.Cinv);

  LALDDestroyVector(&status,&InputData.A->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.A->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.A->history);
  TESTSTATUS( &status );
  LALFree(InputData.A);

  LALDDestroyVector(&status,&InputData.AW->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.AW->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.AW->history);
  TESTSTATUS( &status );
  LALFree(InputData.AW);

  LALDDestroyVector(&status,&InputData.D->directCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.D->recursCoef);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&InputData.D->history);
  TESTSTATUS( &status );
  LALFree(InputData.D);

  /* Free output data */
  LALDDestroyVector(&status,&OutputData.h.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&OutputData.hC.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&OutputData.hR.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.alpha.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.beta.data);
  TESTSTATUS( &status );
  LALZDestroyVector(&status,&OutputData.alphabeta.data);
  TESTSTATUS( &status );
  LALI4DestroyVector(&status,&OutputDQ.data);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();

  return 0;
}

/*******************************************************************************/
#endif
