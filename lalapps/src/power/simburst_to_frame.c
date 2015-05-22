/*
 * simburst_to_frame.c
 * 
 * Copyright 2013 Salvatore Vitale <salvatore.vitale@ligo.org>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>
#include <lal/Units.h>
#include <lal/LALFrStream.h>
#include <lal/LALFrameIO.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

struct options {
	char *output;
	double time_step;
    char *simburst_file;
	char **ifonames;
    char **channames;
	UINT4 nIFO;
	double mdc_gps_start;
    int pad;
	int mdc_duration;
};


static struct options options_defaults(void)
{
	struct options defaults;

	defaults.output = NULL;
	defaults.ifonames=NULL;
        defaults.channames=NULL;
	defaults.nIFO=-1;
	defaults.mdc_gps_start=-1;
	defaults.mdc_duration=-1;
	defaults.simburst_file=NULL;
        defaults.pad=-1;
        defaults.time_step=0.0;
    return defaults;
}

void ParseCharacterOptionString(char *input, char **strings[], UINT4 *n);
static void write_log(SimBurst **injs, TimeSlide * time_slide_table_head,struct options *options,char* fname);


static void print_usage(void)
{
	fprintf(stderr, 
" lalapps_simburst_to_mdc --simburst-file simburst.xml --ifos [IFO1,IFO2,...IFON] [options]\n" \
"\n" \
"Description:\n"\
"Takes a simburst table and a list of ifos and create a frame file  with the\n"\
" signals contained in the burst table. There will be a channel for each IFO. \n"\
"Writes down a log file using cwb style\n"\
"\n"\
"By default the frame will start 10 seconds before the first injection of the\n"\
" xml table and end 10 seconds after the last injection. \n"\
"Each channel will be named as the corresponding IFO:Strain (e.g. H1:Strain) .\n"\
"\n" \
"Options:\n" \
"\n" \
"[--gps-start seconds]: Override the default starting time of the frame file. Note that this may lead to a frame with *fewer* signals than the original xml.\n" \
"[--duration seconds]: Override the default duration of frame file. Note that this may lead to a frame file which contain *fewer* signals than the original xml \n" \
"[--channels [A,B,..,N]: Override the default suffix of the channels to be IFO1:A, IFO2:B, etc. \n" \
"[--pad seconds]: Override the default padding before the first and after the last signals\n"\
);
}

static struct options parse_command_line(int *argc, char **argv[])
{
	struct options options = options_defaults();
	int c;
	int option_index;
    UINT4 nchannels=0;
	struct option long_options[] = {
    {"simburst-file",required_argument,NULL,'0'},
		{"help", no_argument, NULL, 'C'},
		{"gps-start", required_argument,0,1728},
		{"duration", required_argument,0,1729},
    {"channels",required_argument,NULL,1732},
		{"pad", required_argument, NULL, 1731},
    {"ifos",required_argument,NULL,1730},
		{NULL, 0, NULL, 0}
	};
	do switch(c = getopt_long(*argc, *argv, "", long_options, &option_index)) {
    case '0':
    options.simburst_file = optarg;
		break;
	case 'C':
		print_usage();
		exit(0);
	case 'V':
		options.output = optarg;
		break;
	case 1728:
		options.mdc_gps_start=atof(optarg);
		break;
	case 1729:
		options.mdc_duration=  atoi(optarg);
		break;
  case 1730:
    ParseCharacterOptionString(optarg, &(options.ifonames), &(options.nIFO));
    break;
  case 1732:
    ParseCharacterOptionString(optarg, &(options.channames), &nchannels);
    break;
  case 1731:
    options.pad=atof(optarg);
    break;
      
	case 0:
		/* option sets a flag */
		break;
	case -1:
		/* end of arguments */
		break;
	case '?':
		/* unrecognized option */
		print_usage();
		exit(1);
	case ':':
		/* missing argument for an option */
		print_usage();
		exit(1);
	} while(c != -1);
    
    
    /* check some of the input parameters for consistency */
    if (options.simburst_file==NULL){
        fprintf(stderr,"Must provide a sim burst file with --simburst-file\n");
        exit(1);
    }
    if (options.mdc_duration==-1){
        fprintf(stdout,"\n mdc_duration not provided. Will default to (max_t - min_t +2 pad) where max_t and min_t are the max and min trigtimes contained in the simburst table. You can change this behavior with --mdc-duration\n\n");
    }
    if (options.mdc_gps_start==-1){
        fprintf(stdout,"mdc-gps-start time not provided. Will default to (min_t - pad) where min_t is the earliest trigtime contained in the simburst table. You can change this behavior with --mdc-gps-start\n\n");
    }
    if (options.pad==-1){
        fprintf(stdout,"mdc padding length not provided. Will default to 10 seconds. You can change this behavior with --pad\n\n");
    }
    if (options.ifonames==NULL){
        fprintf(stderr,"Must provide a csv list of ifos to write into frame: --ifos [H1,L1,V1] \n");
        exit(1);
    }
    if (options.ifonames && options.channames){
        if (options.nIFO != nchannels){
            fprintf(stderr,"Must provide a channel name for each ifo: --channels [H1:inj,L1:inj,V1:inj]\n");
            exit(1);
        }
    }
    if (! options.channames){
        fprintf(stdout,"channels not provided. Will default to ifoname:Strain for each ifo. You can change this with --channels\n\n");
    }
    return options;
}



int main(int argc, char **argv)
{

	struct options options;
	TimeSlide *time_slide_table_head = NULL;
	SimBurst *injs = NULL;
	SimBurst *inj,*cutinjs = NULL;

	/*
	 * Command line and process params table.
	 */
     
	options = parse_command_line(&argc, &argv);
	
    injs=XLALSimBurstTableFromLIGOLw(options.simburst_file,0,0);
    cutinjs=injs;
    
    time_slide_table_head=XLALTimeSlideTableFromLIGOLw(options.simburst_file);
    
    double xml_gps_start;
    double xml_gps_end;
    UINT4 lost_before=0;
        
	double mdc_gps_start=0;
	int mdc_duration=0;
    double pad=options.pad;
    double mdc_max_time=0.; // that is the time of the last trigger we can write into the MDC , which is at least pad seconds far from the end of the frame
    double mdc_min_time=0.; // that is the time of the first trigger we can write into the MDC , which is at least pad seconds far from the beginning of the frame

    LIGOTimeGPS epoch;
    char channame[256];
    
    INT4 detectorFlags;
	LALFrameH *frame=NULL;
	CHAR fname[256];
    REAL8 trigtime=0.0;
	REAL8 srate=16384.;
	REAL8 deltaT=1./srate;
	REAL8 seglen;
    UINT4 events=0;
    UINT4 i;
    inj=injs;
    
    /* if user did not provide pad, default to 10seconds */
    if (pad==-1)
        pad=10.;
    
    xml_gps_start=inj->time_geocent_gps.gpsSeconds+1.0e-9 * inj->time_geocent_gps.gpsNanoSeconds;
    while (inj){
        xml_gps_end=inj->time_geocent_gps.gpsSeconds+1.0e-9 * inj->time_geocent_gps.gpsNanoSeconds;
        events++;
        inj=inj->next;
    }
    
    if (options.mdc_gps_start==-1)
        mdc_gps_start=floor(xml_gps_start-pad);
    else
        mdc_gps_start=options.mdc_gps_start;
        
    XLALGPSSet(&epoch,mdc_gps_start,0); 
    
    if (options.mdc_duration!=-1)
        mdc_duration=options.mdc_duration;
    else
        // Accomodate all the events in the xml table. Allow 1 pad per side
        mdc_duration=(int) (2.*pad + (xml_gps_end-xml_gps_start));
    
    /* Check mdc_duration is larger than 2*pad */
    if (2.*pad >= (double) mdc_duration){
        fprintf(stderr,"ERROR: The duration of the frame (%d secs) is smaller than twice the pad (%.2f secs). Consider increasing mdc-duration or reducing pad\n",mdc_duration,2.*pad);
        exit(1);
    }    
    
    mdc_max_time=xml_gps_start+mdc_duration-2.0*pad;
    mdc_min_time=mdc_gps_start+1.0*pad;
    
    inj=injs;
    i=0;
    
    /* Now loop over the table, and only keep times smaller than mdc_max_time */
    while(inj){
        trigtime=inj->time_geocent_gps.gpsSeconds+1.0e-9 * inj->time_geocent_gps.gpsNanoSeconds;
        if (trigtime < mdc_min_time){
            /* This event is happening before the start of the frame, or is too close to the beginning of the frame. Skip it moving the start of injs by 1 */
            cutinjs=inj->next;
            inj=inj->next;
            lost_before++;
            continue;
        }
        if (trigtime > mdc_max_time){
            /* This event is happening after the end of the frame, or is too close to the end of the frame. That can happen only if this time is the first time of the xml (otherwise we would not have arrived here, because of the last if of this block) */
            i=0;
            break;
        }

        i++;

        if (inj->next){
            if( inj->next->time_geocent_gps.gpsSeconds+1.0e-9 * inj->next->time_geocent_gps.gpsNanoSeconds<=mdc_max_time){
            /* Next event is still going to be in range, so continue*/
            inj=inj->next;
            continue;
            }
        }
        
        /* Next event is going to be out of range. Set next to NULL*/
        inj->next=NULL;
        break;
    }
    
    if (i==0){
        fprintf(stderr,"ERROR: The frame will not contain any injections. Check the gps-start, duration, and pad\n");
        exit(1);
    }
    if (i<events){
        if (i<events-lost_before)
            fprintf(stdout,"WARNING: The frame will miss the last %d injections of the XML file, to follow duration and padding request!!\n\n",events-lost_before-i);
        if (lost_before>0)
            fprintf(stdout,"WARNING: The frame will miss the first %d injections of the XML file, to follow duration and padding request!!\n\n",lost_before);
    }
    
    
    /* point injs to cutinjs, i.e. to the first row with time >mdc_min_time */
    injs=cutinjs;
    inj=injs;
    
    CHAR frameType[256];
    sprintf(frameType,"%s",inj->waveform);
    
    /* Create a string with the first char of the IFO name to use as prefix for the frame name */
    char IFOs[10]="";
    i=0;
    while(i<options.nIFO){
         sprintf(IFOs,"%s%c",IFOs,options.ifonames[i][0]);
         i++;
    }
    
    snprintf( fname, 256, "%s-%s-%d-%d.gwf", IFOs,frameType, (int) mdc_gps_start, mdc_duration );
    
    /* set detector flags */
	detectorFlags = LAL_GEO_600_DETECTOR_BIT | LAL_LHO_4K_DETECTOR_BIT |
			LAL_LHO_2K_DETECTOR_BIT | LAL_LLO_4K_DETECTOR_BIT |
			LAL_TAMA_300_DETECTOR_BIT | LAL_VIRGO_DETECTOR_BIT;
		
	seglen = mdc_duration*srate;
	frame = XLALFrameNew( &epoch, mdc_duration, "LIGO", 0, 1,detectorFlags );
	
	for (i=0;i<options.nIFO;i++){
		
        if (options.channames)
            sprintf(channame,"%s:%s",options.ifonames[i],options.channames[i]);
        else
            sprintf(channame,"%s:Science",options.ifonames[i]);
        
		REAL8TimeSeries *soft=NULL;
		soft = XLALCreateREAL8TimeSeries(channame,&epoch,0.0,deltaT,&lalStrainUnit,	seglen);
		memset(soft->data->data,0.0,soft->data->length*sizeof(REAL8));
		XLALBurstInjectSignals(soft,inj,time_slide_table_head , NULL);
		/*char foutname[50]="";
        sprintf(foutname,"MDC_create_time_%s",IFOname);
		FILE * fout = fopen(foutname,"w");
        UINT4 j=0;
		for (j=0;j<soft->data->length;j++)
		fprintf(fout,"%lf %10.10e\n", epoch.gpsSeconds+j*deltaT, soft->data->data[j]);
		fclose(fout);*/
		XLALFrameAddREAL8TimeSeriesSimData( frame, soft );
		XLALDestroyREAL8TimeSeries(soft);
	}
    
	XLALFrameWrite( frame, fname);
	XLALFrameFree(frame);
    
    write_log(&injs, time_slide_table_head, &options, fname);
    
    if(options.ifonames)
        XLALFree(options.ifonames);
	return 0;
}

static void write_log(SimBurst **injs, TimeSlide * time_slide_table_head,struct options *options,char* fname){
        
    SimBurst *inj=&(*injs[0]);
    (void) time_slide_table_head;
    char log_filename[256];
    sprintf(log_filename,"%s.log",fname);
    FILE * log_file=fopen(log_filename,"w");
    REAL8 geoc_time;
    REAL8 ra, dec,psi,Fplus,Fcross;
    LALStatus status;
    memset(&status,0,sizeof(LALStatus));
    REAL8 q=0.0;
    REAL8 f=0.0;
    REAL8 hrss=0.0;
    SkyPosition currentEqu, currentGeo;
    LIGOTimeGPS injtime;
    REAL8 gmst,timedelay;
    char IFOnames[5][4]={"GEO","H1","H2","L1","V1"};    
    int nifos=5;
    int i=0;
    
    currentEqu.system = COORDINATESYSTEM_EQUATORIAL;
    currentGeo.system = COORDINATESYSTEM_GEOGRAPHIC;


    while(inj){
        
    
    q=inj->q;
    f=inj->frequency;
    hrss=inj->hrss;
    ra=inj->ra;
    dec=inj->dec;
    psi=inj->psi;
    geoc_time=inj->time_geocent_gps.gpsSeconds + 1.e-9* inj->time_geocent_gps.gpsNanoSeconds;  
    XLALGPSSet(&injtime,inj->time_geocent_gps.gpsSeconds ,inj->time_geocent_gps.gpsNanoSeconds);
    
    currentEqu.latitude = dec;
    currentEqu.longitude = ra;
    double mdc_gps_start=options->mdc_gps_start;
    
    LALEquatorialToGeographic(&status, &currentGeo, &currentEqu, &injtime);
    
    char wf[4];
    if (!strcmp("SineGaussian",inj->waveform))
        strcpy(wf,"SG");
    else if (!strcmp("SineGaussianF",inj->waveform))
        strcpy(wf,"SGF");   
    else if (!strcmp("Gaussian",inj->waveform))
        strcpy(wf,"GA");
    else if (!strcmp("GaussianF",inj->waveform))
        strcpy(wf,"GAF");
    else if (!strcmp("StringCusp",inj->waveform))
        strcpy(wf,"SC");
    else if (!strcmp("BTLWNB",inj->waveform))
        strcpy(wf,"WN");
    else{
        fprintf(stderr,"Uknown waveform value %s. Exiting\n",inj->waveform);
        exit(1);
    }
	
    char WFname[256]; 
    sprintf(WFname,"%s%dQ%dd%d",wf,(int) f,(int) floor(q),(int) (10.*(q-floor(q))));
    fprintf(log_file,"%s\
    -100 \
    -100 \
    %.10e\
    -100\
    -100\
    %.10e\
    %.10e\
    %.10e\
    %d\
    %10.10f\
    %s\
    -100\
    -100\
    -100 ",WFname,hrss, cos(LAL_PI_2 - currentGeo.latitude),currentGeo.longitude,psi,(int) (mdc_gps_start), geoc_time,WFname);
    
    i=0;
    while (i<nifos){
        LALDetector* detector=XLALCalloc(1,sizeof(LALDetector));
        if (!strcmp(IFOnames[i],"H1") || !strcmp(IFOnames[i],"H2"))
            memcpy(detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
        else if (!strcmp(IFOnames[i],"L1"))
             memcpy(detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
        else if (!strcmp(IFOnames[i],"V1"))
            memcpy(detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));
        else if(!strcmp(IFOnames[i],"GEO")) 
            memcpy(detector,&lalCachedDetectors[LALDetectorIndexGEO600DIFF],sizeof(LALDetector));
        gmst=XLALGreenwichMeanSiderealTime(&injtime);
        XLALComputeDetAMResponse(&Fplus, &Fcross,(const REAL4(*)[3]) detector->response, ra, dec, psi, gmst);
        timedelay = XLALTimeDelayFromEarthCenter(detector->location,ra, dec, &injtime);
    
        fprintf(log_file, "%s %10.10f %.10e %.10e ", IFOnames[i], geoc_time+timedelay,Fplus,Fcross);
    
        i++;
        XLALFree(detector);
    }    
    
    fprintf(log_file, "\n");
    inj=inj->next;
    }
    
    fclose(log_file);    

}


void ParseCharacterOptionString(char *input, char **strings[], UINT4 *n)
/* "Stolen" from LALInference.c */
/* parses a character string (passed as one of the options) and decomposes   */
/* it into individual parameter character strings. Input is of the form      */
/*   input   :  "[one,two,three]"                                            */
/* and the resulting output is                                               */
/*   strings :  {"one", "two", "three"}                                      */
/* length of parameter names is for now limited to 512 characters.           */
/* (should 'theoretically' (untested) be able to digest white space as well. */
/* Irrelevant for command line options, though.)                             */
{
  UINT4 i,j,k,l;
  /* perform a very basic well-formedness-check and count number of parameters: */
  i=0; j=0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==',')) ++*n;
    if ((j==1) & (input[i]==']')) {++*n; j=2;}
    ++i;
  }
  if (j!=2) XLAL_ERROR_VOID(XLAL_EINVAL, "Argument vector \"%s\" is not well-formed!", input);
  /* now allocate memory for results: */
  *strings  = (char**)  XLALMalloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*strings)[i] = (char*) XLALMalloc(sizeof(char)*512);
  i=0; j=0;
  k=0; /* string counter    */
  l=0; /* character counter */
  while ((input[i] != '\0') & (j<3)) {
    /* state transitions: */
    if ((j==0) & ((input[i]!='[') & (input[i]!=' '))) j=1;
    if (((j==1)|(j==2)) & (input[i]==',')) {(*strings)[k][l]='\0'; j=2; ++k; l=0;}
    if ((j==1) & (input[i]==' ')) j=2;
    if ((j==1) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & (input[i]==']')) {(*strings)[k][l]='\0'; j=3;}
    if ((j==2) & ((input[i]!=']') & (input[i]!=',') & (input[i]!=' '))) j=1;
    /* actual copying: */
    if (j==1) {
      if (l>=511) {
        XLAL_PRINT_WARNING("Character \"%s\" argument too long!", (*strings)[k]);
      }
      else {
        (*strings)[k][l] = input[i];
        ++l;
      }
    }
    ++i;
  }
}
