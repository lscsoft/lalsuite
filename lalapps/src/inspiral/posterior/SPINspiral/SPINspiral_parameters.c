/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   SPINspiral_parameters.c:   routines to read/write input files, set constants and set true and null parameters
   
   
   Copyright 2007, 2008, 2009 Christian Roever, Marc van der Sluys, Vivien Raymond, Ilya Mandel
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <getopt.h>

#include <lal/LIGOMetadataTables.h>
//#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>

#include <SPINspiral.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif



/**
 * \file SPINspiral_parameters.c
 * \brief Contains routines that handle reading and randomisation of parameters
 */



// ****************************************************************************************************************************************************  
/**
 * \brief Read and parse command-line arguments/options
 */
// ****************************************************************************************************************************************************  
void readCommandLineOptions(int argc, char* argv[], struct runPar *run)
{
  int i = 0;
  int C;
  int nIFO = 0;
  int nChannel = 0;
  int nCache = 0;
  char **networkseq=NULL;
  char **channelseq=NULL;
  run->outputPath=NULL;
  run->triggerMc = 0.0;               // Chirp mass from the command line - zero means no value
  run->triggerEta = 0.0;              // Eta from the command line - zero means no value
  run->triggerTc = 0.0;               // Time of coalescence from the command line - zero means no value
  run->triggerDist = 0.0;             // Distance (in Mpc) from the command line - zero means no value
  run->PSDstart = 0.0;                // GPS start of the PSD - zero means no value
  run->tukey1 = 0.15;
  run->tukey2 = 0.01;//
  for(i=0;i<3;i++) strcpy(run->channelname[i],"");
  
  if(argc > 1) printf("   Parsing %i command-line arguments:\n",argc-1);
  
  
  //Set up struct with long (--) options:
  static struct option long_options[] =
    {
		{"injXMLfile",      required_argument, 0,               0},
		{"injXMLnr",        required_argument, 0,               0},
		{"mChirp",          required_argument, 0,             'm'},
		{"eta",             required_argument, 0,             'e'},
		{"tc",              required_argument, 0,             't'},
		{"dist",            required_argument, 0,             'd'},
		{"template",        required_argument, 0,             'w'},
		{"tukey1",          required_argument, 0,             'y'},
		{"tukey2",          required_argument, 0,             'z'},
		{"nIter",           required_argument, 0,             'n'},
		{"nSkip",           required_argument, 0,             's'},
		{"rseed",           required_argument, 0,             'r'},
		{"network",         required_argument, 0,             'a'},
		{"channel",         required_argument, 0,             'a'},
		{"downsample",      required_argument, 0,             'a'},
		{"beforetc",        required_argument, 0,             'a'},
		{"aftertc",         required_argument, 0,             'a'},
		{"Flow",            required_argument, 0,             'a'},
		{"Fhigh",           required_argument, 0,             'a'},
		{"nPSDsegment",     required_argument, 0,             'a'},
		{"lPSDsegment",     required_argument, 0,             'a'},
		{"PSDstart",        required_argument, 0,             'a'},
		{"outputPath",      required_argument, 0,             'o'},
		{"cache",           required_argument, 0,             'c'},      
		{0, 0, 0, 0}
    };
  
  
  int option_index = 0;
  while( (C = getopt_long(argc, argv, "i:m:e:t:d:n:d:a:o:c:r:w:y:x:",long_options, &option_index)) != -1) {
    switch(C) {
      
      
      // *** Treat (untranslated) long options:
    case 0:
      if(strcmp(long_options[option_index].name,"injXMLfile")==0) {
        run->injXMLfilename=(char*)malloc(strlen(optarg)+1);
        strcpy(run->injXMLfilename,optarg);
        printf("    - reading injection parameters from the XML file %s\n",run->injXMLfilename);
      }
      if(strcmp(long_options[option_index].name,"injXMLnr")==0) {
        run->injXMLnr = atoi(optarg);
        printf("    - using injection %d from the injection XML file\n",run->injXMLnr);
      }
      break; //For case 0: long options
      
      
      // *** Treat the short and translated long options:
    case 'i':
      strcpy(run->mainFilename,optarg);
      printf("    - using main input file %s\n",run->mainFilename);
      break;
      
    case 'm':           
      run->triggerMc = atof(optarg);
      printf("    - From command line, trigger value for mChirp\t\t= %f\n",run->triggerMc);
      break;
      
    case 'e':           
      run->triggerEta = atof(optarg);
      printf("    - From command line, trigger value for eta\t\t\t= %f\n",run->triggerEta);
      break;
      
    case 't':           
      run->triggerTc = atof(optarg);
      printf("    - From command line, trigger value for tc\t\t\t= %f\n",run->triggerTc);
      break;
      
    case 'd':           
      run->triggerDist = atof(optarg);
      printf("    - From command line, trigger value for the distance (Mpc)\t= %f\n",run->triggerDist);
      break;
      
    case 'n':           
      run->nIter = atoi(optarg);
      run->commandSettingsFlag[0] = 1;
      printf("    - From command line, number of iterations\t\t\t= %d\n",run->nIter);
      break;
      
    case 's':           
      run->thinOutput = atoi(optarg);
      run->commandSettingsFlag[1] = 1;
      printf("    - From command line, thin output by\t\t\t\t= %d\n",run->thinOutput);
      break;
      
	case 'r':           
		run->MCMCseed = atoi(optarg);
		run->commandSettingsFlag[16] = 1;
		printf("    - From command line, random seed\t\t\t\t= %d\n",run->MCMCseed);
		break;
		
	case 'w':           
		run->mcmcWaveform = atoi(optarg);
		run->commandSettingsFlag[17] = 1;
		if (run->mcmcWaveform == 1) {
			printf("    - From command line, mcmc template\t\t\t\t= %d\t 1.5PN 12-parameter SINGLE SPIN Apostolatos (Unless specified otherwise in an input file)\n",run->mcmcWaveform);
			run->nMCMCpar=12;
			fprintf(stderr,"        1.5PN 12-parameter SINGLE SPIN Apostolatos not available on command line \n");
			exit(1);
		} else if (run->mcmcWaveform == 2) {
			printf("    - From command line, mcmc template\t\t\t\t= %d\t 12-parameter SINGLE SPIN LAL (Unless specified otherwise in an input file)\n",run->mcmcWaveform);
			run->nMCMCpar=12;
			fprintf(stderr,"        12-parameter SINGLE SPIN LAL  not available on command line \n");
			exit(1);
		} else if (run->mcmcWaveform == 3) {
			printf("    - From command line, mcmc template\t\t\t\t= %d\t 15-parameter DOUBLE SPIN LAL (Unless specified otherwise in an input file)\n",run->mcmcWaveform);
			run->nMCMCpar=15;
		} else if (run->mcmcWaveform == 4) {
			printf("    - From command line, mcmc template\t\t\t\t= %d\t 9-parameter NO SPIN LAL (Unless specified otherwise in an input file)\n",run->mcmcWaveform);
			run->nMCMCpar=9;
		} else {
			fprintf(stderr,"    - unknown waveform chosen as MCMC template: %d.   Available tested waveforms are:\n",run->mcmcWaveform);
			fprintf(stderr,"        1: Apostolatos, simple precession, 12 parameters\n");
			fprintf(stderr,"        2: LAL, single spin, 12 parameters\n");
			fprintf(stderr,"        3: LAL, double spin, 15 parameters\n");
			fprintf(stderr,"        4: LAL, non-spinnig, 9 parameters\n");
			fprintf(stderr,"      Please set mcmcWaveform in %s to one of these values. Input file or default DOUBLE SPIN will be choosen.\n\n",run->mcmcFilename);
		}
		break;
	
	case 'y':           
		run->tukey1 = atof(optarg);
		printf("    - From command line, tukey window alpha_1\t\t\t= %f\n",run->tukey1);
		break;
		
	case 'z':           
		run->tukey2 = atof(optarg);
		printf("    - From command line, tukey window alpha_2\t\t\t= %f\n",run->tukey2);
		break;		
			
    case 'a': //Detector options
      if(strcmp(long_options[option_index].name,"network")==0) {
        parseCharacterOptionString(optarg,&networkseq,&nIFO);
        run->networkSize = nIFO;
        run->commandSettingsFlag[2] = 1;
        printf("    - From command line, network size\t\t\t\t= %d\n",run->networkSize);
        for(i=0;i<run->networkSize;i++) {
          run->selectifos[i] = atoi(networkseq[i]);
          printf("    - From command line, IFO%d\t\t\t\t\t= %d\n",(i+1),run->selectifos[i]);
        }
        run->commandSettingsFlag[3] = 1;
      }
      
      if(strcmp(long_options[option_index].name,"channel")==0) {
        parseCharacterOptionString(optarg,&channelseq,&nChannel);
        if (run->networkSize != nChannel) {printf(" ERROR: number of IFOs %d should be the same as number of channels %d\n",nIFO,nChannel); exit(1);}
        else {
          for(i=0;i<run->networkSize;i++) {
            strcpy(run->channelname[i],channelseq[i]);
            printf("    - From command line, channel %d\t\t\t\t= %s\n",(i+1),run->channelname[i]);
          }
          run->commandSettingsFlag[4] = 1;
        }
      }
      
      
      if(strcmp(long_options[option_index].name,"downsample")==0) {
        run->downsampleFactor = atoi(optarg);
        run->commandSettingsFlag[6] = 1;
        printf("    - From command line, downsample factor\t\t\t= %d\n",run->downsampleFactor);
      }
      if(strcmp(long_options[option_index].name,"beforetc")==0) {
        run->dataBeforeTc = atof(optarg);
        run->commandSettingsFlag[7] = 1;
        printf("    - From command line, before tc\t\t\t\t= %f\n",run->dataBeforeTc);
      }
      if(strcmp(long_options[option_index].name,"aftertc")==0) {
        run->dataAfterTc = atof(optarg);
        run->commandSettingsFlag[8] = 1;
        printf("    - From command line, after tc\t\t\t\t= %f\n",run->dataAfterTc);
      }
      if(strcmp(long_options[option_index].name,"Flow")==0) {
        run->lowFrequencyCut = atof(optarg);
        run->commandSettingsFlag[9] = 1;
        printf("    - From command line, low frequency cut\t\t\t= %f\n",run->lowFrequencyCut);
      }
      if(strcmp(long_options[option_index].name,"Fhigh")==0) {
        run->highFrequencyCut = atof(optarg);
        run->commandSettingsFlag[10] = 1;
        printf("    - From command line, high frequency cut\t\t\t= %f\n",run->highFrequencyCut);
      }
      if(strcmp(long_options[option_index].name,"nPSDsegment")==0) {
        run->PSDsegmentNumber = atoi(optarg);
        run->commandSettingsFlag[11] = 1;
        printf("    - From command line, number of PSD segments\t\t\t= %d\n",run->PSDsegmentNumber);
      }
      if(strcmp(long_options[option_index].name,"lPSDsegment")==0) {
        run->PSDsegmentLength = atof(optarg);
        run->commandSettingsFlag[12] = 1;
        printf("    - From command line, length of PSD segments\t\t\t= %f\n",run->PSDsegmentLength);
      }
      if(strcmp(long_options[option_index].name,"PSDstart")==0) {
        run->PSDstart = atof(optarg);
        run->commandSettingsFlag[13] = 1;
        printf("    - From command line, start of PSD segments\t\t\t= %f\n",run->PSDstart);
      }
      
      break;            
      
    case 'o':           
      run->outputPath=(char*)malloc(strlen(optarg)+1);
      strcpy(run->outputPath,optarg);
      printf("    - From command line, output path\t\t\t\t= %s\n",run->outputPath);
      break;
    case 'c':
      
      parseCharacterOptionString(optarg,&(run->cacheFilename),&nCache);
      if (run->networkSize != nCache) {printf(" ERROR: number of IFOs %d should be the same as number of cache files %d\n",nIFO,nCache); exit(1);}
      else {
        for(i=0;i<run->networkSize;i++) {
          printf("    - From command line, cache file %d\t\t\t\t= %s\n",(i+1),run->cacheFilename[i]);
          readCachefile(run,i);
        }
        run->commandSettingsFlag[15] = 1;
      }
      break;
      
    default:
      //fprintf(stderr,"   Unrecognised option: %d\n",c);  // This notice is already produced by getopt_long()
      fprintf(stderr,USAGE); 
      exit(1);
      break;
      
    } // switch()
  } // while()
  
  
  // Print any remaining command line arguments (the ones that are not options, i.e. lack a starting '-'):
  if(optind < argc) {
    printf("   SPINspiral - unused command-line arguments: ");
    while(optind < argc) printf ("%s ", argv[optind++]);
    printf("\n");
  }
  
} // End void readCommandLineOptions(argc,argv)
// ****************************************************************************************************************************************************  




// ****************************************************************************************************************************************************  
/**
 * \brief parse strings from "[one,two,three]" to {"one", "two", "three"}
 *
 */
// ****************************************************************************************************************************************************  

void parseCharacterOptionString(char *input, char **strings[], int *n)
// ****************************************************************************************************************************************************  
/**
 * \brief Read and parse command-line arguments/options
 */
// ****************************************************************************************************************************************************  
/**
 * Parses a character string (passed as one of the options) and decomposes
 * it into individual parameter character strings. Input is of the form
 * input   :  "[one,two,three]"
 * and the resulting output is
 * strings :  {"one", "two", "three"}
 * length of parameter names is by now limited to 512 characters.
 * (should 'theoretically' (untested) be able to digest white space as well.
 * Irrelevant for command line options, though.)
 */
{
  int i,j,k,l;
  /* perform a very basic well-formedness-check and count number of parameters: */
  i=0; j=0;
  *n = 0;
  while (input[i] != '\0') {
    if ((j==0) & (input[i]=='[')) j=1;
    if ((j==1) & (input[i]==',')) ++*n;
    if ((j==1) & (input[i]==']')) {++*n; j=2;}
    ++i;
  }
  if (j!=2) {fprintf(stderr, " ERROR: argument vector \"%s\" not well-formed!\n", input); exit(1);}
  /* now allocate memory for results: */
  *strings  = (char**)  malloc(sizeof(char*) * (*n));
  for (i=0; i<(*n); ++i) (*strings)[i] = (char*) malloc(sizeof(char)*512);
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
        fprintf(stderr, " WARNING: character argument too long!\n");
        fprintf(stderr, " \"%s\"\n",(*strings)[k]);
      }
      else {
        (*strings)[k][l] = input[i];
        ++l;
      }
    }
    ++i;
  } 
}  // End void parseCharacterOptionString(char *input, char **strings[], int *n)
// ****************************************************************************************************************************************************  






// ****************************************************************************************************************************************************  
/**
 * \brief Read the main input file
 *
 * All parameters that are read here should be(come) members of the runvar struct and lose their global status
 */
// ****************************************************************************************************************************************************  
void readMainInputfile(struct runPar *run)
{
  int i=0;
  char tmpStr[500],*cstatus=NULL;
  FILE *fin;
  char UNUSED *str;
  
	if((fin = fopen(run->mainFilename,"r")) == NULL) {
		fprintf(stderr, "   No main input file: %s, using default values.\n",run->mainFilename);
		
		run->doSNR=1;
		run->doMCMC=1;
		run->doMatch=0;
		run->writeSignal=0;
		run->beVerbose=1;
		
		sprintf(run->mcmcFilename, "SPINspiral.input.mcmc");
		sprintf(run->dataFilename, "SPINspiral.input.data");
		sprintf(run->injectionFilename, "SPINspiral.input.injection");
		sprintf(run->parameterFilename, "SPINspiral.input.parameters");
		sprintf(run->systemFilename, "SPINspiral.input.system");
		
		//exit(1);
		
  } else {
    printf("   Using main input file: %s.\n",run->mainFilename);
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) cstatus = fgets(tmpStr,500,fin);  //Read first 3 lines
  
  //Operation and output:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->doSNR);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->doMCMC);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->doMatch);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->writeSignal);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->beVerbose);
  
  
  //Secondary input files:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin); //Read the empty and comment lines
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%s",run->mcmcFilename);
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%s",run->dataFilename);
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%s",run->injectionFilename);
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%s",run->parameterFilename);
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%s",run->systemFilename);
  
  fclose(fin);
	    }
 
  if(!cstatus) {
      printf("Error reading string\n");
      exit(1);
  }

  
}  //End of readMainInputfile
// ****************************************************************************************************************************************************  










// ****************************************************************************************************************************************************  
/**
 * \brief Read the MCMC input file
 *
 * All parameters that are read here should be(come) members of the runvar struct and lose their global status
 */
// ****************************************************************************************************************************************************  
void readMCMCinputfile(struct runPar *run)
{
  int i=0,istatus=0;
  double tmpdbl=0.0;
  char tmpStr[500],*cstatus=NULL;
  FILE *fin;
  
	if((fin = fopen(run->mcmcFilename,"r")) == NULL) {
		fprintf(stderr, "   No MCMC input file: %s, using default values.\n",run->mcmcFilename);
		
		if(run->commandSettingsFlag[0] == 0) {
			run->nIter = 2000000;
		}
		
		if(run->commandSettingsFlag[1] == 0) {        
			run->thinOutput = 100;
		}
		
		run->thinScreenOutput = 10000;
		
		if(run->commandSettingsFlag[16] == 0) {        
			run->MCMCseed = 0;
		}
		
		run->adaptiveMCMC = 1;
		run->acceptRateTarget = 0.25;
		run->minlogL = 0.0;
		run->blockFrac = 0.012;
		
		//Correlated update proposals:
		run->correlatedUpdates = 0;
		run->corrFrac = 0.7;
		run->nCorr = 1000;
		run->matAccFr = 0.5;
		run->prMatrixInfo = 0;
		
		
		//Annealing:
		run->annealTemp0 = 1.00;
		run->annealNburn = 100000;
		run->annealNburn0 = 100000;
		
		//Parallel tempering:
		run->parallelTempering = 1;
		run->nTemps = 5;
		run->maxTemp = 40.0;
		run->saveHotChains = 100;
		run->prParTempInfo = 1;
		
		//Manual temperature ladder for parallel tempering:
		run->tempLadder[0] = 1.00;
		run->tempLadder[1] = 2.51;
		run->tempLadder[2] = 6.32;
		run->tempLadder[3] = 15.91;
		run->tempLadder[4] = 40.00;
		
		//exit(1);
	} else {
    printf("   Using MCMC input file: %s.\n",run->mcmcFilename);

  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=3;i++) { //Read first 3 lines
    cstatus = fgets(tmpStr,500,fin);
  }
  
  //Basic settings
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  
  cstatus = fgets(tmpStr,500,fin);  
  if(run->commandSettingsFlag[0] == 0) {
    sscanf(tmpStr,"%lg",&tmpdbl);
    run->nIter = (int)tmpdbl;
  }
  
  cstatus = fgets(tmpStr,500,fin);
  if(run->commandSettingsFlag[1] == 0) {        
    sscanf(tmpStr,"%d",&run->thinOutput);
  }
  
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->thinScreenOutput);
	
	cstatus = fgets(tmpStr,500,fin);
	if(run->commandSettingsFlag[16] == 0) {
		sscanf(tmpStr,"%d",&run->MCMCseed);
	}	 

  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->adaptiveMCMC);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->acceptRateTarget);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->minlogL);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->blockFrac);
  
  
  
  //Correlated update proposals:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->correlatedUpdates);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->corrFrac);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lg",&tmpdbl);
  run->nCorr = (int)tmpdbl;
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->matAccFr);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->prMatrixInfo);
  
  
  //Annealing:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->annealTemp0);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lg",&tmpdbl);
  run->annealNburn = (int)tmpdbl;
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lg",&tmpdbl);
  run->annealNburn0 = (int)tmpdbl;
  
  //Parallel tempering:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->parallelTempering);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->nTemps);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->maxTemp);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->saveHotChains);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->prParTempInfo);
  
  //Manual temperature ladder for parallel tempering:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin); //Read the empty and comment line
  for(i=0;i<run->nTemps;i++) istatus = fscanf(fin,"%lf",&run->tempLadder[i]);  //Read the array directly, because sscanf cannot be in a loop...
  
  fclose(fin);
	}
  
  if(!cstatus || !istatus) {                                                                                                           
     printf("Error reading string\n");                                                                                   
     exit(1);                                                                                                             
  }  
} //End of void readMCMCinputfile(struct runPar *run)
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Read the data input file
 */
// ****************************************************************************************************************************************************  
void readDataInputfile(struct runPar *run, struct interferometer ifo[])
{
  int i=0,j=0,istatus=0;
  double lati=0.0,longi=0.0,leftArm=0.0,rightArm=0.0;
  char tmpStr[500],*cstatus=NULL, subdir[500];
  FILE *fin;
  int dump = 0;
  
	if((fin = fopen(run->dataFilename,"r")) == NULL) {
		fprintf(stderr, "   No data file: %s, using default values.\n",run->dataFilename);
		
		
		sprintf(run->datasetName, "default");
		
		if(run->commandSettingsFlag[2] == 0) run->networkSize = 1;
		if(run->commandSettingsFlag[3] == 0) run->selectifos[0]=1;
		
		//Data handling:
		if(run->commandSettingsFlag[6] == 0) run->downsampleFactor = 8;
		if(run->commandSettingsFlag[7] == 0) run->dataBeforeTc = 6.0;
		if(run->commandSettingsFlag[8] == 0) run->dataAfterTc = 2.0;
		if(run->commandSettingsFlag[9] == 0) run->lowFrequencyCut = 40.0;
		if(run->commandSettingsFlag[10] == 0) run->highFrequencyCut = 250.0;
		run->tukeyWin = 0.1;
		
		
		//Read input for PSD estimation:
		if(run->commandSettingsFlag[11] == 0) run->PSDsegmentNumber = 32; 
		if(run->commandSettingsFlag[12] == 0) run->PSDsegmentLength = 8.0;
		
		sprintf(ifo[0].name, "Hanford");
		lati = 46.45;
		longi = -119.41;
		rightArm = 36.80;
		leftArm = 126.80;
		
		ifo[0].lati      = lati     /180.0*pi;
		ifo[0].longi     = longi    /180.0*pi;
		ifo[0].rightArm  = rightArm /180.0*pi;
		ifo[0].leftArm   = leftArm  /180.0*pi;
		
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[0].ch1name, "H1:LSC-STRAIN");}
		else{ strcpy(ifo[0].ch1name,run->channelname[0]);}
		sprintf(subdir, ".");
		sprintf(ifo[0].ch1filepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[0].ch1fileprefix, "H1_RDS_C03_L2");
		sprintf(ifo[0].ch1filesuffix, "-128.gwf");
		ifo[0].ch1filesize = 128;
		ifo[0].ch1fileoffset = 0;
		ifo[0].ch1doubleprecision = 0;
		ifo[0].add2channels = 0;
		
		ifo[0].noiseGPSstart = 0;
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[0].noisechannel, "H1:LSC-STRAIN");}
		else{ strcpy(ifo[0].noisechannel,run->channelname[0]);}  
		sprintf(subdir, ".");
		sprintf(ifo[0].noisefilepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[0].noisefileprefix, "H1_RDS_C03_L2");
		sprintf(ifo[0].noisefilesuffix, "-128.gwf");
		ifo[0].noisefilesize = 128;
		ifo[0].noisefileoffset = 0;
		ifo[0].noisedoubleprecision = 0;  
		
		
		sprintf(ifo[1].name, "Livingston");
		lati = 30.56;
		longi = -90.77;
		rightArm = 108.00;
		leftArm = 198.00;
		
		ifo[1].lati      = lati     /180.0*pi;
		ifo[1].longi     = longi    /180.0*pi;
		ifo[1].rightArm  = rightArm /180.0*pi;
		ifo[1].leftArm   = leftArm  /180.0*pi;
		
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[1].ch1name, "L1:LSC-STRAIN");}
		else{ strcpy(ifo[1].ch1name,run->channelname[1]);}
		sprintf(subdir, ".");
		sprintf(ifo[1].ch1filepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[1].ch1fileprefix, "L1_RDS_C03_L2");
		sprintf(ifo[1].ch1filesuffix, "-128.gwf");
		ifo[1].ch1filesize = 128;
		ifo[1].ch1fileoffset = 0;
		ifo[1].ch1doubleprecision = 0;
		ifo[1].add2channels = 0;
		
		ifo[1].noiseGPSstart = 0;
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[1].noisechannel, "L1:LSC-STRAIN");}
		else{ strcpy(ifo[1].noisechannel,run->channelname[1]);}  
		sprintf(subdir, ".");
		sprintf(ifo[1].noisefilepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[1].noisefileprefix, "L1_RDS_C03_L2");
		sprintf(ifo[1].noisefilesuffix, "-128.gwf");
		ifo[1].noisefilesize = 128;
		ifo[1].noisefileoffset = 0;
		ifo[1].noisedoubleprecision = 0;  
		
		
		sprintf(ifo[2].name, "Pisa");
		lati = 43.63;
		longi = 10.50;
		rightArm = 341.50;
		leftArm = 71.50;
		
		ifo[2].lati      = lati     /180.0*pi;
		ifo[2].longi     = longi    /180.0*pi;
		ifo[2].rightArm  = rightArm /180.0*pi;
		ifo[2].leftArm   = leftArm  /180.0*pi;
		
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[2].ch1name, "V1:STRAIN");}
		else{ strcpy(ifo[2].ch1name,run->channelname[2]);}
		sprintf(subdir, ".");
		sprintf(ifo[2].ch1filepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[2].ch1fileprefix, "V1_RDS_C03_L2");
		sprintf(ifo[2].ch1filesuffix, "-128.gwf");
		ifo[2].ch1filesize = 128;
		ifo[2].ch1fileoffset = 0;
		ifo[2].ch1doubleprecision = 0;
		ifo[2].add2channels = 0;
		
		ifo[2].noiseGPSstart = 0;
		if(run->commandSettingsFlag[4] ==0){ sprintf(ifo[2].noisechannel, "V1:STRAIN");}
		else{ strcpy(ifo[2].noisechannel,run->channelname[2]);}  
		sprintf(subdir, ".");
		sprintf(ifo[2].noisefilepath,"%s%s%s",run->dataDir,"/",subdir);
		sprintf(ifo[2].noisefileprefix, "V1_RDS_C03_L2");
		sprintf(ifo[2].noisefilesuffix, "-128.gwf");
		ifo[2].noisefilesize = 128;
		ifo[2].noisefileoffset = 0;
		ifo[2].noisedoubleprecision = 0;  
		
		
		//exit(1);
	} else {
    printf("   Using data input file: %s.\n",run->dataFilename);
  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(j=1;j<=3;j++) cstatus = fgets(tmpStr,500,fin);  //Read first 3 lines
  cstatus = fgets(run->datasetName,80,fin);  cstatus = fgets(tmpStr,500,fin);  //Read name of the data set used, and then the rest of the line
  
  //Detector network:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  
  cstatus = fgets(tmpStr,500,fin);
  if(run->commandSettingsFlag[2] == 0) sscanf(tmpStr,"%d",&run->networkSize);
  
  for(i=0;i<run->networkSize;i++) {
    if(run->commandSettingsFlag[3] == 0) {
      istatus = fscanf(fin,"%d",&run->selectifos[i]);  //Read the array directly, because sscanf cannot be in a loop...
    } else {
      istatus = fscanf(fin,"%d",&dump);
    }
  }
  cstatus = fgets(tmpStr,500,fin);  //Read the rest of the line
  
  
  //Data handling:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin); 
  if(run->commandSettingsFlag[6] == 0) sscanf(tmpStr,"%d",&run->downsampleFactor);
  cstatus = fgets(tmpStr,500,fin); 
  if(run->commandSettingsFlag[7] == 0) sscanf(tmpStr,"%lf",&run->dataBeforeTc);
  cstatus = fgets(tmpStr,500,fin); 
  if(run->commandSettingsFlag[8] == 0) sscanf(tmpStr,"%lf",&run->dataAfterTc);
  cstatus = fgets(tmpStr,500,fin); 
  if(run->commandSettingsFlag[9] == 0) sscanf(tmpStr,"%lf",&run->lowFrequencyCut);
  cstatus = fgets(tmpStr,500,fin); 
  if(run->commandSettingsFlag[10] == 0) sscanf(tmpStr,"%lf",&run->highFrequencyCut);
  cstatus = fgets(tmpStr,500,fin); sscanf(tmpStr,"%lf",&run->tukeyWin);
  
  
  //Read input for PSD estimation:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment lines
  cstatus = fgets(tmpStr,500,fin);  
  if(run->commandSettingsFlag[11] == 0) sscanf(tmpStr,"%d",&run->PSDsegmentNumber);
  cstatus = fgets(tmpStr,500,fin);  
  if(run->commandSettingsFlag[12] == 0) sscanf(tmpStr,"%lf",&run->PSDsegmentLength);
  
  
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment lines
  for(i=0;i<run->maxIFOdbaseSize;i++){
    cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment lines
    
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",ifo[i].name);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&lati);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&longi);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&rightArm);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&leftArm);
    
    ifo[i].lati      = lati     /180.0*pi;
    ifo[i].longi     = longi    /180.0*pi;
    ifo[i].rightArm  = rightArm /180.0*pi;
    ifo[i].leftArm   = leftArm  /180.0*pi;
    
    cstatus = fgets(tmpStr,500,fin);  //Read the empty line
    
    cstatus = fgets(tmpStr,500,fin);
    if(run->commandSettingsFlag[4] ==0){ sscanf(tmpStr,"%s",ifo[i].ch1name);}
    else{ strcpy(ifo[i].ch1name,run->channelname[i]);}
    //cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",&ifo[i].ch1filepath);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",subdir);
    sprintf(ifo[i].ch1filepath,"%s%s%s",run->dataDir,"/",subdir);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",ifo[i].ch1fileprefix);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",ifo[i].ch1filesuffix);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].ch1filesize);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].ch1fileoffset);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].ch1doubleprecision);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].add2channels);
    
    cstatus = fgets(tmpStr,500,fin);  //Read the empty line
    
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%ld",&ifo[i].noiseGPSstart);
    cstatus = fgets(tmpStr,500,fin);
    if(run->commandSettingsFlag[4] ==0){ sscanf(tmpStr,"%s",ifo[i].noisechannel);}
    else{ strcpy(ifo[i].noisechannel,run->channelname[i]);}  
    //cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",&ifo[i].noisefilepath);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",subdir);
    sprintf(ifo[i].noisefilepath,"%s%s%s",run->dataDir,"/",subdir);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",ifo[i].noisefileprefix);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%s",ifo[i].noisefilesuffix);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].noisefilesize);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].noisefileoffset);
    cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&ifo[i].noisedoubleprecision);
    
  }
  fclose(fin);
	}
  
   if(!cstatus||!istatus) {
       printf("Error reading string\n");
       exit(1);
   }  
}  //End of readDataInputfile
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Read the injection input file
 *
 * All parameters that are read in here should be(come) members of the runvar struct
 */
// ****************************************************************************************************************************************************  
void readInjectionInputfile(struct runPar *run)
{
  int i=0,istatus=0;
  char tmpStr[500],*cstatus=NULL;
  FILE *fin;
  
  // Open injection input file:
	if((fin = fopen(run->injectionFilename,"r")) == NULL) {
		fprintf(stderr, "   No injection input file: %s, using default values.\n",run->injectionFilename);
		
		run->injectSignal = 0;
		run->injectionWaveform = 3;
		run->injectionPNorder = 3.5;
		run->injectionSNR = 0.0;
		run->injRanSeed = 12345;
		run->nInjectPar = 0;
		
		//exit(1);
	} else {
    printf("   Using injection input file: %s.\n",run->injectionFilename);

  
  
  // Read injection input file:
  // use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=2;i++) cstatus = fgets(tmpStr,500,fin);  //Read first 2 lines
  
  // Read general injection data:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->injectSignal);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->injectionWaveform);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->injectionPNorder);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->injectionSNR);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->injRanSeed);
  
  //Get the number of injection parameters from the injectionWaveform
  if(run->injectSignal >= 1) {
    if(run->injectionWaveform==1) {
      if(run->beVerbose>=1) printf("    - using Apostolatos, 1.5-pN, 12-parameter waveform for the software injection.\n"); // Only the 1.5-pN order is available
      run->nInjectPar=12;
    } else if(run->injectionWaveform==2) {
      if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 12-parameter waveform for the software injection.\n",run->injectionPNorder);
      run->nInjectPar=12;
    } else if(run->injectionWaveform==3) {
      if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 15-parameter waveform for the software injection.\n",run->injectionPNorder);
      run->nInjectPar=15;
    } else if(run->injectionWaveform==4) {
      if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, non-spinning waveform for the software injection.\n",run->injectionPNorder);
      run->nInjectPar=9;
	} else if(run->injectionWaveform==5) {
	  if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 15-parameter PhenSpinTaylorRD waveform for the software injection.\n",run->injectionPNorder);
	  run->nInjectPar=15;
	} else if(run->injectionWaveform==9) {
	  if(run->beVerbose>=1) printf("    - using analytic Likelihood.\n");
	  run->nInjectPar=15;		
    } else {
      fprintf(stderr,"    - unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->injectionWaveform);
      fprintf(stderr,"        1: Apostolatos, simple precession, 12 parameters\n");
      fprintf(stderr,"        2: LAL, single spin, 12 parameters\n");
      fprintf(stderr,"        3: LAL, double spin, 15 parameters\n");
      fprintf(stderr,"        4: LAL, non-spinning, 9 parameters\n");
	  fprintf(stderr,"        9: LAL, analytic Likelihood\n");	
      fprintf(stderr,"      Please set injectionWaveform in %s to one of these values.\n\n",run->injectionFilename);
      exit(1);
    }
  }
  
  
  
  // Read injection parameters:
  for(i=1;i<=5;i++) cstatus = fgets(tmpStr,500,fin);  //Read empty and comment lines
  
  for(i=0;i<run->nInjectPar;i++) {
    istatus = fscanf(fin,"%d %d %lf %d %lf %d %lf %lf",&run->injNumber[i],&run->injID[i],&run->injParValOrig[i],&run->injRanPar[i],&run->injSigma[i],&run->injBoundType[i],&run->injBoundLow[i],&run->injBoundUp[i]);
    cstatus = fgets(tmpStr,500,fin);  //Read rest of the line
    
    //printf("%d %d %lf %d %lf %d %lf %lf\n",run->injNumber[i],run->injID[i],run->injParValOrig[i],run->injRanPar[i],run->injSigma[i],run->injBoundType[i],run->injBoundLow[i],run->injBoundUp[i]);
    
    
    // Some sanity checks:
    if(run->injNumber[i] != i+1) {
      fprintf(stderr, "\n\n   ERROR reading injection input file %s:  parameter %d has number %d.\n   Aborting...\n\n",run->injectionFilename,i+1,run->injNumber[i]);
      exit(1);
    }
    
    if(run->parDef[run->injID[i]] != 1) {
      fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d:\n     parameter ID %d is not defined.\n   Aborting...\n\n",
              run->injectionFilename,run->injNumber[i],run->injID[i]);
      exit(1);
    }
    
    // Get reverse parameter ID:
    run->injRevID[run->injID[i]] = i;
    
    // Set this parameter as being used for the injection:
    run->injParUse[run->injID[i]] = 1;
    
    // Get the desired injection boundaries:
    switch (run->injBoundType[i]) {
    case 1 :  // Keep as is: boundLow-boundUp
      break; 
      
    case 2 :  // InjectionValue + BoundLow - InjectionValue + BoundUp
      if(run->injBoundLow[i] > 0.0 || run->injBoundUp[i] < 0.0) {
        fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     for injBoundType = 2, injBoundLow and injBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
                run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]]);
        exit(1);
      }
      run->injBoundLow[i] = run->injParValOrig[i] + run->injBoundLow[i];
      run->injBoundUp[i]  = run->injParValOrig[i] + run->injBoundUp[i];
      break;
      
    case 3 :  // InjectionValue * BoundLow - InjectionValue * BoundUp
      if(run->injBoundLow[i] > 1.0 || run->injBoundUp[i] < 1.0) {
        fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     for injBoundType = 3, injBoundLow and injBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
                run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]]);
        exit(1);
      }
      run->injBoundLow[i] = run->injParValOrig[i] * run->injBoundLow[i];
      run->injBoundUp[i]  = run->injParValOrig[i] * run->injBoundUp[i];
      break;
      
    default :
      fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     %d is not a valid option for injBoundType.\n   Aborting...\n\n",
              run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injBoundType[i]);
      exit(1);
    } //End switch (run->injBoundType[i])
    
    
    
    // More sanity checks:
    
    // Check whether value for injRanPar is valid
    if(run->injRanPar[i] < 0 || run->injRanPar[i] > 2) {
      fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     %d is not a valid option for injRanPar.\n   Aborting...\n\n",
              run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injRanPar[i]);
      exit(1);
    }      
    
    //Check whether the lower boundary < the upper
    if(run->injBoundLow[i] >= run->injBoundUp[i]) {
      fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
              run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]],run->injBoundLow[i],run->injBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower boundary <= injection value <= upper boundary
    if( (run->injParValOrig[i] < run->injBoundLow[i] || run->injParValOrig[i] > run->injBoundUp[i])  &&  run->injRanPar[i] == 0) {
      fprintf(stderr, "\n\n   ERROR reading injection input file %s, parameter %d (%s):\n     the injection value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
              run->injectionFilename,run->injNumber[i],run->parAbrev[run->injID[i]], run->injParValOrig[i], run->injBoundLow[i], run->injBoundUp[i]);
      exit(1);
    }
    
  } //End for (i): injection parameter
  
  
  // Pick random injection values:
  setRandomInjectionParameters(run);    //Copy the injection parameters from injParValOrig to injParVal, and randomise where wanted
  run->geocentricTc = run->injParVal[run->injRevID[11]];    // This value must be overwritten by the 'best' value in readParameterInputfile() which is called next, in the case of no SW injection
  
  fclose(fin);
	}
  
  
  // Print injection parameters and prior ranges to screen:
  if(run->beVerbose>=1) {
    if(run->injectSignal >= 1) {
      printf("\n      %2i software-injection parameters:\n        Nr:  ID: Name:           Injection value:     Obtained:\n",run->nInjectPar);
      for(i=0;i<run->nInjectPar;i++) {
        if(run->injRanPar[i]==0) {
          printf("        %2d  %3i  %-11s     %15.4lf      Taken from the value set in %s\n",i,run->injID[i],run->parAbrev[run->injID[i]],run->injParVal[i],
                 run->injectionFilename);
        } else if(run->injRanPar[i]==1) {
          printf("        %2d  %3i  %-11s     %15.4lf      Drawn randomly from a Gaussian distribution with centre  %lf  and width  %lf  (seed=%i)\n",i,run->injID[i],
                 run->parAbrev[run->injID[i]],run->injParVal[i],run->injParValOrig[i],run->injSigma[i],run->injRanSeed);
        } else if(run->injRanPar[i]==2) {
          printf("        %2d  %3i  %-11s     %15.4lf      Drawn randomly from a uniform distribution  %lf  -  %lf  (seed=%i)\n",i,run->injID[i],run->parAbrev[run->injID[i]]
                 ,run->injParVal[i],run->injBoundLow[i],run->injBoundUp[i],run->injRanSeed);
        }
      }
      printf("\n");
    } else {
      printf("    - not injecting a signal\n");
    }
  }
  
   if(!cstatus||!istatus) {
       printf("Error reading string\n");
       exit(1);
   }  
  
}  //End of readInjectionInputfile
// ****************************************************************************************************************************************************  













// ****************************************************************************************************************************************************  
/**
 * \brief Read the (MCMC) parameter input file
 *
 * All parameters that are read in here should be(come) members of the runvar struct
 */
// ****************************************************************************************************************************************************  
void readParameterInputfile(struct runPar *run)
{
  int i=0,iInj=0,istatus=0;
  char tmpStr[500],*cstatus=NULL;
  FILE *fin;
  int warnings = 0;
  
	if((fin = fopen(run->parameterFilename,"r")) == NULL) {
		fprintf(stderr, "   No parameter input file: %s, using default values.\n",run->parameterFilename);
		
		run->mcmcPNorder = 3.5;
		run->priorSet = 1;
		run->offsetMCMC = 1;
		run->offsetX = 10.0;	
		
		if(run->commandSettingsFlag[17] == 1 && run->mcmcWaveform == 4) {
			
			if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 9-parameter non-spinning waveform as the MCMC template.\n",run->mcmcPNorder);
			run->nMCMCpar=9;
			
			for(i=0;i<run->nMCMCpar;i++) {
				run->parNumber[i] = i+1;
				run->parFix[i] = 0;
			}
			
			run->parID[0] = 65;run->parBestVal[0] = 1.3;		    run->parStartMCMC[0] = 2;run->parSigma[0] = 0.0250;		run->priorType[0] = 13;		run->priorBoundLow[0] = 0.890;		run->priorBoundUp[0] = 1.125;
			run->parID[1] = 62;run->parBestVal[1] = 0.2;		    run->parStartMCMC[1] = 5;run->parSigma[1] = 0.0250;		run->priorType[1] = 11;		run->priorBoundLow[1] = 0.0312;		run->priorBoundUp[1] = 0.25;
			run->parID[2] = 11;run->parBestVal[2] = 0.0;			run->parStartMCMC[2] = 2;run->parSigma[2] = 0.0100;		run->priorType[2] = 12;		run->priorBoundLow[2] = -0.10000;	run->priorBoundUp[2] = 0.10000;
			run->parID[3] = 21;run->parBestVal[3] = 125000.0;		run->parStartMCMC[3] = 2;run->parSigma[3] = 1.0000;		run->priorType[3] = 11;		run->priorBoundLow[3] = 1.0;		run->priorBoundUp[3] = 1000000.0;
			run->parID[4] = 31;run->parBestVal[4] = 0.0;			run->parStartMCMC[4] = 5;run->parSigma[4] = 0.1000;		run->priorType[4] = 21;		run->priorBoundLow[4] = 0.0;		run->priorBoundUp[4] = 6.283185;
			run->parID[5] = 32;run->parBestVal[5] = 0.0;			run->parStartMCMC[5] = 5;run->parSigma[5] = 0.1000;		run->priorType[5] = 11;		run->priorBoundLow[5] = -0.999999;	run->priorBoundUp[5] = 0.999999;
			run->parID[6] = 51;run->parBestVal[6] = 0.0;			run->parStartMCMC[6] = 5;run->parSigma[6] = 0.1000;		run->priorType[6] = 11;		run->priorBoundLow[6] = -0.999999;	run->priorBoundUp[6] = 0.999999;
			run->parID[7] = 41;run->parBestVal[7] = 0.0;			run->parStartMCMC[7] = 5;run->parSigma[7] = 0.1000;		run->priorType[7] = 21;		run->priorBoundLow[7] = 0.0;		run->priorBoundUp[7] = 6.283185;
			run->parID[8] = 52;run->parBestVal[8] = 0.0;			run->parStartMCMC[8] = 5;run->parSigma[8] = 0.1000;		run->priorType[8] = 22;		run->priorBoundLow[8] = 0.0;		run->priorBoundUp[8] = 3.141593;
			
		} else {
			run->mcmcWaveform = 3;
			if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 15-parameter spinning waveform as the MCMC template.\n",run->mcmcPNorder);
			run->nMCMCpar=15;
			
			for(i=0;i<run->nMCMCpar;i++) {
				run->parNumber[i] = i+1;
				run->parFix[i] = 0;
			}
			
			run->parID[0] = 65;run->parBestVal[0] = 1.3;		run->parStartMCMC[0] = 2;run->parSigma[0] = 0.0250;		run->priorType[0] = 13;		run->priorBoundLow[0] = 0.890;		run->priorBoundUp[0] = 1.125;
			run->parID[1] = 62;run->parBestVal[1] = 0.2;		run->parStartMCMC[1] = 5;run->parSigma[1] = 0.0250;		run->priorType[1] = 11;		run->priorBoundLow[1] = 0.0312;		run->priorBoundUp[1] = 0.25;
			run->parID[2] = 11;run->parBestVal[2] = 0.0;			run->parStartMCMC[2] = 2;run->parSigma[2] = 0.0100;		run->priorType[2] = 12;		run->priorBoundLow[2] = -0.10000;	run->priorBoundUp[2] = 0.10000;
			run->parID[3] = 21;run->parBestVal[3] = 125000.0;		run->parStartMCMC[3] = 2;run->parSigma[3] = 1.0000;		run->priorType[3] = 11;		run->priorBoundLow[3] = 1.0;		run->priorBoundUp[3] = 1000000.0;
			run->parID[4] = 31;run->parBestVal[4] = 0.0;			run->parStartMCMC[4] = 5;run->parSigma[4] = 0.1000;		run->priorType[4] = 21;		run->priorBoundLow[4] = 0.0;		run->priorBoundUp[4] = 6.283185;
			run->parID[5] = 32;run->parBestVal[5] = 0.0;			run->parStartMCMC[5] = 5;run->parSigma[5] = 0.1000;		run->priorType[5] = 11;		run->priorBoundLow[5] = -0.999999;	run->priorBoundUp[5] = 0.999999;
			run->parID[6] = 51;run->parBestVal[6] = 0.0;			run->parStartMCMC[6] = 5;run->parSigma[6] = 0.1000;		run->priorType[6] = 11;		run->priorBoundLow[6] = -0.999999;	run->priorBoundUp[6] = 0.999999;
			run->parID[7] = 41;run->parBestVal[7] = 0.0;			run->parStartMCMC[7] = 5;run->parSigma[7] = 0.1000;		run->priorType[7] = 21;		run->priorBoundLow[7] = 0.0;		run->priorBoundUp[7] = 6.283185;
			run->parID[8] = 52;run->parBestVal[8] = 0.0;			run->parStartMCMC[8] = 5;run->parSigma[8] = 0.1000;		run->priorType[8] = 22;		run->priorBoundLow[8] = 0.0;		run->priorBoundUp[8] = 3.141593;
			run->parID[9] = 71;run->parBestVal[9] = 0.000001;		run->parStartMCMC[9] = 5;run->parSigma[9] = 0.1000;		run->priorType[9] = 11;		run->priorBoundLow[9] = 0.000001;	run->priorBoundUp[9] = 0.999999;
			run->parID[10] = 72;run->parBestVal[10] = 0.0;		run->parStartMCMC[10] = 5;run->parSigma[10] = 0.1000;	run->priorType[10] = 11;	run->priorBoundLow[10] = -0.999999;	run->priorBoundUp[10] = 0.999999;
			run->parID[11] = 73;run->parBestVal[11] = 0.0;		run->parStartMCMC[11] = 5;run->parSigma[11] = 0.1000;	run->priorType[11] = 21;	run->priorBoundLow[11] = 0.0;		run->priorBoundUp[11] = 6.283185;
			run->parID[12] = 81;run->parBestVal[12] = 0.000001;	run->parStartMCMC[12] = 5;run->parSigma[12] = 0.1000;	run->priorType[12] = 11;	run->priorBoundLow[12] = 0.000001;	run->priorBoundUp[12] = 0.999999;
			run->parID[13] = 82;run->parBestVal[13] = 0.0;		run->parStartMCMC[13] = 5;run->parSigma[13] = 0.1000;	run->priorType[13] = 11;	run->priorBoundLow[13] = -0.999999;	run->priorBoundUp[13] = 0.999999;
			run->parID[14] = 83;run->parBestVal[14] = 0.0;		run->parStartMCMC[14] = 5;run->parSigma[14] = 0.1000;	run->priorType[14] = 21;	run->priorBoundLow[14] = 0.0;		run->priorBoundUp[14] = 6.283185;
		}
		
		
		
		for(i=0;i<run->nMCMCpar;i++) {
			
			//Check if trigger values from the command line should replace these values:
			if(run->parID[i] == 11 && fabs(run->triggerTc) > 1.e-10) run->parBestVal[i] = run->triggerTc;               // Tc
			if(run->parID[i] == 21 && fabs(run->triggerDist) > 1.e-10) run->parBestVal[i] = pow(run->triggerDist,3.0);  // Distance (d^3)
			if(run->parID[i] == 22 && fabs(run->triggerDist) > 1.e-10) run->parBestVal[i] = log(run->triggerDist);      // Distance log(d)
			if(run->parID[i] == 61 && fabs(run->triggerMc) > 1.e-10) run->parBestVal[i] = run->triggerMc;               // Chirp mass
			if(run->parID[i] == 62 && fabs(run->triggerEta) > 1.e-10) run->parBestVal[i] = run->triggerEta;             // Eta
			if(run->parID[i] == 65 && fabs(run->triggerMc) > 1.e-10) run->parBestVal[i] = pow(run->triggerMc,2*c3rd);
			
			//printf("%d:  %d %d %lf %d %lf %d %lf %lf\n",i,run->parNumber[i],run->parID[i],run->parBestVal[i],run->parStartMCMC[i],run->parSigma[i],
			//run->priorType[i],run->priorBoundLow[i],run->priorBoundUp[i]);
			
			
			if(run->parNumber[i] != i+1) {
				fprintf(stderr, "\n\n   ERROR reading parameter:  parameter %d has number %d.\n   Aborting...\n\n",i+1,run->parNumber[i]);
				exit(1);
			}
			
			if(run->parDef[run->parID[i]] != 1) {
				fprintf(stderr, "\n\n   ERROR reading parameter, parameter %d:\n     parameter ID %d is not defined.\n   Aborting...\n\n",
						run->parNumber[i],run->parID[i]);
				exit(1);
			}
			
			
			
			// Set the reverse parameter ID:
			run->parRevID[run->parID[i]] = i;
			
			// Set this parameter as being used for MCMC:
			run->mcmcParUse[run->parID[i]] = 1;
			
			
			
			// Get the desired boundary conditions:
			switch (run->priorType[i]) {
				case 11 : // General range, BoundLow-BoundUp (as is)
					break;
					
				case 12 : // General range, best value+BoundLow - best value+BoundUp
					if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
						fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 12, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
								run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
						exit(1);
					}
					run->priorBoundLow[i] = run->parBestVal[i] + run->priorBoundLow[i];
					run->priorBoundUp[i]  = run->parBestVal[i] + run->priorBoundUp[i];
					break;
					
				case 13 : // General range, best value*BoundLow - best value*BoundUp
					if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
						fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 13, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
								run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
						exit(1);
					}
					run->priorBoundLow[i] = run->parBestVal[i] * run->priorBoundLow[i];
					run->priorBoundUp[i]  = run->parBestVal[i] * run->priorBoundUp[i];
					break;
					
				case 14 : // General range, injection value+BoundLow - injection value+BoundUp
					if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
						fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 14, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
								run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
						exit(1);
					}
					run->priorBoundLow[i] = run->injParVal[i] + run->priorBoundLow[i];
					run->priorBoundUp[i]  = run->injParVal[i] + run->priorBoundUp[i];
					break;
					
				case 15 : // General range, injection value*BoundLow - injection value*BoundUp
					if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
						fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 15, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
								run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
						exit(1);
					}
					run->priorBoundLow[i] = run->injParVal[i] * run->priorBoundLow[i];
					run->priorBoundUp[i]  = run->injParVal[i] * run->priorBoundUp[i];
					break;
					
				case 21 : // Periodic boundaries, 0-2pi
					run->priorBoundLow[i] = 0.0;
					run->priorBoundUp[i]  = tpi;
					break;
					
				case 22 : // Periodic boundaries, 0-pi
					run->priorBoundLow[i] = 0.0;
					run->priorBoundUp[i]  = pi;
					break;
					
				default :
					fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for priorType.\n   Aborting...\n\n",
							run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorType[i]);
					exit(1);
			} //End switch
			
			
			// Check whether value for fix is valid
			if(run->parFix[i] < 0 || run->parFix[i] > 2) {
				fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parFix.\n   Aborting...\n\n",
						run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parFix[i]);
				exit(1);
			}      
			
			// Check whether value for start is valid
			if(run->parStartMCMC[i] < 1 || run->parStartMCMC[i] > 5) {
				fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parStartMCMC.\n   Aborting...\n\n",
						run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parStartMCMC[i]);
				exit(1);
			}
			if((run->parStartMCMC[i] == 3 || run->parStartMCMC[i] == 4) && run->injectSignal <= 0) {
				fprintf(stdout, "    - no software injection was performed, so I'll use bestValue rather than injectionValue for %s.\n",run->parAbrev[run->parID[i]]);
				run->parStartMCMC[i] -= 2;
			}
			
			
			//Check whether the lower prior boundary < the upper
			if(run->priorBoundLow[i] >= run->priorBoundUp[i]) {
				fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
						run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorBoundLow[i],run->priorBoundUp[i]);
				exit(1);
			}
			
			//Check whether  lower prior boundary <= best value <= upper boundary
			if( (run->parBestVal[i] < run->priorBoundLow[i] || run->parBestVal[i] > run->priorBoundUp[i])  && run->parStartMCMC[i] == 1 ) {
				fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the best value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
						run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]], run->parBestVal[i], run->priorBoundLow[i], run->priorBoundUp[i]);
				exit(1);
			}
			
			//Check whether  lower prior boundary <= INJECTION value <= upper boundary
			iInj = run->injRevID[run->parID[i]];  //Get the index of this parameter in the injection set.  -1 if not available.
			if(iInj >= 0) {
				if(run->injParVal[iInj] < run->priorBoundLow[i] || run->injParVal[iInj] > run->priorBoundUp[i]) {
					fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the injection value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
							run->parameterFilename, run->parNumber[i], run->parAbrev[run->parID[i]], run->injParVal[iInj], run->priorBoundLow[i], run->priorBoundUp[i]);
					exit(1);
				}
			} else {
				if(run->injectSignal != 0) {
					if(warnings==0) fprintf(stderr, "\n");
					printf("    * Warning:  MCMC parameter %i (%s) does not occur in the injection template;  I cannot verify whether the injection value lies within the prior range *\n",
						   run->parNumber[i],run->parAbrev[run->parID[i]]);
					warnings += 1;
				}
			} 
		} //End for (i)
		
		if(run->injectSignal<=0) {
			run->geocentricTc = run->parBestVal[run->parRevID[11]];    // This value overwrites the injection value from readInjectionInputfile(), in the case of no SW injection
			for(i=0;i<run->nMCMCpar;i++) run->injParVal[i] = run->parBestVal[i];   //CHECK Needed to avoid SegFault in the case of t_c
		}
		
		
		//exit(1);
	} else {
    printf("   Using parameter input file: %s.\n",run->parameterFilename);

  
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  
  for(i=1;i<=2;i++) cstatus = fgets(tmpStr,500,fin);  //Read first 2 lines
  
  //Priors:
  cstatus = fgets(tmpStr,500,fin); cstatus = fgets(tmpStr,500,fin);  //Read the empty and comment line
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->mcmcWaveform);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->mcmcPNorder);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->priorSet);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%d",&run->offsetMCMC);
  cstatus = fgets(tmpStr,500,fin);  sscanf(tmpStr,"%lf",&run->offsetX);
  
  if(run->mcmcWaveform==1) {
    if(run->beVerbose>=1) printf("    - using Apostolatos, 1.5-pN, 12-parameter waveform as the MCMC template.\n");  //Only 1.5-pN order is available
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==2) {
    if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 12-parameter waveform as the MCMC template.\n",run->mcmcPNorder);
    run->nMCMCpar=12;
  } else if(run->mcmcWaveform==3) {
    if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 15-parameter waveform as the MCMC template.\n",run->mcmcPNorder);
    run->nMCMCpar=15;
  } else if(run->mcmcWaveform==4) {
    if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, non-spinning waveform as the MCMC template.\n",run->mcmcPNorder);
    run->nMCMCpar=9;
  } else if(run->mcmcWaveform==5) {
	if(run->beVerbose>=1) printf("    - using LAL,%4.1f-pN, 15-parameter PhenSpinTaylorRD waveform as the MCMC template.\n",run->mcmcPNorder);
	run->nMCMCpar=15;	  
  } else if(run->mcmcWaveform==9) {
	if(run->beVerbose>=1) printf("    - using analytic Likelihood.\n");
	run->nMCMCpar=15;			  
  } else {
    fprintf(stderr,"    - unknown waveform chosen as MCMC template: %d.   Available waveforms are:\n",run->mcmcWaveform);
    fprintf(stderr,"        1: Apostolatos, simple precession, 12 parameters\n");
    fprintf(stderr,"        2: LAL, single spin, 12 parameters\n");
    fprintf(stderr,"        3: LAL, double spin, 15 parameters\n");
    fprintf(stderr,"        4: LAL, non-spinnig, 9 parameters\n");
	fprintf(stderr,"        9: LAL, analytic Likelihood\n");
    fprintf(stderr,"      Please set mcmcWaveform in %s to one of these values.\n\n",run->mcmcFilename);
    exit(1);
  }
  
  
  
  
  //Parameters:
  for(i=1;i<=5;i++) cstatus = fgets(tmpStr,500,fin);  //Read empty and comment lines
  
	
  for(i=0;i<run->nMCMCpar;i++) {
    istatus = fscanf(fin,"%d %d %lf %d %d %lf %d %lf %lf",&run->parNumber[i],&run->parID[i],&run->parBestVal[i],&run->parFix[i],&run->parStartMCMC[i],&run->parSigma[i],&run->priorType[i],&run->priorBoundLow[i],&run->priorBoundUp[i]);
    cstatus = fgets(tmpStr,500,fin);  //Read rest of the line
    
    //Check if trigger values from the command line should replace these values:
    if(run->parID[i] == 11 && fabs(run->triggerTc) > 1.e-10) run->parBestVal[i] = run->triggerTc;               // Tc
    if(run->parID[i] == 21 && fabs(run->triggerDist) > 1.e-10) run->parBestVal[i] = pow(run->triggerDist,3.0);  // Distance (d^3)
    if(run->parID[i] == 22 && fabs(run->triggerDist) > 1.e-10) run->parBestVal[i] = log(run->triggerDist);      // Distance log(d)
    if(run->parID[i] == 61 && fabs(run->triggerMc) > 1.e-10) run->parBestVal[i] = run->triggerMc;               // Chirp mass
    if(run->parID[i] == 62 && fabs(run->triggerEta) > 1.e-10) run->parBestVal[i] = run->triggerEta;             // Eta
    
    //printf("%d:  %d %d %lf %d %lf %d %lf %lf\n",i,run->parNumber[i],run->parID[i],run->parBestVal[i],run->parStartMCMC[i],run->parSigma[i],
    //run->priorType[i],run->priorBoundLow[i],run->priorBoundUp[i]);
    
    
    if(run->parNumber[i] != i+1) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s:  parameter %d has number %d.\n   Aborting...\n\n",run->parameterFilename,i+1,run->parNumber[i]);
      exit(1);
    }
    
    if(run->parDef[run->parID[i]] != 1) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d:\n     parameter ID %d is not defined.\n   Aborting...\n\n",
              run->injectionFilename,run->parNumber[i],run->parID[i]);
      exit(1);
    }
    
    
    
    // Set the reverse parameter ID:
    run->parRevID[run->parID[i]] = i;
    
    // Set this parameter as being used for MCMC:
    run->mcmcParUse[run->parID[i]] = 1;
    
    
    
    // Get the desired boundary conditions:
    switch (run->priorType[i]) {
    case 11 : // General range, BoundLow-BoundUp (as is)
      break;
      
    case 12 : // General range, best value+BoundLow - best value+BoundUp
      if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
        fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 12, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
                run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
        exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] + run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] + run->priorBoundUp[i];
      break;
      
    case 13 : // General range, best value*BoundLow - best value*BoundUp
      if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
        fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 13, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
                run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
        exit(1);
      }
      run->priorBoundLow[i] = run->parBestVal[i] * run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->parBestVal[i] * run->priorBoundUp[i];
      break;
      
    case 14 : // General range, injection value+BoundLow - injection value+BoundUp
      if(run->priorBoundLow[i] > 0.0 || run->priorBoundUp[i] < 0.0) {
        fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 14, priorBoundLow and priorBoundUp must be <= 0 and >= 0 respectively.\n   Aborting...\n\n",
                run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
        exit(1);
      }
      run->priorBoundLow[i] = run->injParVal[i] + run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->injParVal[i] + run->priorBoundUp[i];
      break;
      
    case 15 : // General range, injection value*BoundLow - injection value*BoundUp
      if(run->priorBoundLow[i] > 1.0 || run->priorBoundUp[i] < 1.0) {
        fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     for priorType = 15, priorBoundLow and priorBoundUp must be <= 1 and >= 1 respectively.\n   Aborting...\n\n",
                run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]]);
        exit(1);
      }
      run->priorBoundLow[i] = run->injParVal[i] * run->priorBoundLow[i];
      run->priorBoundUp[i]  = run->injParVal[i] * run->priorBoundUp[i];
      break;
      
    case 21 : // Periodic boundaries, 0-2pi
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = tpi;
      break;
      
    case 22 : // Periodic boundaries, 0-pi
      run->priorBoundLow[i] = 0.0;
      run->priorBoundUp[i]  = pi;
      break;
      
    default :
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for priorType.\n   Aborting...\n\n",
              run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorType[i]);
      exit(1);
    } //End switch
    
    
    // Check whether value for fix is valid
    if(run->parFix[i] < 0 || run->parFix[i] > 2) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parFix.\n   Aborting...\n\n",
              run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parFix[i]);
      exit(1);
    }      
    
    // Check whether value for start is valid
    if(run->parStartMCMC[i] < 1 || run->parStartMCMC[i] > 5) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     %d is not a valid option for parStartMCMC.\n   Aborting...\n\n",
              run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->parStartMCMC[i]);
      exit(1);
    }
    if((run->parStartMCMC[i] == 3 || run->parStartMCMC[i] == 4) && run->injectSignal <= 0) {
      fprintf(stdout, "    - no software injection was performed, so I'll use bestValue rather than injectionValue for %s.\n",run->parAbrev[run->parID[i]]);
      run->parStartMCMC[i] -= 2;
    }
    
    
    //Check whether the lower prior boundary < the upper
    if(run->priorBoundLow[i] >= run->priorBoundUp[i]) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the lower boundary of the prior is larger than or equal to the upper boundary (%lf vs. %lf).\n   Aborting...\n\n",
              run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]],run->priorBoundLow[i],run->priorBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower prior boundary <= best value <= upper boundary
    if( (run->parBestVal[i] < run->priorBoundLow[i] || run->parBestVal[i] > run->priorBoundUp[i])  && run->parStartMCMC[i] == 1 ) {
      fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the best value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
              run->parameterFilename,run->parNumber[i],run->parAbrev[run->parID[i]], run->parBestVal[i], run->priorBoundLow[i], run->priorBoundUp[i]);
      exit(1);
    }
    
    //Check whether  lower prior boundary <= INJECTION value <= upper boundary
    iInj = run->injRevID[run->parID[i]];  //Get the index of this parameter in the injection set.  -1 if not available.
    if(iInj >= 0) {
      if(run->injParVal[iInj] < run->priorBoundLow[i] || run->injParVal[iInj] > run->priorBoundUp[i]) {
        fprintf(stderr, "\n\n   ERROR reading parameter input file %s, parameter %d (%s):\n     the injection value (%lf) lies outside the prior range (%lf - %lf).\n   Aborting...\n\n",
                run->parameterFilename, run->parNumber[i], run->parAbrev[run->parID[i]], run->injParVal[iInj], run->priorBoundLow[i], run->priorBoundUp[i]);
        exit(1);
      }
    } else {
      if(run->injectSignal != 0) {
        if(warnings==0) fprintf(stderr, "\n");
        printf("    * Warning:  MCMC parameter %i (%s) does not occur in the injection template;  I cannot verify whether the injection value lies within the prior range *\n",
               run->parNumber[i],run->parAbrev[run->parID[i]]);
        warnings += 1;
      }
    } 
  } //End for (i)
  //if(warnings >0) printf("\n");
  
  
  if(run->injectSignal<=0) {
    run->geocentricTc = run->parBestVal[run->parRevID[11]];    // This value overwrites the injection value from readInjectionInputfile(), in the case of no SW injection
    for(i=0;i<run->nMCMCpar;i++) run->injParVal[i] = run->parBestVal[i];   //CHECK Needed to avoid SegFault in the case of t_c
  }
  
  fclose(fin);
	}
  
  
  //Print MCMC parameters and prior ranges to screen:
  char FixStr[3][99];
  strcpy(FixStr[0],"No, let it free");
  strcpy(FixStr[1],"Yes, to best value");
  strcpy(FixStr[2],"Yes, to injection");
  
  char StartStr[6][99];
  strcpy(StartStr[1],"From the best value");
  strcpy(StartStr[2],"Randomly from a Gaussian around the best value");
  strcpy(StartStr[3],"From the injection value");
  strcpy(StartStr[4],"Randomly from a Gaussian around the injection");
  strcpy(StartStr[5],"Randomly from the prior");
  
  
  //Write parameter choice to screen:
  if(run->beVerbose>=1) {
    printf("\n      %2i MCMC parameters:\n        Nr:  ID: Name:                Best value:     Prior:     min:            max:    Fix parameter?        Start chain:\n",run->nMCMCpar);
    for(i=0;i<run->nMCMCpar;i++) {
      printf("        %2d  %3i  %-11s     %15.4lf     %15.4lf %15.4lf     %-20s  %-45s",i,run->parID[i],run->parAbrev[run->parID[i]],run->parBestVal[i],
             run->priorBoundLow[i],run->priorBoundUp[i],  FixStr[run->parFix[i]],StartStr[run->parStartMCMC[i]]);
      if(run->parStartMCMC[i]==2 || run->parStartMCMC[i]==4 || run->parStartMCMC[i]==5) {
	printf(" (seed=%i)\n",run->MCMCseed);
      } else {
	printf("\n");
      }
     
    }
    printf("\n");
  }
  
   if(!cstatus||!istatus) {
       printf("Error reading string\n");
       exit(1);
   }  
  
}  //End of readParameterInputfile
// ****************************************************************************************************************************************************  








// ****************************************************************************************************************************************************  
/**
 * \brief Read the input file for system (system-dependent) variables, e.g. SPINspiral.input.system
 */
// ****************************************************************************************************************************************************  
void readSystemInputfile(struct runPar *run)
{
  int i=0,istatus=0;
  char tmpStr[500],*cstatus=NULL;
  FILE *fin;
  
	if((fin = fopen(run->systemFilename,"r")) == NULL) {
		fprintf(stderr, "   No system file: %s.\n",run->systemFilename);
		sprintf(run->dataDir,"/");
		//exit(1);
	} else {
    printf("   Using system input file: %s.\n",run->parameterFilename);
  
  //Use and l for floats: %lf, %lg, etc, rather than %f, %g
  for(i=1;i<=3;i++) { //Read first 3 lines
    cstatus = fgets(tmpStr,500,fin);
  }  
  
  //Data directory:
  istatus = fscanf(fin, "%s",run->dataDir);
  
  fclose(fin);
	}
  
   if(!cstatus||!istatus) {
       printf("Error reading string\n");
       exit(1);
   } 
    
}  //End of readSystemInputfile
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Read an XML injection file if specified
 */
// ****************************************************************************************************************************************************  
void readInjectionXML(struct runPar *run)
{
  int NumInj=0,i=0,j=0;
  SimInspiralTable *injTable = NULL;
  double Sx=0.0,Sy=0.0,Sz=0.0,Sxy=0.0,S=0.0;
  
  printf("  Reading injection XML file %s, injection %d.\n",run->injXMLfilename,run->injXMLnr);
  NumInj = SimInspiralTableFromLIGOLw(&injTable,run->injXMLfilename,0,0);
  if(NumInj <= 0) {
    fprintf(stderr,"\n\n   ERROR opening XML file %s.\n   Aborting.\n\n",run->injXMLfilename);
    exit(1);
  }
  if(run->injXMLnr >= NumInj) {
    fprintf(stderr,"\n\n  ERROR: requested injection number %d larger than number of injections (%d) in XML file %s.\n  Aborting.\n\n",run->injXMLnr,NumInj-1,run->injXMLfilename);
    exit(1);
  }
  
  j=0;
  while(j<run->injXMLnr) {j++; injTable = injTable->next;}  // Select injection
  
  i=0;
  for(i=0;i<run->nInjectPar;i++) {
    
    //Time:
    if(run->injID[i] == 11) run->injParVal[i] = injTable->geocent_end_time.gpsSeconds+injTable->geocent_end_time.gpsNanoSeconds*1.0e-9;  // t_c in geocentre
    
    
    //Distance:
    if(run->injID[i] == 21) run->injParVal[i] = pow(injTable->distance,3.0);  // d_L^3
    if(run->injID[i] == 22) run->injParVal[i] = log(injTable->distance);      // log(d_L/Mpc)
    
    
    //Sky location:
    if(run->injID[i] == 31) run->injParVal[i] = injTable->longitude;          // RA
    if(run->injID[i] == 32) run->injParVal[i] = sin(injTable->latitude);      // sin Dec
    
    
    //Orbital/GW phase:
    if(run->injID[i] == 41) run->injParVal[i] = injTable->coa_phase;          // \phi_c
    
    
    //Binary orientation:
    if(run->injID[i] == 51) run->injParVal[i] = cos(injTable->inclination);   // cos(i)
    if(run->injID[i] == 52) run->injParVal[i] = injTable->polarization;       // \psi
    if(run->injID[i] == 53) {
      printf("\n    * You're reading injection data from an XML file, while using the Apostolatos J_0.  This may lead to unexpected and unwanted results *\n");
      run->injParVal[i] = sin(injTable->inclination);   // sin(\theta_J0), should probably not be used
    }
    if(run->injID[i] == 54) {
      printf("\n    * You're reading injection data from an XML file, while using the Apostolatos J_0.  This may lead to unexpected and unwanted results *\n");
      run->injParVal[i] = injTable->polarization;       // \phi_J0, should probably not be used
    }
    
    
    //Mass:
    if(run->injID[i] == 61) run->injParVal[i] = injTable->mchirp;             // Mc
    if(run->injID[i] == 62) run->injParVal[i] = injTable->eta;                // \eta
    if(run->injID[i] == 63) run->injParVal[i] = injTable->mass1;              // M1
    if(run->injID[i] == 64) run->injParVal[i] = injTable->mass2;              // M2
    
    
    //Spin 1:
    Sx = injTable->spin1x;
    Sy = injTable->spin1y;
    Sz = injTable->spin1z;
    Sxy = sqrt(Sx*Sx+Sy*Sy);
    S = sqrt(Sxy*Sxy+Sz*Sz);
    if(run->injID[i] == 71) run->injParVal[i] = S;                            // Spin1 magnitude (a_1)
    if(run->injID[i] == 72) {
      run->injParVal[i] = Sz/S;                                               // cos(\theta_1): angle between S1 and L (at t_c(?)))
      if(fabs(S)<1.0e-9) run->injParVal[i] = 1.0;                             // In case S1=0
    }
    if(run->injID[i] == 73) run->injParVal[i] = atan2(Sy,Sx);                 // \phi_1: phase of spin in orbital plane (at t_c(?))
    
    if(run->injID[i] == 75) run->injParVal[i] = Sx;                           // Spin1_x magnitude
    if(run->injID[i] == 76) run->injParVal[i] = Sy;                           // Spin1_y magnitude
    if(run->injID[i] == 77) run->injParVal[i] = Sz;                           // Spin1_z magnitude
    
    
    //Spin 2:
    Sx = injTable->spin2x;
    Sy = injTable->spin2y;
    Sz = injTable->spin2z;
    Sxy = sqrt(Sx*Sx+Sy*Sy);
    S = sqrt(Sxy*Sxy+Sz*Sz);
    if(run->injID[i] == 81) run->injParVal[i] = S;                            // Spin2 magnitude (a_2)
    if(run->injID[i] == 82) {
      run->injParVal[i] = Sz/S;                                               // cos(\theta_2): angle between S2 and L (at t_c(?)))
      if(fabs(S)<1.0e-9) run->injParVal[i] = 1.0;                             // In case S2=0
    }
    if(run->injID[i] == 83) run->injParVal[i] = atan2(Sy,Sx);                 // \phi_2: phase of spin in orbital plane (at t_c(?))
    
    if(run->injID[i] == 85) run->injParVal[i] = Sx;                           // Spin2_x magnitude
    if(run->injID[i] == 86) run->injParVal[i] = Sy;                           // Spin2_y magnitude
    if(run->injID[i] == 87) run->injParVal[i] = Sz;                           // Spin2_z magnitude
    
    
    //Merger, ringdown, ...:
    //if(run->injID[i] == 91) run->injParVal[i] = injTable->;                   // 
    
    printf("    - %2d  %-11s     %15.4lf\n",run->injNumber[i],run->parAbrev[run->injID[i]],run->injParVal[i]);
    
  } // i (injectPar)
  
  run->geocentricTc = run->injParVal[run->injRevID[11]];    // This value must be overwritten by the 'best' value in readParameterInputfile() which is called next, in the case of no SW injection
  run->lowFrequencyCutInj = injTable->f_lower;  // May be 0.0!
  
  printf("\n");
  
} // End void readInjectionXML()
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Read a Cache file. Returns an array of what is in the cache file.
 *
 */
// ****************************************************************************************************************************************************  
void readCachefile(struct runPar *run, int ifonr)
{
  int i;
  int line=0;
  char tmpStr[2048],*cstatus=NULL;
  FILE *fin;
  
  if((fin = fopen(run->cacheFilename[ifonr],"r")) == NULL) {
    fprintf(stderr, "\n\n   ERROR opening cache file: %s, aborting.\n\n\n",run->cacheFilename[ifonr]);
    exit(1);
  } else {
    printf("   Reading cache file: %s.\n",run->cacheFilename[ifonr]);
  }
  
  while ( ! feof (fin) ) //just to get the number of line. TO CHECK : last line of .cache file always empty ?
    {
      cstatus = fgets(tmpStr , 2048 , fin);
      line++;
    }
  fclose (fin); 
  
  run->nFrame[ifonr] = line - 1;
  
  run->FrameDetector[ifonr]  = (char**)  malloc(sizeof(char*) * (line));
  for (i=0; i<(line); ++i) (run->FrameDetector[ifonr])[i] = (char*) malloc(sizeof(char)*5);
  run->FramePrefix[ifonr]  = (char**)  malloc(sizeof(char*) * (line));
  for (i=0; i<(line); ++i) (run->FramePrefix[ifonr])[i] = (char*) malloc(sizeof(char)*512);
  run->FrameGPSstart[ifonr] = (int*) malloc(sizeof(int)* (line));
  run->FrameLength[ifonr] = (int*) malloc(sizeof(int)* (line));
  run->FrameName[ifonr]  = (char**)  malloc(sizeof(char*) * (line));
  for (i=0; i<(line); ++i) (run->FrameName[ifonr])[i] = (char*) malloc(sizeof(char)*512);
  
  fin = fopen(run->cacheFilename[ifonr],"r");
  for(i=0;i<(line-1);i++) {
    //Read line by line:
    cstatus = fgets(tmpStr,2048,fin); sscanf(tmpStr,"%s %s %d %d %s",run->FrameDetector[ifonr][i],run->FramePrefix[ifonr][i],&(run->FrameGPSstart[ifonr][i]),&(run->FrameLength[ifonr][i]),run->FrameName[ifonr][i]);
    //      printf("%s %s %d %d %s %d %d\n",run->FrameDetector[ifonr][i],run->FramePrefix[ifonr][i],run->FrameGPSstart[ifonr][i],run->FrameLength[ifonr][i],run->FrameName[ifonr][i],i,run->nFrame[ifonr]);
    
    //remove file://localhost at the beginning of the file name if present.
    if(strncmp("file://localhost",run->FrameName[ifonr][i],16)==0)
      {
	run->FrameName[ifonr][i] = &(run->FrameName[ifonr][i][16]);
      }
    
    
  }
  
  fclose(fin);
  
   if(!cstatus) {
       printf("Error reading string\n");
       exit(1);
   } 
    
}  //End of readCachefile
// ****************************************************************************************************************************************************  









// ****************************************************************************************************************************************************  
/**
 * \brief Get random values for the injection parameters.
 */
// ****************************************************************************************************************************************************  
void setRandomInjectionParameters(struct runPar *run)  
{
  int i=0;
  gsl_rng *ran;
  double ranGauss = 0.0, ranUnif=0.0, db=0.0;
  ran = gsl_rng_alloc(gsl_rng_mt19937);  // GSL random-number seed
  
  // Manually select a random seed, *** USE ONLY FOR TESTING ***
  if(1==1 && run->injRanSeed == 0) {
    printf("\n  *** SELECTING RANDOM INJECTION SEED ***  This should only be done while testing!!! setRandomInjectionParameters() \n\n");
    run->injRanSeed = 0;
    setSeed(&run->injRanSeed);
    printf("   Picking seed from the system clock for random injection parameters: %d\n", run->injRanSeed);
  }
  
  gsl_rng_set(ran, run->injRanSeed);     // Set seed
  
  for(i=0;i<run->nInjectPar;i++) {
    ranGauss = gsl_ran_gaussian(ran,run->injSigma[i]);                                    //Make sure you always draw the same number of random variables
    ranUnif = gsl_rng_uniform(ran);                                                       //Make sure you always draw the same number of random variables
    if(run->injRanPar[i]==0) {                  
      run->injParVal[i] = run->injParValOrig[i];                                          //Keep the suggested value
    } else if(run->injRanPar[i]==1) {                                                     
      run->injParVal[i] = run->injParValOrig[i] + ranGauss;                               //Draw random number from Gaussian with width ranGauss
      run->injParVal[i] = max(run->injParVal[i],run->injBoundLow[i]);                     //Stick to the boundary, rather than redrawing to keep number of random numbers constant
      run->injParVal[i] = min(run->injParVal[i],run->injBoundUp[i]);
    } else if(run->injRanPar[i]==2) {
      db = run->injBoundUp[i]-run->injBoundLow[i];                                        //Width of the range
      run->injParVal[i] = run->injBoundLow[i] + ranUnif*db;                               //Draw random number from uniform range with width db
    }
  }
  
  gsl_rng_free(ran);
} // End setRandomInjectionParameters()
// ****************************************************************************************************************************************************  





// ****************************************************************************************************************************************************  
/**
 * \brief Set parameters names in the hardcoded parameter database
 */
// ****************************************************************************************************************************************************  
void setParameterNames(struct runPar * run)
{
  
  //Set 01: time
  strcpy(run->parAbrev[11], "t_c");
  strcpy(run->parAbrv[11], "t_c");
  run->parDef[11] = 1;
  strcpy(run->parAbrev[12], "t_40");
  strcpy(run->parAbrv[12], "t_40");
  run->parDef[12] = 1;
  
  //Set 02: distance
  strcpy(run->parAbrev[21], "d^3");
  strcpy(run->parAbrv[21], "d^3");
  run->parDef[21] = 1;
  strcpy(run->parAbrev[22], "log(d)");
  strcpy(run->parAbrv[22], "logD");
  run->parDef[22] = 1;
  
  //Set 03: sky position
  strcpy(run->parAbrev[31], "R.A.");
  strcpy(run->parAbrv[31], "RA");
  run->parDef[31] = 1;
  strcpy(run->parAbrev[32], "sin(dec)");
  strcpy(run->parAbrv[32], "sdec");
  run->parDef[32] = 1;
  
  //Set 04: phase
  strcpy(run->parAbrev[41], "phi_orb");
  strcpy(run->parAbrv[41], "phio");
  run->parDef[41] = 1;
  
  //Set 05: orientation
  strcpy(run->parAbrev[51], "cos(i)");
  strcpy(run->parAbrv[51], "cosi");
  run->parDef[51] = 1;
  strcpy(run->parAbrev[52], "psi");
  strcpy(run->parAbrv[52], "psi");
  run->parDef[52] = 1;
  strcpy(run->parAbrev[53], "sin_th_J0");
  strcpy(run->parAbrv[53], "thJ0");
  run->parDef[53] = 1;
  strcpy(run->parAbrev[54], "phi_J0");
  strcpy(run->parAbrv[54], "phJ0");
  run->parDef[54] = 1;
  
  //Set 06: mass
  strcpy(run->parAbrev[61], "Mc");
  strcpy(run->parAbrv[61], "Mc");
  run->parDef[61] = 1;
  strcpy(run->parAbrev[62], "eta");
  strcpy(run->parAbrv[62], "eta");
  run->parDef[62] = 1;
  strcpy(run->parAbrev[63], "M1");
  strcpy(run->parAbrv[63], "M1");
  run->parDef[63] = 1;
  strcpy(run->parAbrev[64], "M2");
  strcpy(run->parAbrv[64], "M2");
  run->parDef[64] = 1;
  strcpy(run->parAbrev[65], "Mc_16");
  strcpy(run->parAbrv[65], "Mc_16");
  run->parDef[65] = 1;
	
	
  //Set 07: spin1
  strcpy(run->parAbrev[71], "a_spin1");
  strcpy(run->parAbrv[71], "asp1");
  run->parDef[71] = 1;
  strcpy(run->parAbrev[72], "cs_th_sp1");
  strcpy(run->parAbrv[72], "ths1");
  run->parDef[72] = 1;
  strcpy(run->parAbrev[73], "phi_spin1");
  strcpy(run->parAbrv[73], "phs1");
  run->parDef[73] = 1;
  
  strcpy(run->parAbrev[75], "S1_x");
  strcpy(run->parAbrv[75], "S1x");
  run->parDef[75] = 1;
  strcpy(run->parAbrev[76], "S1_y");
  strcpy(run->parAbrv[76], "s1y");
  run->parDef[76] = 1;
  strcpy(run->parAbrev[77], "S1_z");
  strcpy(run->parAbrv[77], "S1z");
  run->parDef[77] = 1;
  
  //Set 08: spin2
  strcpy(run->parAbrev[81], "a_spin2");
  strcpy(run->parAbrv[81], "asp2");
  run->parDef[81] = 1;
  strcpy(run->parAbrev[82], "cs_th_sp2");
  strcpy(run->parAbrv[82], "ths2");
  run->parDef[82] = 1;
  strcpy(run->parAbrev[83], "phi_spin2");
  strcpy(run->parAbrv[83], "phs2");
  run->parDef[83] = 1;
  
  strcpy(run->parAbrev[85], "S2_x");
  strcpy(run->parAbrv[85], "S2x");
  run->parDef[85] = 1;
  strcpy(run->parAbrev[86], "S2_y");
  strcpy(run->parAbrv[86], "s2y");
  run->parDef[86] = 1;
  strcpy(run->parAbrev[87], "S2_z");
  strcpy(run->parAbrv[87], "S2z");
  run->parDef[87] = 1;
  
  //Set 09: merger, ringdown
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
	//Set 10: analytic
	strcpy(run->parAbrev[185], "x_1");
	strcpy(run->parAbrv[185], "x1");
	run->parDef[185] = 1;		
	strcpy(run->parAbrev[186], "x_2");
	strcpy(run->parAbrv[186], "x2");
    run->parDef[186] = 1;	
	strcpy(run->parAbrev[187], "x_3");
	strcpy(run->parAbrv[187], "x3");
	run->parDef[187] = 1;		
	strcpy(run->parAbrev[188], "x_4");
	strcpy(run->parAbrv[188], "x4");
	run->parDef[188] = 1;	
	strcpy(run->parAbrev[189], "x_5");
	strcpy(run->parAbrv[189], "x5");
	run->parDef[189] = 1;		
	strcpy(run->parAbrev[190], "x_6");
	strcpy(run->parAbrv[190], "x6");
	run->parDef[190] = 1;	
	strcpy(run->parAbrev[191], "x_7");
	strcpy(run->parAbrv[191], "x7");
	run->parDef[191] = 1;		
	strcpy(run->parAbrev[192], "x_8");
	strcpy(run->parAbrv[192], "x8");
	run->parDef[192] = 1;	
	strcpy(run->parAbrev[193], "x_9");
	strcpy(run->parAbrv[193], "x9");
	run->parDef[193] = 1;		
	strcpy(run->parAbrev[194], "x_10");
	strcpy(run->parAbrv[194], "x10");
	run->parDef[194] = 1;	
	strcpy(run->parAbrev[195], "x_11");
	strcpy(run->parAbrv[195], "x11");
	run->parDef[195] = 1;		
	strcpy(run->parAbrev[196], "x_12");
	strcpy(run->parAbrv[196], "x12");
	run->parDef[196] = 1;	
	strcpy(run->parAbrev[197], "x_13");
	strcpy(run->parAbrv[197], "x13");
	run->parDef[197] = 1;		
	strcpy(run->parAbrev[198], "x_14");
	strcpy(run->parAbrv[198], "x14");
	run->parDef[198] = 1;	
	strcpy(run->parAbrev[199], "x_15");
	strcpy(run->parAbrv[199], "x15");
	run->parDef[199] = 1;		
	
  //Set:
  //strcpy(run->parAbrev[], "");
  //run->parDef[] = 1;
  
} // End setParameterNames()
// ****************************************************************************************************************************************************  













// ****************************************************************************************************************************************************  
/**
 * \brief Copy some of the elements of the struct runPar to the struct MCMCvariables
 */
// ****************************************************************************************************************************************************  
void copyRun2MCMC(struct runPar run, struct MCMCvariables *mcmc)
{
  int i=0,j=0;
  
  mcmc->beVerbose = run.beVerbose;                      // Print long stretches of output
  
  mcmc->maxnPar = run.maxnPar;                          // Absolute maximum number of mcmc/template parameters allowed
  mcmc->parDBn = run.parDBn;                            // Number of elements in hardcoded parameter database
  mcmc->nMCMCpar = run.nMCMCpar;                        // Number of mcmc/template parameters
  mcmc->nInjectPar = run.nInjectPar;                    // Number of mcmc/template parameters
  
  mcmc->mcmcWaveform = run.mcmcWaveform;                // Waveform used as MCMC template
  mcmc->mcmcPNorder = run.mcmcPNorder;                  // pN order used for MCMC template
  mcmc->injectSignal = run.injectSignal;                // Inject a signal in the data or not
  mcmc->injectionWaveform = run.injectionWaveform;      // Waveform used as injection template
  mcmc->injectionPNorder = run.injectionPNorder;        // pN order used for injection template
  mcmc->networkSize = run.networkSize;                  // Network size
  mcmc->seed = run.MCMCseed;                            // MCMC seed
  mcmc->blockFrac = run.blockFrac;                      // Fraction of non-correlated updates that is a block update
  mcmc->corrFrac = run.corrFrac;                        // Fraction of MCMC updates that used the correlation matrix
  mcmc->matAccFr = run.matAccFr;                        // Fraction of elements on the diagonal that must 'improve' in order to accept a new covariance matrix
  mcmc->offsetMCMC = run.offsetMCMC;                    // Start MCMC offset (i.e., not from injection values) or not
  mcmc->offsetX = run.offsetX;                          // Start offset chains from a Gaussian distribution offsetX times wider than parSigma
  
  mcmc->nIter = run.nIter;                              // Number of MCMC iterations to compute
  mcmc->thinOutput = run.thinOutput;                    // Save every thiOutput-th MCMC iteration to file
  mcmc->thinScreenOutput = run.thinScreenOutput;        // Save every thiOutput-th MCMC iteration to screen
  mcmc->adaptiveMCMC = run.adaptiveMCMC;                // Use adaptive MCMC
  mcmc->acceptRateTarget = run.acceptRateTarget;        // Target acceptance rate for MCMC
  mcmc->minlogL = run.minlogL;                          // Minimum value for the log Likelihood to accept a jump. We used 0 for a long time, this number shouldn't be positive!
  
  
  mcmc->correlatedUpdates = run.correlatedUpdates;      // Switch to do correlated update proposals
  mcmc->nCorr = run.nCorr;                              // Number of iterations for which the covariance matrix is calculated
  mcmc->prMatrixInfo = run.prMatrixInfo;                // Print information to screen on proposed matrix updates
  
  mcmc->annealTemp0 = run.annealTemp0;                  // Starting temperature of the annealed chain
  mcmc->annealNburn = run.annealNburn;                  // Number of iterations for the annealing burn-in phase
  mcmc->annealNburn0 = run.annealNburn0;                // Number of iterations during which temp=annealTemp0
  
  mcmc->parallelTempering = run.parallelTempering;      // Switch to use parallel tempering
  mcmc->maxTemp = run.maxTemp;                          // Maximum temperature in automatic parallel-tempering ladder
  mcmc->saveHotChains = run.saveHotChains;              // Save hot (T>1) parallel-tempering chains
  mcmc->prParTempInfo = run.prParTempInfo;              // Print information on the temperature chains
  
  
  mcmc->chTemp = max(mcmc->annealTemp0,1.0);            // Current temperature
  
  mcmc->geocentricTc = run.geocentricTc;                // Geocentric Tc
  mcmc->baseTime = (double)((floor)(mcmc->geocentricTc/100.0)*100);  //'Base' time, gets rid of the first 6-7 digits of GPS time
  
  
  for(i=0;i<mcmc->maxnPar;i++) {
    mcmc->injParVal[i] = run.injParVal[i];
    mcmc->injID[i] = run.injID[i];
    
    mcmc->parNumber[i] = run.parNumber[i];
    mcmc->parID[i] = run.parID[i];
    mcmc->parBestVal[i] = run.parBestVal[i];
    mcmc->parFix[i] = run.parFix[i];
    mcmc->parStartMCMC[i] = run.parStartMCMC[i];
    mcmc->parSigma[i] = run.parSigma[i];
    
    mcmc->priorBoundLow[i] = run.priorBoundLow[i];
    mcmc->priorBoundUp[i] = run.priorBoundUp[i];
    mcmc->priorType[i] = run.priorType[i];
    
    //mcmc->[i] = run.[i];
  }
  
  //Parameter database:
  for(i=0;i<mcmc->parDBn;i++) {
    mcmc->injRevID[i] = run.injRevID[i];
    mcmc->parRevID[i] = run.parRevID[i];
    mcmc->mcmcParUse[i] = run.mcmcParUse[i];
    mcmc->injParUse[i]  = run.injParUse[i];
    
    
    for(j=0;j<99;j++) {
      mcmc->parName[i][j] = run.parName[i][j];
      mcmc->parAbrev[i][j] = run.parAbrev[i][j];
      mcmc->parAbrv[i][j] = run.parAbrv[i][j];
    }
  }
  
  mcmc->nTemps = run.nTemps;                            // Size of temperature ladder
  for(i=0;i<mcmc->nTemps;i++) mcmc->tempLadder[i] = run.tempLadder[i];
  
} // End copyRun2MCMC()
// ****************************************************************************************************************************************************  












// ****************************************************************************************************************************************************  
/**
 * \brief Set global variables
 *
 * \todo This routine should eventually contain mathematical and (astro)physical constants only.
 */
// ****************************************************************************************************************************************************  
void setConstants()
{
  // Mathematical constants:
  pi   = 3.141592653589793;   // pi
  tpi  = 6.283185307179586;   // 2 pi
  mtpi = 6.283185307179586e6; // Large multiple of 2 pi (2 megapi)
  
  // Define some physical constants:
  G    = 6.67259e-11;         // 6.674215e-11; */ /* gravity constant (SI)
  c    = 299792458.0;         // speed of light (m/s)
  
  Ms   = 1.9889194662e30;     // solar mass (kg)
  Mpc  = 3.08568025e22;       // metres in a Mpc  (LAL: 3.0856775807e22)
  Mpcs = 1.029272137e14;      // seconds in a Mpc  (Mpc/c)
} // End setConstants
// ****************************************************************************************************************************************************  








// ****************************************************************************************************************************************************  
/**
 * \brief Returns the 'injection values' to the parameter set par
 *
 * \todo Remove par->mc etc. struct elements
 */
// ****************************************************************************************************************************************************  
void getInjectionParameters(struct parSet *par, int nInjectionPar, double *injParVal)
{
  int i=0;
  for(i=0;i<nInjectionPar;i++) {
    par->par[i] = injParVal[i];
  }
  par->nPar = nInjectionPar;
  
} // End getInjectionParameters()
// ****************************************************************************************************************************************************  




// ****************************************************************************************************************************************************  
/**
 * \brief Returns the 'best-guess values' to the parameter set par
 */
// ****************************************************************************************************************************************************  
void getStartParameters(struct parSet *par, struct runPar run)  //Set the parameters to the starting values for the MCMC chain
{
  int i=0;
  for(i=0;i<run.nMCMCpar;i++) {
    par->par[i]      = run.parBestVal[i];
  }
  par->nPar = run.nMCMCpar;
} // End getStartParameters
// ****************************************************************************************************************************************************  



// ****************************************************************************************************************************************************  
/**
 * \brief Allocate memory for the vectors in the struct parSet
 */
// ****************************************************************************************************************************************************  
void allocParset(struct parSet *par, int networkSize)
{
  par->loctc    = NULL;
  par->localti  = NULL;
  par->locazi   = NULL;
  
  par->loctc    = (double*)calloc(networkSize,sizeof(double));
  par->localti  = (double*)calloc(networkSize,sizeof(double));
  par->locazi   = (double*)calloc(networkSize,sizeof(double));
} // End allocParset
// ****************************************************************************************************************************************************  




// ****************************************************************************************************************************************************  
/**
 * \brief Deallocate memory for the vectors in the struct parSet
 */
// ****************************************************************************************************************************************************  
void freeParset(struct parSet *par)
{
  free(par->loctc);         par->loctc        = NULL;
  free(par->localti);       par->localti      = NULL;
  free(par->locazi);        par->locazi       = NULL;
} // End freeParset
// ****************************************************************************************************************************************************  







