/***********************************************\
* stripAdc.c                                    *
* 02/18/02                                      *
*                                               *
* Routine to strip ADC channels.  Outputs       *
* merged frames.  Input are in a config file.   *
* There is only one option available:           *
*                                               *
* F - copy full adc data                        *
*                                               *
\***********************************************/
/* $Id$ */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <FrameL.h>
  
#define TRUE             1
#define NCHAR           31
#define NCHARMAX       100
#define NADC           100
#define NFILES          50
#define COMPRESSION     -1

int main(int argc, char *argv[])
{
   struct FrFile     *iFile, *oFile;
   struct FrameH     *frame;
   struct FrHistory  *history;

   short int         longitudeD;
   short int         longitudeM;
   float             longitudeS;
   short int         latitudeD;
   short int         latitudeM;
   float             latitudeS;
   float             elevation;
   float             armXazimuth;
   float             armYazimuth;
   float             armLength;

   struct in {
    unsigned int     startTime;
    unsigned int     delT;
    unsigned int     framesPerFile;
    char             outFrameFile[NCHARMAX];
    char            *inFrameFile[NFILES];
    char            *channel[NADC];
    char            *operation[NADC];
    unsigned int     nFrameFiles;
    unsigned int     nChannels;
   } in;

   char             *allChannels;
   char              c;
   char            **channel;
   char              command[200];
   unsigned int      elapsedTime;
   short int         ENDRDS;
   time_t            endTime;
   unsigned int      frameCtr, newFrameCtr;
   short int         grepCtr;
   unsigned int      i;
   unsigned int      inCtr;
   unsigned int      inFileCtr;
   unsigned int      nChannels;
   unsigned int      nFrameFiles;
   char            **operation;
   unsigned int      outFileCtr;
   time_t            startTime;
   time_t            startTimeFrames;
   char              stringTemp[NCHARMAX];
   struct tm        *tmEndTime;
   struct tm        *tmFrames;

   FILE             *cfp;
   FILE             *ifp1, *ifp3;
   FILE             *ofp2;
   char              configFile[60];
   char              inFile3[60];
   char              outFile1[60], outFile2[60];
   char              iFileName[100];

   startTime = time(NULL);

/*** read input from file ***/

   if (argc != 2)
   {
    printf("Invalid number of arguments.\n");
    exit(-1);
   }
   strcpy(configFile, argv[1]);
   if (access(configFile, F_OK) != 0)
   {
    printf("Configuration file %s does not exist.\n", configFile);
    exit(-1);
   }
   
   if ((cfp = fopen(configFile, "r")) == NULL)
   {
      printf("Error opening input file %s.\n", configFile);
      exit(-1);
   }

   for (i = 0; i < NFILES; i++)
    in.inFrameFile[i] = malloc(NCHARMAX*sizeof(char));
   for (i = 0; i < NADC; i++)
   {
    in.channel[i] = malloc(NCHAR*sizeof(char));
    in.operation[i] = malloc(2*sizeof(char));
   }

   in.nFrameFiles = 0;
   in.nChannels   = 0;
   inCtr          = 0;
   fseek(cfp,0,0);
   while ((c = fgetc(cfp)) != EOF)
   {
    if (c == '=')
    {
     inCtr++;
     switch (inCtr)
     {
      case 1 :  fscanf(cfp,"%d", &(in.startTime));
                break;
      case 2 :  fscanf(cfp,"%d", &(in.delT));
                break;
      case 3 :  fscanf(cfp,"%d", &(in.framesPerFile));
                break;
      case 4 :  i = 0;
                while (TRUE)
                { 
                 fscanf(cfp,"%s", stringTemp);
                 if (strcmp(stringTemp,"input") == 0) break;
                 strcpy(in.outFrameFile,stringTemp);
                 i++;
                }
                if (i > 1)
                {
                 printf("Invalid number of output file prefixes.\n");
                 exit(-1);
                }
                break;
      case 5 :  i = 0;
                while (TRUE)
                {
                 fscanf(cfp,"%s", stringTemp);
                 if (strcmp(stringTemp,"channels") == 0) break;
                 strcpy(in.inFrameFile[i],stringTemp);
                 i++;
                }
                in.nFrameFiles = i;
                break;
      case 6 :  i = 0;
                while (fscanf(cfp,"%s %s", in.channel[i], in.operation[i]) != EOF)
                {
                 if (strcmp(in.operation[i],"F") != 0 &&
                     strcmp(in.operation[i],"f") != 0)
                 {
                  printf("Invalid operation code:  %s\n", in.operation[i]);
                  exit(-1);
                 }
                 i++;
                }
                in.nChannels = i;
                break;
     }
    }
   }
   fclose(cfp);

   printf("%d\n", in.startTime);
   printf("%d\n", in.delT);
   printf("%d\n", in.framesPerFile);
   printf("%s\n", in.outFrameFile);

/*** open output log file ***/

   startTimeFrames = time(NULL);
   tmFrames = gmtime(&startTimeFrames);
   sprintf(outFile2,"rdslog_%02d-%02d-%02d.txt", tmFrames->tm_mon+1, tmFrames->tm_mday,\
                                                 tmFrames->tm_year-100);
   if ((ofp2 = fopen(outFile2,"a")) == NULL)
   {
    printf("Error opening log file %s.\n", outFile2);
    exit(-1);
   }

   allChannels = NULL;
   allChannels = malloc(NADC*NCHAR*sizeof(char));
   if (allChannels == NULL)
   {
    printf("Cannot allocate allChannels.\n");
    fprintf(ofp2,"Cannot allocate allChannels.\n");
    exit(-1);
   }

   for (i = 0; i < in.nFrameFiles; i++)
   {
    if (i == 0)
    {
     sprintf(command,"ls %s > framefiletemp.ffl", in.inFrameFile[i]);
    }
    else
    {
     sprintf(command,"ls %s >> framefiletemp.ffl", in.inFrameFile[i]);
    }
    system(command);
   }
   if (access("channelnametemp.txt", F_OK) == 0)
    system("rm channelnametemp.txt");
   for (i = 0; i < in.nChannels; i++)
   {
    grepCtr = 0;
    sprintf(command,"grep '%s' adc_lho.txt", in.channel[i]);
    ifp1 = popen(command, "r");
    while (fscanf(ifp1,"%s", stringTemp) != EOF)
    {
     grepCtr++;
     strcat(allChannels, " ");
     strcat(allChannels, stringTemp);
     strcat(stringTemp, " ");
     if (strcmp(in.operation[i],"f") == 0) strcpy(in.operation[i],"F");
     strcat(stringTemp, in.operation[i]);
     if (access("channelnametemp.txt", F_OK) == 0)
     {
      sprintf(command,"echo %s >> channelnametemp.txt", stringTemp);
     }
     else
     {
      sprintf(command,"echo %s > channelnametemp.txt", stringTemp);
     }
     system(command);
    }
    pclose(ifp1);     
    if (grepCtr == 0)
    {
     printf("Invalid channel name %s\n", in.channel[i]);
     exit(-1);
    }
    free(in.channel[i]);
    free(in.operation[i]);
   }
   printf("%s\n", allChannels);

   nChannels = 0;
   while (access("channelnametemp.txt", F_OK) != 0) {}
   ifp1 = popen("wc -l channelnametemp.txt", "r");
   fscanf(ifp1,"%s", stringTemp);
   pclose(ifp1);
   nChannels = atoi(stringTemp);
   nFrameFiles = 0;
   ifp1 = popen("wc -l framefiletemp.ffl", "r");
   fscanf(ifp1,"%s", stringTemp);
   pclose(ifp1);
   nFrameFiles = atoi(stringTemp);
   printf("number of channels :  %d\n", nChannels);
   printf("number of files    :  %d\n", nFrameFiles);

   ifp1 = popen("cat channelnametemp.txt", "r");
   channel = malloc(nChannels*sizeof(char *));
   operation = malloc(nChannels*sizeof(char *));
   for (i = 0; i < nChannels; i++)
   {
    channel[i] = malloc(NCHAR*sizeof(char));
    operation[i] = malloc(2*sizeof(char));
    fscanf(ifp1,"%s %s", channel[i], operation[i]);
   }
   pclose(ifp1);


/*** open the input frame files ***/

  strcpy(inFile3, "framefiletemp.ffl");
  if ((ifp3 = fopen(inFile3, "r")) == NULL)
  {
   printf("Error opening input file %s.\n", inFile3);
   exit(-1);
  }

  longitudeD  = 0;
  longitudeM  = 0;
  longitudeS  = 0.0;
  latitudeD   = 0;
  latitudeM   = 0;
  latitudeS   = 0.0;
  elevation   = 0.0;
  armXazimuth = 0.0;
  armYazimuth = 0.0;
  armLength   = 0.0;

/*** read all the frames ***/
  
  startTimeFrames = time(NULL);
  tmFrames = gmtime(&startTimeFrames);
  printf("Start reduced data set:  UTC %sReading all frames...\n",\
                                            asctime(tmFrames));
  fprintf(ofp2,"Start reduced data set:  UTC %sReading all frames...\n",\
                                            asctime(tmFrames));
  frameCtr    = 0;
  newFrameCtr = 0;
  outFileCtr  = 0;
  inFileCtr   = 0;
  ENDRDS      = 0;
  fseek(ifp3,0,0);
  while (fscanf(ifp3,"%s",iFileName) != EOF)
  {
   iFile = NULL;
   if ((iFile = FrFileINew(iFileName)) == NULL)
   {
    printf("Error opening input file %s!\n", iFileName);
    fprintf(ofp2,"Error opening input file %s!\n", iFileName);
    break;
   }
   iFile->compress = 1;
   inFileCtr++;
   if (inFileCtr == 1)
   {
    frame = NULL;
    if ((frame = FrameRead(iFile)) != NULL)
    {
     longitudeD  = frame->detectProc->longitudeD;
     longitudeM  = frame->detectProc->longitudeM;
     longitudeS  = frame->detectProc->longitudeS;
     latitudeD   = frame->detectProc->latitudeD;
     latitudeM   = frame->detectProc->latitudeM;
     latitudeS   = frame->detectProc->latitudeS;
     elevation   = frame->detectProc->elevation;
     armXazimuth = frame->detectProc->armXazimuth;
     armYazimuth = frame->detectProc->armYazimuth;
     armLength   = frame->detectProc->armLength;

     history = NULL;
     history = FrHistoryAdd(NULL,frame->history->comment);
     if (history == NULL)
     {
      printf("Error reading history.\n");
      fprintf(ofp2,"Error reading history.\n");
      FrameFree(frame);
      FrFileIEnd(iFile);
      break;
     }
     history->time = frame->history->time;
     strcpy(history->name,frame->history->name);

     FrameFree(frame);
    }
    else
    {
     printf("Error reading frame.\n");
     fprintf(ofp2,"Error reading frame.\n");
     FrFileIEnd(iFile);
     break;
    }
   }

   frame = NULL;
   if ((frame = FrameReadTAdc(iFile, 0, allChannels)) != NULL)
   {
    frameCtr++;
    if (frame->GTimeS < in.startTime)
    {
     FrameFree(frame);
     FrFileIEnd(iFile);
     continue;
    }
    elapsedTime = frame->GTimeS - in.startTime;
    if (elapsedTime >= in.delT)
    {
     FrDetectorFree(frame->detectProc);
     FrHistoryFree(frame->history);
     FrameFree(frame);
     FrFileIEnd(iFile);
     ENDRDS = 1;
     break;
    }

    frame->detectProc = NULL;
    frame->detectProc = FrDetectorNew("LIGO_1");
    if (frame->detectProc == NULL)
    {
     printf("Error creating detectProc.\n");
     fprintf(ofp2,"Error creating detectProc.\n");
     FrameFree(frame);
     FrFileIEnd(iFile);
     break;
    }
    frame->detectProc->longitudeD  = longitudeD;
    frame->detectProc->longitudeM  = longitudeM;
    frame->detectProc->longitudeS  = longitudeS;
    frame->detectProc->latitudeD   = latitudeD;
    frame->detectProc->latitudeM   = latitudeM;
    frame->detectProc->latitudeS   = latitudeS;
    frame->detectProc->elevation   = elevation;
    frame->detectProc->armXazimuth = armXazimuth;
    frame->detectProc->armYazimuth = armYazimuth;
    frame->detectProc->armLength   = armLength;

    frame->history = NULL;
    frame->history = FrHistoryAdd(frame,history->comment);
    if (frame->history == NULL)
    {
     printf("Error adding history.\n");
     fprintf(ofp2,"Error adding history.\n");
     FrDetectorFree(frame->detectProc);
     FrameFree(frame);
     FrFileIEnd(iFile);
     break;
    }
    frame->history->time = history->time;
    strcpy(frame->history->name,history->name);
  
    if (newFrameCtr == 0)
    {
     sprintf(outFile1,"%s-%d-16.gwf", in.outFrameFile, frame->GTimeS);
     oFile = NULL;
     if ((oFile = FrFileONewH(outFile1,COMPRESSION, "stripAdc2.c")) == NULL)
     {
      printf("Error opening output file %s!\n", outFile1);
      fprintf(ofp2,"Error opening output file %s!\n", outFile1);
      FrDetectorFree(frame->detectProc);
      FrHistoryFree(frame->history);
      FrameFree(frame);
      FrFileIEnd(iFile);
      break;
     }
    }
    if (FrameWrite(frame,oFile) != FR_OK)
    {
     printf("%s:  Error during frame write!\n", outFile1);
     fprintf(ofp2,"%s:  Error during frame write!\n", outFile1);
     FrDetectorFree(frame->detectProc);
     FrHistoryFree(frame->history);
     FrameFree(frame);
     FrFileOEnd(oFile);
     FrFileIEnd(iFile);
     outFileCtr++;
     newFrameCtr = 0;
     break;
    }
    else
    {
     newFrameCtr++;
    }
    if (newFrameCtr == in.framesPerFile)
    {
     FrFileOEnd(oFile);
     outFileCtr++;
     newFrameCtr = 0;
    }
   
    printf("%u processed\n", frame->GTimeS);
    fprintf(ofp2,"%u processed\n", frame->GTimeS);
    fflush(ofp2);
    FrDetectorFree(frame->detectProc);
    FrHistoryFree(frame->history);
    FrameFree(frame);
   } /* frame loop */
   else
   {
    printf("Error reading frame %s.\n", iFileName);
    fprintf(ofp2,"Error reading frame %s.\n", iFileName);
    FrFileIEnd(iFile);
    break;
   }
   FrFileIEnd(iFile);
/* if (ENDRDS == TRUE) break; */
  } /* input frame file loop */
  
/*** close files ***/

   if (newFrameCtr > 0 && newFrameCtr < in.framesPerFile)
   {
    FrFileOEnd(oFile);
    outFileCtr++;
    printf("%d number of frames saved -- less than %d\n", newFrameCtr, in.framesPerFile);
   }
   fclose(ifp3);

/*** end program ***/

   endTime = time(NULL);
   printf("total processing time :  %u seconds\n",\
                                          (unsigned int)difftime(endTime,startTime));
   fprintf(ofp2,"total processing time :  %u seconds\n",\
                                          (unsigned int)difftime(endTime,startTime));
   tmEndTime = gmtime(&endTime);
   printf("End reduced data set:  UTC %s", asctime(tmEndTime));
   fprintf(ofp2,"End reduced data set:  UTC %s", asctime(tmEndTime));
   fclose(ofp2);
   return(0);
}
