/*
*  Copyright (C) 2007 Craig Robinson , Thomas Cokelaer
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

/*
 [PURPOSE]
	 convert ascii file from Bankefficiency code into an xml file

 [INPUT]
	 an ascii file provideed by BankEfficiency code
	 which name is "Trigger.dat"
 
 optional:
	 a prototype called BE_Proto.xml given by BankEfficiency code 
	 for non standalone jobs
 	
 [OUTPUT]
	 the equivalent xml file called "Trigger.xml"
 
 [USAGE]
	 BEAscii2Xml 
	 no optino or argument

 [AUTHOR]
 	Thomas Cokelaer April 2004

[Compilation]
        gcc BEAscii2Xml.c -o BEAscii2Xml  -I. -I/home/cokelaer/share/include -I../lalapps -I/software/geopptools/include -I../
	gcc BEAscii2Xml.c -o BEAscii2Xml -I/<include file>
*/



#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralHeaders.h>
#include "BankEfficiency.h"

#define BEASCII2XML_INPUT1 "Trigger.dat"
#define BEASCII2XML_INPUT2 "BE_Proto.xml"
#define BEASCII2XML_OUTPUT "Trigger.xml"
#define BEASCII2XML_BANK   "BE_Bank.xml"


void
BEAscii2XmlHelp(void);



int
main (int argc, char **argv ) 
{
  UINT4 line=0, i;
UINT8  id=0;
  UINT4 start = 0;
  ResultIn trigger;
  REAL4 tau0, tau3, tau0I, tau3I, psi0, psi3, phaseI, coaTime;
  FILE *input1, *input2, *bank;
  FILE *output;      
  SnglInspiralTable     *inputData = NULL;
  INT4 numFileTriggers = 0;

  MetadataTable         templateBank;




  char sbuf[512];
  
  if (argc>1) {
    if (strcmp(argv[1], "--help")==0) {
      BEAscii2XmlHelp();
    }  
    if (strcmp(argv[1], "-h")==0) { 
      BEAscii2XmlHelp();
    }  
  }
  /* Main program starts here */
  /* First we open the file containing the ascii results */
  fprintf(stderr,"opening the xml output data file -- %s", BEASCII2XML_OUTPUT);
  output = fopen(BEASCII2XML_OUTPUT,"w");
  fprintf(stderr,"done\n");

  fprintf(stderr,"opening the xml input data file -- %s", BEASCII2XML_INPUT1);
  if  ( (input1  = fopen(BEASCII2XML_INPUT1,"r"))==NULL) {
    fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT1);
    exit(0);
  }
  fprintf(stderr,"done\n");

  fprintf(stderr,"opening the xml prototype (argument of BankEfficiency code) -- ");
  if  ( (input2  = fopen(BEASCII2XML_INPUT2,"r"))==NULL)
    {
      fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT2);
      fprintf(stderr,"the xml file will not contains parameters information\n");
      PRINT_LIGOLW_XML_HEADER(output);
      fprintf(stderr,"creating the header file -- done\n");
    }
  else 
    {
      /* read prototype and save in outputfile */
      fprintf(stderr,"parsing the prototype  -- ");
      while(fgets(sbuf,1024,input2) !=NULL)
	fputs(sbuf, output);
      fprintf(stderr," done\n");

      
    }
  

  /* insert the template bank here */
  if  ( (bank  = fopen(BEASCII2XML_BANK,"r"))==NULL)
    {
      fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_BANK);
      fprintf(stderr,"the xml file will not contains the bank table\n");
    }
  else 
    {
      /* read prototype and save in outputfile */
      fprintf(stderr,"parsing the bank  -- ");
      fprintf( stdout, "reading triggers from file: %s\n", BEASCII2XML_BANK );
      numFileTriggers = 
        LALSnglInspiralTableFromLIGOLw( &inputData,BEASCII2XML_BANK , 0, -1 );
      fprintf(stderr," done %d\n", numFileTriggers);      
       myfprintf(output, LIGOLW_XML_SNGL_INSPIRAL );
       while(inputData)
	{
	  /*	  id = inputData->event_id->id;*/
	  
	  fprintf(output, SNGL_INSPIRAL_ROW, 
		  inputData->ifo,
		  inputData->search,
		  inputData->channel,
		  inputData->end_time.gpsSeconds,
		  inputData->end_time.gpsNanoSeconds,
		  inputData->end_time_gmst,
		  inputData->impulse_time.gpsSeconds,
		  inputData->impulse_time.gpsNanoSeconds,
		  inputData->template_duration,
		  inputData->event_duration,
		  inputData->amplitude,
		  inputData->eff_distance,
		  inputData->coa_phase,
		  inputData->mass1,
		  inputData->mass2,
		  inputData->mchirp,
		  inputData->mtotal,
		  inputData->eta,
		  inputData->tau0,
		  inputData->tau2,
		  inputData->tau3,
		  inputData->tau4,
		  inputData->tau5,
		  inputData->ttotal,
		  inputData->psi0,
		  inputData->psi3,
		  inputData->alpha,
		  inputData->alpha1,
		  inputData->alpha2,
		  inputData->alpha3,
		  inputData->alpha4,
		  inputData->alpha5,
		  inputData->alpha6,
		  inputData->beta,
		  inputData->f_final,
		  inputData->snr,
		  inputData->chisq,
		  inputData->chisq_dof,
		  inputData->sigmasq,
		  id);
	  inputData = inputData->next;
	  fprintf(output, "\n");

	}
       myfprintf(output, LIGOLW_XML_TABLE_FOOTER );

    }



  PRINT_LIGOLW_XML_BANKEFFICIENCY(output);
  fprintf(stderr,"done\n");
  /* read ascii input and save in xml format */
  fprintf(stderr,"reading the ascii file -- and saving xml file");
  do  
    {
      fscanf(input1,BANKEFFICIENCY_PARAMS_ROW_SPACE,
	     &trigger.psi0_triggerU,
	     &trigger.psi3_triggerU,
	     &trigger.psi0_trigger,
	     &trigger.psi3_trigger,
	     &psi0, 
	     &psi3,
	     &tau0, 
	     &tau3, 
	     &tau0I, 
	     &tau3I,
	     &trigger.fend_triggerU, 
	     &trigger.fend_trigger, 
	     &trigger.fend_inject,
	     &trigger.mass1_inject,
	     &trigger.mass2_inject,
	     &trigger.rho_finalU,
	     &trigger.phaseU,
	     &phaseI,
	     &trigger.alphaFU, 
	     &trigger.layerU,
	     &trigger.binU,
	     &trigger.rho_final,
	     &trigger.snrAtCoaTime, 
	     &phaseI,
	     &trigger.phase,
	     &trigger.alphaF, 
	     &trigger.layer,
	     &trigger.bin, 
	     &coaTime);

     if (start==0){
	      start+=1;
      }
      else 
      {
	      fprintf(output,",\n");
      }
      fprintf(output, BANKEFFICIENCY_PARAMS_ROW,
	      trigger.psi0_triggerU,
	      trigger.psi3_triggerU,
	      trigger.psi0_trigger,
	      trigger.psi3_trigger,
	      psi0, 
	      psi3,
	      tau0,
	      tau3,
	      tau0I,
	      tau3I,
	      trigger.fend_triggerU, 
	      trigger.fend_trigger, 
	      trigger.fend_inject,
	      trigger.mass1_inject,
	      trigger.mass2_inject,
	      trigger.rho_finalU,
	      phaseI,
	      trigger.phaseU,
	      trigger.alphaFU, 
	      trigger.layerU,
	      trigger.binU,
	      trigger.rho_final,
	      trigger.snrAtCoaTime,
	      phaseI,
	      trigger.phase,
	      trigger.alphaF, 
	      trigger.layer,
	      trigger.bin, 
	      coaTime);
         line++;
    }
   while(!feof(input1));

    fprintf(stderr,"read %d lines...done\n", line);
  PRINT_LIGOLW_XML_TABLE_FOOTER(output);
  PRINT_LIGOLW_XML_FOOTER(output);

  fclose(output);
    fprintf(stderr,"closing xml file\n");

  return 0;
}


void BEAscii2XmlHelp(void)
{
  fprintf(stderr, "BEAscii2Xml help:\n");
  fprintf(stderr, "=================\n");
  fprintf(stderr, "[PURPOSE]\n\t%s\n \t%s\n \t%s\n \t%s\n", 
		  "That code reads a file in ascii format generated by BankEfficiency",
		  "code and generate the appropriate xml files. It is used if one forget",
		  "to use the appropriate option within BankEfficiency code (--print-result-xml)",
		  "or when the bankEfficiency code has been used within the condor ",
		  "submit file. In that case indeed several ascii files have been ",
		  "generated but no xml files.\n");
  
  fprintf(stderr, "[INPUT/OUTPUT]\n\t%s\n \t%s\n \t%s\n \t%s %s %s %s\n", 
		  "That code do not need any input argument but needs an input file",
		  "in ascii format (the user has to be sure it is in the appropriate",
		  "format) and then it writes an xml file into an output file \n", 
		  "the input and output file name are ", BEASCII2XML_INPUT1, "and", BEASCII2XML_OUTPUT);

  exit(1);
  

}
