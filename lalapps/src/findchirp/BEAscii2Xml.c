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

#include <stdio.h>
#include <stdlib.h>

#include <lal/LIGOLwXMLHeaders.h>
#include "BankEfficiency.h"

#define BEASCII2XML_INPUT1 "Trigger.dat"
#define BEASCII2XML_INPUT2 "BE_Proto.xml"
#define BEASCII2XML_OUTPUT "Trigger.xml"



#define BANKEFFICIENCY_PARAMS_INPUT \
"         %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %d %d %d"




int
main (int argc, char **argv ) 
{
  UINT4 start = 0;
  ResultIn trigger;
  FILE *input1, *input2;
  FILE *output;
  char sbuf[512];
  
  output = fopen(BEASCII2XML_OUTPUT,"w");

  if  ( (input1  = fopen(BEASCII2XML_INPUT1,"r"))==NULL){
    fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT1);
    exit(0);
  }

  if  ( (input2  = fopen(BEASCII2XML_INPUT2,"r"))==NULL){
    fprintf(stderr,"error while opening input file %s\n",BEASCII2XML_INPUT2);
    fprintf(stderr,"the xml file will not contains parameters information\n");
    fprintf(output,"%s", LIGOLW_XML_HEADER);
  }
  else 
    {
      /* read prototype and save in outputfile */
      while(fgets(sbuf,512,input2) !=NULL)
	fputs(sbuf, output);
      
    }

    fprintf(output,"%s", LIGOLW_XML_BANKEFFICIENCY);

  /* read ascii input and save in xml format */
    do  
    {
      fscanf(input1,BANKEFFICIENCY_PARAMS_INPUT,
	     &trigger.psi0_trigger,
	     &trigger.psi3_trigger,
	     &trigger.psi0_triggerC,
	     &trigger.psi3_triggerC,
	     &trigger.psi0_inject, 
	     &trigger.psi3_inject,
	     &trigger.fend_trigger, 
	     &trigger.fend_triggerC, 
	     &trigger.fend_inject,
	     &trigger.totalMass_trigger,
	     &trigger.eta_trigger,
	     &trigger.totalMass_triggerC,
	     &trigger.eta_triggerC,
	     &trigger.mass1_inject,
	     &trigger.mass2_inject,
	     &trigger.rho_final,
	     &trigger.phase,
	     &trigger.alpha,
	     &trigger.alpha_f, 
	     &trigger.layer,
	     &trigger.bin,
	     &trigger.rho_finalC,
	     &trigger.phaseC,
	     &trigger.alphaC,
	     &trigger.alpha_fC, 
	     &trigger.layerC,
	     &trigger.binC, &trigger.coaTime


	     );
      if (start==0){
	      start+=1;
      }
      else 
      {
	      fprintf(output,"\n");
      }
      fprintf(output, BANKEFFICIENCY_PARAMS_ROW,
	      trigger.psi0_trigger,
	      trigger.psi3_trigger,
	      trigger.psi0_triggerC,
	      trigger.psi3_triggerC,
	      trigger.psi0_inject, 
	      trigger.psi3_inject,
	      trigger.fend_trigger, 
	      trigger.fend_triggerC, 
	      trigger.fend_inject,
	      trigger.totalMass_trigger,
	      trigger.eta_trigger,
	      trigger.totalMass_triggerC,
	      trigger.eta_triggerC,
	      trigger.mass1_inject,
	      trigger.mass2_inject,
	      trigger.rho_final,
	      trigger.phase,
	      trigger.alpha,
	      trigger.alpha_f, 
	      trigger.layer,
	      trigger.bin,
	      trigger.rho_finalC,
	      trigger.phaseC,
	      trigger.alphaC,
	      trigger.alpha_fC, 
	      trigger.layerC,
	      trigger.binC, trigger.coaTime);
    }
   while(!feof(input1));
  fprintf(output,"%s", LIGOLW_XML_TABLE_FOOTER);
  fprintf(output,"%s", LIGOLW_XML_FOOTER);


  return 0;
}
