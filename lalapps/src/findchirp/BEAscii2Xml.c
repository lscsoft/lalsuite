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
	gcc BEAscii2Xml.c -o BEAscii2Xml -I/<include file>
*/

#include <stdio.h>
#include <stdlib.h>

#include <lal/LIGOLwXMLHeaders.h>


#define BEASCII2XML_INPUT1 "Trigger.dat"
#define BEASCII2XML_INPUT2 "BE_Proto.xml"

#define BEASCII2XML_OUTPUT "Trigger.xml"


#define BANKEFFICIENCY_PARAMS_OUTPUT \
"         %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,"

#define BANKEFFICIENCY_PARAMS_INPUT \
"         %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d "


#define LIGOLW_XML_BANKEFFICIENCY \
"   <Table Name=\"bankefficiencygroup:bankefficiency:table\">\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0T\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3T\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fI\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:totalMassT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:etaT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass1I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass2I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:overlap\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:phase\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alpha\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alpha_f\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:layerT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:bin\" Type=\"real_4\"/>\n" \
"      <Stream Name=\"bankefficiencygroup:bankefficiency:table\" Type=\"Local\" Delimiter=\",\">\n"

typedef struct{
  REAL4 psi0_trigger;
  REAL4 psi3_trigger;
  REAL4 psi0_inject;
  REAL4 psi3_inject;
  REAL4 fend_trigger;
  REAL4 fend_inject;
  REAL4 mass1_inject;
  REAL4 mass2_inject;
  REAL4 totalMass_trigger;
  REAL4 eta_trigger;
  INT4 layer;
  REAL4 rho_final;
  REAL4 alpha;
  REAL4 alpha_f;
  INT4 bin;
  REAL4 phase;
  UINT4 ntrial;
} ResultIn;


int
main (int argc, char **argv ) 
{
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
      /*lire le proto et sauver dans ouput*/
      while(fgets(sbuf,512,input2) !=NULL)
	fputs(sbuf, output);
      
    }

    fprintf(output,"%s", LIGOLW_XML_BANKEFFICIENCY);

  /* lit ascii input and save in xml format */
  while(!feof(input1))
    {
      fscanf(input1,BANKEFFICIENCY_PARAMS_INPUT,
	     &trigger.psi0_trigger,
	      &trigger.psi3_trigger,
	      &trigger.psi0_inject, 
	      &trigger.psi3_inject,
	      &trigger.fend_trigger, 
	      &trigger.fend_inject,
	      &trigger.totalMass_trigger,
	      &trigger.eta_trigger,
	      &trigger.mass1_inject,
	      &trigger.mass2_inject,
	      &trigger.rho_final,
	      &trigger.phase,
	      &trigger.alpha,
	      &trigger.alpha_f, 
	      &trigger.layer,
	      &trigger.bin);
      fprintf(output, BANKEFFICIENCY_PARAMS_OUTPUT"\n",
	      trigger.psi0_trigger,
	      trigger.psi3_trigger,
	      trigger.psi0_inject, 
	      trigger.psi3_inject,
	      trigger.fend_trigger, 
	      trigger.fend_inject,
	      trigger.totalMass_trigger,
	      trigger.eta_trigger,
	      trigger.mass1_inject,
	      trigger.mass2_inject,
	      trigger.rho_final,
	      trigger.phase,
	      trigger.alpha,
	      trigger.alpha_f, 
	      trigger.layer,
	      trigger.bin);
    }
  fprintf(output,"%s", LIGOLW_XML_TABLE_FOOTER);
  fprintf(output,"%s", LIGOLW_XML_FOOTER);


  return 0;
}
