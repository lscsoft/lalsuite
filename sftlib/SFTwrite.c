/* $Id$ */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "SFTReferenceLibrary.h"

/* some local prototypes */
FILE *openfile(const char* name);
void printerror(int err);


/* function definitions */

FILE *openfile(const char* name) {
  FILE *fp=fopen(name, "w");
  if (!fp) {
    if (errno)
      perror("Unable to open file for writing");
    fprintf(stderr,"Unable to open file %s for writing\n", name);
    exit(SFTENULLFP);
  }
  return fp;
}

void printerror(int err) {
  if (err) {
    fprintf(stderr, "WriteSFT failed with return value %d\n", err);
    exit(err);
  }
  return;
}

void dosystem(const char *command) {
  int returnval=system(command);
 if (returnval){
    fprintf(stderr, "Unable to execute the command: %s (exit status %d)\n", command, returnval);
    exit(returnval);
  }
}

/* define some detectors for convenience */
#define DET1 "H1"
#define DET2 "L1"
#define DET3 "X3"

int main(void) {
  FILE *fp;
  float data[]={1.0, 0.0, 2.0, -1.0, 3.0, -2.0, 4.0, -3.0};
  
  printerror(WriteSFT(fp=openfile("SFT-test1"), 12345,          6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test1", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test2"), 12345+60,       6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test2", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test3"), 12345+60+60,    6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test3", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test4"), 12345+60+60+60, 6789, 50, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test4", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test5"), 12345+60+60+60, 6789, 60, 1100, sizeof(data)/(2*sizeof(float)),  DET1, "test5", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test6"), 12345+60+60+60, 6789, 60, 1000, sizeof(data)/(2*sizeof(float))-1,DET1, "test6", data)); fclose(fp);
  printerror(WriteSFT(fp=openfile("SFT-test7"), 12345+60+60+60, 6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET2, "test7", data)); fclose(fp);

  dosystem("cat SFT-test[123] > SFT-good");
  /* GPS times not increasing */
  dosystem("cat SFT-test[123] SFT-test3 > SFT-bad1");
  /* Different time baselines */
  dosystem("cat SFT-test[124] > SFT-bad2");
  /* Different first frequency indexes */
  dosystem("cat SFT-test[125] > SFT-bad3");
  /* Different numbers of samples */
  dosystem("cat SFT-test[126] > SFT-bad4");
  /* Different detectors */
  dosystem("cat SFT-test[127] > SFT-bad5");
  
  printf("To test SFTs, do for example:\n"
	 "./SFTvalidate SFT-good SFT-test[1234567]\n"
	 "./SFTvalidate SFT-bad1\n"
	 "./SFTvalidate SFT-bad2\n"
	 "./SFTvalidate SFT-bad3\n"
	 "./SFTvalidate SFT-bad4\n"
	 "./SFTvalidate SFT-bad5\n"
	 "(checking exit status after each command) or you can also try\n"
	 "./SFTdumpheader SFT-good SFT-test[1234567]\n"
	 "./SFTdumpheader SFT-bad1\n"
	 "./SFTdumpheader SFT-bad2\n"
	 "./SFTdumpheader SFT-bad3\n"
	 "./SFTdumpheader SFT-bad4\n"
	 "./SFTdumpheader SFT-bad5\n"
	 "or you can also replace SFTdumpheader with SFTdumpall.\n");

  return 0;
}
