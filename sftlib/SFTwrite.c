/* $Id$ */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "SFTReferenceLibrary.h"

/* some local prototypes */
FILE *openfile(const char* name);
void printerror(int err);
void dosystem(const char *command);
void modify_bytes(const char *filename, int byte_offset, const char *new_value, int nbytes);
int isbigendian(void);
void modify_checksum(const char *filename, unsigned long long newchecksum);
void swap2(char *location);
void swap4(char *location);
void swap8(char *location);


/* function definitions */
int isbigendian(void){
  short i=0x0100;
  char *tmp=(char *)&i;
  return *tmp;
}

/* swap 2, 4 or 8 bytes.  Point to low address */
void swap2(char *location){
  char tmp=*location;
  *location=*(location+1);
  *(location+1)=tmp;
  return;
}
  
void swap4(char *location){
  char tmp=*location;
  *location=*(location+3);
  *(location+3)=tmp;
  swap2(location+1);
  return;
}
                                                                                                                  
void swap8(char *location){
  char tmp=*location;
  *location=*(location+7);
  *(location+7)=tmp;
  tmp=*(location+1);
  *(location+1)=*(location+6);
  *(location+6)=tmp;
  swap4(location+2);
  return;
}

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

void modify_bytes(const char *filename, int byte_offset, const char *new_value, int nbytes) {
  FILE *fp=fopen(filename, "rb+");
  if (!new_value) {
    fprintf(stderr, "modify_byes() called with null pointer for new data!\n");
    exit(SFTENULLFP);
  }
  if (!fp) {
    if (errno)
      perror("Unable to open file for writing");
    fprintf(stderr,"Unable to open file %s for reading/writing\n", filename);
    exit(SFTENULLFP);
  }
  if (fseek(fp, byte_offset, SEEK_SET)) {
    if (errno)
      perror("Failed fseek()");
    fprintf(stderr,"Failed seek in file %s\n", filename);
    exit(SFTESEEK);
  }
  if (1 != fwrite((const void *)new_value, nbytes, 1, fp)) {
    if (errno)
      perror("Failed fwrite()");
    fprintf(stderr,"Failed write to file %s\n", filename);
    exit(SFTESEEK);
  }
  fclose(fp);

  return;
}

void modify_checksum(const char *filename, unsigned long long newchecksum) {
  if (isbigendian())
    swap8((char *)&newchecksum);
  modify_bytes(filename, 32, (char *)&newchecksum, 8);
  return;
}


/* define some detectors for convenience */
#define DET1 "H1"
#define DET2 "L1"

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
  /* Inconsistent checksum (corrupted comment, change test1 to Test1) */
  dosystem("cat SFT-test1 > SFT-bad6");
  modify_bytes("SFT-bad6", 48, "T", 1);
  /* Nonexistent detector (change H1 to I1) checksum */
  dosystem("cat SFT-test1 > SFT-bad7");
  modify_bytes("SFT-bad7", 40, "I", 1);
  modify_checksum("SFT-bad7", 2594301065140926588ULL);
  /* SFT with comment containing hidden data */
  dosystem("cat SFT-test1 > SFT-bad8");
  modify_bytes("SFT-bad8", 49, "", 1);
  modify_checksum("SFT-bad8", 9546122005026447583ULL);
  /* SFT with comment no terminating null character containing hidden data */
  printerror(WriteSFT(fp=openfile("SFT-bad9"), 12345,          6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test567", data)); fclose(fp);
  modify_bytes("SFT-bad9", 55, "8", 1);
  modify_checksum("SFT-bad9", 8104828364422822013ULL);

  printf("To test SFTs, do for example:\n"
	 "./SFTvalidate SFT-good SFT-test[1234567]\n"
	 "./SFTvalidate SFT-bad1\n"
	 "./SFTvalidate SFT-bad2\n"
	 "./SFTvalidate SFT-bad3\n"
	 "./SFTvalidate SFT-bad4\n"
	 "./SFTvalidate SFT-bad5\n"
	 "./SFTvalidate SFT-bad6\n"
	 "./SFTvalidate SFT-bad7\n"
	 "./SFTvalidate SFT-bad8\n"
	 "./SFTvalidate SFT-bad9\n"
	 "(checking exit status after each command) or you can also try\n"
	 "./SFTdumpheader SFT-good SFT-test[1234567]\n"
	 "./SFTdumpheader SFT-bad1\n"
	 "./SFTdumpheader SFT-bad2\n"
	 "./SFTdumpheader SFT-bad3\n"
	 "./SFTdumpheader SFT-bad4\n"
	 "./SFTdumpheader SFT-bad5\n"
	 "./SFTdumpheader SFT-bad6\n"
	 "./SFTdumpheader SFT-bad7\n"
	 "./SFTdumpheader SFT-bad8\n"
	 "./SFTdumpheader SFT-bad9\n"
	 "or you can also replace SFTdumpheader with SFTdumpall.\n");

  return 0;
}
