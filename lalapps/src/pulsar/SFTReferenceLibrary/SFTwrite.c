/*
 *  Copyright (C) 2004 Bruce Allen
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

/**
 * \author Bruce Allen
 * \file
 * \brief
 * Makes a test SFT (called SFT-test)
 *
 * You can do this on little-endian and big-endian machines to generate
 * both flavors.  This produces a set of good and a set of bad SFTs.  The
 * good SFTs are:
 * SFT-test[1234567] and SFT-good
 * and the bad SFTs are
 * SFT-bad[123456789] and SFT-bad1[0-4]
 */

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
void modify_checksum(const char *filename, unsigned long long newchecksumle, unsigned long long newchecksumbe);

/* function definitions */
int isbigendian(void){
  short i=0x0100;
  char *tmp=(char *)&i;
  return *tmp;
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

void modify_checksum(const char *filename, unsigned long long newchecksumle, unsigned long long newchecksumbe) {
  if (isbigendian())
    modify_bytes(filename, 32, (char *)&newchecksumbe, 8);
  else
    modify_bytes(filename, 32, (char *)&newchecksumle, 8);
  return;
}


/* define some detectors for convenience */
#define DET1 "H1"
#define DET2 "L1"

int main(void) {
  FILE *fp;
  float data[]={1.0, 0.0, 2.0, -1.0, 3.0, -2.0, 4.0, -3.0};
  int tmp;
  double dtmp;

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
  modify_checksum("SFT-bad7", 2594301065140926588ULL, 1040429337927860515ULL);
  /* SFT with comment containing hidden data */
  dosystem("cat SFT-test1 > SFT-bad8");
  modify_bytes("SFT-bad8", 49, "", 1);
  modify_checksum("SFT-bad8", 9546122005026447583ULL, 12540893872916896128ULL);
  /* SFT with comment no terminating null character containing hidden data */
  printerror(WriteSFT(fp=openfile("SFT-bad9"), 12345,          6789, 60, 1000, sizeof(data)/(2*sizeof(float)),  DET1, "test567", data)); fclose(fp);
  modify_bytes("SFT-bad9", 55, "8", 1);
  modify_checksum("SFT-bad9", 8104828364422822013ULL, 6488197905334075682ULL);
  /* GPS nsec too big/small */
  tmp=1000000000;
  dosystem("cat SFT-test1 > SFT-bad10");
  modify_bytes("SFT-bad10", 12, (char *)&tmp, 4);
  modify_checksum("SFT-bad10",14553276594141125530ULL, 13919267759677155442ULL);
  tmp=-1;
  dosystem("cat SFT-test1 > SFT-bad11");
  modify_bytes("SFT-bad11", 12, (char *)&tmp, 4);
  modify_checksum("SFT-bad11",9905878389211410491ULL, 12381252028459463649ULL);
  /* first frequency index negative */
  tmp=-1;
  dosystem("cat SFT-test1 > SFT-bad12");
  modify_bytes("SFT-bad12", 24, (char *)&tmp, 4);
  modify_checksum("SFT-bad12",5605290681458522985ULL, 7307622666058252853ULL);
  /* number of samples negative */
  tmp=-1;
  dosystem("cat SFT-test1 > SFT-bad13");
  modify_bytes("SFT-bad13", 28, (char *)&tmp, 4);
  modify_checksum("SFT-bad13", 14214363586458317283ULL, 6130423751169108856ULL);
  /* time base not positive */
  dtmp=-60.0;
  dosystem("cat SFT-test1 > SFT-bad14");
  modify_bytes("SFT-bad14", 16, (char *)&dtmp, 8);
  modify_checksum("SFT-bad14", 9944972421627148413ULL, 15720824585133081082ULL);

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
	 "./SFTvalidate SFT-bad10\n"
	 "./SFTvalidate SFT-bad11\n"
	 "./SFTvalidate SFT-bad12\n"
	 "./SFTvalidate SFT-bad13\n"
	 "./SFTvalidate SFT-bad14\n"
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
	 "./SFTdumpheader SFT-bad10\n"
	 "./SFTdumpheader SFT-bad11\n"
	 "./SFTdumpheader SFT-bad12\n"
	 "./SFTdumpheader SFT-bad13\n"
	 "./SFTdumpheader SFT-bad14\n"
	 "or you can also replace SFTdumpheader with SFTdumpall.\n");

  return 0;
}
