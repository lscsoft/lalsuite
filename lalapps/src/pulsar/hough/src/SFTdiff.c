/* program for comparing two sfts in a frequency band */

#include <lal/SFTfileIO.h>


/* Default parameters. */

NRCSID (SFTCLEANC, "$Id$");


INT4 lalDebugLevel=0;


#define MAXFILENAMELENGTH 512
/* defaults chosen for L1 */

#define FILE1 "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs/CAL_SFT.734359206"
#define FILE2 "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean/CLEAN_SFT.734359206"
#define STARTFREQ 150.0
#define BANDFREQ 1.0



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){ 
  static  LALStatus   status;  /* LALStatus pointer */ 
  SFTtype     *sft1, *sft2;
  REAL8       startFreq, bandFreq; 
  /*  CHAR file1[MAXFILENAMELENGTH], file2[MAXFILENAMELENGTH]; */
  CHAR *file1, *file2;
  INT4 arg;

  /* set defaults */
  file1 = FILE1;
  file2 = FILE2;
  startFreq = STARTFREQ;
  bandFreq = BANDFREQ;

  /********************************************************/  
  /* Parse argument list.  i stores the current position. */
  /********************************************************/
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      }
    }  
    /* parse first sft filename */
    else if ( !strcmp( argv[arg], "-A" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        file1 = argv[arg++];
      } 
    }  
    /* parse second sft filename */
    else if ( !strcmp( argv[arg], "-B" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        file2 = argv[arg++];
      }
    }  
    /* parse start freq */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        startFreq = atof(argv[arg++]);
      }
    }  
    /* parse bandwidth  */
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        bandFreq = atof(argv[arg++]);
      }
    }  
    /* Unrecognized option. */
    else {
      printf("unknown argument\n");
      printf("options are: \n");
      printf("-d LALdebuglevel -A firstsftfile -B secondsftfile -f startfrequency -b freqband\n");
      exit(0);
    }
  } 

  /* End of argument parsing loop. */
  /******************************************************************/   


  
  sft1 = NULL;
  sft2 = NULL;
  LALReadSFTfile (&status, &sft1, startFreq, startFreq + bandFreq, file1);
  REPORTSTATUS( &status);

  LALReadSFTfile (&status, &sft2, startFreq, startFreq + bandFreq, file2);
  /*REPORTSTATUS( &status);*/
  {
    UINT4 j, nBins;
    REAL8 diff;
    COMPLEX8 *data1, *data2;
    nBins = sft1->data->length;

    for (j=0; j<nBins; j++)
      {
	data1 = sft1->data->data + j;
	data2 = sft2->data->data + j;
	diff = (data1->re - data2->re)*(data1->re - data2->re) + (data1->im - data2->im)*(data1->im - data2->im);
	printf("%1.3e\n", sqrt(diff));
      }
  }

  LALDestroySFTtype (&status, &sft1);
  /*REPORTSTATUS( &status);*/

  LALDestroySFTtype (&status, &sft2);
  /*REPORTSTATUS( &status);*/

  LALCheckMemoryLeaks(); 
  /*REPORTSTATUS( &status);*/

  return status.statusCode;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













