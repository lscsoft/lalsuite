/*
*  Copyright (C) 2012 Matt Pitkin
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

/* Many of these functions are (with occasional modification) taken directly 
 * from the TEMPO2 software package http://www.atnf.csiro.au/research/pulsar/tempo2/
 * written by George Hobbs and Russell Edwards */

/* create an ascii text ephemeris file of time delays to either convert from
 * TT to TDB (as in TEMPO) or TT to TCB/Teph (from Irwin and Fukushima, 1999)
 * as used in TEMPO2 */

#include "create_time_correction_ephemeris.h"

int verbose = 0;

int main(int argc, char **argv){
  inputParams_t inputParams;
  
  FILE *fp = NULL;
  
  int i=0, ne=0;
  long double mjdtt=0.;
  double deltaT=0.;
  
  struct tm beg, fin;
  
  char outputfile[MAXFNAME];
  
  /* get input arguments */
  get_input_args(&inputParams, argc, argv);
  
  /* create output file name */
  beg = *XLALGPSToUTC( &beg, (int)inputParams.startT );
  fin = *XLALGPSToUTC( &fin, (int)inputParams.endT );
  
  /* make sure the year in the filename only covers full years */
  if ( beg.tm_yday != 0 && beg.tm_hour != 0 && beg.tm_min != 0 && beg.tm_sec != 0 ) beg.tm_year++;
  if ( fin.tm_yday != 0 && fin.tm_hour != 0 && fin.tm_min != 0 && fin.tm_sec != 0 ) fin.tm_year--;
  
  if( inputParams.et == TT2TCB ){
    if( !sprintf(outputfile, "%s/te405_%d-%d.dat", 
      inputParams.outputpath, beg.tm_year+1900, fin.tm_year+1900) ){
      fprintf(stderr, "Error... problem creating output file name!\n");
      exit(1);
    }
  }
  else if( inputParams.et == TT2TDB ){
    if( !sprintf(outputfile, "%s/tdb_%d-%d.dat", 
      inputParams.outputpath, beg.tm_year+1900, fin.tm_year+1900) ){
      fprintf(stderr, "Error... problem creating output file name!\n");
      exit(1);
    }
  }
  
  /* open the output file */
  if( (fp = fopen(outputfile, "w")) == NULL ){
    fprintf(stderr, "Error... can't open output file: %s!\n", outputfile);
    exit(1);
  }
  
  ne = floor((inputParams.endT-inputParams.startT)/inputParams.interval);
  
  if( verbose ){
    fprintf(stdout, "Creating ephemeris file from %lf to %lf, with %.2lf second\
 time steps\n", inputParams.startT,  inputParams.endT, inputParams.interval);
  }
  
  /* output header info: start time, end time, time step, no. of entries */
  fprintf(fp, "%lf\t%lf\t%lf\t%d\n", inputParams.startT, inputParams.endT,
    inputParams.interval, ne);
  
  for( i=0; i<ne; i++){
    /* convert GPS to equivalent TT value (in MJD) */
    mjdtt = MJDEPOCH + (inputParams.startT + (double)i*inputParams.interval +
      GPS2TT)/DAYSTOSEC;
    
    /* these basically calculate the einstein delay */
      
    /* using TEMPO2/TT to TCB/Teph file */
    if( inputParams.et == TT2TCB ) deltaT = IF_deltaT(mjdtt);
    /* using TEMPO/TT to TDB file */
    else if( inputParams.et == TT2TDB ){
      deltaT = FB_deltaT(mjdtt, inputParams.ephemfile);
    }
    else{
      fprintf(stderr, "Error... unknown ephemeris type!\n");
      exit(1);
    }
    
    /* print output */
    fprintf(fp, "%.16lf\n", deltaT);
  }
  
  IFTE_close_file();
  
  fclose(fp);
  
  return 0;  
}

void get_input_args(inputParams_t *inputParams, int argc, char *argv[]){
  struct option long_options[] =
  {
    { "help",          no_argument,       0, 'h' },
    { "verbose",       no_argument, &verbose, 1 },
    { "ephem-type",    required_argument, 0, 't'},
    { "output-path",   required_argument, 0, 'o'},
    { "start",         required_argument, 0, 's'},
    { "interval",      required_argument, 0, 'i'},
    { "end",           required_argument, 0, 'e'},
    { 0, 0, 0, 0 }
  };

  char args[] = "ht:o:s:i:e:";
  char *program = argv[0];

  char *tempo2path;
  
  if(argc == 1){
    fprintf(stderr, "Error... no input parameters given! Try again!!!\n");
    exit(0);
  }

  if( (tempo2path = getenv("TEMPO2")) == NULL ){
    fprintf(stderr, "Error... TEMPO2 environment variable not set!\n");
    exit(1);
  }
  
  /* set defaults */
  inputParams->interval = 60.; /* default to 60 seconds between */

  inputParams->startT = 883180814.;  /* default start to 1 Jan 2008 00:00 UTC */
  inputParams->endT = 1072569615.;   /* default end to 1 Jan 2014 00:00 UTC */

  inputParams->outputpath = NULL;
  
  /* parse input arguments */
  while ( 1 ){
    int option_index = 0;
    int c;

    c = getopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch( c ){
      case 0:
        if( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error passing option %s with argument %s\n",
            long_options[option_index].name, optarg);
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 't':
        inputParams->ephemtype = strdup(optarg);
        break;
      case 'o':
        inputParams->outputpath = strdup(optarg);
        break;
      case 's':
        inputParams->startT = atof(optarg);
        break;
      case 'e':
        inputParams->endT = atof(optarg);
        break;
      case 'i':
        inputParams->interval = atof(optarg);
        break;
      case '?':
        fprintf(stderr, "Unknown error while parsing options\n");
      default:
        fprintf(stderr, "Unknown error while parsing options\n");
    }
  }
  
  /* set the ephemeris file */
  if( !strcmp(inputParams->ephemtype, "TEMPO") ||
      !strcmp(inputParams->ephemtype, "TDB") ){ /* using TEMPO/TDB file */
    inputParams->et = TT2TDB;
    sprintf(inputParams->ephemfile, "%s%s", tempo2path, TT2TDB_FILE);
  
    if( verbose ){
      fprintf(stdout, "Using TEMPO-style TT to TDB conversion file: %s\n",
        inputParams->ephemfile); 
    }
  }
  /* using TEMPO2/TCB/Teph file */
  else if( !strcmp(inputParams->ephemtype, "TEMPO2") ||
           !strcmp(inputParams->ephemtype, "TCB") ){
    inputParams->et = TT2TCB;
    sprintf(inputParams->ephemfile, "%s%s", tempo2path, IFTEPH_FILE);
  
    if( verbose ){
      fprintf(stdout, "Using TEMPO2-style TT to te405 conversion file: %s\n",
        inputParams->ephemfile); 
    }
    
    /* open and initialise file */
    IFTE_init( inputParams->ephemfile );
  }
  else{
    fprintf(stderr, "Error... unrecognised ephemeris type %s given.\n",
      inputParams->ephemtype);
    exit(1);
  }
}


/* read in the TDB-TT ephemeris file and output the time delay. Copied from
 the TEMPO2 tt2tdb.C file */
double FB_deltaT(long double mjd_tt, char fname[MAXFNAME]){
  double ctatv;    /* output TDB-TDT */
  long double tdt; /* jd1 + jd2  (Julian date) */
  static int tdbnrl=-1;
  int nr,j,k,np;
  double jda,jdb,tdbd1,tdbd2,t[2];
  int tdbdt,tdbncf;
  double dna,dt1,temp,pc[18],tc,twot;
  static double buf[16];

  /* Set up the TDB-TDT ephemeris file for reading */

  /* Now do the calculations */
  open_file(fname); /* Open a Fortran made file for reading in C */
  tdbd1 = read_double();
  tdbd2 = read_double();
  tdbdt = read_int();
  tdbncf = read_int();
  jda = read_double(); /* use jda as a dummy variable here */
  jda = read_double(); /* use jda as a dummy variable here */
  jda = read_double(); /* use jda as a dummy variable here */
  jda = read_double(); /* use jda as a dummy variable here */
  jda = read_double(); /* use jda as a dummy variable here */
              
  /* Use the corrected TT time and convert to Julian date */
  tdt = mjd_tt + 2400000.5; 
  if (tdt - (int)tdt >= 0.5) {
    jda = (int)tdt + 0.5;
    jdb = tdt - (int)tdt - 0.5;
  }
  else {
    jda = (int)tdt - 0.5;
    jdb = tdt - (int)tdt + 0.5;
  }
  nr = (int)((jda-tdbd1)/tdbdt)+2; 
  if (nr < 1 || tdt > tdbd2){
    printf("ERROR [CLK4]: Date %.10f out of range of TDB-TDT \
table\n",(double)tdt);
    exit(1);
  }
  if (nr!=tdbnrl) {
    tdbnrl = nr;
    /* MUST JUMP TO RECORD NR IN THE FILE */
    for (j=0;j<nr-1;j++)
      for (k=0;k<tdbncf;k++)
        buf[k] = read_double();
  }
  t[0] = ((jda-((nr-2)*tdbdt+tdbd1))+jdb)/tdbdt; /* Fraction within record */
  t[1] = 1.0; /* Unused */
              
  /* Interpolation: call interp(buf,t,tdbncf,1,  1,   1,   ctatv) */
  np = 2;
  twot = 0.0; 
              
  pc[0] = 1.0; pc[1]=0.0;
              
  dna = 1.0;
  dt1 = (int)(t[0]);
  temp = dna * t[0];
  tc = 2.0*(fortran_mod(temp,1.0)+dt1)-1.0;
              
  if (tc != pc[1]){
    np = 2;
    pc[1] = tc;
    twot = tc+tc;         
  } 
  if (np<tdbncf) {
    for (k=np+1;k<=tdbncf;k++) pc[k-1] = twot*pc[k-2]-pc[k-3];
    np = tdbncf;
  }
  ctatv = 0.0;
  for (j=tdbncf;j>=1;j--) ctatv = ctatv +pc[j-1]*buf[j-1];
  
  close_file();

  return ctatv;
}


double IF_deltaT( long double mjd_tt ){
  return IFTE_DeltaT(2400000.0+(int)mjd_tt, 0.5+(mjd_tt-(int)mjd_tt))*DAYSTOSEC;
}

/******************* FORTRAN FILE READING FUNCTIONS ************************/
/* taken from TEMPO2 read_fortran.h */

FILE *c_fileptr;
int swapByte;

int open_file(char fname[MAXFNAME]){
  int i;
  
  /* NOTE MUST PUT BACK TO 0 FOR SOLARIS!!!!! */
  union Convert{
    char asChar[sizeof(int)];
    int asInt;
  } intcnv;
  
  intcnv.asInt = 1;
  if (intcnv.asChar[0]==1) swapByte = 1;
  else swapByte = 0;
  
  /* Look for blank character in filename */
  for (i=0;fname[i];i++) {
    if (fname[i]==' '){
      fname[i]='\0';
      break;
    }
  }
  if ((c_fileptr = fopen(fname,"rb")) == NULL){
    fprintf(stderr, "Error: Unable to open filename %s\n", fname);
    exit(1);
  }
  return 0;
}

void close_file(void){
  fclose(c_fileptr);
}

int read_int(void){
  int i;
  
  union Convert{
    char asChar[sizeof(int)];
    int asInt;
  } intcnv;

  if (swapByte==1){
    for (i=sizeof(int)-1;i>=0;i--)
      intcnv.asChar[i]=fgetc(c_fileptr);
  }
  else {
    for (i=0;i<(int)sizeof(int);i++)
      intcnv.asChar[i]=fgetc(c_fileptr);
  }
  return intcnv.asInt;
}

double read_double(void){
  int i;

  union Convert {
    char asChar[sizeof(double)];
    double asDouble;
  } dblcnv;

  if (swapByte==1) {
    for (i=sizeof(double)-1;i>=0;i--)
      dblcnv.asChar[i]=fgetc(c_fileptr);
  }
  else{
    for (i=0;i<(int)sizeof(double);i++)
      dblcnv.asChar[i]=fgetc(c_fileptr);
  }
  return dblcnv.asDouble;
}

/* taken from TEMPO2 temp2Util.C */
double fortran_mod(double a, double p){
  double ret;
  ret = a - (int)(a/p)*p;
  return ret;
}

/**** FUNCTIONS FROM TEMPO2 ifteph.C file ***********/

void IFTE_init(const char fname[MAXFNAME]) {
  FILE *f;
  char buf[1024];
  int ncon;
  double double_in;

  if( (f = fopen(fname, "r")) == NULL ){
    fprintf(stderr, "Error... opening time ephemeris file '%s'\n",
            fname);
    exit(1);
  }

  /* read in header info */
  fread(buf, 1, 252, f); /* read CHARACTER*6 TTL(14,3) */
  fread(buf, 1, 12, f); /* read CHARACTER*6 CNAM(2) */
  fread(&ifte.startJD, 1, 8, f);
  fread(&ifte.endJD, 1, 8, f);
  fread(&ifte.stepJD, 1, 8, f);
  fread(&ncon, 1, 4, f);
  
  ifte.swap_endian = (ncon!=2); /* check for endianness */

  if (ifte.swap_endian)  IFTswapInt(&ncon);
  if (ncon!=2){ /* check that we can decode the file */
    fprintf(stderr, "Cannot understand format of time ephemeris file '%s'!\n",
            fname);
    exit(1);
  }
  if (ifte.swap_endian) {
    IFTswapDouble(&ifte.startJD); 
    IFTswapDouble(&ifte.endJD); 
    IFTswapDouble(&ifte.stepJD); 
  }

  fread(ifte.ipt, 8, 3, f);
  if (ifte.swap_endian) IFTswapInts(&ifte.ipt[0][0], 6); 

  /* figure out the record length */
  ifte.reclen = 4 * 2*(ifte.ipt[1][0]-1 + 3*ifte.ipt[1][1]*ifte.ipt[1][2]);
  
  /* get the constants from record "2" */
  fseek(f, ifte.reclen, SEEK_SET);
  fread(&double_in, 8, 1, f);
  if (ifte.swap_endian) IFTswapDouble(&double_in); 
  ifte.ephver = (int)floor(double_in);
  fread(&ifte.L_C, 8, 1, f);
  if (ifte.swap_endian) IFTswapDouble(&ifte.L_C); 

  ifte.f = f;
  ifte.iinfo.np = 2;
  ifte.iinfo.nv = 3;
  ifte.iinfo.pc[0] = 1.0;
  ifte.iinfo.pc[1] = 0.0;
  ifte.iinfo.vc[1] = 1.0;
  ifte.irec = -1;
  /* Note: file is not closed as it is used by other routines */
}

void IFTE_close_file(void){
  if (ifte.f != NULL)
    fclose(ifte.f);
}

/* functions for swapping endianess */
void IFTswap4(char *word) {
  char tmp;
  tmp = word[0]; word[0] = word[3]; word[3] = tmp;
  tmp = word[1]; word[1] = word[2]; word[2] = tmp;
} 

void IFTswapInt(int *word) {
  IFTswap4((char *)word);
}

void IFTswapInts(int *word, int n) {
  int i;
  for (i=0; i < n; i++) IFTswap4((char *)(word+i));
}

void IFTswap8(char *dword) {
  char tmp;
  tmp = dword[0]; dword[0] = dword[7]; dword[7] = tmp;
  tmp = dword[1]; dword[1] = dword[6]; dword[6] = tmp;
  tmp = dword[2]; dword[2] = dword[5]; dword[5] = tmp;
  tmp = dword[3]; dword[3] = dword[4]; dword[4] = tmp;
}

void IFTswapDouble(double *dbl) {
  IFTswap8((char *)dbl);
}

void IFTswap8N(char *dwords, int n){
  char tmp;
  int i;
  for (i=0; i < n; i++){
    tmp = dwords[0]; dwords[0] = dwords[7]; dwords[7] = tmp;
    tmp = dwords[1]; dwords[1] = dwords[6]; dwords[6] = tmp;
    tmp = dwords[2]; dwords[2] = dwords[5]; dwords[5] = tmp;
    tmp = dwords[3]; dwords[3] = dwords[4]; dwords[4] = tmp;
    dwords += 8;
  }
}

void IFTswapDoubles(double *dbl, int N) {
  IFTswap8N((char *)dbl, N);
}

/* general purpose value-getter */
void IFTE_get_Vals(double JDeph0, double JDeph1, int kind,
                   double *res){
  double whole0, whole1, frac0, frac1;
  int irec;
  double t[2];
  int ncoeff = ifte.reclen/8;
  size_t nread;
  
  whole0 = floor(JDeph0-0.5);
  frac0 = JDeph0-0.5-whole0;
  whole1 = floor(JDeph1);
  frac1 = JDeph1-whole1;
  whole0 += whole1 + 0.5;
  frac0 += frac1;
  whole1 = floor(frac0);
  frac1 = frac0-whole1;
  whole0 += whole1;

  JDeph0 = whole0;
  JDeph1 = frac1;

  if (JDeph0 < ifte.startJD) {
    fprintf (stderr, "Error: Requested JD=%lf is less than start JD=%lf\n",
             JDeph0, ifte.startJD);
    exit(1);
  }

  /* CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL */
  irec = (int)floor((JDeph0-ifte.startJD)/ifte.stepJD)+2;

  if (JDeph0 == ifte.endJD) irec--;

  t[0] = (JDeph0-(ifte.startJD+ifte.stepJD*(irec-2))+JDeph1)/ifte.stepJD;
  t[1] = ifte.stepJD;

  if (irec != ifte.irec){
    if (fseek(ifte.f, ifte.reclen*irec, SEEK_SET) < 0) {
      fprintf(stderr, "Error reading time ephemeris");
      exit(1);
    }
    nread = fread(ifte.buf, 8, ncoeff, ifte.f);
    if ((int)nread < ncoeff){
      fprintf(stderr, "Error reading time ephemeris: Only read %zd coefficients,\
 wanted %d!\n", nread, ncoeff);
      exit(1);
    }
    if (ifte.swap_endian) IFTswapDoubles(ifte.buf, ncoeff); 
  }
  
  /*  INTERPOLATE time ephemeris */
  if (kind==1)
    IFTEinterp(&ifte.iinfo, ifte.buf+ifte.ipt[0][0]-1, t, 
               ifte.ipt[0][1], 1, ifte.ipt[0][2], 2, res);
  else
    IFTEinterp(&ifte.iinfo, ifte.buf+ifte.ipt[1][0]-1, t, 
               ifte.ipt[1][1], 3, ifte.ipt[1][2], 2, res);
}

/* convenience interfaces */
void IFTE_get_DeltaT_DeltaTDot(double Teph0, double Teph1,
                               double *DeltaT, double *DeltaTDot) {
  double res[2];
  IFTE_get_Vals(Teph0, Teph1, 1, res);
  *DeltaT = res[0];
  *DeltaTDot = res[1];
}

double IFTE_DeltaT(double Teph0, double Teph1) {
  double DeltaT, DeltaTDot;
  IFTE_get_DeltaT_DeltaTDot(Teph0, Teph1, &DeltaT, &DeltaTDot);
  return DeltaT;
}

/*  The following routine is borrowed from the JPL ephemeris C code */
/*****************************************************************************
*        *****    jpl planetary and lunar ephemerides    *****     C ver.1.2 *
******************************************************************************
*                                                                            *
*  This program was written in standard fortran-77 and it was manually       *
*  translated to C language by Piotr A. Dybczynski (dybol@phys.amu.edu.pl),  *
*  subsequently revised heavily by Bill J Gray (pluto@gwi.net).              *
*                                                                            *
******************************************************************************/

/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**      iinfo   stores certain chunks of interpolation info,  in hopes      **
**              that if you call again with similar parameters,  the        **
**              function won't have to re-compute all coefficients/data.    **
**                                                                          **
**       coef   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**    posvel   interpolated quantities requested.  dimension                **
**              expected is posvel[ncm*ifl], double precision.              **
**                                                                          **
*****************************************************************************/
static void IFTEinterp( struct IFTE_interpolation_info *iinfo,
                 const double coef[], const double t[2], const int ncf,
                 const int ncm, const int na, const int ifl, double posvel[]){
  double dna, dt1, temp, tc, vfac, temp1;
  double *pc_ptr;
  int l, i, j; 

  /*  entry point. get correct sub-interval number for this set
      of coefficients and then get normalized chebyshev time
      within that subinterval.                                             */

  dna = (double)na;
  modf( t[0], &dt1);
  temp = dna * t[0];
  l = (int)(temp - dt1);

  /*  tc is the normalized chebyshev time (-1 <= tc <= 1)    */
  tc = 2.0 * (modf( temp, &temp1) + dt1) - 1.0;

  /*  check to see whether chebyshev time has changed,
      and compute new polynomial values if it has.
      (the element iinfo->pc[1] is the value of t1[tc] and hence
      contains the value of tc on the previous call.)     */

  if(tc != iinfo->pc[1]){
    iinfo->np = 2;
    iinfo->nv = 3;
    iinfo->pc[1] = tc;
    iinfo->twot = tc+tc;
  }

  /*  be sure that at least 'ncf' polynomials have been evaluated
      and are stored in the array 'iinfo->pc'.    */

  if(iinfo->np < ncf){
    pc_ptr = iinfo->pc + iinfo->np;

    for(i=ncf - iinfo->np; i; i--, pc_ptr++)
      *pc_ptr = iinfo->twot * pc_ptr[-1] - pc_ptr[-2];
    iinfo->np=ncf;
  }

  /*  interpolate to get position for each component  */
  for( i = 0; i < ncm; ++i){  /* ncm is a number of coordinates */
    const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
    const double *pc_ptr2 = iinfo->pc + ncf;

    posvel[i]=0.0;
    for( j = ncf; j; j--) posvel[i] += (*--pc_ptr2) * (*--coeff_ptr);
  }

  if(ifl <= 1) return;

  /*  if velocity interpolation is wanted, be sure enough
      derivative polynomials have been generated and stored.    */

  vfac=(dna+dna)/t[1];
  iinfo->vc[2] = iinfo->twot + iinfo->twot;
  if( iinfo->nv < ncf){
    double *vc_ptr = iinfo->vc + iinfo->nv;
    const double *pc_ptr2 = iinfo->pc + iinfo->nv - 1;

    for( i = ncf - iinfo->nv; i; i--, vc_ptr++, pc_ptr2++)
      *vc_ptr = iinfo->twot * vc_ptr[-1] + *pc_ptr2 + *pc_ptr2 - vc_ptr[-2];
    iinfo->nv = ncf;
  }

  /*  interpolate to get velocity for each component    */
  for( i = 0; i < ncm; ++i){
    double tval = 0.;
    const double *coeff_ptr = coef + ncf * (i + l * ncm + 1);
    const double *vc_ptr = iinfo->vc + ncf;

    for( j = ncf; j; j--) tval += (*--vc_ptr) * (*--coeff_ptr);
    posvel[ i + ncm] = tval * vfac;
  }
  
  return;
}
