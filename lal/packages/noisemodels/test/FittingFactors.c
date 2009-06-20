/*
*  Copyright (C) 2007 Stas Babak
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

/* ****************************************************************
 * Author: Babak, S., B. Krishnan
 * ************************************************************** */
#include "lal/FFUtils.h"

/* --- version information --- */
NRCSID( FITTINGFACTORSC, "$Id$");
RCSID(  "$Id$");

#define CVS_ID_STRING_C      "$Id$"
#define CVS_REVISION_C      "$Revision$"




InspiralCoarseBankIn          bankIn;        /* bank input parameters */
SnglInspiralTable            *tmplt  = NULL; /* linked list of templates */
CHAR  			 ifo[3];	     /* two char ifo code (need to make it up for LISA) */

int
main (INT4 argc, CHAR **argv)
{

  LALStatus	status = blank_status;

 /* --- Variables ---*/
  INT4 	        i;
  INT4          j;
  INT4	 	ntrials = 0;    /* number of trial (signals) */
  INT4		length = 0; 	/* duration of the template/signal */
  INT4		psdSize = 0;  	/*size of PSD vector */
  REAL4         match = 0;      /* overlap */
  REAL4		fitfac = 0;	/* fitting factor */

/* --- signal related --- */
  REAL4Vector                   signal;
  RandomSignalIn 	randIn;
  REAL4			normS;


  /* --- Bank related structures --- */

  INT4                          numTmplts    = 0;
  InspiralTemplate             	*bankHead     = NULL;
  InspiralTemplate             	*bankCurrent  = NULL;
  InspiralTemplateNode		*tmpltCurrent = NULL;
  InspiralTemplateNode	        *tmpltHead    = NULL;
  InspiralTemplate              bestTemplate;
  REAL4				normT;

 /* -- others --*/
  RealFFTPlan 		*fwdp = NULL;
  RealFFTPlan 		*revp = NULL;

/* set up inital debugging values */

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );


/* call the argument parser and check function */

  arg_parse_check( argc, argv);

/* Find the maximum length of the template */

 UINT4 checkLength = 0;
 GetMaximumTemplateSize(&status, &bankIn , &checkLength);


/* Compute noise psd */

 UINT psdLength = checkLength/2 + 1;
 bankIn.shf.deltaF  = ...;
 bankIn.shf.f0 = ...;
 bankIn.shf.data = NULL;
 LALDCreateVector( &status, &(bankIn.shf.data), psdLength);
 /** add option to read psd from the ascii file */
 GetNoisePSD(&status, bankIn.shf, ifo, df);


/***     Generate Bank    ******/

 LALInspiralBankGenreation( &status, bankIn, &tmplt, &numTmplts);

 if(numTmplts == 0, .....);


/*** Now we have read or generate signal ***/

/*** !!! Create structure *signal*  and allocate memories  ***/

if (readSignal){

 /* Here we read signal(s) from the file */

}else{

 /* Here we generate signal using lal functions */

}

/*** !!!! Make sure signal is NOT longer than max template length (checkLength) */


LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);


/**** Finding fitting factors ****/

/*  Create loop over signals */

for (i=0; i<ntrials; i++){

 /** Allocate memory for fSignal */
 CreateFSignal(&status, &fwdp, signal, &fSignal);

 fitfac = 0.0;
 /** loop over templates */
 for (j=0; j<numTmplts; j++){

   /** allocate memory for fTemplate */
   ...;
   CreateFTemplate(&status, &fwdp, bankCurrent, &fTemplate);

   ComputeOverlap(&status, &revp, &fSignal, &fTemplate, bankIn.shf, &match);

   if (fitfac < match){
       /* record template parameters into 'bestTemplate' */
       fitfac = match;
   }
 }
 if (fifac > 0.0 ){

   /* record: best template, signal index, fitting factor */

 }

}

/*** Deallocate all memories ***/


 /* Clean up */

  LALDestroyRealFFTPlan(&status,&fwdp);
  LALDestroyRealFFTPlan(&status,&revp);




  return(0);

}/* end of main */





/************************************************************************
*************************************************************************/


/*** This will be modified according to what we need */


#define USAGE\
  "lalapps_BankEfficiencyNew [options]\n\n"\
"  --help		display this message\n"\
"  --verbose		print progress info\n"\
"  --version		printversion and exit\n"\
"  --debug-level  LEVEL 	set LAL debug level to LEVEL\n"\
"  --user-tag STRING		set the user tag to STRING\n"\
"  --comment STRING		set process table comment to STRING\n"\
"  --bank-file FILE             read template bank from the xml FILE\n"\
"  --psd-file FILE             	read PSD from the data FILE\n"\
"  --sample-rate F              filter data at F Hz\n"\
"  --num-trials N		performe filtering of N random signals\n"\
"  --seed S			set seed for random number generator to S\n"\
"  --duration T			set the maximum length of template to T, it should correspond to the freq. resolution df in PSD!\n"\
"  --faithfulness		compute overlap between signal and template with the same parameters\n"\
"  --fast-simulation		compute overlap with templates in vicinity of the injection\n"\
"  --inj-mass1 M1		force first mass of injected signal to be M1 (comes together with M2)\n"\
"  --inj-mass2 M2		force second mass of injected signal to be M2 (comes together with M1)\n"\
"  --inj-tau0 T0		force injection to have tua0=T0 value (comes together with tau3)\n"\
"  --inj-tau3 T3		force injection to have tua3=T3 value (comes together with tau0)\n"\
"  --inj-psi0 P0		force injection to have psi0=P0 value (comes together with psi3)\n"\
"  --inj-psi3 P3		force injection to have psi3=P3 value (comes together with psi3)\n"\
"  --inj-alpha A		force injection to have alpha=A >= 0 value (BCV1)\n"\
"  --alpha-constraint C 	perform BCV1 filtering with alpha constraint if C>0, otherwise without \n"\
"  --inj-phase	PHASE		forces all signals to have initial phase = PHASE >= 0\n"\
"  --signal-approx APPR		set approximant model for injection: APPR should be one of these: TaylorT1, TaylorT2, TaylorT3, TaylorF2, PadeT1, EOB, BCV, SpinTaylorT3\n"\
"  --signal-order N		set PN order of injections:  onePointFivePN, twoPN, twoPointFivePN, threePN, threePointFivePN \n"\
"  --signal-min-mass MMIN	set min mass of injected signal to MMIN \n"\
"  --signal-max-mass MMAX	set max mass of injected signal to MMAX \n"\
"  --max-total-mass M		set maximum total mass of random signal to M\n"\
"  --signal-max-spin SMAX	set max. spin of injection to SMAX (not implemented yet)\n"\
"  --signal-min-spin SMIN	set min. spin of injection to SMIN (not implemented yet)\n"\
"  --signal-psi0-min PSI0MIN 	set min. psi0 for BCV injections to PSI0MIN \n"\
"  --signal-psi0-max PSI0MAX	set max psi0 for BCV iinjections to PSI0MAX \n"\
"  --signal-psi3-min PSI3MIN 	set min psi3 for BCV injections to PSI3MIN \n"\
"  --signal-psi3-max PSI3MAX	set max psi3 for BCV injections to PSI3MAX \n"\
"  --signal-beta-min BETAMIN	set min beta for BCV2 injections to BETAMIN \n"\
"  --signal-beta-max BETAMAX	set max beta for BCV2 injections to BETAMAX \n"\
"  --signal-fl   Fs		set the lower cut off frequency of signal to Fs\n"\
"  --signal-ffinal Ff		set the final frequency min(F_isco,Ff)\n"\
"  --template-approx APPR	set template model to be APPR: TaylorT1, TaylorT2, TaylorT3,TaylorF2, PadeT1, EOB, SpinTaylorT3 \n3"\
"  --template-order N		set PN order of the template to N (options the same as for signal)\n"\
"  --template-fl Ft      	set the lower cut off frequnecy of template to Ft\n"\
"  --template-ffinal Ff		set the final frequency of templates to min(Ff, ISCO, bank-defined)\n"\
"\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
  {"verbose",		`		no_argument, 	&vrbflg, 	1},
  {"disable-alpha-constraint",		no_argument,    &alphConst,	0},
  {"faithfulness",			no_argument,    &faithfulness	1},
  {"fast-simulation",			no_argument,    &fastsim	1},
    /* these options require argument */
  {"debug-level",			required_argument, 0,		'a'},
  {"user-tag",				required_argument, 0,		'A'},
  {"comment",				required_argument, 0,		'b'},
  {"bank-file",				required_argument, 0,		'B'},
  {"duration",				required_argument, 0,		'd'},
  {"psd-file",				required_argument, 0,		'D'},
  {"sample-rate",			required_argument, 0, 		'c'},
  {"num-trials",			required_argument, 0,		'C'},
  {"seed",				required_argument, 0,		'd'},
  {"inj-mass1",				required_argument, 0, 		'E'},
  {"inj-mass2",				required_argument, 0, 		'f'},
  {"inj-tau0", 				required_argument, 0, 		'F'},
  {"inj-tau3",				required_argument, 0,		'g'},
  {"inj-psi0",				required_argument, 0, 		'G'},
  {"inj-psi3",				required_argument, 0, 		'h'},
  {"inj-alpha",				required_argumnet, 0,		'H'},
  {"inj-phase",				required_argument, 0,		'I'},
  {"signal-approx",			required_argument, 0,		'k'},
  {"signal-order",			required_argument, 0,		'K'},
  {"signal-min-mass",			required_argument, 0,		'l'},
  {"signal-max-mass",			required_argument, 0,		'L'},
  {"max-total-mass",			required_argument, 0,		'm'},
  {"signal-max-spin",			required_argument, 0,		'M'},
  {"signal-min-spin",			required_argument, 0,		'n'},
  {"signal-psi0-min",			required_argument, 0,		'N'},
  {"signal-psi0-max",			required_argument, 0,		'o'},
  {"signal-psi3-min",			required_argument, 0,		'O'},
  {"signal-psi3-max",			required_argument, 0,		'p'},
  {"signal-beta-min",			required_argument, 0,		'P'},
  {"signal-beta-max",			required_argument, 0,		'q'},
  {"signal-fl",				required_argument, 0,		'Q'},
  {"signal-ffinal",			required_argument, 0,		'r'},
  {"template-approx",			required_argument, 0,		'R'},
  {"template-order",			required_argument, 0,		's'},
  {"template-fl",			required_argument, 0,		'S'},
  {"template-ffinal",			required_argument, 0,		't'},
  {"version",				no_argument,	   0,		'V'},
  {0, 0, 0, 0}
  };
  int c;

  LALStatus             status = blank_status;

  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {

    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv,
          "-A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:V"
	  "a:b:c:d:f:g:h:j:k:l:m:n:o:p:q:r:s:t:",
	  long_options, &options_index);
/* detect the end of the options */
    if ( c == - 1 )
    {
      break;
    }
    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;
      case 'a':
	set_debug_level(optarg);
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
      case 'A':
      /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;
      case 'b':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;
      case 'B':
        optarg_len = strlen( optarg ) + 1;
        bankFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( bankFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
      case 'd':
	duration = (REAL4) atof(optarg);
	if ( duration <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "duration must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, duration );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", duration );
        break;

      case 'D':
        optarg_len = strlen( optarg ) + 1;
        psdFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( psdFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'c':
        sampleRate = (REAL4) atof( optarg );
        if ( sampleRate < 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "rate must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, sampleRate );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", sampleRate );
        break;
      case 'C':
        numInjections = (INT4) atoi(optarg);
        if( numInjections <= 0){
	   fprintf(stderr, "invalid argument to --%s:\n"
			   "number of injections must be >0: "
			   "(%d specified)]\n",
		long_options[option_index].name, numInjections);
	   exit(1);
	}
	ADD_PROCESS_PARAM("int", "%d", numInjections);
	break;
      case 'd':
        useed = (INT4) atoi(optarg);
        ADD_PROCESS_PARAM("int", "%d", useed);
	break;
      case 'E':
	injMass1 = (REAL4) atof(optarg);
	if ( injMass1 <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "mass1 must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, injMass1 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", injMass1 );
	break;
      case 'f':
       	injMass2 = (REAL4) atof(optarg);
	if ( injMass2 <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "mass2 must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, injMass2 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", injMass2 );
	break;
      case 'F':
	injTau0 = (REAL4) atof(optarg);
	if ( injTau0 <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "tau0 must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, injTau0 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", injTau0 );
	break;
      case 'g':
	injTau3 = (REAL4) atof(optarg);
	if ( injTau3 <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "tau3 must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, injTau3 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", injTau3 );
	break;
     case 'G':
	injPsi0 = (REAL4) atof(optarg);
	if ( injPsi0 <= 0.0)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "psi0 must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, injPsi0 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", injPsi0 );
	break;
     case 'h':
	injPsi3 = (REAL4) atof(optarg);
        ADD_PROCESS_PARAM( "float", "%e", injPsi3 );
	break;
     case 'H':
	injAlpha = (REAL4) atof(optarg);
        ADD_PROCESS_PARAM( "float", "%e", injAlpha );
	break;
     case 'I':
	injPhase = (REAL4) atof(optarg);
        ADD_PROCESS_PARAM( "float", "%e", injPhase );
	break;
    case 'k':
        /* create storage for the approximant name */
        optarg_len = strlen( optarg ) + 1;
        approximantNameSignal = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( approximantNameSignal, optarg, optarg_len );
        if ( ! strcmp( "TaylorT1", optarg ) )
        {
          signalModel = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg ) )
        {
          signalModel = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg ) )
        {
           signalModel = TaylorT3;
        }
        else if ( ! strcmp( "TaylorF2", optarg ) )
        {
          signalModel = TaylorF2;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          signalModel  = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          signalModel  = EOB;
        }
        else if ( ! strcmp( "SpinTaylorT3", optarg ) )
        {
          signalModel = SpinTaylorT3;
        }
        else if ( ! strcmp( "BCV", optarg ) )
        {
          signalModel = BCV;
        }
       /* else if ( ! strcmp( "BCVSpin", optarg ) )
        {
          approximant = BCVSpin;
        }*/
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown model specified: "
              "%s (use --help for model list)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveSignalApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
    case 'K':
        if ( ! strcmp( "onePointFivePN", optarg ) )
        {
          signalOrder = onePointFivePN;
        }
	else if  ( ! strcmp( "twoPN", optarg ) )
        {
          signalOrder = twoPN;
        }
	else if  ( ! strcmp( "twoPointFivePN", optarg ) )
        {
          signalOrder = twoPointFivePN;
        }
	else if  ( ! strcmp( "threePN", optarg ) )
        {
          signalOrder = threePN;
        }
	else if  ( ! strcmp( "threePointFivePN", optarg ) )
        {
          signalOrder = threePointFivePN;
        }
	else{
	  fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s see help for possible order\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        haveSignalOrder = 1;
        break;
    case 'l':
        signalMinMass = atof(optarg);
	if ( signalMinMass <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "min mass must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMinMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'L':
        signalMaxMass = atof(optarg);
	if ( signalMaxMass <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "max mass must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMaxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case  'm':
        signalMaxTotMass = atof(optarg);
	if ( signalMaxTotMass <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "max total mass must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMaxTotMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'M':
        signalMaxSpin = atof(optarg);
	if ( signalMaxSpin <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "max spin must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMaxSpin );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'n':
        signalMinSpin = atof(optarg);
	if ( signalMinSpin <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "min spin must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMinSpin );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'o':
        signalMaxPsi0 =  atof(optarg);
	if ( signalMaxPsi0 < 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Psi0 of injections must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMaxPsi0 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'N':
        signalMinPsi0 =  atof(optarg);
	if ( signalMinPsi0 < 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Psi0 of injections must be positive: "
              "(%f specified)\n",
              long_options[option_index].name, signalMinPsi0 );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'p':
        signalMaxPsi3 =  atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'O':
        signalMinPsi3 =  atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'P':
        signalMinBeta =  atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'q':
        signalMaxBeta =  atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    case 'Q':
        signalfLow = (REAL4) atof( optarg );
        if ( signalfLow < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, signalfLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", signalfLow );
        break;
    case 'r':
        signalfFinal = (REAL4) atof( optarg );
        if ( signalfFinal <= signalfLow )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "upper frequency cutoff is less than lower: "
              "(%f is less than %f)\n",
              long_options[option_index].name, signalfFinal, signalfLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", signalfFinal );
        break;
    case 'R':
        /* create storage for the approximant name */
        optarg_len = strlen( optarg ) + 1;
        approximantNameTmplt = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( approximantNameTmplt, optarg, optarg_len );
        if ( ! strcmp( "TaylorT1", optarg ) )
        {
           tmpltModel = TaylorT1;
        }
        else if ( ! strcmp( "TaylorT2", optarg ) )
        {
           tmpltModel = TaylorT2;
        }
        else if ( ! strcmp( "TaylorT3", optarg ) )
        {
           tmpltModel = TaylorT3;
        }
        else if ( ! strcmp( "TaylorF2", optarg ) )
        {
           tmpltModel = TaylorF2;
        }
        else if ( ! strcmp( "PadeT1", optarg ) )
        {
          tmpltModel  = PadeT1;
        }
        else if ( ! strcmp( "EOB", optarg ) )
        {
          tmpltModel  = EOB;
        }
        else if ( ! strcmp( "SpinTaylorT3", optarg ) )
        {
          tmpltModel = SpinTaylorT3;
        }
        else if ( ! strcmp( "BCV", optarg ) )
        {
          tmpltModel = BCV;
        }
        else if ( ! strcmp( "BCVSpin", optarg ) )
        {
          tmpltModel = BCVSpin;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown template model specified: "
              "%s (use --help for model list)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveTmpltApprox = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
    case 's':
        if ( ! strcmp( "onePointFivePN", optarg ) )
        {
          tmpltOrder = onePointFivePN;
        }
	else if  ( ! strcmp( "twoPN", optarg ) )
        {
          tmpltOrder = twoPN;
        }
	else if  ( ! strcmp( "twoPointFivePN", optarg ) )
        {
          tmpltOrder = twoPointFivePN;
        }
	else if  ( ! strcmp( "threePN", optarg ) )
        {
          tmpltOrder = threePN;
        }
	else if  ( ! strcmp( "threePointFivePN", optarg ) )
        {
          tmpltOrder = threePointFivePN;
        }
	else{
	  fprintf( stderr, "invalid argument to --%s:\n"
              "unknown order specified: "
              "%s see help for possible order\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
	haveTmpltOrder = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
     case 'S':
        tmpltfLow = (REAL4) atof( optarg );
        if ( tmpltfLow < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "low frequency cutoff of tmplt is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, tmpltfLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", tmpltfLow );
        break;
    case 't':
        tmpltFinal = (REAL4) atof( optarg );
        if ( tmpltfFinal <= tmpltfLow )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "upper frequency cutoff of tmplt is less than lower: "
              "(%f is less than %f)\n",
              long_options[option_index].name, tmpltfFinal, tmpltfLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", tmpltfFinal );
        break;
    case 'V':
        fprintf( stdout, "Bank eficiency code\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

    case '?':
        exit( 1 );
        break;

    default:
        fprintf( stderr, "unknown error while parsing options (%d)\n", c );
        exit( 1 );
    } /* endo switch */

  } /* end of while loop */

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  if ( alphConst == 0 )
  {
    snprintf( procparams.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-alpha-constraint" );
    snprintf( procparams.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( procparams.processParamsTable->value,
        LIGOMETA_VALUE_MAX, " " );
  }
  if ( faithfulness == 1 )
  {
    snprintf( procparams.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--faithfullness" );
    snprintf( procparams.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( procparams.processParamsTable->value,
        LIGOMETA_VALUE_MAX, " " );
  }
  if ( fastsim == 1 )
  {
    snprintf( procparams.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( procparams.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--fast-simulation" );
    snprintf( procparams.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( procparams.processParamsTable->value,
        LIGOMETA_VALUE_MAX, " " );
  }

/*  Check presence of compalsory arguments */


  if ( ! bankFileName )
  {
    fprintf( stderr, "--bank-file must be specified\n" );
    exit( 1 );
  }
  if ( ! psdFileName )
  {
    fprintf( stderr, "--psd-file must be specified\n" );
    exit( 1 );
  }
  if (sampleRate <= 0.0){
    fprintf( stderr, "--sample-rate must be specified\n" );
    exit( 1 );
  }
  if (numInjections <= 0.0){
    fprintf( stderr, "--num-trials must be specified\n" );
    exit( 1 );
  }
  if ( haveTmpltOrder = 0 || haveTmpltApprox = 0){
    fprintf( stderr, "--template-order and --template-approx  must be specified\n" );
    exit(1);
  }
  if ( haveSignalOrder = 0 || haveSignalApprox = 0){
    fprintf( stderr, "--signal-order and --signal-approx  must be specified\n" );
    exit(1);
  }
  if (useed < 0 ){
     fprintf( stderr, "--seed  must be specified \n");
     exit(1);
  }
  if (signalfLow == 0.0){
     fprintf(stderr, "--signal-fl must be specified\n" );
     exit(1);
  }
  if (tmpltfLow == 0.0){
     fprintf(stderr, "--template-fl must be specified\n" );
     exit(1);
  }
  if (duration == 0.0){
     fprintf(stderr, "--duration must be specified\n" );
     exit(1);
  }


 return 0;
}



