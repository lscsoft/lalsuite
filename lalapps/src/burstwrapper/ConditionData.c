#include <burstdso.h>
#include "Translation.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>

  /*****************************************/
  /* prepare output symbols */
  /* For each output statement in the algorithms, we
     have to allocate an entry in the symbols table
     for the stand-alone datacondAPI, so that the data
     can be sent back to the analysis code. In the 
     output statement, the comment field must comtain
     the datatype of the output data in square brackets,
     for instance:
     output(cavfaccplx,_,_,H2:CAL-CAV_FAC,H2 cavity factor [COMPLEX8TimeSeries]);
     Note that this is not necessary if the output type
     is a single precision time series. */
  /*****************************************/
int OutputSymbols(char *algorithms, 
		  int *Nsymbols,
		  datacond_symbol_type **symbols) {

  /*
    Arguments:
    algorithms: string with datacondAPI actions
    Nsymbols: number of datacondAPI symbols
    symbols: datacondAPI symbols chain

    Returns 0 if OK, 1 on error.
  */

  /* NOTE: allocates name in symbol table */
  
  char *p0, *p1;  /* pointers for string manipulation */
  char *buf;      /* local copy of algorithms string */

  buf = (char *)calloc(1+strlen(algorithms),sizeof(char));


  /* copy first line of algorithms into buf
  p0 = algorithms;
  p1 = strchr(algorithms,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  /* loop as long as we have actions */
  while(strlen(buf) || buf[0]=='\n') {

#ifdef DEBUGBURST
    fprintf(stderr,"Algo: %s\n",buf);
#endif

    /* check if we have an output statement */
    if(strstr(buf,"output(")) {

      char *r, *s, *tmp, *name; /* pointers for string manipulation */

      tmp = buf;

      /* find fourth argument of comma separated list: */
      r = strchr(tmp,',');
      if(!r) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }

      r = strchr(r+1,',');
      if(!r) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }

      r = strchr(r+1,',');
      if(!r) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }

      r++;

      s = strchr(r,',');
      if(!s) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }

      *s = 0;

      /* this is fourth argument (name): */
      name = (char *)calloc(1 + strlen(r), sizeof(char));
      strcpy(name, r);

      /* remove spaces */
      {
	int i,j;
	for(i=0;i<strlen(name);i++) {
	  if(name[i] == 32) {
	    for(j=i;j<strlen(name);j++) {
	      name[j] = name[j+1];
	    }
	  }
	}
      }

      /* look at comment (fifth argument) */
      s = s+1;
      r = strchr(s,')');
      if(!r) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }
      *r = 0;

      /* add to datacond symbols */
      *symbols = (datacond_symbol_type *)realloc(*symbols, (1 + *Nsymbols) * sizeof(datacond_symbol_type));
      (*symbols + *Nsymbols)->s_direction = DATACOND_SYMBOL_OUTPUT; /* output */
      /* choose translation function according to datatype */
      if(strstr(s,"COMPLEX8FrequencySeries")) {
	(*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8FrequencySeries;
      } else {
	if(strstr(s, "COMPLEX8TimeSeries")) {
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8TimeSeries;
	} else {
	  /* default is REAL4 time series */
	  (*symbols + *Nsymbols)->s_translator = TranslateREAL4TimeSeries;
	}
      }

      (*symbols + *Nsymbols)->s_symbol_name = name;
      (*symbols + *Nsymbols)->s_aux_data = NULL;
      (*symbols + *Nsymbols)->s_user_data = NULL;
      (*Nsymbols)++;

    }

    /* update line being processed */
    if(p1) {
      p0 = p1+1;
      p1 = strchr(p0,'\n');
      bzero(buf,strlen(buf));
      if(p1) {
	memcpy(buf,p0,p1-p0+1);
      } else {
	strcpy(buf,p0);
      }
    } else {
      break;
    }

  }

  free(buf);

  return 0;
}


/*****************************************/
/* condition data by running the datacondAPI */
/*****************************************/
int ConditionData(int Nsymbols,
		  datacond_symbol_type *symbols,
		  char *algorithms) {

  /*
    Arguments:
    Nsymbols: number of symbols
    symbols: datacondAPI symbols
    algorithms: datacondAPI actions as a string
  */

  int	error = 0;
  char* aliases = "";
  char* error_message = (char*)NULL;

  /*
   * Make call into algorithm section
   */
  error = DatacondCaller( algorithms,
			  aliases,
			  symbols,
			  Nsymbols,
			  &error_message );

  /* check return status: */
  if ( error != DATACOND_OK )
    {
      switch( error )
	{
	case DATACOND_PARSE_ALIASES_FAILURE:
	  fprintf( stderr, "FAILED: %s\n", aliases );
	  break;
	case DATACOND_PARSE_FAILURE:
	  fprintf( stderr, "FAILED: %s\n", algorithms );
	  break;
	default:
	  fprintf( stderr, "FAILED: %s\n", error_message );
	  break;
	}
      free( error_message );
      return 1;
    }

  /* by definition, symbols[0] is the GW channel and must contain data */
  if ( symbols[0].s_user_data == NULL )
    {
      fprintf( stderr, "Produced no results\n" );
      return 1;
    }


  return 0;
}


