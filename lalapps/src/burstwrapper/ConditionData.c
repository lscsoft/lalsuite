#include <burstdso.h>
#include "Translation.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>

int OutputSymbols(char *algorithms, 
		  int *Nsymbols,
		  datacond_symbol_type **symbols) {

  /* NOTE: allocates name in symbol table */
  
  char *p0, *p1;
  char *buf;

  buf = (char *)calloc(1+strlen(algorithms),sizeof(char));

  p0 = algorithms;
  p1 = strchr(algorithms,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  while(strlen(buf) || buf[0]=='\n') {

    printf("Algo: %s\n",buf);

    if(strstr(buf,"output(")) {

      char *r, *s, *tmp, *name;

      tmp = buf;

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

      /* look at comment */
      s = s+1;
      r = strchr(s,')');
      if(!r) { fprintf(stderr,"Malformed output: %s\n",buf); return 1; }
      *r = 0;

      /* add to datacond */
      *symbols = (datacond_symbol_type *)realloc(*symbols, (1 + *Nsymbols) * sizeof(datacond_symbol_type));
      (*symbols + *Nsymbols)->s_direction = DATACOND_SYMBOL_OUTPUT;
      
      if(strstr(s,"COMPLEX8FrequencySeries")) {
	(*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8FrequencySeries;
      } else {
	if(strstr(s, "COMPLEX8TimeSeries")) {
	  (*symbols + *Nsymbols)->s_translator = TranslateCOMPLEX8TimeSeries;
	} else {
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



int ConditionData(int Nsymbols,
		  datacond_symbol_type *symbols,
		  char *algorithms) {

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

  /* by definition, symbols[0] is the GW channel */
  if ( symbols[0].s_user_data == NULL )
    {
      fprintf( stderr, "Produced no results\n" );
      return 1;
    }


  return 0;
}


