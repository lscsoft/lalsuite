#include "burstdso.h"
#include "datacondAPI/DatacondCaller.h"
#include "datacondAPI/TransFunc.h"

/*****************************************/
/* get non-frame data */

/* loop over lines in rFiles, add to datacond callchain & algorithm */
/* format is:
   file push/pass datacond_symbol
   pass option not supported yet
   file can be a URL */
/*****************************************/

int getNonFrameData(char *rFiles, 
		    char **algorithms, 
		    int *Nsymbols,
		    datacond_symbol_type **symbols) {

  /*
    Arguments:
    rFiles: string with responseFiles
    algorithms: modifiable string with algorithms (not used by function!)
    Nsymbols: number of symbols
    symbols: datacond symbols
  */

  char *p0, *p1;  /* pointers for string manipulation */
  char *buf;      /* local copy of rfiles */

  buf = (char *)calloc(1+strlen(rFiles),sizeof(char));

  /* process line by line */

  /* 1st line */
  p0 = rFiles;
  p1 = strchr(rFiles,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  /* loop over non-empty lines */
  while(strlen(buf)) {

    char file[1024], pp[64], name[256]; /* fields */

    if(buf[0] != '\n') {

#ifdef DEBUGBURST
      fprintf(stderr,"nonframe: %s\n",buf);
#endif

      /* FORMAT: file push/pass name */
      sscanf(buf,"%s\t%s\t%s",file,pp,name);

      if(!strcmp(pp,"push")) { /* send data to datacondAPI */

	char *alias, *fname;

	/* permanent memory for alias */
	alias = (char *)calloc(1 + strlen(name), sizeof(char));
	strcpy(alias,name);

	/* permanent memory for file name */
	fname = (char *)calloc(1 + strlen(file), sizeof(char));
	strcpy(fname, file);

	/* sending file to datacond by creating new symbol */
	*symbols = (datacond_symbol_type *)realloc(*symbols, (1 + *Nsymbols) * sizeof(datacond_symbol_type));
	(*symbols + *Nsymbols)->s_direction = DATACOND_SYMBOL_INPUT;
	(*symbols + *Nsymbols)->s_translator = TranslateILwdFile;
	(*symbols + *Nsymbols)->s_symbol_name = alias;
	(*symbols + *Nsymbols)->s_aux_data = NULL;
	(*symbols + *Nsymbols)->s_user_data = fname;
	(*Nsymbols)++;
	
      } else {
	fprintf(stderr,"option %s not supported\n",pp);
	return 1;
      }
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
