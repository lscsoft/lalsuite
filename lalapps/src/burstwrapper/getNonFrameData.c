#include "burstdso.h"
#include "datacondAPI/DatacondCaller.h"
#include "datacondAPI/TransFunc.h"

int getNonFrameData(char *rFiles, 
		    char **algorithms, 
		    int *Nsymbols,
		    datacond_symbol_type **symbols) {

  char *p0, *p1;
  char *buf;

  buf = (char *)calloc(1+strlen(rFiles),sizeof(char));

  /* process line by line */

  p0 = rFiles;
  p1 = strchr(rFiles,'\n');
  if(p1) {
    memcpy(buf,p0,p1-p0+1);
  } else {
    strcpy(buf,p0);
  }

  while(strlen(buf)) {

    char file[1024], pp[64], name[256];

    if(buf[0] != '\n') {

#ifdef DEBUGBURST
      fprintf(stderr,"nonframe: %s\n",buf);
#endif

      /* FORMAT: file push/pass name */
      sscanf(buf,"%s\t%s\t%s",file,pp,name);

      if(!strcmp(pp,"push")) {

	char *alias, *fname;

	alias = (char *)calloc(1 + strlen(name), sizeof(char));
	strcpy(alias,name);

	fname = (char *)calloc(1 + strlen(file), sizeof(char));
	strcpy(fname, file);

	/* sending file to datacond */
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
