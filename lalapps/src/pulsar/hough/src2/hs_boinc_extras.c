/* Extras for building an Einstein@Home BOINC App from HierarchicalSearch
   Bernd Machenschalk for Einstein@Home
   $Id$
*/

/* TODO:
   - error handling
   - signal handling
   - checkpointing
   - boinc_init(), boinc finish
*/

#include "hs_boinc_extras.h"


#define MAX_PATH_LEN 256


/* compare strings s1 and s2 up to the length of s1 (w/o the '\0') and set l to the length */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))



/* a local structure to keep information about the output files */
static struct {
  char **outfiles = NULL;        /* the names  of the output files of the program */
  int  noutfiles  = 0;           /* the number of the output files of the program */
  char resultfile[MAX_PATH_LEN]; /* the name of the file / zip archive to return */
} boinc_output;



/* this registers a new output file to be zipped into the archive that is returned
   to the server as a result file */
void register_output_file(char*filename) {
  boinc_output.names = realloc(boinc_output.names,(boinc_output.noutfiles+1)*sizeof(char*));
  if (boinc_output.names == NULL)
    LogPrintf (LOG_FATAL, "ERROR: Can't allocate output filename '%s'\n", filename);
  boinc_output.names[boinc_output.noutfiles] = malloc(strlen(filename));
  if (boinc_output.names[boinc_output.noutfiles] == NULL)
    LogPrintf (LOG_FATAL, "ERROR: Can't allocate output filename '%s'\n", filename);
  strcpy(boinc_output.names[boinc_output.noutfiles],filename);
  boinc_output.noutfiles++;
}



/* to be inmplemented... */
void set_checkpoint(char*filename,double rac,double dec, long tpl_count, long tpl_total);
void get_checkpoint(char*filename);
void show_progress(double rac, double dec, long tpl_count, long tpl_total);



/* check if given file is a zip-file by checking the zip-magic header 'PK\003\044'
 * RETURN: 1 = true, 0 = false, -1 = error
 */
static int is_zipped ( const char *fname ) {
  FILE *fp;
  CHAR zip_magic[] = {'P', 'K', 3, 4 };
  CHAR file_header[4];

  if ( (fp = fopen( fname, "rb")) == NULL ) {
    LogPrintf (LOG_CRITICAL, "Failed to open '%s' for reading.\n", fname);
    return -1;
  }
  if ( 4 != fread ( file_header, sizeof(CHAR), 4, fp ) ) {
    LogPrintf (LOG_CRITICAL, "Failed to read first 4 bytes from '%s'.\n", fname);
    return -1;
  }
  fclose(fp);

  if ( memcmp ( file_header, zip_magic, 4 ) )
    return 0;	/* false: no zip-file */
  else
    return 1;	/* yep, found magic zip-header */
} /* is_zipped() */



/* unzip a file in-place if it is zipped */
/* TODO: error handling */
static int unzip_if_necessary(char*filename) {
  char zipname[MAX_PATH_LEN];
  int zipped;
  zipped = is_zipped (filename);
  if (zipped<0) {
    return(-1);
  } else if (zipped) { 
    strncpy(zipname,filename,sizeof(zipname));
    strncat(zipname,".zip",sizeof(zipname));
    boinc_delete_file(zipname);
    boinc_rename(filename,zipname);
    boinc_zip(UNZIP_IT,zipname,filename);
    boinc_delete_file(zipname);
  }
  return(0);
}




/* the main function of the BOINC App
*/
main (int argc, char*argv[]) {
  int i,j,l;
  int rargc = argc;
  char **rargv = NULL;
  char tempstr[MAX_PATH_LEN];

  /* PATCH THE COMMAND LINE

     The actual parsing of the command line will be left to the
     main_hierarchical_search(). However, command line arguments
     that can be identified as filenames must be boinc_resolved
     before passing them to the main function.
     We will also look if the files are possibly zipped and unzip
     them as needed.
  */
  rargv = (char**)malloc(argc*sizeof(char*));
  rargv[0] = argv[0];
  /* TODO: ephermis files */
  /* TODO: unzippig */
  for (i=1; i<argc; i++) {
    
    /* config file */
    if (argv[i][0] == '@') {
      rargv[i] = (char*)malloc(MAX_PATH_LEN);
      rargv[i][0] = '@';
      if (boinc_resolve_filename(argv[i]+1,rargv[i]+1,255)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[i]+1);
      }

    /* skygrid file */
    } else if (MATCH_START("--skyGridFile=",argv[i],l)) {
      rargv[i] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[i],argv[i],l);
      if (boinc_resolve_filename(argv[i]+l,rargv[i]+l,MAX_PATH_LEN-l)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve skygrid file '%s'\n", argv[i]+1);
      }

    /* file to return (zip archive) */
    } else if (MATCH_START("--BOINCresfile=",argv[i],l)) {
      if (boinc_resolve_filename(argv[i]+l,boinc_output.resultfile,sizeof(boinc_output.resultfile))) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve skygrid file '%s'\n", argv[i]+1);
      }
      rargc--; /* this argument is not passed to the main worker function */

    /* ephermis files/directory */
    } else if (MATCH_START("--ephemDir=",argv[i],l)) {
      rargv[i] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[i],argv[i],l);
      if (boinc_resolve_filename(argv[i]+l,rargv[i]+l,MAX_PATH_LEN-l)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve skygrid file '%s'\n", argv[i]+1);
      }

    /* SFT files (no unzipping, but dealing with multiple files separated by ';' */
    } else if (0 == strncmp("--DataFiles",argv[i],11)) {
      rargv[i] = (char*)malloc(1024);
      /* copy & skip the "[1|2]=" characters, too */
      strncpy(rargv[i],argv[i],13);
      appc = rargv[i][13];
      startc = argv[i]+13;
      /* look for multiple paths separated by ';' */
      while(endc=srchr(startc,';')) {
	*endc = '\0';
	if (boinc_resolve_filename(startc,appc,255)) {
	  LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", tempstr);
	}
	/* append a ';' */
	appc = appc + strlen(appc) + 1;
	*(appc-1) = ';';
	*appc = '\0';
	/* skip the ';' */
	startc = endc+1;
      }
      /* handle last or only filename */
      if (boinc_resolve_filename(startc,appc,255)) {
	LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve input file '%s'\n", tempstr);
      }

    /* any other argument */
     } else 
      rargv[i] = argv[i];
  } /* for all command line arguments */

  if (!boinc_output.resultfile) {
      LogPrintf (LOG_ERROR, "ERROR: no result file has been specified");

  res=main_hierarchical_search(rargc,rargv);
  if (res) {
    LogPrintf (LOG_ERROR, "ERROR: main worker returned with error '%d'\n",res);
  }

  if(boinc_output.noutfiles == 0)
    LogPrintf (LOG_ERROR, "ERROR: no output file has been specified");
  for(i=0;i<boinc_output.noutfiles;i++)
    if ( boinc_zip(ZIP_IT, boinc_output.resultfile, boinc_output.outfiles[i]) ) {
      LogPrintf (LOG_NORMAL, "WARNING: Can't zip output file '%s'\n", boinc_output.outfiles[i]);
    }
}
