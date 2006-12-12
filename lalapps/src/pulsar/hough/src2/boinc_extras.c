#define MAX_PATH_LEN 256

/* compare strings s1 and s2 up to the length of s1 (w/o the '\0') and set l to the length */
#define MATCH_START(s1,s2,l) (0 == strncmp(s1,s2,(l=strlen(s1))-1))

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

/* the main function of the BOINC App
*/

/* TODO: boinc_init, boinc_finish */
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
  /* TODO: zipping of output file(s) */
  for (i=1; i<argc; i++) {
    
    /* config file */
    if (argv[i][0] == '@') {
      rargv[i] = (char*)malloc(MAX_PATH_LEN);
      rargv[i][0] = '@';
      if (boinc_resolve_filename(argv[i]+1,rargv[i]+1,255)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve config file '%s'\n", argv[i]+1);
      }

    /* skygrid file */
    } else if (MATCH_START("--skyGridFile=",argv[i],l) {
      rargv[i] = (char*)malloc(MAX_PATH_LEN);
      strncpy(rargv[i],argv[i],l);
      if (boinc_resolve_filename(argv[i]+l,rargv[i]+l,MAX_PATH_LEN-l)) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve skygrid file '%s'\n", argv[i]+1);
      }

    /* file to return (zip archive) */
    } else if (MATCH_START("--BOINCresfile=",argv[i],l) {
      if (boinc_resolve_filename(argv[i]+l,boinc_output.resultfile,sizeof(boinc_output.resultfile))) {
        LogPrintf (LOG_NORMAL, "WARNING: Can't boinc-resolve skygrid file '%s'\n", argv[i]+1);
      }
      rargc--; /* this argument is not passed to the main worker function */

    /* ephermis files/directory */
    } else if (MATCH_START("--ephemDir=",argv[i],l) {
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
    } else {
      rargv[i] = argv[i];
    }
  }
}
