#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#define MAXSTR 1024
#define TABSIZE 8
#define BADEXIT 1
#define OKEXIT	0
#define PHILE "file="
#define HASHDEFINE "#define"

extern FILE *ptrLALErrorFile;
extern int lineNum;

typedef void *(func) (void *p);

typedef struct 
tagLALEnviroment
{  
    char *OnFlag     ; /* String that turns the environment on          */
    char *OffFlag    ; /* String that turns the environment off         */
    char *closer     ; /* closeing brace ie ">"                         */
    int   On         ; /* 1 when in evironment, 0 when not              */
    FILE *InFilePtr  ; /* file ptr to input source                      */
    FILE *OutFilePtr ; /* file ptr  for output                          */
    char *fileName   ; /* file Name for output                          */
    char *sourceFile ; /* source file being parsed                      */
    char *allCaps    ; /* string file.c -> FILEC                        */
    char *dfltFile   ; /* string file.H -> fileH.tex                    */
    char *suffix     ; /* string added to filename ie ".tex"            */
    char *errCodePrfx; /* error codes prefix. Should agree with allCaps */
    int  cnvtnVioltn ; /* 1 if some convention is violated, 0 otherwise */
    int  lineCount   ; /* numbr of lines in environment (counting ends) */
    func *Preamble   ; /* a function that prepares top of output        */
    func *PostFix    ; /* a function that prepares end of output        */
} LALEnvironment;

FILE *OpenAFile(char * file, char *readwriteappend , int timestamp  );
int   CloseAFile(FILE *filePtr, int timeStamp );

int LalDocGetLine(char *line, int max, FILE *fpin);
int CheckEnvironment(char *line, LALEnvironment *Env);
int FindNameOfOutPutFile( char *line , LALEnvironment *LaTeXEnv );
int FindDefaultFileName( LALEnvironment *Env );
int ParseCommandLine( int argcount , char **argvector , char *source , 
                      FILE **ptrptrInput);
int LALDocErr( char *message , char *sourceCodeFileName,
               int laldocLineNumber, char *laldocSourceCodeFile,
               int fatal );
int CleanExit( FILE *ptrInput , FILE *ptrOutput , 
               LALEnvironment *LaTeXEnv    ,
               LALEnvironment *VerbatimEnv ,
               LALEnvironment *ErrTableEnv );
int CheckForEnvConflicts( LALEnvironment *LaTeXEnv    ,
                          LALEnvironment *VerbatimEnv ,
                          LALEnvironment *ErrTableEnv );

/* Functions specific to the LaTeX environment */
int SetLaTeXFlags(LALEnvironment *LaTeXEnv, char *source , FILE *in);
int WriteLaTeX(char *line, LALEnvironment *Env);
void *TopOfLaTeXEnv( void *recastToEnv );
void *EndOfLaTeXEnv( void *recastToEnv );

/* Functions specific to the Verbatim environment */
int SetVerbatimFlags( LALEnvironment *VerbatimEnv, char *sourcefile , FILE *in);
int WriteVerbatim(char *line , LALEnvironment *Env );
void *TopOfVerbatimEnv( void *recastToEnv );
void *EndOfVerbatimEnv( void *recastToEnv );
void TabAlign(char *out, char *in, int n, int tabsize);

/* Functions specific to the Error Table environment */
int SetErrTableFlags( LALEnvironment *ErrTableEnv, char *sourcefile , FILE *in);
int WriteErrTable(char *line , LALEnvironment *Env );
void *TopOfErrTableEnv( void *recastToEnv );
void *EndOfErrTableEnv( void *recastToEnv );
int LALDocConstructErrorTable(  LALEnvironment *Env  );

