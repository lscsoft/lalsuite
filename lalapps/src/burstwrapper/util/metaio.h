#ifndef METAIO_H
#define METAIO_H

#include <setjmp.h>

typedef int   METAIO_INT_4S;
typedef short METAIO_INT_2S;
typedef long long METAIO_INT_8S;
typedef float METAIO_REAL_4;
typedef double METAIO_REAL_8;
typedef char  METAIO_CHAR;
typedef unsigned char METAIO_CHAR_U;

/*-- This must be kept consistent with "const TypeText" in metaio.c !! --*/
enum METAIO_Type {
    METAIO_TYPE_ILWD_CHAR,
    METAIO_TYPE_ILWD_CHAR_U,
    METAIO_TYPE_INT_4S,
    METAIO_TYPE_INT_2S,
    METAIO_TYPE_INT_8S,
    METAIO_TYPE_LSTRING,
    METAIO_TYPE_REAL_4,
    METAIO_TYPE_REAL_8,
    METAIO_TYPE_UNKNOWN
};

/* A structure to hold variable-sized C-style strings */
struct MetaioString {
    METAIO_CHAR* data;      /* A pointer to the data, null terminated */
    size_t       len;       /* The length of the string */
    size_t       datasize;  /* The length of the memory pointed to by data */
};

typedef struct MetaioString* MetaioStringPtr;

/* A structure to hold variable-sized byte strings */
struct MetaioStringU {
    METAIO_CHAR_U* data;    /* A pointer to the data, which may include null */
    size_t         len;       /* The length of the string */
    size_t         datasize;  /* The length of the memory pointed to by data */
};

typedef struct MetaioStringU* MetaioStringUPtr;

struct MetaioColumn {
    char* name;
    int   data_type;
};

struct MetaioStream {
    char* name;
    char* type;
    char  delimiter;
};

struct MetaioRowElement {
    struct MetaioColumn* col;
    int valid;
    union MetaioRowEltData {
	METAIO_REAL_4        real_4;
	METAIO_REAL_8        real_8;
	METAIO_INT_4S        int_4s;
	METAIO_INT_2S        int_2s;
	METAIO_INT_8S        int_8s;
	struct MetaioString  lstring;
	struct MetaioStringU ilwd_char;
	struct MetaioStringU ilwd_char_u;
    } data;
};

#define METAIOMAXCOLS 100

struct MetaioTable {
    char*                   name;
    char*                   comment;
    struct MetaioColumn     col[METAIOMAXCOLS];
    struct MetaioRowElement elt[METAIOMAXCOLS];
    int                     numcols;
    int                     maxcols;
    struct MetaioStream     stream;
};

struct MetaioLigo_lw {
    char*        name;
    char*        comment;
    struct MetaioTable table;
};

struct MetaioFileRecord {
    char*  name;
    FILE*  fp;
    size_t lineno;
    int nrows;
    char mode;
    char* tablename;
};

typedef struct MetaioFileRecord* MetaioFile;

struct MetaioParseEnvironment {
    struct MetaioFileRecord fileRec;
    MetaioFile              file;
    struct MetaioString     buffer;
    int                     token;
    int                     mierrno;
    struct MetaioString     mierrmsg;
    jmp_buf                 jmp_env;
    struct MetaioLigo_lw    ligo_lw;
};

typedef struct MetaioParseEnvironment* MetaioParseEnv;

/*
  Returns the name of the data type as a string.
  Usually only useful for debugging.
*/
extern
const char* const MetaioTypeText(const int data_type);

/*
  This function parses the contents of filename up to the beginning
  of the first row of data. It fills the parse environment structure
  'env' with information describing the format that MetaioGetRow() will
  expect the rows to be in. It must be called before using MetaioGetRow()

  Returns 0 if successful, non-zero otherwise.
*/
extern
int MetaioOpen(MetaioParseEnv const env, const char* const filename);

/*
  This function is like MetaioOpen, but it allows the user to specify the
  name of the table to be read from a multi-table file.  If a blank table name
  is specified, it reads the first table from the file (i.e. equivalent to
  MetaioOpen).  The table name is case-insensitive.  An error occurs if the
  specified table does not exist in the file.

  Returns 0 if successful, non-zero otherwise.
*/
extern
int MetaioOpenTable(MetaioParseEnv const env, const char* const filename,
		    const char* const tablename);

/*
  Parse the next row in the input file and insert it into 'env',
  overwriting the current row. Row element i can be accessed
  via env->ligo_lw.table.elt[i].

  Returns 1 if a row was obtained, 0 if no row was obtained and
  a negative number if an error was encountered.
*/
extern
int MetaioGetRow(MetaioParseEnv const env);

/*
  Finish off parsing the file (looking for closing tags and so on), close
  the file and free resources owned by 'env'. After calling this, accessing
  any part of 'env' is undefined.

  Returns 0 if successful, or a negative number if an error was encountered.
*/
extern
int MetaioClose(MetaioParseEnv const env);

/*
  Immediately close the file, without trying to parse the rest of it,
  and free resources owned by 'env'. After calling this, accessing
  any part of 'env' is undefined.

  Returns 0 if successful, or a negative number if an error was encountered.
*/
extern
int MetaioAbort(MetaioParseEnv const env);

/*
  Look up a the name of the column with the given index, ignoring any "prefix"
  string(s) delimited by colons.  For instance, if the file contains a column
  named "processgroup:process:program", this routine returns a pointer to the
  first letter of "program".

  Returns a pointer to the column name, or if the column number is invalid.
*/
extern
char *MetaioColumnName( const MetaioParseEnv env, const int icol );

/*
  Find a column by name.  The name comparison is case-insensitive and ignores
  any "prefix" strings, delimited by colons, in the column name in the file.

  Returns the index of the column with the specified name, or -1 if no such
  column is found.
*/
extern
int MetaioFindColumn( const MetaioParseEnv env, const char *name );

/*
  Compares element values in a way which depends on the particular data type.

  Returns 0 if the values are equal; returns -1 if the first value is "less"
  than the second value; returns 1 if the first value is "greater" than the
  second value.
*/
extern
int MetaioCompareElements( struct MetaioRowElement *elt1,
			   struct MetaioRowElement *elt2 );

/*
  Print the value of a row element to the given file (which may be a standard
  stream such as stdout, stderr) as a string.
  
  Returns the number of characters printed.
*/
extern
int MetaioFprintElement(FILE* const f, const struct MetaioRowElement* elt);

/*
  Opens a file for writing, and writes the LIGO_LW header.
  Returns 0 if successful, nonzero if there was an error creating the file.
*/
extern
int MetaioCreate( const MetaioParseEnv env, const char* const filename );

/*
  Copies column definitions, etc., from one metaio environment to another.
  Returns 0 if successful, nonzero if there was an error.
*/
extern
int MetaioCopyEnv( const MetaioParseEnv dest, const MetaioParseEnv source );

/*
  Copies row contents from one metaio stream to another.
  Returns 0 if successful, nonzero if there was an error.
*/
extern
int MetaioCopyRow( const MetaioParseEnv dest, const MetaioParseEnv source );

/*
  Writes out the current row.
  Returns 0 if successful, nonzero if there was an error.
*/
extern
int MetaioPutRow( const MetaioParseEnv env );

#endif /* METAIO_H */

