#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <assert.h>

#include "metaio.h"

#if 0
#define DEBUG
#endif

#ifdef  DEBUG
#define DEBUGMSG1(s)    printf(s##"\n")
#define DEBUGMSG2(s,d)  printf(s, (d))
#else
#define DEBUGMSG1(s)
#define DEBUGMSG2(s,d)
#endif

enum Token {
    CLOSE_COLUMN,
    CLOSE_COMMENT,
    CLOSE_LIGO_LW,
    CLOSE_STREAM,
    CLOSE_TABLE,
    COLUMN,
    COMMENT,
    DELIMITER,
    LIGO_LW,
    NAME,
    STREAM,
    TABLE,
    TYPE,
    GREATER_THAN,
    EQUALS,
    FORWARD_SLASH,
    UNKNOWN,
    END_OF_FILE
};

static
const char* const TokenText[UNKNOWN] = {
    "</COLUMN",
    "</COMMENT",
    "</LIGO_LW",
    "</STREAM",
    "</TABLE",
    "<COLUMN",
    "<COMMENT",
    "DELIMITER",
    "<LIGO_LW",
    "NAME",
    "<STREAM",
    "<TABLE",
    "TYPE",
    ">",
    "=",
    "/"
};

/*-- This must be kept consistent with "enum METAIO_Type" in metaio.h !! --*/
static
const char* const TypeText[METAIO_TYPE_UNKNOWN] = {
    "ILWD:CHAR",
    "ILWD:CHAR_U",
    "INT_4S",
    "INT_2S",
    "INT_8S",
    "LSTRING",
    "REAL_4",
    "REAL_8"
};

static
void string_resize(MetaioStringPtr const str, const size_t len)
{
    const size_t ResizeIncrement = 512;
    const size_t new_datasize = ((len+1)/ResizeIncrement + 1)*ResizeIncrement;

    if (new_datasize != str->datasize)
    {
	str->datasize = new_datasize;
	str->data = realloc(str->data, str->datasize);
    }
}

static
void append_char(MetaioStringPtr const str, const int c)
{
    str->len += 1;

    if (str->len > str->datasize)
    {
	string_resize(str, str->len);
    }

    str->data[str->len - 1] = c;
    str->data[str->len] = '\0';
}

static
void append_cstr(MetaioStringPtr const str, const char* const s)
{
    const size_t s_len = strlen(s);
    const size_t new_len = str->len + s_len;

    if (new_len > str->datasize)
    {
	string_resize(str, new_len);
    }

    str->len = new_len;
    strcpy(&str->data[str->len - s_len], s);
}

static
void stringu_resize(MetaioStringUPtr const str, const size_t len)
{
    const size_t ResizeIncrement = 512;
    const size_t new_datasize = ((len+1)/ResizeIncrement + 1)*ResizeIncrement;

    if (new_datasize != str->datasize)
    {
	str->datasize = new_datasize;
	str->data = realloc(str->data, str->datasize);
    }
}

static
void append_char_u(MetaioStringUPtr const str, const unsigned int c)
{
    str->len += 1;

    if (str->len > str->datasize)
    {
	stringu_resize(str, str->len);
    }

    str->data[str->len - 1] = c;
    str->data[str->len] = '\0';
}

static
void parse_error(MetaioParseEnv const env, const int code,
		 const char* const format, ...)
{
#define ERRBUFSIZE 256

    char errbuf[ERRBUFSIZE + 1];
    va_list args;

    assert(code != 0);

    snprintf(errbuf, ERRBUFSIZE, "error in %s line %d: ",
	     env->file->name, env->file->lineno);
    
    append_cstr(&(env->mierrmsg), errbuf);

    va_start(args, format);
    vsnprintf(errbuf, ERRBUFSIZE, format, args);
    va_end(args);
    
    append_cstr(&(env->mierrmsg), errbuf);

    env->mierrno = code;
    
    DEBUGMSG2("Buffer is <%s>\n", env->buffer.data);

    longjmp(env->jmp_env, code);

#undef ERRBUFSIZE
}

static
int init_parse_env(MetaioParseEnv const env, const char* const filename,
		    const char* mode )
{
    /* Error handling must be initiliased first */
    env->mierrno = 0;

    env->mierrmsg.data = 0;
    env->mierrmsg.len = 0;
    env->mierrmsg.datasize = 0;

    env->file = &(env->fileRec);
    env->file->name = strdup(filename);
    env->file->fp = 0;
    env->file->lineno = 1;
    env->file->nrows = 0;
    env->file->mode = mode[0];
    env->file->tablename = NULL;

    env->token = UNKNOWN;

    env->buffer.data = 0;
    env->buffer.len = 0;
    env->buffer.datasize = 0;

    env->ligo_lw.name = 0;
    env->ligo_lw.comment = 0;

    env->ligo_lw.table.name = 0;
    env->ligo_lw.table.comment = 0;
    env->ligo_lw.table.maxcols = METAIOMAXCOLS;
    env->ligo_lw.table.numcols = 0;
    
    env->ligo_lw.table.stream.name = 0;
    env->ligo_lw.table.stream.type = 0;
    env->ligo_lw.table.stream.delimiter = '\0';

    /*-- Now try to open the file --*/
    switch ( mode[0] )
    {
    case 'r':
      if ((env->file->fp = fopen(filename, "r")) == 0)
	{
	  parse_error(env, -1, "file not found");
	  return 1;
	}
      break;

    case 'w':
      if ((env->file->fp = fopen(filename, "w")) == 0) { return 1; }

    }

    env->file->nrows = 0;

    return 0;
}

static
int destroy_parse_env(MetaioParseEnv const env)
{
    int ret = 0;
    int i = 0;

    for (i = 0; i < env->ligo_lw.table.numcols; i++)
    {
	/* Delete the row data */
	const int type = env->ligo_lw.table.elt[i].col->data_type;

	if (type == METAIO_TYPE_LSTRING
	    && env->ligo_lw.table.elt[i].data.lstring.data != 0)
	{
	    free(env->ligo_lw.table.elt[i].data.lstring.data);
	    env->ligo_lw.table.elt[i].data.lstring.data = 0;
	    env->ligo_lw.table.elt[i].data.lstring.len = 0;
	    env->ligo_lw.table.elt[i].data.lstring.datasize = 0;
	}
	else if (type == METAIO_TYPE_ILWD_CHAR
	    && env->ligo_lw.table.elt[i].data.ilwd_char.data != 0)
	{
	    free(env->ligo_lw.table.elt[i].data.ilwd_char.data);
	    env->ligo_lw.table.elt[i].data.ilwd_char.data = 0;
	    env->ligo_lw.table.elt[i].data.ilwd_char.len = 0;
	    env->ligo_lw.table.elt[i].data.ilwd_char.datasize = 0;
	}
	else if (type == METAIO_TYPE_ILWD_CHAR_U
	    && env->ligo_lw.table.elt[i].data.ilwd_char_u.data != 0)
	{
	    free(env->ligo_lw.table.elt[i].data.ilwd_char_u.data);
	    env->ligo_lw.table.elt[i].data.ilwd_char_u.data = 0;
	    env->ligo_lw.table.elt[i].data.ilwd_char_u.len = 0;
	    env->ligo_lw.table.elt[i].data.ilwd_char_u.datasize = 0;
	}
	
	env->ligo_lw.table.elt[i].col = 0;

	/* Delete the column data */
	if (env->ligo_lw.table.col[i].name != 0)
	{
	    free(env->ligo_lw.table.col[i].name);
	    env->ligo_lw.table.col[i].name = 0;
	}
    }

    /* Delete the stream */
    if (env->ligo_lw.table.stream.name != 0)
    {
	free(env->ligo_lw.table.stream.name);
	env->ligo_lw.table.stream.name = 0;
    }

    if (env->ligo_lw.table.stream.type != 0)
    {
	free(env->ligo_lw.table.stream.type);
	env->ligo_lw.table.stream.type = 0;
    }

    /* Delete the table */
    if (env->ligo_lw.table.name != 0)
    {
	free(env->ligo_lw.table.name);
	env->ligo_lw.table.name = 0;
    }

    if (env->ligo_lw.table.comment != 0)
    {
	free(env->ligo_lw.table.comment);
	env->ligo_lw.table.comment = 0;
    }

    /* Delete the ligo_lw */
    if (env->ligo_lw.name != 0)
    {
	free(env->ligo_lw.name);
	env->ligo_lw.name = 0;
    }

    if (env->ligo_lw.comment != 0)
    {
	free(env->ligo_lw.comment);
	env->ligo_lw.comment = 0;
    }
    
    /* Delete the buffer */
    if (env->buffer.data != 0 && env->buffer.datasize != 0)
    {
	free(env->buffer.data);
	env->buffer.data = 0;
	env->buffer.len = 0;
	env->buffer.datasize = 0;
    }

    /* Delete the file */
    if (env->file->name != 0)
    {
	free(env->file->name);
	env->file->name = 0;
    }

    if (env->file->fp != 0)
    {
	ret = fclose(env->file->fp);
	env->file->fp = 0;
    }

    env->file->nrows = 0;

    return ret;
}

static
int is_in(const char* s, const int c)
{
    while (*s != '\0')
    {
	if (c == *s++)
	{
	    return 1;
	}
    }

    return 0;
}

static
int is_token_char(const int c)
{
    /* Only these characters are allowed in token names */
    return (isalnum(c) || (c == '_') || (c == ':'));
}

/*
  Return the single character corresponding to a character entity string.
*/
static
int character_entity(const char* const s)
{
    int c = 0;

    /* Note that character entities are ALWAYS in lower case */
    if (strcmp(s, "&gt;") == 0)
    {
	c = '>';
    }
    else if (strcmp(s, "&lt;") == 0)
    {
	c = '<';
    }
    else if (strcmp(s, "&amp;") == 0)
    {
	c = '&';
    }
    else if (strcmp(s, "&nbsp;") == 0)
    {
	c = ' ';
    }
    else
    {
	c = '\0';
    }
    
    return c;
}

/*
  A filter for fgetc to guarantee we always count lines correctly
*/
static
int get_char(MetaioFile const f)
{
    const int c = fgetc(f->fp);

    if (c == '\n')
    {
	f->lineno++;
    }

    return c;
}

/*
  A filter for ungetc to guarantee we always count lines correctly
*/
static
int unget_char(MetaioFile const f, const int c)
{
    if (c == '\n')
    {
	f->lineno--;
    }

    return ungetc(c, f->fp);
}

/*
  Return the next non-whitespace character.

  Since this gets used so much I'm inlining the calls to get_char
*/
static
int skip_whitespace(MetaioFile const f)
{
    FILE* const fp = f->fp;
    int c = fgetc(fp);

    while (isspace(c))
    {
	if (c == '\n')
	{
	    f->lineno++;
	}
	c = fgetc(fp);
    }

    return c;
}

static
int fscanf_lstring(MetaioFile const f, MetaioStringPtr const s,
		   const char* const terminator, const char* const specials)
{
    const size_t start_len = s->len;
    int c = get_char(f);

    while ((!is_in(terminator, c)) && (c != EOF))
    {
	if (c == '&')
	{
	    const size_t pos = s->len; /* The position of the '&' in s */
	    append_char(s, c);
	    
	    fscanf_lstring(f, s, ";", "" );

	    /* Consume the terminal semi-colon */
	    append_char(s, get_char(f));

	    /* Identify the character entity */
	    c = character_entity(&(s->data[pos]));
	    s->len -= strlen(&(s->data[pos]));
	}
	else if (c == '\\')
	{
	  /* Check whether the next character is one of the special characters
	     which needs to be escaped.  If so, ignore the backslash; if not,
	     then the backslash is real data, so keep it.
	  */
	  c = get_char(f);
	  if ( ! is_in(specials,c) ) {
	    unget_char(f, c);
	    c = '\\';
	  }

	}
	append_char(s, c);
	c = get_char(f);
    }
  
    /* Check that there was a terminator */
    if (c == EOF)
    {

    }

    /* Unconsume the character that caused us to halt */
    unget_char(f, c);

    return s->len - start_len;
}

/*
  Scan the file for a ilwd:char string until the quote character is reached.
  The quote character us not consumed. An ilwd:char staring consists of a
  combination of the following:

  * Unsigned bytes represented by 3 octal digits eg. \123

  * Non-space character eg. xyz

  * Backslash-escaped characters eg. backslash \\, space \ , comma \,
    and quotes \" or \'.

  The string may be any combination of the above. Spaces are treated as
  separators between two bytes and are ignored, unless escaped using a
  backslash. eg. '\057 a\ b  \"x' would (in ASCII) be translated to 'Wa b"x'.

  Like scanf(), this function returns the number of objects assigned,
  which in this case means the number of Byte objects, *not* the number
  of bytes - that value is contained in the 'len' field of the Byte object.

  If one or more valid bytes is read, 1 is returned.
  
  If no valid bytes are read, 0 is returned.
  
  If EOF is reached before any valid bytes are read, EOF is returned.

  :NOTE: 
  
  A slight flaw in this code is that it will accept octal
  strings that are shorter than 3 characters eg. \10.
*/
static
int fscanf_ilwd_char(MetaioFile const f, MetaioStringUPtr const b,
		     const char* const enders, const char* const specials)
{
    unsigned int val = 0;
    int count = 1;
    int ret   = 0;

    b->len = 0;

    /*-- Get first character which is not a space --*/
    do {
      val = fgetc(f->fp);
    } while ( val == ' ' );

    while ((count == 1) && (val <= 256) && ! is_in(enders,val) && (val != EOF))
    {
	if (val == '&')
	{
	    int i=0;
	    char ch[2]=" ", tmp[8]="&";

	    /* NOTE: no protection against missing terminal semicolon */
	    do {
		i++;
	        ch[0] = get_char(f);
	        strcat( tmp, ch );
	    } while ( ch[0] != ';' && i < 7 );

	    /* Identify the character entity */
	    val = character_entity(tmp);
	}
	else if (val == '\\')
	{
	    val = get_char(f);
	    if (isdigit(val))
	    {
		/*
		  We got something of the form \[digit][digit]...
		  read it as an octal.
		*/
		unget_char(f, (int) val);
		count = fscanf(f->fp, "%3o", &val);
	    }
	    else if ( is_in(specials,val) )
	    {
	      /* We got a properly-escaped character (such as space) */
	      count = 1;	      
	    }
	    else
	    {
	      /*-- This character should not have been escaped --*/
	      return -1;
	    }
	}
	append_char_u(b, val);

	/*-- Get next character, skipping any spaces. --*/
	do {
	  val = fgetc(f->fp);
	} while ( val == ' ' );
    }
    
    /* Unget the character that caused the loop to finish */
    unget_char(f, (int) val);

    /*
      A val of > 256 means that a valid octal string or character was read
      but it was too large to fit in a byte. This invalidates the whole
      string, so we abandon any data we got.
    */
    
    if (val > 256)
    {
	b->len = 0;
	ret = 0;
    }
    else 
    {	
	if (b->len == 0)
	{
	    ret = ((count == EOF) ? (EOF) : (0));
	}
	else
	{
	    ret = 1;
	}
    }

    return ret;
}

int MetaioFprintElement(FILE* const f,
			const struct MetaioRowElement* const elt)
{
    const int data_type = elt->col->data_type;
    int ret = 0;

    switch (data_type)
    {
    case METAIO_TYPE_INT_4S:
	ret = fprintf(f, "%d", elt->data.int_4s);
	break;
    case METAIO_TYPE_INT_2S:
	ret = fprintf(f, "%hd", elt->data.int_2s);
	break;
    case METAIO_TYPE_INT_8S:
	ret = fprintf(f, "%lld", elt->data.int_8s);
	break;
    case METAIO_TYPE_REAL_4:
	ret = fprintf(f, "%.8g", elt->data.real_4);
	break;
    case METAIO_TYPE_REAL_8:
	ret = fprintf(f, "%.17g", elt->data.real_8);
	break;
    case METAIO_TYPE_LSTRING:
	ret = fprintf(f, "%s", elt->data.lstring.data);
	break;
    case METAIO_TYPE_ILWD_CHAR:
	{
	    size_t i = 0;
	    for (i = 0; i < elt->data.ilwd_char.len; i++)
	    {
		ret += fprintf(f, "%c", elt->data.ilwd_char.data[i]);
	    }
	}
	break;
    case METAIO_TYPE_ILWD_CHAR_U:
	{
	    size_t i = 0;
	    for (i = 0; i < elt->data.ilwd_char_u.len; i++)
	    {
		ret += fprintf(f, "\\%03o", elt->data.ilwd_char_u.data[i]);
	    }
	}
	break;
    default:
	ret = fprintf(f, "{unknown}");
	break;
    }

    return ret;
}

static
int match_token_string(const char* const s)
{
    int i = 0;
    for (i = 0; i < UNKNOWN; i++)
    {
	if (strcasecmp(TokenText[i], s) == 0)
	{
	    return i;
	}
    }

    return UNKNOWN;
}

static
int match_type_string(const char* const s)
{
    char stemp[64];

    int i = 0;
    for (i = 0; i < METAIO_TYPE_UNKNOWN; i++)
    {
	if (strcasecmp(TypeText[i], s) == 0)
	{
	    return i;
	}
    }

    /*-- If we get here, then the type did not match any known types --*/
    /*-- Special case: if the Stream is empty, then LDAS apparently writes the
      type as "char" or "char_u" rather than "ilwd:char" or "ilwd:char_u" --*/
    strcpy( stemp, "ilwd:" );
    strcat( stemp, s );
    for (i = 0; i < METAIO_TYPE_UNKNOWN; i++)
    {
	if (strcasecmp(TypeText[i], stemp) == 0)
	{
	    return i;
	}
    }

    return METAIO_TYPE_UNKNOWN;
}

static
const char* const token_text(const int token)
{
    if (token >= 0 && token < UNKNOWN)
    {
	return TokenText[token];
    }
    else if (token == END_OF_FILE)
    {
	return "END_OF_FILE";
    }
    else
    {
	return "UNKNOWN";
    }
}

const char* const MetaioTypeText(const int type)
{
    if (type >= 0 && type < METAIO_TYPE_UNKNOWN)
    {
	return TypeText[type];
    }
    else
    {
	return "TYPE_UNKNOWN";
    }
}

static
int read_token_string(MetaioParseEnv const env, MetaioStringPtr const s)
{
    const size_t start_len = s->len;
    int c = get_char(env->file);

    while (is_token_char(c) && (c != EOF))
    {
	append_char(s, c);
	c = get_char(env->file);
    }

    /* Unconsume the character that caused us to halt */
    unget_char(env->file, c);

    return s->len - start_len;
}

static
void read_literal_or_token(MetaioParseEnv const env, MetaioStringPtr const s)
{
    const int c = skip_whitespace(env->file);

    if (c == '\"' || c == '\'')
    {
	char quote[2] = "\0";
	quote[0] = c;
	fscanf_lstring(env->file, s, quote, "");

	/* Consume the terminating quote */
	get_char(env->file);
    }
    else
    {
	append_char(s, c);
	read_token_string(env, s);
    }
}

static
void get_next_token(MetaioParseEnv const env)
{
    const int c = skip_whitespace(env->file);

    env->buffer.len = 0;
    append_char(&(env->buffer), c);

    switch(c)
    {
    case '<':
	append_char(&(env->buffer), skip_whitespace(env->file));
	read_token_string(env, &(env->buffer));
	env->token = match_token_string(env->buffer.data);
	break;
    case '>':
	env->token = GREATER_THAN;
	break;
    case '=':
	env->token = EQUALS;
	break;
    case '/':
	env->token = FORWARD_SLASH;
	break;
    case EOF:
	env->token = END_OF_FILE;
	break;
    default:
	read_token_string(env, &(env->buffer));
	env->token = match_token_string(env->buffer.data);
	break;
    }
}

static
void unexpected_token(MetaioParseEnv const env, const int expected)
{
    parse_error(env, -1, "got %s when expecting %s\n",
		token_text(env->token), token_text(expected));
}

static
void match(MetaioParseEnv const env, const int token)
{
    if (env->token == token)
    {
	get_next_token(env);
    }
    else
    {
	unexpected_token(env, token);
    }
}

static
void name(MetaioParseEnv const env, char** const s)
{
    switch(env->token)
    {
    case NAME:
	match(env, NAME);
	env->buffer.len = 0;
	read_literal_or_token(env, &(env->buffer));
	*s = strdup(env->buffer.data);
	match(env, EQUALS);
	break;
    default:
	*s = strdup("");
	break;
    }
}

static
void type(MetaioParseEnv const env, char** const s)
{
  switch(env->token)
    {
    case TYPE:

      match(env, TYPE);
      env->buffer.len = 0;
      read_literal_or_token(env, &(env->buffer));
      *s = strdup(env->buffer.data);
      match(env, EQUALS);
      break;

    default:
      *s = strdup("Local");
      break;
    }

  return;
}

static
void data_type(MetaioParseEnv const env, int* const type)
{
    match(env, TYPE);

    env->buffer.len = 0;
    read_literal_or_token(env, &(env->buffer));

    *type = match_type_string(env->buffer.data);
    if (*type == METAIO_TYPE_UNKNOWN)
    {
	parse_error(env, -1, "unknown data type \"%s\"\n", env->buffer.data);
    }

    match(env, EQUALS);
}

static
void delimiter(MetaioParseEnv const env, char* const c)
{
  switch(env->token)
    {
    case DELIMITER:

      match(env, DELIMITER);

      env->buffer.len = 0;
      read_literal_or_token(env, &(env->buffer));

      if (env->buffer.len > 1)
	{
	  parse_error(env, -1, 
		    "invalid delimiter \"%s\", must be a single character\n", 
		    env->buffer.data);
	}
    
      *c = env->buffer.data[0];
    
      /*
	Some delimiter characters are banned, including most things that
	can appear in a row element.
      */
      if (*c == '\0' || is_in("<\"\'\\0123456789+-.Ee", *c))
	{
	  parse_error(env, -1, 
		      "character <%c> is invalid as a delimiter\n", *c);
	}
	
      match(env, EQUALS);

      break;

    default:
      *c = ',';
      break;
    }

  return;
}

static
void comment(MetaioParseEnv const env, char** const s)
{
    switch(env->token)
    {
    case COMMENT:
	match(env, COMMENT);
	env->buffer.len = 0;
	fscanf_lstring(env->file, &(env->buffer), "<", "");
	*s = strdup(env->buffer.data);
	match(env, GREATER_THAN);
	match(env, CLOSE_COMMENT);
	match(env, GREATER_THAN);
	break;
    default:
	*s = strdup("");
	break;
    }
}

static
void column_attr(MetaioParseEnv const env)
{
    const int colnum = env->ligo_lw.table.numcols - 1;

    name(env, &(env->ligo_lw.table.col[colnum].name));
    DEBUGMSG2("NAME = \"%s\"\n", env->ligo_lw.table.col[colnum].name);

    data_type(env, &(env->ligo_lw.table.col[colnum].data_type));
    DEBUGMSG2("DATA TYPE = \"%s\"\n",
	      MetaioTypeText(env->ligo_lw.table.col[colnum].data_type));

    env->ligo_lw.table.elt[colnum].col = &(env->ligo_lw.table.col[colnum]);
    env->ligo_lw.table.elt[colnum].data.int_4s = 0;
    env->ligo_lw.table.elt[colnum].data.int_2s = 0;
    env->ligo_lw.table.elt[colnum].data.int_8s = 0LL;
    env->ligo_lw.table.elt[colnum].data.real_4 = 0.0;
    env->ligo_lw.table.elt[colnum].data.real_8 = 0.0;
    env->ligo_lw.table.elt[colnum].data.lstring.data = 0;
    env->ligo_lw.table.elt[colnum].data.lstring.len = 0;
    env->ligo_lw.table.elt[colnum].data.lstring.datasize = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char.data = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char.len = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char.datasize = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char_u.data = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char_u.len = 0;
    env->ligo_lw.table.elt[colnum].data.ilwd_char_u.datasize = 0;

    match(env, FORWARD_SLASH);
    match(env, GREATER_THAN);
}

static
void column(MetaioParseEnv const env)
{
    if (env->ligo_lw.table.numcols < env->ligo_lw.table.maxcols)
    {
	env->ligo_lw.table.numcols++;
    }
    else
    {
	parse_error(env, -1, "number of columns exceeds %d\n",
		    env->ligo_lw.table.maxcols);
    }

    match(env, COLUMN);
    column_attr(env);
}

static
void match_numeric(MetaioParseEnv const env,
		   struct MetaioRowElement* const elt)
{
    const int data_type = elt->col->data_type;
    int count = 0;

    /*
      Need to skip whitespace with this function, since I don't
      want to fscanf() to eat newlines without incrementing the
      line number
    */
    const int c = skip_whitespace(env->file);
    unget_char(env->file, c);

    switch(data_type)
    {
    case METAIO_TYPE_REAL_4:
        elt->data.real_4 = 0.0;
	count = fscanf(env->file->fp, "%f", &(elt->data.real_4));
	break;
    case METAIO_TYPE_REAL_8:
        elt->data.real_8 = 0.0;
	count = fscanf(env->file->fp, "%lf", &(elt->data.real_8));
	break;
    case METAIO_TYPE_INT_4S:
        elt->data.int_4s = 0;
	count = fscanf(env->file->fp, "%d", &(elt->data.int_4s));
	break;
    case METAIO_TYPE_INT_2S:
        elt->data.int_2s = 0;
	count = fscanf(env->file->fp, "%hd", &(elt->data.int_2s));
	break;
    case METAIO_TYPE_INT_8S:
        elt->data.int_8s = 0LL;
	count = fscanf(env->file->fp, "%lld", &(elt->data.int_8s));
	break;
    default:
	/* The row element has a non-numeric data type */
	fprintf(stderr, "BUG at %s line %d\n", __FILE__, __LINE__);
	count = 0;
	break;
    }

    /*
      A zero count indicates that there was nothing which looked like
      a numeric value.  So just take this to be a null value.  (NOT an error.)

    */ 
    if (count > 0) { elt->valid = 1; }
}

static
void match_lstring(MetaioParseEnv const env,
		   struct MetaioRowElement* const elt)
{
    char delim = env->ligo_lw.table.stream.delimiter;
    char *cptr;
    int c = skip_whitespace(env->file);

    char enders[3];
    char specials[2];

    enders[0] = delim;
    enders[1] = '<';
    enders[2] = '\0';
    specials[0] = delim;
    specials[1] = '\0';
    
    if (c == '\"')
    {
	elt->data.lstring.len = 0;
	/*-- Read everything up to the next non-escaped delimiter --*/
	fscanf_lstring(env->file, &(elt->data.lstring), enders, specials);
	elt->valid = 1;

	/*-- Ignore everything from the terminating quote onward --*/
	cptr = strrchr( elt->data.lstring.data, c );
	if ( cptr == 0 ) {
	  parse_error(env, -1, "unmatched quote when reading lstring\n" );
	} else {
	  *cptr = '\0';
	  elt->data.lstring.len = (size_t) (cptr - elt->data.lstring.data);
	}

    }
    else if ( is_in(enders,c) )
    {
	/* Whitespace between two delimiters maps to the null string */
	unget_char(env->file, c);
	elt->data.lstring.len = 0;
	/* Make sure there is some memory allocated for the string */
	if ( elt->data.lstring.datasize == 0 ) {
	  string_resize( &(elt->data.lstring), 1 );
	}
	strcpy(elt->data.lstring.data, "");
    }
    else
    {
	/* An lstring *must* be quoted. May as well unget what we read. */
	unget_char(env->file, c);
	parse_error(env, -1, "missing quote when reading lstring\n");
    }
}

static
void match_ilwd_char_u(MetaioParseEnv const env,
		       struct MetaioRowElement* const elt)
{
    char delim = env->ligo_lw.table.stream.delimiter;
    unsigned char *cptr;
    int c = skip_whitespace(env->file);
    
    char enders[3];
    char specials[4];

    enders[0] = delim;
    enders[1] = '<';
    enders[2] = '\0';
    specials[0] = delim;
    specials[1] = '\\';
    specials[2] = ' ';
    specials[3] = '\0';

    if (c == '\"')
    {
	const int quote = c;
	int count = 0;
	
	count = fscanf_ilwd_char(env->file, &(elt->data.ilwd_char_u),
				 enders, specials);
	if ( count < 0 ) {
	    parse_error(env, -1, "Error parsing ilwd:char_u element\n");
	}
	elt->valid = 1;

	/*-- Ignore everything from the terminating quote onward --*/
	for ( cptr = elt->data.ilwd_char_u.data+elt->data.ilwd_char_u.len-1;
	      cptr >= elt->data.ilwd_char_u.data; cptr-- ) {
	  if ( *cptr == quote ) { break; }
	}

	if ( cptr < elt->data.ilwd_char_u.data ) {
	  parse_error(env, -1, "unmatched quote when reading ilwd:char_u\n" );
	} else {
	  *cptr = '\0';
	  elt->data.ilwd_char_u.len =
	    (size_t) (cptr - elt->data.ilwd_char_u.data);
	}

    }
    else if ( is_in(enders,c) )
    {
	/* Whitespace between two delimiters maps to an empty byte string */
	unget_char(env->file, c);
	elt->data.ilwd_char_u.len = 0;
    }
    else
    {
	/* An ilwd:char_u *must* be quoted. May as well unget what we read. */
	unget_char(env->file, c);
	parse_error(env, -1, "missing quote when reading ilwd:char_u\n");
    }
}


static
void match_ilwd_char(MetaioParseEnv const env,
		     struct MetaioRowElement* const elt)
{
    char delim = env->ligo_lw.table.stream.delimiter;
    unsigned char *cptr;
    int c = skip_whitespace(env->file);
    
    char enders[3];
    char specials[4];

    enders[0] = delim;
    enders[1] = '<';
    enders[2] = '\0';
    specials[0] = delim;
    specials[1] = '\\';
    specials[2] = ' ';
    specials[3] = '\0';

    if (c == '\"')
    {
	const int quote = c;
	int count = 0;

	count = fscanf_ilwd_char(env->file, &(elt->data.ilwd_char),
				 enders, specials);
	if ( count < 0 ) {
	    parse_error(env, -1, "Error parsing ilwd:char element\n");
	}
	elt->valid = 1;

	/*-- Ignore everything from the terminating quote onward --*/
	for ( cptr = elt->data.ilwd_char.data+elt->data.ilwd_char.len-1;
	      cptr >= elt->data.ilwd_char.data; cptr-- ) {
	  if ( *cptr == quote ) { break; }
	}

	if ( cptr < elt->data.ilwd_char.data ) {
	  parse_error(env, -1, "unmatched quote when reading ilwd:char\n" );
	} else {
	  *cptr = '\0';
	  elt->data.ilwd_char.len =
	    (size_t) (cptr - elt->data.ilwd_char.data);
	}

    }
    else if ( is_in(enders,c) )
    {
	/* Whitespace between two delimiters maps to an empty byte string */
	unget_char(env->file, c);
	elt->data.ilwd_char.len = 0;
    }
    else
    {
	/* An ilwd:char *must* be quoted. May as well unget what we read. */
	unget_char(env->file, c);
	parse_error(env, -1, "missing quote when reading ilwd:char\n");
    }
}

static
void match_delimiter(MetaioParseEnv const env)
{
    const int c = skip_whitespace(env->file);

    if (c != env->ligo_lw.table.stream.delimiter)
    {
	unget_char(env->file, c);
	parse_error(env, -1, "expected delimiter character \'%c\'\n",
		    env->ligo_lw.table.stream.delimiter);
    }
}

static
void row_element(MetaioParseEnv const env, struct MetaioRowElement* const elt)
{
    const int data_type = elt->col->data_type;

    /*-- Mark this item as invalid, until we successfully parse it --*/
    elt->valid = 0;

    switch (data_type)
    {
    case METAIO_TYPE_REAL_4:
    case METAIO_TYPE_REAL_8:
    case METAIO_TYPE_INT_4S:
    case METAIO_TYPE_INT_2S:
    case METAIO_TYPE_INT_8S:
	match_numeric(env, elt);
	break;
    case METAIO_TYPE_LSTRING:
	match_lstring(env, elt);
	break;
    case METAIO_TYPE_ILWD_CHAR:
	match_ilwd_char(env, elt);
	break;
    case METAIO_TYPE_ILWD_CHAR_U:
	match_ilwd_char_u(env, elt);
	break;
    default:
	parse_error(env, -1, "row element has unknown type\n");
	break;
    }
}

static
int row(MetaioParseEnv const env)
{
    /* Peek ahead to the next non-whitespace character */
    const int c = skip_whitespace(env->file);

    unget_char(env->file, c);

    if (c != '<')
    {
	const int numcols = env->ligo_lw.table.numcols;
	int col = 0;

	/* Process the whole row */
	row_element(env, &(env->ligo_lw.table.elt[0]));

	for (col = 1; col < numcols; col++)
	{
	    match_delimiter(env);
	    row_element(env, &(env->ligo_lw.table.elt[col]));
	}

	return 1;
    }
    else
    {
	return 0;
    }
}

static
void stream_attr(MetaioParseEnv const env)
{
    delimiter(env, &(env->ligo_lw.table.stream.delimiter));
    DEBUGMSG2("DELIMITER = \"%c\"\n", env->ligo_lw.table.stream.delimiter);

    name(env, &(env->ligo_lw.table.stream.name));
    DEBUGMSG2("NAME = \"%s\"\n", env->ligo_lw.table.stream.name);

    type(env, &(env->ligo_lw.table.stream.type));
    DEBUGMSG2("TYPE = \"%s\"\n", env->ligo_lw.table.stream.type);

    /*
      This is a special case. We would normally do a match to
      GREATER_THAN but that would also trigger an advance to the
      next token. Since we expect to be reading data from the
      stream object after this, the next data read shouldn't be tokenized.
      Instead we match the GREATER_THAN after we return from reading
      the stream

      It is a little annoying that the stream_attr function needs to
      be clairvoyant but in practise the initial parsing will only read
      to here anyway. Reading the rows will be done on an individual basis
      by the user.
      
      Ideally what we would do is unget all the text read by match(), but
      unfortunately ungetting is only guaranteed for one character.
    */

    /*-- If the stream is empty, the end of the tag will be "/>" --*/
    if (env->token == FORWARD_SLASH) {
      get_next_token(env);
    }

    if (env->token != GREATER_THAN)
    {
	unexpected_token(env, env->token);
    }
}

static
void stream(MetaioParseEnv const env)
{
    match(env, STREAM);
    stream_attr(env);
}

static
void table_attr(MetaioParseEnv const env)
{
    name(env, &(env->ligo_lw.table.name));
    DEBUGMSG2("NAME = \"%s\"\n", env->ligo_lw.table.name); 

    match(env, GREATER_THAN);
}

static
void table_body(MetaioParseEnv const env)
{
    comment(env, &(env->ligo_lw.table.comment));
    DEBUGMSG2("COMMENT = \"%s\"\n", env->ligo_lw.table.comment); 

    while (env->token == COLUMN)
    {
	column(env);
    }

    stream(env);
}

static
void table(MetaioParseEnv const env)
{
    int len, cmp;
    char *namePtr;
    char *colonPtr;
    char *cptr;

    while( 1 ) {
        match(env, TABLE);
        table_attr(env);

	/*-- If no particular table name was specified, use this table --*/
	if ( env->file->tablename == NULL ) break;
	if ( strlen(env->file->tablename) == 0 ) break;

	/*-- Get a pointer to the table name --*/
	namePtr = env->ligo_lw.table.name;

	/*-- Discard the ":table" postfix, if present --*/
	colonPtr = NULL;
	len = strlen(namePtr);
	if ( len > 6 ) {
	  if ( strcasecmp( namePtr+len-6, ":table" ) == 0 ) {
	    /*-- Temporarily replace the colon in ":table" with a null --*/
	    colonPtr = namePtr+len-6;
	    *colonPtr = '\0';
	  }
	}

	/*-- Skip over any colon-delimited prefixes --*/
	cptr = strrchr( namePtr, ':' );
	if ( cptr != NULL ) {
	  namePtr = cptr + 1;
	}

	/*-- Now do a case-insensitive comparison --*/
	cmp = strcasecmp( namePtr, env->file->tablename );

	/*-- Restore the colon in ":table", if we nulled it --*/
	if ( colonPtr ) {
	  *colonPtr = ':';
	}

	/*-- If name matches, then we have the table we want --*/
	if ( cmp == 0 ) break;

	/*-- If we get here, then table name doesn't match --*/
	/*-- Skip to the end of the stream for this table --*/
	while ( env->token != CLOSE_STREAM ) {
	  get_next_token(env);
	}
	get_next_token(env);
	match(env, GREATER_THAN);

	/* Close the <TABLE> */
	match(env, CLOSE_TABLE);
	match(env, GREATER_THAN);
    }

    table_body(env);
}

static
void ligo_lw_attr(MetaioParseEnv const env)
{
    name(env, &(env->ligo_lw.name));
    DEBUGMSG2("NAME = \"%s\"\n", env->ligo_lw.name); 

    match(env, GREATER_THAN);
}

static
void ligo_lw_body(MetaioParseEnv const env)
{
    comment(env, &(env->ligo_lw.comment));
    DEBUGMSG2("COMMENT = \"%s\"\n", env->ligo_lw.comment); 

    table(env);
}

static
void ligo_lw(MetaioParseEnv const env)
{
    match(env, LIGO_LW);
    ligo_lw_attr(env);
    ligo_lw_body(env);
}

static
void leading_junk(MetaioParseEnv const env)
{
    get_next_token(env);

    while (env->token != LIGO_LW && env->token != END_OF_FILE)
    {
	get_next_token(env);
    }
}


static
void putstring( FILE *fp, const char *string, const char *specials )
{
  char *work;
  size_t slen;
  char *sptr;
  char schar;

  work = (char *) string;
  while ( *work != '\0' ) {
    /*-- Find first character in the specials list --*/
    slen = strcspn( work, specials );
    sptr = work + slen;

    /*-- Temporarily replace the special character with a null --*/
    schar = *sptr;
    *sptr = '\0';

    /*-- Output everything up to the special character --*/
    if ( slen > 0 ) fprintf( fp, "%s", work );

    /*-- If special character is a null, then we're done with this string --*/
    if ( schar == '\0' ) break;

    /*-- Restore the special character --*/
    *sptr = schar;

    /*-- Now output the special character, suitably escaped --*/
    switch ( schar )
    {
    case '>':  fputs( "&gt;", fp ); break;
    case '<':  fputs( "&lt;", fp ); break;
    case '&':  fputs( "&amp;", fp ); break;
    default:  fprintf( fp, "\\%c", schar );
    }

    /*-- Update pointer --*/
    work = sptr + 1;
  }

  return;
}


static
void putbinary( FILE *fp, const METAIO_CHAR_U* buf, const size_t len )
{
  int i;

  for ( i=0; i<len; i++ ) {
    fprintf( fp, "\\%03o", buf[i] );
  }

  return;
}


static
void putelement( FILE *fp, struct MetaioRowElement *elt )
{
  /*-- If value is null, just return --*/
  if ( elt->valid == 0 ) return;

  switch ( elt->col->data_type )
  {
  case METAIO_TYPE_INT_4S:
    fprintf( fp, "%d", elt->data.int_4s );
    break;
  case METAIO_TYPE_INT_2S:
    fprintf( fp, "%hd", elt->data.int_2s );
    break;
  case METAIO_TYPE_INT_8S:
    fprintf( fp, "%lld", elt->data.int_8s );
    break;
  case METAIO_TYPE_REAL_4:
    fprintf( fp, "%.8g", elt->data.real_4 );
    break;
  case METAIO_TYPE_REAL_8:
    fprintf( fp, "%.17g", elt->data.real_8 );
    break;
  case METAIO_TYPE_LSTRING:
    fputc( '"', fp );
    putstring( fp, elt->data.lstring.data, "," );
    fputc( '"', fp );
    break;
  case METAIO_TYPE_ILWD_CHAR:
    fputc( '"', fp );
    /*-- Check whether value contains any nulls --*/
    if ( strlen(elt->data.ilwd_char.data) == elt->data.ilwd_char.len ) {
      /*-- No nulls, so interpret this as a string --*/
      putstring( fp, elt->data.ilwd_char.data, " ,\\" );
    } else {
      putbinary( fp, elt->data.ilwd_char.data, elt->data.ilwd_char.len );
    }
    fputc( '"', fp );
    break;
  case METAIO_TYPE_ILWD_CHAR_U:
    fputc( '"', fp );
    putbinary( fp, elt->data.ilwd_char_u.data, elt->data.ilwd_char_u.len );
    fputc( '"', fp );
    break;
  }

  return;
}


static
void putheader( const MetaioParseEnv env )
{
  FILE *fp = env->file->fp;
  struct MetaioTable *table = &(env->ligo_lw.table);
  int icol;

  fprintf( fp, "\n   <Table Name=\"%s\">", table->name );

  /*-- Comment (if any) --*/
  if ( table->comment ) {
    fprintf( fp, "\n      <Comment>" );
    putstring( fp, table->comment, "><&" );
    fprintf( fp, "</Comment>" );
  }

  /*-- Column elements --*/
  for ( icol = 0; icol < table->numcols; icol++ ) {
    fprintf( fp, "\n      <Column Name=\"%s\" Type=\"%s\"/>",
	     table->col[icol].name, TypeText[table->col[icol].data_type] );
  }

  /*-- Stream element tag --*/
  fprintf( fp, "\n      <Stream Name=\"%s\" Type=\"%s\" Delimiter=\"%c\">",
	   table->name, "Local", ',' );

  return;
}


int MetaioOpen(MetaioParseEnv const env, const char* const filename)
{
    int ret = 0;

    ret = setjmp(env->jmp_env);
    
    if (ret == 0)
    {
	init_parse_env(env, filename, "r");
	leading_junk(env);
	ligo_lw(env);
    }

    return ret;
}

int MetaioOpenTable(MetaioParseEnv const env, const char* const filename,
		    const char* const tablename)
{
    int ret = 0;

    ret = setjmp(env->jmp_env);
    
    if (ret == 0)
    {
	init_parse_env(env, filename, "r");
	if ( tablename != NULL ) {
	  env->file->tablename = strdup(tablename);
	}
	leading_junk(env);
	ligo_lw(env);
    }

    return ret;
}

int MetaioGetRow(MetaioParseEnv const env)
{
    int ret = 0;

    /* Reset the position to which errors jump */
    ret = setjmp(env->jmp_env);
    
    if (ret == 0)
    {
	if (row(env) != 0)
	{
	    const int c = skip_whitespace(env->file);
	    if (c != env->ligo_lw.table.stream.delimiter)
	    {
		unget_char(env->file, c);
	    }
	    /*-- Increment the count of the number of rows --*/
	    env->file->nrows++;
	    return 1;  /*-- Success --*/
	}
	else
	{
	    return 0;  /*-- End of file --*/
	}
    } else {
        /* We longjmp'ed to here */
        return -1;  /*-- Parsing error --*/
    }
}

int MetaioClose(MetaioParseEnv const env)
{
    int ret = 0;

    /* If file has already been closed, just return */
    if ( env->file->fp == 0 ) {
	return 0;
    }

    /* Reset the position to which errors jump */
    ret = setjmp(env->jmp_env);

    /* Try to parse the rest of the file (input files only).  If there is an
       error, then longjmp back here and just close the file. */
    if (ret == 0 && env->file->fp != 0 && env->file->mode == 'r' )
    {
	/* Consume any remaining rows */
	if (row(env) != 0)
	{
	    const int c = skip_whitespace(env->file);
	    if (c != env->ligo_lw.table.stream.delimiter)
	    {
		unget_char(env->file, c);
	    }
	}

	/*-- Skip any tokens found until we get to the CLOSE_TABLE.  The
	  actual set of tokens encountered will vary depending on whether the
	  stream was empty (in which case there is no separate </STREAM> tag)
	  or non-empty. --*/
	while ( env->token != CLOSE_TABLE ) {
	  get_next_token(env);
	}
	
	/* Close the <TABLE> */
	match(env, CLOSE_TABLE);
	match(env, GREATER_THAN);
	
	/* Close the <LIGO_LW> */
	match(env, CLOSE_LIGO_LW);
	match(env, GREATER_THAN);
	match(env, END_OF_FILE);
    }

    /*-- Handle output file --*/
    if (ret == 0 && env->file->fp != 0 && env->file->mode == 'w' ) {
      /*-- If header was never written out, write it out now --*/
      if ( env->file->nrows == 0 ) putheader( env );

      /*-- Write out the file trailer stuff before closing --*/
      fprintf( env->file->fp, "\n      </Stream>\n   </Table>\n</LIGO_LW>" );
    }

    return destroy_parse_env(env);
}


int MetaioAbort(MetaioParseEnv const env)
{
    return destroy_parse_env(env);
}


char *MetaioColumnName( const MetaioParseEnv env, const int icol )
/*--
  Written 31 Jan 2001 by Peter Shawhan.
  Returns a pointer to the column name, ignoring any "prefix" string(s)
  delimited by colons.  For instance, if the file contains a column named
  "processgroup:process:program", this routine returns a pointer to the
  first letter of "program".  If an invalid column number is passed,
  this routine returns 0.
--*/
{
  char *colnamePtr, *colonPtr;

  if ( icol < 0 || icol > env->ligo_lw.table.numcols-1 ) { return 0; }

  colnamePtr =  env->ligo_lw.table.col[icol].name;
  colonPtr = strrchr( colnamePtr, ':' );
  if ( colonPtr == 0 ) {
    return colnamePtr;
  } else {
    return (colonPtr+1);
  }
}


int MetaioFindColumn( const MetaioParseEnv env, const char *name )
/*--
  Written 31 Jan 2001 by Peter Shawhan.
  Returns the index of the column with the specified name, or -1 if no such
  column is found.  The name comparison is case-insensitive and ignores any
  "prefix" strings, delimited by colons, in the column name in the file.
--*/
{
  int icol;
  char *cptr;

  for ( icol = 0; icol < env->ligo_lw.table.numcols; icol++ ) {
    cptr = MetaioColumnName( env, icol );
    if ( cptr == 0 ) { return -1; }
    if ( strcasecmp( name, cptr ) == 0 ) {
      return icol;
    }
  }

  /*-- Column was not found --*/
  return -1;
}


int MetaioCompareElements( struct MetaioRowElement *elt1,
			   struct MetaioRowElement *elt2 )
/*--
  Written 1 Feb 2001 by Peter Shawhan.
  Does an appropriate comparison of the element values, which depends on the
  particular data type.  Returns 0 if the values are equal; returns -1 if
  the first value is "less" than the second value; returns 1 if the first
  value is "greater" than the second value.  If the two elements cannot be
  compared (e.g. they have different types), returns 2.
--*/
{
  int retval;
  int type1, type2;
  METAIO_INT_8S ival1, ival2;
  METAIO_REAL_4 rval1, rval2;
  METAIO_REAL_8 dval1, dval2;
  size_t size1, size2, minsize;
  METAIO_CHAR_U *ptr1, *ptr2;

  type1 = elt1->col->data_type;
  type2 = elt2->col->data_type;

  /* lstring and ilwd:char can be compared!  Map type to METAIO_TYPE_LSTRING.*/
  /* int_2s, int_4s, and int_8s can be compared!  Map type to
     METAIO_TYPE_INT_8S for all of them to do the comparison. */
  /* If one item is REAL_4 and other is REAL_8, do comparison as REAL_8. */
  switch ( type1 ) {
  case METAIO_TYPE_LSTRING:
    ptr1 = (METAIO_CHAR_U*) elt1->data.lstring.data;
    size1 = elt1->data.lstring.len;
    break;
  case METAIO_TYPE_ILWD_CHAR:
    ptr1 = elt1->data.ilwd_char.data;
    size1 = elt1->data.ilwd_char.len;
    type1 = METAIO_TYPE_LSTRING;
    break;
  case METAIO_TYPE_ILWD_CHAR_U:
    ptr1 = elt1->data.ilwd_char_u.data;
    size1 = elt1->data.ilwd_char_u.len;
    break;
  case METAIO_TYPE_INT_4S:
    ival1 = (METAIO_INT_8S) elt1->data.int_4s;
    type1 = METAIO_TYPE_INT_8S;
    break;
  case METAIO_TYPE_INT_2S:
    ival1 = (METAIO_INT_8S) elt1->data.int_2s;
    type1 = METAIO_TYPE_INT_8S;
    break;
  case METAIO_TYPE_INT_8S:
    ival1 = elt1->data.int_8s;
    break;
  case METAIO_TYPE_REAL_4:
    if ( type2 == METAIO_TYPE_REAL_8 ) {
      dval1 = (METAIO_REAL_8) elt1->data.real_4;
      type1 = METAIO_TYPE_REAL_8;
    } else {
      rval1 = elt1->data.real_4;
    }
    break;
  case METAIO_TYPE_REAL_8:
    dval1 = elt1->data.real_8;
    break;
  }

  switch ( type2 ) {
  case METAIO_TYPE_LSTRING:
    ptr2 = (METAIO_CHAR_U*) elt2->data.lstring.data;
    size2 = elt2->data.lstring.len;
    break;
  case METAIO_TYPE_ILWD_CHAR:
    ptr2 = elt2->data.ilwd_char.data;
    size2 = elt2->data.ilwd_char.len;
    type2 = METAIO_TYPE_LSTRING;
    break;
  case METAIO_TYPE_ILWD_CHAR_U:
    ptr2 = elt2->data.ilwd_char_u.data;
    size2 = elt2->data.ilwd_char_u.len;
    break;
  case METAIO_TYPE_INT_4S:
    ival2 = (METAIO_INT_8S) elt2->data.int_4s;
    type2 = METAIO_TYPE_INT_8S;
    break;
  case METAIO_TYPE_INT_2S:
    ival2 = (METAIO_INT_8S) elt2->data.int_2s;
    type2 = METAIO_TYPE_INT_8S;
    break;
  case METAIO_TYPE_INT_8S:
    ival2 = elt2->data.int_8s;
    break;
  case METAIO_TYPE_REAL_4:
    if ( type1 == METAIO_TYPE_REAL_8 ) {
      dval2 = (METAIO_REAL_8) elt2->data.real_4;
      type2 = METAIO_TYPE_REAL_8;
    } else {
      rval2 = elt2->data.real_4;
    }
    break;
  case METAIO_TYPE_REAL_8:
    dval2 = elt2->data.real_8;
    break;
  }

  if ( type1 != type2 ) { return 2; }

  /*-- Check whether one or both values are null.  If one is null but the other
    is not, then take the non-null one to be "greater" than the null one. --*/
  if ( elt1->valid == 0 ) {
    if ( elt2->valid == 0 ) { return 0; } else { return -1; }
  } else if ( elt2->valid == 0 ) {
    return 1;
  }

  /*-- Now compare numerical or lexical values --*/
  switch ( type1 ) {
  case METAIO_TYPE_REAL_4:
    if ( rval1 > rval2 ) {
      retval = 1;
    } else if ( rval1 < rval2 ) {
      retval = -1;
    } else {
      retval = 0;
    }
    break;
  case METAIO_TYPE_REAL_8:
    if ( dval1 > dval2 ) {
      retval = 1;
    } else if ( dval1 < dval2 ) {
      retval = -1;
    } else {
      retval = 0;
    }
    break;
  case METAIO_TYPE_INT_8S:
    if ( ival1 > ival2 ) {
      retval = 1;
    } else if ( ival1 < ival2 ) {
      retval = -1;
    } else {
      retval = 0;
    }
    break;
  case METAIO_TYPE_LSTRING:
  case METAIO_TYPE_ILWD_CHAR_U:

    /*-- For lstrings, we want to ignore trailing spaces --*/
    while ( size1 > 0 && *(ptr1+size1-1) == ' ' ) { size1--; }
    while ( size2 > 0 && *(ptr2+size2-1) == ' ' ) { size2--; }

    minsize = ( size1<size2 ? size1 : size2 );
    
    /*-- Note: memcmp compares arrays of UNsigned bytes --*/
    retval = memcmp( ptr1, ptr2, minsize );

    /*-- If arrays are the same up to the point where the shorter one ends,
      then consider the longer array to be "greater" --*/
    if ( retval == 0 ) {
      if ( size1 > size2 ) {
	retval = 1;
      } else if ( size1 < size2 ) {
	retval = -1;
      }
    }
    break;
  }

  return retval;
}


int MetaioCreate( const MetaioParseEnv env, const char* const filename )
/*--
  Written 15 Jul 2002 by Peter Shawhan.
  Opens a file for writing, and writes the LIGO_LW header.
  Returns 0 if successful, nonzero if there was an error creating the file.
--*/
{
  int status;

  status = init_parse_env( env, filename, "w" );
  if ( status ) { return status; }

  /*-- Write out the LIGO_LW header --*/
  /* JS
  fprintf( env->file->fp, "<?xml version='1.0' encoding='utf-8' ?>
<!DOCTYPE LIGO_LW [
<!ELEMENT LIGO_LW ((LIGO_LW|Comment|Param|Table|Array|Stream|IGWDFrame|AdcData|AdcInterval|Time|Detector)*)>
<!ATTLIST LIGO_LW
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED>

<!ELEMENT Comment (#PCDATA)>

<!ELEMENT Param (#PCDATA|Comment)*>
<!ATTLIST Param 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Start CDATA #IMPLIED
          Scale CDATA #IMPLIED
          Unit CDATA #IMPLIED
          DataUnit CDATA #IMPLIED>

<!ELEMENT Table (Comment?,Column*,Stream?)>
<!ATTLIST Table 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED>

<!ELEMENT Column EMPTY>
<!ATTLIST Column
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Unit CDATA #IMPLIED>

<!ELEMENT Array (Dim*,Stream?)>
<!ATTLIST Array 
          Name CDATA #IMPLIED
          Type CDATA #IMPLIED
          Unit CDATA #IMPLIED>

<!ELEMENT Dim (#PCDATA)>
<!ATTLIST Dim 
          Name  CDATA #IMPLIED
          Unit CDATA #IMPLIED
          Start CDATA #IMPLIED
          Scale CDATA #IMPLIED>

<!ELEMENT Stream (#PCDATA)>
<!ATTLIST Stream 
          Name      CDATA #IMPLIED
          Type      (Remote|Local) \"Local\"
          Delimiter CDATA \",\"
          Encoding  CDATA #IMPLIED
          Content   CDATA #IMPLIED>

<!ELEMENT IGWDFrame ((Comment|Param|Time|Detector|AdcData|LIGO_LW|Stream?|Array|IGWDFrame)*)>
<!ATTLIST IGWDFrame 
          Name CDATA #IMPLIED>

<!ELEMENT Detector ((Comment|Param|LIGO_LW)*)>
<!ATTLIST Detector 
          Name CDATA #IMPLIED>

<!ELEMENT AdcData ((AdcData|Comment|Param|Time|LIGO_LW|Array)*)>
<!ATTLIST AdcData 
          Name CDATA #IMPLIED>

<!ELEMENT AdcInterval ((AdcData|Comment|Time)*)>
<!ATTLIST AdcInterval 
          Name CDATA #IMPLIED
          StartTime CDATA #IMPLIED
          DeltaT CDATA #IMPLIED>

<!ELEMENT Time (#PCDATA)>
<!ATTLIST Time 
          Name CDATA #IMPLIED
          Type (GPS|Unix|ISO-8601) \"ISO-8601\">
]>

<LIGO_LW Name=\"ligo:metaio:file\">" );
  */

  return 0;
}


int MetaioCopyEnv( const MetaioParseEnv dest, const MetaioParseEnv source )
/*--
  Copies column definitions, etc., from one metaio environment to another.
  Returns 0 if successful, nonzero if there was an error.
  Written 15 Jul 2002 by Peter Shawhan.
--*/
{
  int icol;
  struct MetaioRowElement* elt;

  dest->ligo_lw.table.name = malloc( strlen(source->ligo_lw.table.name) + 1 );
  strcpy( dest->ligo_lw.table.name, source->ligo_lw.table.name );

  dest->ligo_lw.table.comment =
    malloc( strlen(source->ligo_lw.table.comment) + 1 );
  strcpy( dest->ligo_lw.table.comment, source->ligo_lw.table.comment );

  dest->ligo_lw.table.numcols = source->ligo_lw.table.numcols;

  for ( icol = 0; icol < dest->ligo_lw.table.numcols; icol++ ) {

    dest->ligo_lw.table.col[icol].name =
      malloc( strlen(source->ligo_lw.table.col[icol].name) + 1 );
    strcpy( dest->ligo_lw.table.col[icol].name,
	    source->ligo_lw.table.col[icol].name );

    dest->ligo_lw.table.col[icol].data_type = 
      source->ligo_lw.table.col[icol].data_type;

    elt = &(dest->ligo_lw.table.elt[icol]);

    elt->col = &(dest->ligo_lw.table.col[icol]);
    elt->data.int_4s = 0;
    elt->data.int_2s = 0;
    elt->data.int_8s = 0LL;
    elt->data.real_4 = 0.0;
    elt->data.real_8 = 0.0;
    elt->data.lstring.data = 0;
    elt->data.lstring.len = 0;
    elt->data.lstring.datasize = 0;
    elt->data.ilwd_char.data = 0;
    elt->data.ilwd_char.len = 0;
    elt->data.ilwd_char.datasize = 0;
    elt->data.ilwd_char_u.data = 0;
    elt->data.ilwd_char_u.len = 0;
    elt->data.ilwd_char_u.datasize = 0;
  }

  return 0;
}


int MetaioCopyRow( const MetaioParseEnv dest, const MetaioParseEnv source )
/*--
  Copies row contents from one metaio stream to another.
  Returns 0 if successful, nonzero if there was an error.
  Written 15 Jul 2002 by Peter Shawhan.
--*/
{
  int icol;
  struct MetaioRowElement *selt, *delt;

  for ( icol = 0; icol < dest->ligo_lw.table.numcols; icol++ ) {
    /*-- Get pointers to source and dest elements, for convenience --*/
    selt = &(source->ligo_lw.table.elt[icol]);
    delt = &(dest->ligo_lw.table.elt[icol]);

    delt->valid = selt->valid;

    switch ( dest->ligo_lw.table.col[icol].data_type )
    {
    case METAIO_TYPE_ILWD_CHAR:
      if ( delt->data.ilwd_char.datasize < selt->data.ilwd_char.datasize ) {
	   stringu_resize( &(delt->data.ilwd_char),
			   selt->data.ilwd_char.datasize );
      }
      memcpy( delt->data.ilwd_char.data, selt->data.ilwd_char.data,
	      selt->data.ilwd_char.len );
      delt->data.ilwd_char.len = selt->data.ilwd_char.len;
      break;

    case METAIO_TYPE_ILWD_CHAR_U:
      if ( delt->data.ilwd_char_u.datasize < selt->data.ilwd_char_u.datasize ) {
	   stringu_resize( &(delt->data.ilwd_char_u),
			   selt->data.ilwd_char_u.datasize );
      }
      memcpy( delt->data.ilwd_char_u.data, selt->data.ilwd_char_u.data,
	      selt->data.ilwd_char_u.len );
      delt->data.ilwd_char_u.len = selt->data.ilwd_char_u.len;
      break;

    case METAIO_TYPE_LSTRING:
      if ( delt->data.lstring.datasize < selt->data.lstring.datasize ) {
	   string_resize( &(delt->data.lstring), selt->data.lstring.datasize );
      }
      memcpy( delt->data.lstring.data, selt->data.lstring.data,
	      selt->data.lstring.len );
      delt->data.lstring.len = selt->data.lstring.len;
      break;

    default:
      delt->data = selt->data;
    }

  }

  return 0;
}


int MetaioPutRow( const MetaioParseEnv env )
/*--
  Writes out the current row.
  Returns 0 if successful, nonzero if there was an error.
  Written 15 Jul 2002 by Peter Shawhan.
--*/
{
  FILE *fp = env->file->fp;
  struct MetaioTable *table = &(env->ligo_lw.table);
  int icol;

  /*-- If file is not open for writing, just return --*/
  if ( fp == NULL ) return 0;
  if ( env->file->mode != 'w' ) return 1;

  /*-- If we have not yet written out any rows, write the table header --*/
  if ( env->file->nrows == 0 ) putheader(env);

  /*-- Write out the data for this row --*/
  if ( env->file->nrows > 0 ) fputs( ",", fp );
  for ( icol = 0; icol < table->numcols; icol++ ) {
    if ( icol == 0 ) {
      fputs( "\n         ", fp );
    } else {
      fputs( ",", fp );
    }

    putelement( fp, &(table->elt[icol]) );
  }

  /*-- Increment the count of the number of rows written out --*/
  env->file->nrows++;

  return 0;
}
