#ifndef _SISTR_H
#define _SISTR_H

/*==========================================================================
SIStr.h - Header file for client API to stream a waveform to the awg front end
Written 28 Sep - 4 Oct 2001 by Peter Shawhan
See comments about compiling and linking at the beginning of the file SIStr.c
==========================================================================*/

#ifdef  __cplusplus
extern "C" {
#endif

/*------ Global variable ---------------------------------------------------*/
#ifdef _SISTR_LIBRARY
int SIStr_debug = 0;
#else
extern int SIStr_debug;
#endif

/*------ Compile-time parameters -------------------------------------------*/

#define SIStr_MAGICVAL 12345678
#define SIStr_MAXCHANNAMELENGTH 64
#define SIStr_MAXCHANLISTSIZE 65536

/* Target "lead time" for sending waveform data, in NANOseconds */
#define SIStr_LEADTIME 13000000000LL

/* Block size in "epoch" units (1/16 sec).  Allowed values are 1 through 16 */
#define SIStr_BLOCKSIZE 16

#define SIStr_MAXBUFSIZE 65536
#define SIStr_MAXBUFS 8

/* Max number of times to try sending same data */
#define SIStr_MAXTRIES 5

/*------ Status codes ------------------------------------------------------*/

#define SIStr_OK        0

/* Status codes returned from RPC calls */
#define SIStr_WFULL     1  /* Data accepted but front-end buffer is now full */
#define SIStr_WDUP      2  /* Data accepted but time was duplicated */
#define SIStr_WGAP      3  /* Data accepted but was not contiguous with prev */
#define SIStr_EBADSLOT -1  /* This awg slot is not set up for stream data */
#define SIStr_EBADDATA -2  /* Invalid data block or size */
#define SIStr_EPAST    -3  /* Time is already past */
#define SIStr_EFUTURE  -4  /* Time is unreasonably far in the future */
#define SIStr_ECONN    -5  /* RPC connection failed */

/* Error codes internal to SIStr */
#define SIStr_EBADARG   -101  /* Bad function argument */
#define SIStr_EBADSTR   -102  /* Stream is not correctly open for writing */
#define SIStr_EBADRATE  -103  /* Invalid sampling rate */
#define SIStr_EGAP      -104  /* Gap detected in stream data */
#define SIStr_EUNINIT   -105  /* Stream was not properly initialized */
#define SIStr_EMALLOC   -107  /* Error allocating memory */
#define SIStr_EOTHER    -108  /* Other error */
#define SIStr_EABORTED  -110  /* Attempted to append data to a stream which had
				 been aborted */
#define SIStr_EBADSTART -111  /* Attempted to set start time to an unreasonable
				 value */
#define SIStr_EINTERNAL -112  /* Unexpected internal error */
#define SIStr_EBUFSIZE  -113  /* Tried to create too large a data buffer */
#define SIStr_ETIMEOUT  -114  /* Timeout while trying to send data */
#define SIStr_ELISTERR  -115  /* Error retrieving channel list */
#define SIStr_ELISTNONE -116  /* Channel list is empty */
#define SIStr_ELISTSIZE -117  /* Channel list is too large */
#define SIStr_EBADCHAN  -118  /* Not a valid excitation channel */
#define SIStr_EDIFFRATE -119  /* Specified rate differs from actual rate */
#define SIStr_ESETSLOT  -120  /* Error setting up an awg slot for channel */
#define SIStr_ESETTP    -121  /* Error setting up test point */
#define SIStr_ESETCOMP  -122  /* Error setting up awgStream component */
#define SIStr_ECLRCOMP  -123  /* Error clearing awgStream component */
#define SIStr_ECLRTP    -124  /* Error clearing test point */
#define SIStr_ECLRSLOT  -125  /* Error freeing awg slot */
#define SIStr_ECLRBOTH  -126  /* Errors clearing test point AND freeing slot */

/*------ Structure definitions ---------------------------------------------*/

typedef struct tagSIStrBuf
{
  int gpstime;  /* Start time of this block (integer number of GPS seconds) */
  int epoch;    /* Start time of this block (epoch counter) */
  int iblock;   /* Block number in sequence, starting with 1 */
  int size;     /* Size of the data array */
  int ndata;    /* Number of values added to the data array so far */
  struct tagSIStrBuf *next;  /* Pointer to next buffer in linked list */
  float *data;
} SIStrBuf;

typedef struct tagSIStream
{
  int magic;    /* "Magic number" to allow sanity checks */
  int id;       /* Internal identifier for this injection stream */
  char channel[SIStr_MAXCHANNAMELENGTH];   /* Channel name */
  int samprate; /* Sampling rate in Hz */
  double starttime; /* Start time of waveform (double-precision GPS seconds) */
  int slot;      /* awg slot number */
  int tp;        /* Flag to indicate whether test point has been set up */
  int comp;      /* Flag to indicate whether awgStream component is set up */
  int blocksize; /* Block size in "epoch" units (1/16 second time intervals)
		    e.g. a block 0.25 seconds long has a blocksize of 4.
		    Allowed values are 1 through 16. */
  int nblocks;  /* Total number of blocks buffered and/or sent so far */
  int curgps;   /* GPS time (integer number of seconds) of current buffer
		   (or next buffer to be created) */
  int curepoch; /* Epoch counter of current/next buffer */

  int sentgps;    /* GPS time (integer seconds) and epoch of last buffer */
  int sentepoch;  /*  sent to the front end */

  int nbufs;    /* Current number of buffers resident in memory */
  SIStrBuf *curbuf;   /* Pointer to current buffer.  A buffer is not created
			 until there is some data to put into it; therefore, it
			 is possible for a stream to have NO current buffer
			 (in which case curbuf==NULL) at various times. */
  SIStrBuf *firstbuf; /* Pointer to first buffer in linked list */
  SIStrBuf *lastbuf;  /* Pointer to last buffer in linked list.  This is the
			 current buffer if there is one; if not, it is the last
			 buffer which was filled. */
  long long lastsend; /* Time that last data was sent to front end, in
			 GPS nanoseconds */
  long long minwait;  /* Enforced minimum on the time interval between RPC
			 calls to transfer data to the front end, in
		         nanoseconds */
  int aborted;        /* Flag to indicate if stream has been aborted */
} SIStream;

/*------ Function prototypes -----------------------------------------------*/

void SIStrAppInfo( char *info );
int SIStrOpen( SIStream *sis, char *channel, int samprate, double starttime );
int SIStrAppend( SIStream *sis, float newdata[], int ndata, float scale );
int SIStrBlank( SIStream *sis, double duration );
int SIStrFlush( SIStream *sis );
int SIStrClose( SIStream *sis );
int SIStrAbort( SIStream *sis );
char * SIStrErrorMsg( int status );

#ifdef  __cplusplus
}
#endif

#endif /* _SISTR_H */
