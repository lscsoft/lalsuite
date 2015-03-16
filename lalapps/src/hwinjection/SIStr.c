/*==========================================================================
SIStr.c - Client API library for streaming a waveform to the awg front end
           (The name 'SIStr' comes from Signal Injection Stream.)
Written 28 Sep - 4 Oct 2001 by Peter Shawhan

This file provides a simple interface to let a client program send a waveform
to the LIGO GDS "arbitrary waveform generator" (awg), which is then injected
onto one of the GDS excitation channels.  The basic skeleton of a client
program is as follows:


#include <stdlib.h>
#include <stdio.h>
#include "SIStr.h"
  SIStream sis;

  SIStrAppInfo( "myclient burst_waveform_2" );

  status = SIStrOpen( &sis, channel, samprate, starttime );
  if ( status != SIStr_OK ) {
    fprintf( stderr, "Error opening stream: %s\n", SIStrErrorMsg(status) );
    return 2;
  }

  while ( there is waveform data to send ) {
    status = SIStrAppend( &sis, data_array, ndata, scale );
    if ( status != SIStr_OK ) {
      fprintf( stderr, "Error streaming data: %s\n", SIStrErrorMsg(status) );
      break;
    }
  }

  status = SIStrClose( &sis );
  if ( status != SIStr_OK ) {
    fprintf( stderr, "Error closing stream: %s\n", SIStrErrorMsg(status) );
    return 2;
  }


Note how we check for an error after each call to a SIStr function, but if
SIStrAppend returns an error, then we do NOT exit out of the program;
we just break out of the loop so we can still call SIStrClose to clean up.

In the call to SIStrOpen, you must pass the name of a valid excitation channel,
and a sampling rate that matches the actual sampling rate for that channel.
You can pass an explicit start time (in GPS seconds, as a double-precision
value), which must be between the current time, and the current time plus
one day.  Or you can pass starttime=0.0, in which case the waveform injection
will start "immediately" (actually, after a delay of several seconds).

In the call to SIStrAppend, data_array is a pointer to an array of single-
precision floating-point values, and ndata is the number of values.  Feel free
to append a single value at a time if you want (i.e. ndata=1), but note that
if you have the value in a scalar variable, then data_array must be a POINTER
to that variable; you can't pass the value directly.  Waveform data passed to
SIStrAppend is assembled into fixed-size buffers which, once filled, are sent
to the front-end awg server at appropriate times using the remote procedure
call (RPC) mechanism.  All the buffering and timing is taken care of
internally by the SIStr library functions, so you do not have to worry about
this as long as you can calculate the waveform and pass it to SIStrAppend at
a rate faster than real-time.

There are a few other functions which you might use in certain circumstances:

  SIStrBlank( &sis, duration )
    Appends zeros for specified number of seconds, up to one day.
    (duration is a double.)

  SIStrFlush( &sis )
    Fills rest of current buffer with zeros, sends all local buffers to front
    end, and sleeps until after the last part of the waveform has actually been
    injected by the front end.  SIStrClose calls this function, so normally
    you do not have to call it explicitly.

  SIStrAbort( &sis )
    Clears all local buffers (but cannot clear buffers which have already been
    sent to the front end).  Also marks the stream as "aborted" so that any
    attempt to append more waveform data to it will fail.  You must still
    be sure to call SIStrClose on this stream to clean up properly.


To compile the object file (SIStr.o):
  cc -c SIStr.c -I$GDS_DIR/src/util -I$GDS_DIR/src/awg -I$GDS_DIR/src/rmem
where GDS_DIR on red      = /opt/CDS/d/gds/diag
              on london   = /opt/LLO/c/gds/diag

To compile and link a client program called 'mystreamer.c' :
  cc mystreamer.c SIStr.o -I$GDS_DIR/src/util -I$GDS_DIR/src/awg \
        -L$GDS_DIR/lib -lawg -ltestpoint -lrt -o mystreamer

To run the client program: $GDS_DIR/lib must be in your LD_LIBRARY_PATH
so that you can load the libawg.so and libtestpoint.so shared objects.

==========================================================================*/

/* The following #define affects whether an "extern" is used in SIStr.h */
#define _SISTR_LIBRARY

/* The following is a debugging tool to fake the behavior of the front end */
#if 0
#define DISABLEAWG
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

/* GDS includes */
#include <tconv.h>   /* For TAInow() */
#include <dtt/awgapi.h>
#include <dtt/testpoint.h>

#include "SIStr.h"
/*====== Prototypes for static local functions =============================*/
static int SIStrInit( SIStream *sis );
static int SIStrCreateBuf( SIStream *sis );
static int SIStrSend( SIStream *sis, int flushflag );
static int SIStrCleanup( SIStream *sis );
static int SIStrLog( SIStream *sis, char *condition, int value );


/*====== Global variables for internal use =================================*/
int SIStr_counter = 0;
char SIStr_appInfo[256] = "";


/*==========================================================================*/
void SIStrAppInfo( char *info )
{
  /* Copy the information string to an internal array */
  strncpy( SIStr_appInfo, info, sizeof(SIStr_appInfo)-1 );
  SIStr_appInfo[sizeof(SIStr_appInfo)-1] = '\0';

  return;
}


/*==========================================================================*/
int SIStrOpen( SIStream *sis, char *channel, int samprate, double starttime )
     /* Initialize a SIStream structure.  Returns SIStr_OK if successful,
	or some other value if there is an error */
{
  int status;
  tainsec_t now_ns;
  double frac;
  int firstsamp, isamp;
  int len;
  char chanlist[SIStr_MAXCHANLISTSIZE];
  char *channelpadded;
  char *cptr;
  int truerate;
  char* errMsg;
  char cmdString[64];

#ifdef DISABLEAWG
  fprintf( stderr, "NOTE: AWG INTERFACE IS DISABLED\n" );
#endif

  if (SIStr_debug) { printf( "In SIStrOpen\n" ); }

  /* Initialize the SIStream structure */
  status = SIStrInit( sis );
  if (SIStr_debug) { printf( "  SIStrInit returned %d\n", status ); }
  if ( status != SIStr_OK ) { return status; }

  /* Make sure the other arguments have valid values */
  if ( channel == NULL ) { return SIStr_EBADARG; }
  if ( strlen(channel) == 0 ) { return SIStr_EBADARG; }
  if ( strlen(channel) > SIStr_MAXCHANNAMELENGTH-1 ) { return SIStr_EBADARG; }
  if ( samprate <= 0 || samprate > 65536 ) { return SIStr_EBADRATE; }

  now_ns = TAInow();

  /* If 0.0 was passed as the start time, pick time several seconds from now */
  if ( starttime == 0.0 ) {
    starttime = (double) ( (now_ns + SIStr_LEADTIME) / 1000000000LL ) + 4.0;
    if (SIStr_debug) {
      printf( "  Current time is %d.%09d\n",
	      (int) (now_ns/1000000000LL), (int) (now_ns%1000000000LL) );
      printf( "  Assigning start time = %.9f\n", starttime );
    }
  } else {
    /* Start time was explicitly specified; make sure it is reasonable */
    if ( starttime < 600000000.0 ) { return SIStr_EBADSTART; }
    if ( starttime > 1800000000.0 ) { return SIStr_EBADSTART; }

    /* Make sure the start time is at least 2 seconds in the future */
    if ( starttime < (double) now_ns * 0.000000001 + 2.0 ) {
      double temp_sistr = (double)SIStr_LEADTIME/1e9;
      fprintf(stderr,"SIStr_LEADTIME = %d\n",temp_sistr);
      /* cethrane
      fprintf(stderr,"starttime = %d \n",starttime);
      double temp_ns = now_ns/1e9;
      fprintf(stderr,"now_ns = %d \n",temp_ns);
      */
      return SIStr_EBADSTART;
    }

    /* Make sure the start time is not more than about a day in the future */
    if ( starttime > (double) now_ns * 0.000000001 + 25.0*3600.0 ) {
      return SIStr_EBADSTART;
    }

  }

  /* Store channel name, sampling rate, and start time */
  strcpy( sis->channel, channel );
  sis->samprate = samprate;
  sis->starttime = starttime;

  /* Calculate the starting block info (integer GPS time and epoch number) */
  sis->curgps = (int) starttime;
  frac = starttime - (double) sis->curgps;
  sis->curepoch = (int) (16.0 * frac);
  /* Do a sanity check */
  if ( sis->curepoch < 0 || sis->curepoch > 15 ) { return SIStr_EINTERNAL; }

  /* Figure out what precise sample (counting from 0) corresponds to the
     start time.  This might be in the middle of an epoch */
  firstsamp = (int) ( (double)samprate * (frac - sis->curepoch/16.0) + 0.5 );
  /* Do a sanity check */
  if ( firstsamp < 0 || firstsamp >= samprate ) { return SIStr_EINTERNAL; }

  if (SIStr_debug) {
    printf( "  Waveform starts at GPS=%d, epoch=%d, sample=%d\n",
	    sis->curgps, sis->curepoch, firstsamp );
  }

  /* If the waveform starts in the middle of an epoch, create the first buffer
     and pad it with zeros up to the point where the waveform should start */
  if ( firstsamp > 0 ) {
    status = SIStrCreateBuf( sis );
    if (SIStr_debug) { printf( "  SIStrCreateBuf returned %d\n", status ); }
    if ( status != SIStr_OK ) { return status; }
    for ( isamp = 0; isamp < firstsamp; isamp++ ) {
      sis->curbuf->data[isamp] = 0.0;
    }
    sis->curbuf->ndata = firstsamp;
  }

  /*-----------------------------------*/
  /* Get the list of valid channel names from the front end */
  /* Put a space at the beginning of the list */
  chanlist[0] = ' ';
#ifdef DISABLEAWG
  strcpy( chanlist+1, "H1:DUMMY 256 H1:DUMMY2 2048 H1:DUMMY3 16384" );
  len = strlen( "H1:DUMMY 256" );
#else
  len = awgGetChannelNames( chanlist+1, SIStr_MAXCHANLISTSIZE-2, 1 );
#endif
  if (SIStr_debug) { printf( "  awgGetChannelNames returned %d\n", len ); }

  /* Check whether we got the list successfully */
  if ( len < 0 ) { return SIStr_ELISTERR; }
  if ( len == 0 ) { return SIStr_ELISTNONE; }

  /* Check whether list was truncated */
  if ( len >= SIStr_MAXCHANLISTSIZE-2 ) { return SIStr_ELISTSIZE; }

  if (SIStr_debug>=2) { printf( "  Channel list: %s\n", chanlist ); }

#ifndef DISABLEAWG
  /* Make sure the specified channel is in the list, and look up its
     sampling rate.  Surround the channel name with spaces to ensure
     an exact match */
  len = strlen(channel);
  channelpadded = (char *) malloc( len+3 );
  if ( channelpadded == NULL ) { return SIStr_EMALLOC; }
  strcpy( channelpadded, " " );
  strcat( channelpadded, channel );
  strcat( channelpadded, " " );
  cptr = strstr( chanlist, channelpadded );
  free( channelpadded );
  if ( cptr == NULL ) { return SIStr_EBADCHAN; }

  truerate = atol( cptr+len+2 );
  if (SIStr_debug) {
    printf( "  Channel %s has true sampling rate %d\n", channel, truerate );
  }

  printf ("truerate = %d, samprate = %d\n", truerate, samprate );
  /* Make sure the specified sampling rate matches the true sampling rate */
  if ( samprate != truerate ) { return SIStr_EDIFFRATE; }
#endif

  /*-----------------------------------*/
  /* Set up an awg slot for this channel */
#ifdef DISABLEAWG
  sis->slot = 1;
#else
  sis->slot = awgSetChannel( channel );
#endif
  if (SIStr_debug) { printf( "  awgSetChannel returned %d\n", sis->slot ); }
  if ( sis->slot < 0 ) {
    fprintf( stderr, "Error code from awgSetChannel: %d\n", sis->slot );
    return SIStr_ESETSLOT;
  }

#ifdef DISABLEAWG
  status = 0;
#else
  /*hack*/
  /*
  printf( "About to call tpCommand show\n" );
  cptr = tpCommand( "show" );
  if ( cptr == NULL ) {
    printf( "tpCommand show returned NULL\n" );
  } else {
    printf( "tpCommand show :\n%s", cptr );
  }
  */

  status = tpRequestName( channel, -1, 0, 0 );
#endif
  if (SIStr_debug) { printf( "  tpRequestName returned %d\n", status ); }
  if ( status < 0 ) {
    fprintf( stderr, "Error code from tpRequestName: %d\n", status );
    status = SIStrCleanup( sis );
    if (SIStr_debug) { printf( "  SIStrCleanup returned %d\n", status ); }
    return SIStr_ESETTP;
  }
  sis->tp = 1;

  sprintf( cmdString, "set %d stream 1.0", sis->slot );
#ifdef DISABLEAWG
  errMsg = "Dummy command";
#else
  errMsg = awgCommand( cmdString );
#endif
  if ( SIStr_debug ) {
    if ( errMsg != NULL ) {
      printf( "  awgCommand returned %s\n", errMsg );
    } else {
      printf( "  awgCommand returned NULL (i.e. no error message)\n", errMsg );
    }
  }
  if( strncmp(errMsg, "error:", 6) == 0 ) {
    fprintf( stderr, "Error message from awgCommand: %s\n", errMsg );
    status = SIStrCleanup( sis );
    if (SIStr_debug) { printf( "  SIStrCleanup returned %d\n", status ); }
    return SIStr_ESETCOMP;
  }
  sis->comp = 1;

#ifndef DISABLEAWG
  /* Append a line to the central log */
  /* If this fails, cancel the injection! We don't want to allow injections
   without a believable log file! */
  status = SIStrLog( sis, "start", 0 );
  if ( status != 0 ) {
    fprintf( stderr, "Aborting because SIStrLog call failed\n" );
  }
#endif

  return SIStr_OK;
}


/*==========================================================================*/
static int SIStrInit( SIStream *sis )
     /* Initialize a SIStream structure.  Returns SIStr_OK if successful,
	or some other value if there is an error */
{
  if (SIStr_debug) { printf( "  In SIStrInit\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }

  /* Set the "magic number" element, to indicate that this structure
    has been properly initialized */
  sis->magic = SIStr_MAGICVAL;

  /* Increment the counter to get an internal identifier */
  SIStr_counter++;
  sis->id = SIStr_counter;

  /* Initialize some other things */
  memset( sis->channel, 0, SIStr_MAXCHANNAMELENGTH );
  sis->samprate = 0;
  sis->starttime = 0.0;
  sis->slot = 0;
  sis->tp = 0;
  sis->comp = 0;
  sis->blocksize = SIStr_BLOCKSIZE;   /* Hard-coded for now */
  sis->nblocks = 0;
  sis->curgps = 0;
  sis->curepoch = 0;
  sis->sentgps = 0;
  sis->sentepoch = 0;
  sis->lastsend = 0LL;
  sis->minwait = (long long) ( sis->blocksize * (1000000000/32) );
  sis->aborted = 0;

  /* Clear the buffer pointers */
  sis->nbufs = 0;
  sis->curbuf = NULL;
  sis->firstbuf = NULL;
  sis->lastbuf = NULL;

  return SIStr_OK;
}


/*==========================================================================*/
int SIStrAppend( SIStream *sis, float newdata[], int ndata, float scale )
     /* Append some waveform data to the stream; when we have a complete
	block, try to send it */
{
  SIStrBuf *curbuf;
  int status;
  int idata;

  if (SIStr_debug>=2) { printf( "In SIStrAppend\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Make sure this stream has not been aborted */
  if ( sis->aborted ) { return SIStr_EABORTED; }

  /* Make sure other arguments are valid */
  if ( newdata == NULL ) { return SIStr_EBADARG; }
  if ( ndata < 0 ) { return SIStr_EBADARG; }

  /* If there is nothing to do, just return */
  if ( ndata == 0 ) { return SIStr_OK; }

  /* If there is no current buffer, create one */
  if ( sis->curbuf == NULL ) {
    status = SIStrCreateBuf( sis );
    if (SIStr_debug) { printf( "  SIStrCreateBuf returned %d\n", status ); }
    if ( status != SIStr_OK ) { return status; }
  }

  /* Get pointer to the current buffer into local variable,
    for convenience */
  curbuf = sis->curbuf;

  /*--- Now loop over input data points and add them to buffers ---*/
  for ( idata=0; idata<ndata; idata++ ) {

    /* Copy the data value, applying the scale factor */
    curbuf->data[curbuf->ndata] = scale * newdata[idata];
    curbuf->ndata++;

    /*-- If this fills the current buffer, close it out and try to send it --*/
    if ( curbuf->ndata == curbuf->size ) {

      if (SIStr_debug) { printf( "  Time to close out this buffer\n" ); }

      /* The stream now has no current buffer */
      sis->curbuf = NULL;

      /* Increment the "current buffer time" for this stream */
      sis->curepoch += sis->blocksize;
      while ( sis->curepoch > 15 ) {
	sis->curgps += 1;
	sis->curepoch -= 16;
      }

      /* Try to send the buffered data to the front end */
      status = SIStrSend( sis, 0 );
      if (SIStr_debug) { printf( "  SIStrSend returned %d\n", status ); }
      if ( status != SIStr_OK ) {
	/* Abort this stream */
	SIStrAbort( sis );
	return status;
      }

      /* If there is more data to be added, create a new buffer */
      if ( idata < ndata-1 ) {
	if (SIStr_debug) { printf( "  Need to create a new buffer\n" ); }
	status = SIStrCreateBuf( sis );
	if (SIStr_debug) {
	  printf( "  SIStrCreateBuf returned %d\n", status );
	}
	if ( status != SIStr_OK ) { return status; }
      }

      /* Update the local variable which points to the current buffer */
      curbuf = sis->curbuf;

    }

  }

  return SIStr_OK;
}


/*==========================================================================*/
static int SIStrCreateBuf( SIStream *sis )
     /* Allocates memory for a new buffer and sets it up as the current
	buffer */
{
  int size;
  SIStrBuf *curbuf;

  if (SIStr_debug) { printf( "  In SIStrCreateBuf\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Make sure there is not a current buffer */
  if ( sis->curbuf != NULL ) { return SIStr_EOTHER; }

  /* Calculate the size of the data array, and make sure it doesn't exceed
     the maximum */
  size = sis->samprate * sis->blocksize / 16;
  if (SIStr_debug) { printf( "    Calculated buffer size is %d\n", size ); }
  if ( size > SIStr_MAXBUFSIZE ) { return SIStr_EBUFSIZE; }

  /* Allocate memory for the buffer object */
  curbuf = (SIStrBuf *) malloc( sizeof(SIStrBuf) );
  if ( curbuf == NULL ) { return SIStr_EMALLOC; }

  /* Allocate memory for the data array */
  curbuf->data = (float *) calloc( size, sizeof(float) );
  if ( curbuf->data == NULL ) {
    free( curbuf );
    return SIStr_EMALLOC;
  }

  /* Update the SIStream structure to know about this buffer */
  sis->nblocks++;
  sis->nbufs++;
  sis->curbuf = curbuf;
  if ( sis->lastbuf != NULL ) {
    sis->lastbuf->next = curbuf;
  } else {
    sis->firstbuf = curbuf;
  }
  sis->lastbuf = curbuf;

  if (SIStr_debug) {
    printf( "    SIStream now has nblocks=%d, nbufs=%d\n",
	    sis->nblocks, sis->nbufs );
  }

  /* Initialize the new buffer */
  curbuf->gpstime = sis->curgps;
  curbuf->epoch = sis->curepoch;
  curbuf->iblock = sis->nblocks;
  curbuf->size = size;
  curbuf->ndata = 0;
  curbuf->next = NULL;

  if (SIStr_debug) {
    printf( "    New buffer has GPS=%d, epoch=%d\n",
	    curbuf->gpstime, curbuf->epoch );
  }

  return SIStr_OK;
}


/*==========================================================================*/
static int SIStrSend( SIStream *sis, int flushflag )
     /* Tries to send the buffered data to the front end.  Normally, this
	function will only send data which is within the "lead time" limit;
	any other data will remain buffered, and this function will return.
	However, if flushflag is nonzero, then this function will be sure to
	send all buffered data, sleeping as necessary to send it at an
	appropriate rate and with the appropriate lead time. */
{
  tainsec_t now_ns, buftime_ns, delta_ns;
  int status;
  int ntries;
  SIStrBuf *buf;
  int i;
  struct timespec wait;  /* On Sun, see /usr/include/sys/time_impl.h */

  if (SIStr_debug) { printf( "  In SIStrSend at %lld\n", TAInow() ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Make sure this stream has not been aborted */
  if ( sis->aborted ) { return SIStr_EABORTED; }

  /*--- If flushflag is zero, we will send AT MOST one buffer.  If flushflag
    is nonzero, we will loop over all buffers and try to send them all. ---*/

  do {

    /* Get pointer to first buffer into a local variable */
    buf = sis->firstbuf;

    /* Make sure this exists and is not the current buffer */
    if ( buf == NULL || buf == sis->curbuf ) { break; }

    /* Calculate the start time of the buffer in GPS NANOseconds */
    buftime_ns = (tainsec_t) buf->gpstime * 1000000000LL
      + (tainsec_t) buf->epoch * 62500000LL;

    /* Look up the current time */
    now_ns = TAInow();

    /* If the time of this buffer is too far in the future, just return --
       UNLESS we have reached the maximum allowed number of local buffers
       (and are not flushing), in which case sleep until it is appropriate
       to send this buffer */
    if ( buftime_ns > (now_ns + SIStr_LEADTIME) ) {
      if ( sis->nbufs < SIStr_MAXBUFS && flushflag == 0 ) {
	if (SIStr_debug) {
	  printf( "    Buffer time is beyond target lead time; returning\n" );
	}
	break;
      }

      /* Figure out how long to sleep for */
      delta_ns = buftime_ns - now_ns - SIStr_LEADTIME;
      wait.tv_sec = (time_t) ( delta_ns / 1000000000LL );
      wait.tv_nsec = (long) ( delta_ns % 1000000000LL );

      /* Sleep, then refresh the current time */
      if (SIStr_debug) {
	printf( "    Sleeping %d.%09d seconds to reach target lead time\n",
	       wait.tv_sec, wait.tv_nsec );
      }
      nanosleep( &wait, NULL );
      now_ns = TAInow();
      if (SIStr_debug) { printf( "    Time is now %lld\n", now_ns ); }
    }

    /* Impose a minimum time between data transfers.  This ensures that we
       do not try to send data at more than twice real-time speed.  If the
       time elapsed since the last transfer is less than this minimum, just
       return, UNLESS we have reached the maximum allowed number of local
       buffers (and are not flushing), in which case sleep until the minimum
       time has elapsed */
    if ( now_ns < (tainsec_t) (sis->lastsend + sis->minwait) ) {
      if ( sis->nbufs < SIStr_MAXBUFS && flushflag == 0 ) {
	if (SIStr_debug) {
	  printf( "    Previous buffer was sent very recently; returning\n" );
	}
	break;
      }

      /* Figure out how long to sleep for */
      delta_ns = (tainsec_t) (sis->lastsend + sis->minwait) - now_ns;
      wait.tv_sec = (time_t) ( delta_ns / 1000000000LL );
      wait.tv_nsec = (long) ( delta_ns % 1000000000LL );

      /* Sleep, then refresh the current time */
      if (SIStr_debug) {
	printf( "    Sleeping %d.%09d seconds to honor minimum wait time\n",
	       wait.tv_sec, wait.tv_nsec );
      }
      nanosleep( &wait, NULL );
      now_ns = TAInow();
      if (SIStr_debug) { printf( "    Time is now %lld\n", now_ns ); }
    }

    /*- At this point, we believe that it is appropriate to send the buffer -*/

    /* Call the routine in awgapi.c which sends a block block via RPC */
#ifdef DISABLEAWG
    if (SIStr_debug) {
      printf( "    Dummy send (%d,%d):", buf->gpstime, buf->epoch );
      for ( i=0; i<buf->ndata; i++ ) { printf( " %f", buf->data[i] ); }
      printf( "\n" );
    }
    status = 0;
#else
    if (SIStr_debug) {
      printf( "    Send (%10d,%2d)     at %lld\n",
	      buf->gpstime, buf->epoch, TAInow() );
      printf( "    Ndata=%d; data0=%f; data1=%f\n",
	      buf->ndata, buf->data[0], buf->data[1]);
    }
    status = awgSendWaveform( sis->slot,
			      (taisec_t) buf->gpstime, buf->epoch,
			      buf->data, buf->ndata );
    if (SIStr_debug) {
      printf( "    awgSendWaveform returned at %lld\n", TAInow() );
    }
#endif
    if (SIStr_debug) {
      printf( "    awgSendWaveform returned %d\n", status );
    }

    /* Record the timestamp for the data in this buffer, plus the time at
       which it was sent to the front end (whether successful or not) */
    sis->sentgps = buf->gpstime;
    sis->sentepoch = buf->epoch;
    sis->lastsend = (long long) now_ns;

    /* If there was a hard error, return with that error code */
    if ( status < 0 ) { return status; }

    /* If we get here, then we know that the data was accepted, so update the
       pointers in the SIStream structure to forget about this buffer, and
       then delete the buffer itself */
    sis->nbufs--;
    sis->firstbuf = buf->next;
    if ( sis->firstbuf == NULL ) { sis->lastbuf = NULL; }
    if ( buf->data != NULL ) { free( buf->data ); }
    free( buf );

    if (SIStr_debug) {
      printf( "    SIStream now has nblocks=%d, nbufs=%d\n",
	      sis->nblocks, sis->nbufs );
    }

  } while ( flushflag != 0 && sis->firstbuf != NULL );

  return SIStr_OK;
}


/*==========================================================================*/
int SIStrBlank( SIStream *sis, double duration )
     /* Append zeros to a stream for the specified amount of time
	(expressed in double-precision GPS seconds) */
{
  int status;
  float *zeros;
  int nsamp, nsend;

  if (SIStr_debug) { printf( "In SIStrBlank with duration=%f\n", duration ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Make sure this stream has not been aborted */
  if ( sis->aborted ) { return SIStr_EABORTED; }

  /* Check the duration argument for special cases */
  if ( duration < 0.0 ) { return SIStr_EBADARG; }
  if ( duration == 0.0 ) { return SIStr_OK; }

  /* Disallow durations of longer than about a day */
  if ( duration > 25.0*3600.0 ) { return SIStr_EBADARG; }

  /* Figure out how many samples this duration is.  Assuming a duration of
     no more than a day, and a sampling rate of no more than 16384, this
     always fits into a 32-bit signed integer. */
  nsamp = (int) ( duration * sis->samprate + 0.5 );
  if (SIStr_debug) {
    printf( "  Appending %d zeros to stream\n", nsamp );
  }

  /* Create an array of zeros, so that we can send them in blocks for
     reasonable efficiency.  Note that calloc is guaranteed to initialize
     its contents to zero. */
  zeros = (float *) calloc( 128, sizeof(float) );
  if ( zeros == NULL ) { return SIStr_EMALLOC; }

  /* Now call SIStrAppend repeatedly to send the correct number of zeros */
  while ( nsamp > 0 ) {
    nsend = ( nsamp>128 ? 128 : nsamp );
    status = SIStrAppend( sis, zeros, nsend, 1.0 );
    if (SIStr_debug>=2) { printf( "    SIStrAppend returned %d\n", status ); }
    if ( status != SIStr_OK ) {
      free( zeros );
      return status;
    }
    nsamp -= nsend;
  }

  free( zeros );
  return SIStr_OK;
}


/*==========================================================================*/
int SIStrFlush( SIStream *sis )
     /* Flushes out locally buffered data to front end, then sleeps until
	the front end buffers are expected to have drained as well */
{
  int status;
  int retval = SIStr_OK;
  float *zeros;
  int nsamp;
  tainsec_t now_ns, draintime_ns, delta_ns;
  struct timespec wait;  /* On Sun, see /usr/include/sys/time_impl.h */

  if (SIStr_debug) { printf( "--In SIStrFlush\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* If there is a partially-full current buffer, fill the rest with zeros */
  if ( sis->curbuf != NULL ) {

    /* Figure out how many zeros are needed to fill up this buffer */
    nsamp = sis->curbuf->size - sis->curbuf->ndata;
    if (SIStr_debug) {
      printf( "    Appending %d zeros to fill last buffer\n", nsamp );
    }

    /* Create an array of zeros of the appropriate size.  Note that calloc is
       guaranteed to initialize its contents to zero. */
    zeros = (float *) calloc( nsamp, sizeof(float) );
    if ( zeros == NULL ) { return SIStr_EMALLOC; }

    /* Send the array of zeros */
    status = SIStrAppend( sis, zeros, nsamp, 1.0 );
    if (SIStr_debug>=2) { printf( "    SIStrAppend returned %d\n", status ); }

    /* Free the array and check the return status */
    free( zeros );
    if ( status < 0 ) { return status; }

    /* Double-check that there is now no current buffer */
    if ( sis->curbuf != NULL ) { return SIStr_EINTERNAL; }
  }

  /* Flush the local buffers */
  status = SIStrSend( sis, 1 );
  if (SIStr_debug) { printf( "    SIStrSend/flush returned %d\n", status ); }
  /* If we got an error, set the return value but don't return yet, since
     we still need to sleep to make sure the front end drains.  In any case,
     if the return code just says that the stream was aborted, don't consider
     that an error. */
  if ( status != SIStr_OK && status != SIStr_EABORTED ) {
    SIStrAbort( sis );
    retval = status;
  }

  /* Sleep until the front-end buffers are sure to have drained
     (one second after ENDING time of last buffer sent to the front end) */
  now_ns = TAInow();
  if (SIStr_debug) { printf( "    Time is now %lld\n", now_ns ); }
  draintime_ns = (tainsec_t) sis->sentgps * 1000000000LL
    + (tainsec_t) (sis->sentepoch + sis->blocksize) * 62500000LL
    + 1000000000LL;
  delta_ns = draintime_ns - now_ns;

  if ( delta_ns > (tainsec_t)0 ) {
    wait.tv_sec = (time_t) ( delta_ns / 1000000000LL );
    wait.tv_nsec = (long) ( delta_ns % 1000000000LL );
    if (SIStr_debug) {
      printf( "    Sleeping %d.%09d seconds to let front-end buffers drain\n",
	     wait.tv_sec, wait.tv_nsec );
    }
    nanosleep( &wait, NULL );
  }

  return retval;
}


/*==========================================================================*/
int SIStrClose( SIStream *sis )
     /* Close a stream */
{
  int status;
  int retval = SIStr_OK;

  if (SIStr_debug) { printf( "In SIStrClose\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Flush all buffered data in the stream */
  status = SIStrFlush( sis );
  if (SIStr_debug) { printf( "  SIStrFlush returned %d\n", status ); }
  /* If an error occurs, set the return value but DO NOT return, since
     it is essential that we call SIStrCleanup to free the awg slot, etc. */
  if ( status != SIStr_OK ) { retval = status; }

  /* Clean up */
  status = SIStrCleanup( sis );
  if (SIStr_debug) { printf( "  SIStrCleanup returned %d\n", status ); }
  if ( status != SIStr_OK ) { return status; }

  return retval;
}


/*==========================================================================*/
static int SIStrCleanup( SIStream *sis )
     /* Cleanly disconnect from the awg.  Always try to do all three
	operations, even if there is an error in one or more of them.
	An error from awgClearWaveforms is not really important, so
        allow it to be overwritten by any other error; but keep track
	of whether errors occur in both of the other operations. */
{
  int status;
  int retval = SIStr_OK;

  if (SIStr_debug) { printf( "  In SIStrCleanup\n" ); }

  if ( sis->comp ) {
#ifdef DISABLEAWG
    status = 0;
#else
    status = awgClearWaveforms( sis->slot );
#endif
    if (SIStr_debug) {
      printf( "    awgClearWaveforms returned %d\n", status );
    }
    if ( status < 0 ) {
      fprintf( stderr, "Error code from awgClearWaveforms: %d\n", status );
      retval = SIStr_ECLRCOMP;
    }
    sis->comp = 0;
  }

  if ( sis->tp ) {
#ifdef DISABLEAWG
    status = 0;
#else
    status = tpClearName( sis->channel );
#endif
    if (SIStr_debug) { printf( "    tpClearName returned %d\n", status ); }
    if ( status < 0 ) {
      fprintf( stderr, "Error code from tpClearName: %d\n", status );
      retval = SIStr_ECLRTP;
    }
    sis->tp = 0;
  }

  if ( sis->slot > 0 ) {
#ifdef DISABLEAWG
    status = 0;
#else
    status = awgRemoveChannel( sis->slot );
#endif
    if (SIStr_debug) {
      printf( "    awgRemoveChannel returned %d\n", status );
    }
    if ( status < 0 ) {
      fprintf( stderr, "Error code from awgRemoveChannel: %d\n", status );
      if ( retval == SIStr_ECLRTP ) {
	retval = SIStr_ECLRBOTH;
      } else {
	retval = SIStr_ECLRSLOT;
      }
    }
    sis->slot = 0;
  }

#ifndef DISABLEAWG
  /* If an error occurred, log it */
  if ( retval != SIStr_OK ) SIStrLog( sis, "error", retval );

  /* Append a line to the central log */
  SIStrLog( sis, "stop", 0 );
#endif

  return retval;
}


/*==========================================================================*/
int SIStrAbort( SIStream *sis )
     /* Deletes all data that was buffered to be sent to the front end */
{
  SIStrBuf *buf;

  if (SIStr_debug) { printf( "In SIStrAbort\n" ); }

  /* Make sure a valid SIStream pointer was passed */
  if ( sis == NULL ) { return SIStr_EBADARG; }
  if ( sis->magic != SIStr_MAGICVAL ) { return SIStr_EUNINIT; }

  /* Clear all buffers */
  for ( buf = sis->firstbuf; buf != NULL; buf = buf->next ) {
    if ( buf->data != NULL ) { free( buf->data ); }
    free( buf );
  }

  sis->nbufs = 0;
  sis->curbuf = NULL;
  sis->firstbuf = NULL;
  sis->lastbuf = NULL;

  if (SIStr_debug) {
    printf( "    SIStream now has nblocks=%d, nbufs=%d\n",
	    sis->nblocks, sis->nbufs );
  }

  /* Mark this stream as aborted, so that any attempt to send more data on
     it will fail */
  sis->aborted = 1;

  return SIStr_OK;
}


/*==========================================================================*/
static int SIStrLog( SIStream *sis, char *condition, int value )
{
  char *chptr;
  char cwd[256];
  char infoline[512];
  char hostname[64];
  int status;
  FILE *file;
  tainsec_t now_ns;

  /*------ Beginning of code ------*/
  if ( sis == NULL ) { return -1; }
  if ( sis->magic != SIStr_MAGICVAL ) { return -1; }

  /*------ Construct the message ------*/

  chptr = getenv( "HOST" );
  if ( chptr == NULL ) { return -1; }
  strcpy( hostname, chptr );

  if ( strcmp(condition,"start") == 0 ) {
    sprintf( infoline, "%17.6f %s start %s %d %d %s", sis->starttime,
	     sis->channel, hostname, getpid(), sis->id, SIStr_appInfo );
  } else if ( strcmp(condition,"error") == 0 ) {
    sprintf( infoline, "%17.6f %s error %s %d %d %d",
	     (double) sis->curgps + (double) sis->curepoch / 16.0,
	     sis->channel, hostname, getpid(), sis->id, value );
  } else {
    now_ns = TAInow();
    sprintf( infoline, "%10d.%06d %s %-5s %s %d %d",
	     (int) (now_ns/1000000000LL), (int) (now_ns%1000000000LL) / 1000,
             sis->channel, condition, hostname, getpid(), sis->id );
  }

  /*------ Write out the message ------*/

  /*-- Look up the current working directory --*/
  chptr = getcwd( cwd, sizeof(cwd) );

  /*-- Change to the directory where the central log file is kept --*/
#ifdef BLIND
  status = chdir( "/tmp/signal_injection_log" );
#else
  status = chdir( "/tmp/signal_injection_log" );
#endif
  if ( status != 0 ) {
    /*-- Try Livingston --*/
#ifdef BLIND
    status = chdir( "/tmp/signal_injection_log" );
#else
    status = chdir( "/tmp/signal_injection_log" );
#endif
    if ( status != 0 ) {
      /*-- Failure! --*/
      if (chptr) chdir( cwd );
      return -1;
    }
  }

  /*-- Open the file for appending --*/
  file = fopen( "injection_log.txt", "a" );
  if (chptr) chdir( cwd );
  if ( file == NULL ) { return -1; }

  /*-- Write the line of information --*/
  fprintf( file, "%s\n", infoline );

  /*-- Close the file and return --*/
  fclose( file );
  return 0;
}


/*==========================================================================*/
char * SIStrErrorMsg( int status )
     /* Returns a pointer to a char array containing the error message for
	the given status code */
{
  switch ( status ) {
  case SIStr_OK:
    return "OK";
  case SIStr_WFULL:
    return "Warning: Data accepted but front-end buffer is now full";
  case SIStr_WDUP:
    return "Warning: Data accepted but time was duplicated";
  case SIStr_WGAP:
    return "Warning: Data accepted but block was not contiguous";
  case SIStr_EBADSLOT:
    return "Specified awg slot is not set up to accept data stream";
  case SIStr_EBADDATA:
    return "Invalid data block or size";
  case SIStr_EPAST:
    return "Block time is already past";
  case SIStr_EFUTURE:
    return "Block time is too far in the future";
  case SIStr_ECONN:
    return "RPC connection failed";
  case SIStr_EBADARG:
    return "Bad function argument";
  case SIStr_EBADSTR:
    return "Stream is not correctly open for writing";
  case SIStr_EBADRATE:
    return "Invalid sampling rate";
  case SIStr_EUNINIT:
    return "Stream was not properly initialized";
  case SIStr_EMALLOC:
    return "Error allocating memory";
  case SIStr_EOTHER:
    return "Other error";
  case SIStr_EABORTED:
    return "Stream was aborted";
  case SIStr_EBADSTART:
    return "Unreasonable start time";
  case SIStr_EINTERNAL:
    return "Unexpected internal error";
  case SIStr_EBUFSIZE:
    return "Tried to create too large a data buffer";
  case SIStr_ETIMEOUT:
    return "Timeout while trying to send data to front end";
  case SIStr_ELISTERR:
    return "Error retrieving channel list from front end";
  case SIStr_ELISTNONE:
    return "Front end returned an empty channel list";
  case SIStr_ELISTSIZE:
    return "Front-end channel list is too large for local variable";
  case SIStr_EBADCHAN:
    return "Not a valid excitation channel";
  case SIStr_EDIFFRATE:
    return "Specified sampling rate differs from actual rate";
  case SIStr_ESETSLOT:
    return "Error setting up an awg slot for the channel";
  case SIStr_ESETTP:
    return "Error setting up test point for the channel";
  case SIStr_ESETCOMP:
    return "Error setting up awgStream component";
  case SIStr_ECLRCOMP:
    return "Error clearing awgStream component";
  case SIStr_ECLRTP:
    return "Error clearing test point";
  case SIStr_ECLRSLOT:
    return "Error freeing awg slot";
  case SIStr_ECLRBOTH:
    return "Error clearing test point AND freeing awg slot";
  }
}
