/* lisatools/io/readxml.h --- Copyright (c) 2006 Michele Vallisneri

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE. */

/* --- data structures --- */

/* Represents one column of double-precision "data" for observable
   "Name", with the specified "TimeOffset", "Length", and "Cadence" */

typedef struct {
    char *Name;
    
    double TimeOffset;
    double Cadence;
    
    long Length;
    
    double *data;
} DataColumn;

/* Represents a time series with "Records" columns (each represented
   by a DataColumn struct, pointed to by the pointers in "Data"),
   with the specified "Name" (usually a comma-separated list of the 
   names of the columns), "TimeOffset", "Cadence", "Length". */

typedef struct {
    char *Name;
    char *FileName;
    
    double TimeOffset;
    double Cadence;
    
    long Length;
    long Records;
    
    DataColumn **Data;
} TimeSeries;

/* Represents a linked list of LISA sampled sources of given position
   in the sky ("EclipticLatitude", "EclipticLongitude", "Polarization"),
   encoded by the hp and hc data contained in the linked TimeSeries;
   NextSource points to the next source in the list, or is zero for
   the last source in the list. */
   
typedef struct LISASource {
    char *Name;
    
    double EclipticLatitude; 
    double EclipticLongitude;                 
    double Polarization;   
    
    TimeSeries *TimeSeries;
    
    struct LISASource *NextSource;
} LISASource;

/* --- data-handling functions --- */

/* Parses XML file "filename" and returns a TimeSeries struct
   corresponding to the FIRST TDIData/TDIObservable element; within
   the structure, the array of pointers DataColumn points to the
   individual struct DataColumn */

TimeSeries *getTDIdata(char *filename);

/* Parses XML file "filename" and returns a linked list of 
   LISASource structs, each of which represents a SampledPlaneWave
   object, and stores hp and hc data in a TimeSeries object */ 

LISASource *getLISASources(char *filename);

/* Frees all memory taken up by the parsed LISA sampled sources */

void freeLISASources(LISASource *lisasource);

/* Frees all memory taken up by the parsed TimeSeries (e.g., TDIData) */

void freeTimeSeries(TimeSeries *timeseries);

/* --- Utility functions --- */

/* Makes a copy of the string "orig", stripping all whitespace-type
   characters; the returned string is a private copy that must be
   deallocated with free(). */

char *stripcopy(const char *orig);

/* Returns the i'th string (i=0,...) in a comma-separated list, stripping
   all whitespace-type characters (' ','\n','\r'); the returned string
   is a private copy that must be deallocated with free(). */

char *splitcopy(const char *orig,int i);

/* Check endianness of platform */

#define BIGENDIAN      0
#define LITTLEENDIAN   1

int testbyteorder(void);
