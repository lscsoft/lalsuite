/*
*  Copyright (C) 2007 Reinhard Prix
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* lisatools/io/readxml.c --- Copyright (c) 2006 Michele Vallisneri

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "ezxml.h"

#include "readxml.h"


char *stripcopy(const char *orig) {
    int pos = 0, len = 0, i;
    char *ret;
    
    /* Strip space-like characters at the beginning */

    while(orig[pos] == ' ' || orig[pos] == '\n' || orig[pos] == '\r')
        pos++;

    /* Walk until the end of the string */
        
    while(orig[pos+len] != 0)
        len++;

    /* Strip space-like characters at the end */

    len--;
    while(orig[pos+len] == ' ' || orig[pos+len] == '\n' || orig[pos+len] == '\r') len--;
    len++;

    /* Copy the string */
    
    ret = malloc( (len+1) * sizeof(char) );
    assert(ret != 0);
    
    for(i=0;i<len;i++) ret[i] = orig[pos+i];
    ret[len] = 0;
    
    return ret;
}


char *splitcopy(const char *orig,int i) {
    int pos = 0, len = 0, comma = 0, j;
    char *ret;
    
    while(comma < i && orig[pos]) {
        if(!orig[pos]) return 0;

        if(orig[pos] == ',') comma++;

        pos++;
    }

    /* Strip space-like characters at the beginning */

    while(orig[pos] == ' ' || orig[pos] == '\n' || orig[pos] == '\r')
        pos++;

    /* Walk until the end of the string */
    
    while(orig[pos+len] != 0 && orig[pos+len] != ',')
        len++;

    /* Strip space-like characters at the end */

    len--;
    while(orig[pos+len] == ' ' || orig[pos+len] == '\n' || orig[pos+len] == '\r') len--;
    len++;

    /* Copy the string */

    ret = malloc( (len+1) * sizeof(char) );
    assert(ret != 0);
        
    for(j=0;j<len;j++) ret[j] = orig[pos+j];
    ret[len] = 0;
    
    return ret;
}


static char *basename(const char *xml, const char *bin) {
    int pos, i, lenbin;
    char *ret;
    
    pos = strlen(xml) - 1;
    while(pos >= 0 && xml[pos] != '/') pos--;

    if(pos < 0) return 0;
    pos++;
    
    lenbin = strlen(bin);
    
    ret = malloc( (lenbin + pos + 1) * sizeof(char) );
    assert(ret != 0);
        
    for(i=0;i<pos;i++) ret[i] = xml[i];
    for(i=0;i<lenbin;i++) ret[i+pos] = bin[i];
    ret[lenbin + pos + 1] = 0;

    return ret;
}


/* Endianness conversion: a bit rough, but should work */
/* From http://www.codeproject.com/cpp/endianness.asp */

#define BIGENDIAN      0
#define LITTLEENDIAN   1

int testbyteorder()
{
   short int word = 0x0001;
   char *byte = (char *) &word;
   return(byte[0] ? LITTLEENDIAN : BIGENDIAN);
}

/* This will work as long as doubles are 8 bytes */
static void convertendianness(double *val) {
    double outval;
    unsigned char *cval = (unsigned char*)&outval;

    cval[0] = ((unsigned char *)val)[7];
    cval[1] = ((unsigned char *)val)[6];
    cval[2] = ((unsigned char *)val)[5];
    cval[3] = ((unsigned char *)val)[4];
    cval[4] = ((unsigned char *)val)[3];
    cval[5] = ((unsigned char *)val)[2];
    cval[6] = ((unsigned char *)val)[1];
    cval[7] = ((unsigned char *)val)[0];

    *val = outval;
}

static TimeSeries *dotimeseries(ezxml_t series,char *xmlname) {
    ezxml_t param, array, dim, stream;
    const char *name=NULL, *timeoffset=NULL, *cadence=NULL, *length=NULL, *records=NULL, *encoding=NULL, *binaryfile=NULL;

    TimeSeries *timeseries;
    double *buffer;
    
    char *pathbinfile;
    FILE *binfile;

    int i,j;

    name = ezxml_attr(series,"Name");
    
    for(param = ezxml_child(series,"Param"); param; param = param->next) {
        if(!strcmp(ezxml_attr(param,"Name"),"TimeOffset")) {
            timeoffset = ezxml_txt(param);
            assert(timeoffset);
            }

        if(!strcmp(ezxml_attr(param,"Name"),"Cadence")) {
            cadence = ezxml_txt(param);
            assert(cadence);
            }
    }
    
    array = ezxml_child(series,"Array");
    assert(array);
    /* Expected Array here */
    
    for(dim = ezxml_child(array,"Dim"); dim; dim = dim->next) {
        if(!strcmp(ezxml_attr(dim,"Name"),"Length")) {
            length = ezxml_txt(dim);
            assert(length);
        }
            
        if(!strcmp(ezxml_attr(dim,"Name"),"Records")) {
            records = ezxml_txt(dim);
            assert(records);
        }
    }
    
    stream = ezxml_child(array,"Stream");
    assert(stream);
    
    ezxml_attr(stream,"Type");
    encoding = ezxml_attr(stream,"Encoding");

    binaryfile = ezxml_txt(stream);
    assert(binaryfile);
    
    /* Now create the TimeSeries structure */
    
    timeseries = malloc(sizeof(TimeSeries));
    assert(timeseries != 0);
    
    timeseries->Name = stripcopy(name);
    timeseries->FileName = stripcopy(binaryfile);
  
    if(timeoffset) {
        timeseries->TimeOffset = strtod(timeoffset,NULL);
    } else {
        timeseries->TimeOffset = 0.0;
    }
    
    timeseries->Cadence = strtod(cadence,NULL);

    timeseries->Length = strtol(length,NULL,10);
    timeseries->Records = strtol(records,NULL,10);
        
    /* fprintf(stderr,"Allocating %ld bytes for read buffer...\n",timeseries->Length * timeseries->Records * sizeof(double)); */
    buffer = malloc(timeseries->Length * timeseries->Records * sizeof(double));
    assert(buffer != 0);

    pathbinfile = basename(xmlname,timeseries->FileName);

    binfile = fopen(pathbinfile,"r");

    if(binfile == 0) {             
      /* fprintf(stderr,"...can't find %s, trying in the working directory...\n",pathbinfile); */
        binfile = fopen(timeseries->FileName,"r");
        assert(binfile != 0);
    }
    
    free(pathbinfile);
    
    size_t len = timeseries->Length * timeseries->Records;
    if ( len != fread(buffer,sizeof(double),len,binfile) )
      fprintf (stderr, "\nFailed to read %zd entries from 'binfile'\n", len );

    /* Do the encoding switch if necessary */

    if( (strstr(encoding,"LittleEndian") && testbyteorder() == BIGENDIAN) ||
        (strstr(encoding,"BigEndian")    && testbyteorder() == LITTLEENDIAN) ) {
        fprintf(stderr,"Converting endianness of binary data!\n");
        
        for(i=0;i<timeseries->Length * timeseries->Records;i++)
            convertendianness(&buffer[i]);
    }
        
    fclose(binfile);

    timeseries->Data = malloc(timeseries->Records * sizeof(DataColumn *));
    assert(timeseries->Data != 0);
    
    for(i=0;i<timeseries->Records;i++) {
        timeseries->Data[i] = malloc(timeseries->Records * sizeof(DataColumn));
        assert(timeseries->Data[i] != 0);
        
        timeseries->Data[i]->Name = splitcopy(timeseries->Name,i);
        
        timeseries->Data[i]->TimeOffset = timeseries->TimeOffset;
        timeseries->Data[i]->Cadence = timeseries->Cadence;

        timeseries->Data[i]->Length = timeseries->Length;
        
        timeseries->Data[i]->data = malloc(timeseries->Length * sizeof(double));
        assert(timeseries->Data[i]->data != 0);

        for(j=0;j<timeseries->Length;j++) {
            timeseries->Data[i]->data[j] = buffer[(j*timeseries->Records) + i];
        }
    }
    
    free(buffer);
    
    return timeseries;
}


void freeTimeSeries(TimeSeries *timeseries) {
    int i;

    for(i=0;i<timeseries->Records;i++) {
        free(timeseries->Data[i]->data);
        free(timeseries->Data[i]->Name);
        
        free(timeseries->Data[i]);
    }
    
    free(timeseries->Data);
    
    free(timeseries->FileName);
    free(timeseries->Name);
    
    free(timeseries);
}


void freeLISASources(LISASource *lisasource) {
    LISASource *next, *current;

    next = lisasource;

    while(next) {
        current = next;
        next = current->NextSource;
        
        free(current->Name);
        free(current->TimeSeries);
    }
}


LISASource *getLISASources(char *filename) {
    ezxml_t tree, section, source, param, series;
    const char *type=NULL, *name=NULL, *elat=NULL, *elon=NULL, *pol=NULL;

    LISASource *first = 0, *current = 0;
    /* TimeSeries *timeseries; */ /* unused */

    tree = ezxml_parse_file(filename);
    assert(tree);

    for(section = ezxml_child(tree, "XSIL"); section; section = section->next) {
        type = ezxml_attr(section, "Type");
        
        if(!strcmp(type,"SourceData")) {
            for(source = ezxml_child(section,"XSIL"); source; source = source->next) { 
                assert(!strcmp(ezxml_attr(source,"Type"),"SampledPlaneWave"));
                /* Expected SampledPlaneWave here... */

                name = ezxml_attr(source,"Name");

                for(param = ezxml_child(source,"Param"); param; param = param->next) {
                    if(!strcmp(ezxml_attr(param,"Name"),"EclipticLatitude")) {
                        elat = ezxml_txt(param);
                        assert(elat);
                    }
            
                    if(!strcmp(ezxml_attr(param,"Name"),"EclipticLongitude")) {
                        elon = ezxml_txt(param);
                        assert(elon);
                    }
                        
                    if(!strcmp(ezxml_attr(param,"Name"),"Polarization")) {
                        pol = ezxml_txt(param);
                        assert(pol);
                    }
                }

                if(!first) {
                    first = malloc(sizeof(LISASource));
                    assert(first != 0);

                    current = first;
                } else {
                    current->NextSource = malloc(sizeof(LISASource));
                    assert(current->NextSource != 0);

                    current = current->NextSource;
                }

                current->Name = stripcopy(name);

                current->EclipticLatitude  = strtod(elat,NULL); 
                current->EclipticLongitude = strtod(elon,NULL);                 
                current->Polarization      = strtod(pol,NULL);
                        
                series = ezxml_child(source,"XSIL");
                assert(!strcmp(ezxml_attr(series,"Type"),"TimeSeries"));
                /* Expected TimeSeries here... */

                current->TimeSeries = dotimeseries(series,filename);

                current->NextSource = 0;
            }
        }
    }

    ezxml_free(tree);
    
    return first;
}


/* For the moment, support only the first TDIObservable/TimeSeries section per file */

TimeSeries *getTDIdata(char *filename) {
    ezxml_t tree, section, obs, series;
    const char *type;
    TimeSeries *timeseries = 0;

    tree = ezxml_parse_file(filename);
    if ( tree == NULL ) {
      fprintf (stderr, "\ngetTDIdata(): call to 'ezxml_parse_file() failed!\n\n");
      return NULL;
    }
    
    for(section = ezxml_child(tree, "XSIL"); section; section = section->next) {
        type = ezxml_attr(section, "Type");
        
        if(!strcmp(type,"TDIData")) {
            for(obs = ezxml_child(section,"XSIL"); obs; obs = obs->next) {
                assert(!strcmp(ezxml_attr(obs,"Type"),"TDIObservable"));
            
                series = ezxml_child(obs,"XSIL");
                if(!strcmp(ezxml_attr(series,"Type"),"TimeSeries")) {
                    timeseries = dotimeseries(series,filename);
                    return timeseries;
                }
            }
        }
    }

    ezxml_free(tree);

    return 0;
}

