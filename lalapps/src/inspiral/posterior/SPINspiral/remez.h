/* 
   
   SPINspiral:                parameter estimation on binary inspirals detected by LIGO, including spins of the binary members
   include/remez.h:           3rd-party routines
   
   
   Copyright 1995, 1998  Jake Janovetz (janovetz@uiuc.edu)
   
   
   This file is part of SPINspiral.
   
   SPINspiral is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   SPINspiral is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with SPINspiral.  If not, see <http://www.gnu.org/licenses/>.
   
*/


//Copied from remez.c, extracted from http://www.janovetz.com/jake/remez/remez-19980711.zip

/**************************************************************************
 * Parks-McClellan algorithm for FIR filter design (C version)
 *************************************************************************/


#ifndef __REMEZ_H__
#define __REMEZ_H__

#define BANDPASS       1
#define DIFFERENTIATOR 2
#define HILBERT        3

#define NEGATIVE       0
#define POSITIVE       1

#define Pi             3.1415926535897932
#define Pi2            6.2831853071795865

#define GRIDDENSITY    16
#define MAXITERATIONS  40

/* Function prototype for remez() - the only function that should need be
 * called from external code
 */
void remez(double h[], int numtaps,
           int numband, double bands[], double des[], double weight[],
           int type);


void CreateDenseGrid(int r, int numtaps, int numband, double bands[],
                     double des[], double weight[], int *gridsize,
                     double Grid[], double D[], double W[],
                     int symmetry);
void InitialGuess(int r, int Ext[], int gridsize);
void CalcParms(int r, int Ext[], double Grid[], double D[], double W[],
	       double ad[], double x[], double y[]);
double ComputeA(double freq, int r, double ad[], double x[], double y[]);
void CalcError(int r, double ad[], double x[], double y[],
               int gridsize, double Grid[],
               double D[], double W[], double E[]);
void Search(int r, int Ext[],
            int gridsize, double E[]);
void FreqSample(int N, double A[], double h[], int symm);
short isDone(int r, int Ext[], double E[]);


#endif /* __REMEZ_H__ */

