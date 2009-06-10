/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Julien Sylvestre
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

#include <lal/LPC.h>
#include <lal/LALStatusMacros.h>

NRCSID (LPCTESTC,"$Id$");


#define CHKST if(status.statusCode != 0) {REPORTSTATUS (&status); return -1;}

INT4 lalDebugLevel = 2;

int main() {

  static LALStatus status;
  REAL4Vector x;
  REAL4Vector a;
  REAL4 data[10] = {1.0,0.9,1.1,0.5,1.2,0.8,0.9,1.4,1.0,0.9};
  REAL4 matlab[4] = {1.0,-0.6861568,-0.3781715,0.218191};
  REAL4 eps = 1e-3;
  UINT4 p = 3, i;

  UINT4 m=5;
  REAL4 ptest1[5] = {1.0,1.0,1.0,1.0,1.0};
  REAL4 ptest2[5] = {1.0,2.0,3.0,4.0,5.0};
  REAL4 ptest3[5] = {0.1,2.0,1.0,2.0,2.0};
  REAL4 ptest1m[5] = {1.0,1.0,1.0,1.0,1.0};
  REAL4 ptest2m[5] = {1.0000,0.8000,0.6000,0.4000,0.2000};
  REAL4 ptest3m[5] = {0.1000,0.0585,0.0595,0.0676,0.0033};
  REAL4Vector ptest;

  /* first test polystab */
  ptest.length = m;
  ptest.data = (REAL4 *)ptest1;

  LALPolystab(&status, &ptest);
  CHKST

  for(i=0;i<m;i++) {
    printf("%g\t%g\n",ptest1[i],ptest1m[i]);
  }

  printf("*********************************************\n");
  ptest.data = (REAL4 *)ptest2;

  LALPolystab(&status, &ptest);
  CHKST

  for(i=0;i<m;i++) {
    printf("%g\t%g\n",ptest2[i],ptest2m[i]);
  }

  printf("*********************************************\n");

  ptest.data = (REAL4 *)ptest3;

  LALPolystab(&status, &ptest);
  CHKST

  for(i=0;i<m;i++) {
    printf("%g\t%g\n",ptest3[i],ptest3m[i]);
  }

  printf("*********************************************\n");

  x.length = 10;
  x.data = (REAL4 *)(&data);

  a.length = 0;
  a.data = NULL;

  LALLPC(&status, &a, &x, p);
  CHKST

  for(i=0;i<a.length;i++) {
    printf("%u\t%g\t%g\n",i,a.data[i],matlab[i]);

    if(fabs(matlab[i]-a.data[i])/fabs(matlab[i]) > eps) {
      fprintf(stderr,"Results don't match matlab results.\n");
      return -1;
    }
  }

  LALFree(a.data);

  return 0;

}
