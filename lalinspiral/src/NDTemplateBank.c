/*
*  Copyright (C) 2007 Chad Hanna, Jolien Creighton, Benjamin Owen, Reinhard Prix
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

/**
 * \author Hanna, C. R.
 * \file
 * \ingroup TemplateBankGeneration_h
 *
 * \brief This module handles template bank generation for up to searches with
 * \f$<=\f$ 12 dimensional parameter spaces.
 *
 * ### Prototypes ###
 *
 * <tt>LALNDTemplateBank()</tt>
 *
 * ### Description ###
 *
 * This module tiles up to a 12 dimensional space when given a metric and a
 * function that determines the search region.
 *
 * ### Algorithm ###
 *
 * The algorithm first draws a rectilinear box in the primed coordinates
 * which includes the distorted box, then steps through along the directions
 * of the primed coordinates.  At each point it tests if the point lies
 * within the distorted box. If the point is inside the distorted box, the
 * algorithm adds a template to the linked list. If not, it continues.
 *
 * ### Uses ###
 *
 * \code
 * LALCalloc()
 * LALFree()
 * \endcode
 *
 * ### Notes ###
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/AVFactories.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALMalloc.h>
#include <lal/MatrixUtils.h>
#include <lal/TemplateBankGeneration.h>

static REAL4 DotProduct(REAL4 *EV, REAL4 *DX){
  INT2 loop = 0;
  REAL4 dot = 0.0;

  for (loop = 0; loop < 12; loop++){
    dot  += EV[loop] * DX[loop];
    }
  return dot;
  }



void

LALNDTemplateBank( LALStatus *status,
                   NDTemplateBankInput *input,
                   NDTemplateBankFunctionPtrs *functionPtrs,
		   NDTemplateBankOutput **output)

{

  INT2 Dimension = input->dimension;
  REAL4Array  *metric =           NULL; /* parameter-space metric 	    */
  REAL4Array  *inverse = 	  NULL;
  REAL4 det = 0;
  UINT4Vector *metricDimensions = NULL; /* contains the dimension of metric */
  REAL4Vector *eigenval =         NULL; /* eigenvalues of metric  	    */
  NDTemplateBankOutput *bank = 	  NULL; /* the template bank 	 	    */
  NDTemplateBankOutput *first =   NULL;
  INT2 dimLoop = 0;
  INT2 loop = 0;
  INT2 testFlag = 0;
  REAL4 dxLoop[12] = 		{0};
  REAL4 coordinatedxs[12] = 	{1,1,1,1,1,1,1,1,1,1,1,1};
  REAL4 EV[12][12] =    	{{0}};
  REAL4 EVinv[12][12] = 	{{0}};
  REAL4 minX[12] = 		{0};
  REAL4 maxX[12] = 		{0};
  REAL4 temp = 0;
  INT4 cnt = 0;
  /* INT4 bccFlag = 0; */

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /* Initialize unused portions of input arrays */
  for(dimLoop = Dimension; dimLoop < 12; dimLoop++){
    input->minCoordinates[dimLoop] = 0.0;
    input->maxCoordinates[dimLoop] = 0.0;
    input->minParameters[dimLoop] = 0.0;
    input->maxParameters[dimLoop] = 0.0;
    }


  /* Set up the metric stuff */
  LALSCreateVector( status->statusPtr, &eigenval, (UINT4) Dimension);
  LALU4CreateVector( status->statusPtr, &metricDimensions, (UINT4) 2 );
  metricDimensions->data[1] = metricDimensions->data[0] = Dimension;
  LALSCreateArray( status->statusPtr, &metric, metricDimensions );
  LALSCreateArray( status->statusPtr, &inverse, metricDimensions );


  /* Call Metric Function */
  functionPtrs->metric(status->statusPtr, input, metric);

  /* Begin all the diagonalization business */
  LALSSymmetricEigenVectors( status->statusPtr, eigenval, metric );
  /*  Extract eigenvectors and determine displacements*/
  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
    for (loop = 0; loop < Dimension; loop++){
      EV[dimLoop][loop] = metric->data[loop*Dimension + dimLoop];
      }
    }

/*
  if (Dimension == 3){
    coordinatedxs[0] = 1.3333333 * sqrt(2.0*input->mm/(eigenval->data[0]));
    coordinatedxs[1] = 1.3333333 * sqrt(2.0*input->mm/(eigenval->data[1]));
    coordinatedxs[2] = 0.6666667 * sqrt(2.0*input->mm/(eigenval->data[2]));
    }
*/

  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
    coordinatedxs[dimLoop] = sqrt(2.0*input->mm/((REAL4) Dimension * eigenval->data[dimLoop]));
    }

  printf("\nCoordinatedxs = %f,%f,%f\n", coordinatedxs[0], coordinatedxs[1], coordinatedxs[2]);

  /* invert the transformation matrix */
  TRY(LALSMatrixInverse(status->statusPtr, &det, metric, inverse), status);


  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
    for (loop = 0; loop < Dimension; loop++){
      EVinv[dimLoop][loop] = inverse->data[loop*Dimension + dimLoop];
      }
    }



  printf("\nmin and max coordinates before transformation\n");
  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
     printf("MinX[%i] = %f\tMaxX[%i] = %f\n", dimLoop, input->minCoordinates[dimLoop], dimLoop, input->maxCoordinates[dimLoop]);
     }

  printf("\nmin and max coordinates after transformation\n");

  /* transform initial coordinates */
  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
    minX[dimLoop] = DotProduct(EVinv[dimLoop], input->minCoordinates);
    maxX[dimLoop] = DotProduct(EVinv[dimLoop], input->maxCoordinates);
    printf("MinX[%i] = %f\tMaxX[%i] = %f\n", dimLoop, minX[dimLoop], dimLoop, maxX[dimLoop]);
    }


  /*switch around the bounds according to sign */
  for (dimLoop = 0; dimLoop < Dimension; dimLoop++){
    if (minX[dimLoop] > maxX[dimLoop]){
      temp = minX[dimLoop];
      minX[dimLoop] = maxX[dimLoop];
      maxX[dimLoop] = temp;
      }
    }


  /* allocate first template bank node */
  first = bank = (NDTemplateBankOutput *) LALCalloc(1, sizeof(NDTemplateBankOutput));
  for(dxLoop[11] =  minX[11];
      dxLoop[11] <= maxX[11];
      dxLoop[11] += coordinatedxs[11]){
   for(dxLoop[10] =  minX[10];
       dxLoop[10] <= maxX[10];
       dxLoop[10] += coordinatedxs[10]){
    for(dxLoop[9] =  minX[9];
        dxLoop[9] <= maxX[9];
        dxLoop[9] += coordinatedxs[9]){
     for(dxLoop[8] =  minX[8];
         dxLoop[8] <= maxX[8];
         dxLoop[8] += coordinatedxs[8]){
      for(dxLoop[7] =  minX[7];
          dxLoop[7] <= maxX[7];
          dxLoop[7] += coordinatedxs[7]){
       for(dxLoop[6] =  minX[6];
           dxLoop[6] <= maxX[6];
           dxLoop[6] += coordinatedxs[6]){
        for(dxLoop[5] =  minX[5];
            dxLoop[5] <= maxX[5];
            dxLoop[5] += coordinatedxs[5]){
         for(dxLoop[4] =  minX[4];
             dxLoop[4] <= maxX[4];
             dxLoop[4] += coordinatedxs[4]){
          for(dxLoop[3] =  minX[3];
              dxLoop[3] <= maxX[3];
              dxLoop[3] += coordinatedxs[3]){
           for(dxLoop[2] =  minX[2];
               dxLoop[2] <= maxX[2];
               dxLoop[2] += coordinatedxs[2]){
            for(dxLoop[1] =  minX[1];
                dxLoop[1] <= maxX[1];
                dxLoop[1] += coordinatedxs[1]){
             for(dxLoop[0] =  minX[0];
                 dxLoop[0] <= maxX[0];
                 dxLoop[0] += coordinatedxs[0]){
               /* test spot */
               for(dimLoop=0; dimLoop < Dimension; dimLoop++){
                 bank->coordinateVals[dimLoop] = DotProduct(EV[dimLoop], dxLoop);
               }
               functionPtrs->test(status->statusPtr, input, bank, &testFlag);
               if (testFlag){
                 bank = bank->next = (NDTemplateBankOutput *) LALCalloc(1, sizeof(NDTemplateBankOutput));
               }

               /* go dx behind */
               for(dimLoop=0; dimLoop < Dimension; dimLoop++){
                 bank->coordinateVals[dimLoop] = DotProduct(EV[dimLoop], dxLoop) -  DotProduct(EV[dimLoop], coordinatedxs);
                 /* set the other values to loop values */
                 for(loop=0; loop < dimLoop; loop++){
 		   bank->coordinateVals[loop] = DotProduct(EV[loop], dxLoop);
                 }
                 /* test the next spot */
                 functionPtrs->test(status->statusPtr, input, bank, &testFlag);
                 /* if it fails test the halfway point */
                 if (!testFlag){
                   bank->coordinateVals[dimLoop] = DotProduct(EV[dimLoop], dxLoop) - 0.5 * DotProduct(EV[dimLoop], coordinatedxs);
		   functionPtrs->test(status->statusPtr, input, bank, &testFlag);
                   /* if the halfway point passes add it */
                   if (testFlag){
                     cnt++;
                     bank = bank->next = (NDTemplateBankOutput *) LALCalloc(1, sizeof(NDTemplateBankOutput));
                   }
                 }
               }
               /* go dx ahead */
               for(dimLoop=0; dimLoop < Dimension; dimLoop++){
                 bank->coordinateVals[dimLoop] = DotProduct(EV[dimLoop], dxLoop) + DotProduct(EV[dimLoop], coordinatedxs);
                 /* set the other values to loop values */
                 for(loop=0; loop < dimLoop; loop++){
                   bank->coordinateVals[loop] = DotProduct(EV[loop], dxLoop);
                 }
                 /* test the next spot */
                 functionPtrs->test(status->statusPtr, input, bank, &testFlag);
                 /* if it fails test the halfway point */
                 if (!testFlag){
                   bank->coordinateVals[dimLoop] = DotProduct(EV[dimLoop], dxLoop) + 0.5 * DotProduct(EV[dimLoop], coordinatedxs);
                   functionPtrs->test(status->statusPtr, input, bank, &testFlag);
                   /* if the halfway point passes add it */
                   if (testFlag){
                     cnt++;
                     bank = bank->next = (NDTemplateBankOutput *) LALCalloc(1, sizeof(NDTemplateBankOutput));
                   }
                 }
               }

             }
            }
 	   }
          }
         }
        }
       }
      }
     }
    }
   }
  }


  printf("\nExtra templates added = %i\n", cnt);
  /* setting the output to the first template bank tile */
  bank->next = NULL;
  *output = first;

  /* if(testFlag) LALFree(bank) */;
  DETATCHSTATUSPTR(status);
  RETURN(status);

}
