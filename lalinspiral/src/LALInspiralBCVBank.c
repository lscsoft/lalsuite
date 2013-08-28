/*
*  Copyright (C) 2007 Thomas Cokelaer
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

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * \defgroup LALInspiralBCVBank_c Module LALInspiralBCVBank.c
 * \ingroup LALInspiralBank_h
 * \brief BCV template bank
 */
/*@{*/

/**
 * \brief Lay a flat grid of BCV templates in the user specified range
 * of the parameters \f$(\psi_0, \psi_3)\f$ in \c coarseIn structure
 * \author Cokelaer, T
 *
 * \heading{Description}
 * Given the range of the parameters \f$(\psi_0, \psi_3),\f$
 * the number of templates in the \c fCut direction,
 * \e minimalMatch, noise spectral density, upper and
 * lower frequency cutoffs (all in the input structure \c coarseIn)
 * this routine outputs the list of templates in the BCV bank
 * for the parameters \f$(\psi_0, \psi_3, f_{\mathrm{cut}}).\f$
 *
 * \heading{Algorithm}
 * A flat signal manifold is assumed and templates are laid
 * uniform in the three dimensions.  See below for an explanation
 * of how templates are chosen in the \c fcut direction.
 */
void
LALInspiralCreateBCVBank (
    LALStatus            *status,	/**< LAL status pointer */
    InspiralTemplateList **list,	/**< [out] an array containing the template bank parameters. */
    INT4                 *nlist,	/**< [out] the number of templates in bank */
    InspiralCoarseBankIn coarseIn	/**< UNDOCUMENTED */
    )

{
  INT4 	j = 0;
  INT4 	nlistOld = 0;
  /*do we really need static declaration here ? */
  static InspiralBankParams 	bankParams;
  static InspiralMetric 	metric;
  static InspiralTemplate 	params;
  static CreateVectorSequenceIn in;
  static REAL4VectorSequence 	*tempList = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( coarseIn.psi0Min > 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi0Max > coarseIn.psi0Min, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi3Min < 0., status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.psi3Max > coarseIn.psi3Min, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
  ASSERT( coarseIn.LowGM < coarseIn.HighGM, status,
      LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );

  /* populate the param structure so as to call
   * ComputeMetricBCV */
  params.fLower 	= coarseIn.fLower;
  params.fCutoff 	= coarseIn.fUpper;
  params.alpha 		= coarseIn.alpha;
  /* Get the BCV metric in psi0psi3 plane. */
  LALInspiralComputeMetricBCV( status->statusPtr,
      &metric, &coarseIn.shf, &params );
  CHECKSTATUSPTR( status );

  /* print the metric if lalinfo is on*/
  if ( lalDebugLevel & LALINFO )
  {
    REAL8 dx0 = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/metric.g00 );
    REAL8 dx1 = sqrt( 2.L * (1.L-coarseIn.mmCoarse)/metric.g11 );
    LALPrintError( "G00=%e G01=%e G11=%e\n",
        metric.G00, metric.G01, metric.G11 );
    LALPrintError( "g00=%e g11=%e theta=%e\n",
        metric.g00, metric.g11, metric.theta );
    LALPrintError( "dp0=%e dp1=%e\n", dx0, dx1 );
  }

  /* We have the metric, which is constant. Now we need to place
   * the templates in the parameter space which is define as follows by
   * the psi0 and psi3 range:
   * */
  bankParams.metric 		= &metric;
  bankParams.minimalMatch 	= coarseIn.mmCoarse;
  bankParams.x0Min 		= coarseIn.psi0Min;
  bankParams.x0Max 		= coarseIn.psi0Max;
  bankParams.x1Min 		= coarseIn.psi3Min;
  bankParams.x1Max 		= coarseIn.psi3Max;

  /* Let us define a temporary list of templates. */
  in.length 		= 1;
  in.vectorLength 	= 2;
  LALSCreateVectorSequence( status->statusPtr, &tempList, &in );
  CHECKSTATUSPTR( status );

  /* First we place templates in the psi0/psi3 plane.
   *
   * Historically we had two template banks. If gridSpacing is set to
   * S2BCV, then, the code generates the bank used during S2. This bank
   * uses a non-oriented square placement in psi0/psi3 plane and the
   * fcut dimension placement is done BCVRegularFcutBank or BCVFCutBank
   * function.
   * If gridSpacing is not S3BCV, then we have the choice between a
   * square or an hexagonal placement, oriented or not and fcut is placed
   * using BankFcutS3S4.
   * */
  if (coarseIn.gridSpacing  != S2BCV)
  {
    LALInspiralCreateFlatBankS3S4( status->statusPtr, tempList, &bankParams , coarseIn);
    CHECKSTATUSPTR( status );
  }
  else
  {
    LALInspiralCreateFlatBank( status->statusPtr, tempList, &bankParams);
    CHECKSTATUSPTR( status );
  }

  *nlist = tempList->length;
  *list = (InspiralTemplateList *)
    LALCalloc( *nlist, sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT (status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
  }
  /* We populate the output.
   * */
  for ( j = 0; j < *nlist; ++j )
  {
    (*list)[j].params.psi0 				= (REAL8) tempList->data[2*j];
    (*list)[j].params.psi3 				= (REAL8) tempList->data[2*j+1];
    (*list)[j].params.fLower 			= params.fLower;
    (*list)[j].params.nStartPad 		= 0;
    (*list)[j].params.nEndPad 			= 0;
    (*list)[j].params.tSampling			= coarseIn.tSampling;
    (*list)[j].params.distance 			=  1.;
    (*list)[j].params.signalAmplitude	= 1.;
    (*list)[j].params.approximant		= BCV;
    (*list)[j].params.massChoice		= psi0Andpsi3;
    (*list)[j].params.order				= LAL_PNORDER_TWO;
    (*list)[j].metric 					= metric;
    (*list)[j].params.alpha 			= coarseIn.alpha;
  }
  nlistOld = *nlist;

  /* Once the psi0/psi3 plane is populated, for each coordinate, we
   * populate along the fcut dimension. Again, for historical reason,
   * we call one of the following functions which slightly differs
   * from each other (see documentation).
   * */

  /* If coarseIn.lowGM == - 1 then LowGM is  unphysical. Hence, we use a
   * Regular grid in cutoff frequency which is independant of LowGM or
   * HighGM and which lays between Flower and Fsampling/2. If
   * coarseIn.lowGM != -1 then, we populate between two frequencyies
   * defines by low and high GM. (i.e lowGM = 6 means fLSO).
   *  */
  if (coarseIn.gridSpacing != S2BCV)
    {
      LALInspiralBCVBankFcutS3S4( status->statusPtr,
				list, nlist, coarseIn);
      CHECKSTATUSPTR( status );
    }
  else  if (coarseIn.LowGM  == -1)
  {
	  LALInspiralBCVRegularFcutBank( status->statusPtr,
	      list, nlist, coarseIn);
	  CHECKSTATUSPTR( status );
  }
  else
  {
	  LALInspiralBCVFcutBank( status->statusPtr,
	      list, nlist, coarseIn);
	  CHECKSTATUSPTR( status );
  }

  if ( lalDebugLevel & LALINFO )
  {
    LALPrintError(
        "LALInspiralBCVBank: template numbers before %d and after %d calling LALInspiralBCVBank\n",
        nlistOld, *nlist );
  }

  LALSDestroyVectorSequence( status->statusPtr, &tempList );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/**
 * The code expects <tt>list-\>vectorLength=2</tt> and allocates just the
 * requisite amount of memory to \c list and returns the number
 * of grid points in <tt>list-\>length</tt>.
 * The data points <tt>list-\>data[2j]</tt>, <tt>j=1,2, ... list-\>length</tt>,
 * contain the \f$x_0\f$-coordinates of the grid and data points <tt>list-\>data[2j+1]</tt>,
 * contain the \f$x_1\f$-coordinates of the grid.
 *
 * \heading{Description}
 * Given the \c metric and the \c minimalMatch this routine calls
 * <tt>bank/LALInspiralUpdateParams</tt> to get the spacings in user coordinates (which are
 * not necessarily the eigen-directions) and lays a uniform grid of templates in
 * the range specified in (<tt>bankParams-\>x0Min</tt>, <tt>bankParams-\>x0Max</tt>) and
 * (<tt>bankParams-\>x1Min</tt>, <tt>bankParams-\>x1Max</tt>).
 *
 * \heading{Algorithm}
 * The algorithm to lay templates is as follows: Given the increments \f$Dx_0\f$ and
 * \f$Dx_1\f$ found from calling <tt>bank/LALInspiralUpdateParams</tt> lay a rectangular
 * grid in the space of \f$(x_0, x_1).\f$
 *
 * <tt>
 * <ul>
 * <li> \f$x_1 = x_1^{\mathrm{min}}\f$
 * <li> do while (\f$x_1 \le x_1^{\mathrm{max}}\f$)<br>
 * <ul>
 * <li> \f$x_0 = x_0^{\mathrm{min}}\f$
 * <li> do while (\f$x_0 \le x_0^{\mathrm{max}}\f$)<br>
 * <ul>
 * <li> Add (\f$x_0, x_1\f$) to list
 * <li> numTemplates++
 * <li> Increment \f$x_0: \; x_0 = x_0 + Dx_0\f$
 * </ul>
 * <li> Increment \f$x_1: \; x_1 = x_1 + Dx_1\f$
 * </ul>
 * </ul>
 * </tt>
 */
void
LALInspiralCreateFlatBank (
    LALStatus            *status,	/**< LAL status pointer */
    REAL4VectorSequence  *list,		/**< [out] an array containing the template bank parameters */
    InspiralBankParams   *bankParams	/**< [in] It is necessary and sufficient to input
                                         * the eigenvalues of the metric and the angle between the \f$x_0\f$ axis and the
                                         * semi-major axis of the ambiguity ellipse, that is,
                                         * <tt>bankParams.metric.g00, bankParams.metric.g11, bankParams.metric.theta</tt>,
                                         * the minimal match, <tt>bankParams.minimalMatch</tt> and the range of the two
                                         * coordinates over which templates must be chosen:
                                         * (<tt>bankParams-\>x0Min</tt>, <tt>bankParams-\>x0Max</tt>) and
                                         * (<tt>bankParams-\>x1Min</tt>, <tt>bankParams-\>x1Max</tt>)
                                         */
    )
{
  InspiralMetric *metric;
  REAL8 minimalMatch;
  REAL8 x0, x1;
  UINT4 nlist = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* From the knowledge of the metric and the minimal match find the */
  /* constant increments bankParams->dx0 and bankParmams->dx1        */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;
  LALInspiralUpdateParams( status->statusPtr,
      bankParams, *metric, minimalMatch );
  CHECKSTATUSPTR( status );

  /* Construct the template bank */
  for (x1 = bankParams->x1Min; x1 <= bankParams->x1Max; x1 += bankParams->dx1)
  {
    for (x0 = bankParams->x0Min; x0 <= bankParams->x0Max; x0 += bankParams->dx0)
    {
      UINT4 ndx = 2 * nlist;
      list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
      if ( !list->data )
      {
        ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
      }
      list->data[ndx] = x0;
      list->data[ndx + 1] = x1;
      ++nlist;
    }
  }

  list->length = nlist;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/**
 * \brief Given a grid of templates with distinct values of \f$(\psi_0, \psi_3)\f$
 * this routine returns a new grid in which every template has \c numFcutTemplates
 * partners differing from one another in the ending frequency <tt>fendBCV</tt>.
 *
 * A call to this function should be preceeded by a call to LALInspiralCreateFlatBank()
 * or a similar function, that gives a grid in \f$(\psi_0, \psi_3)\f$ space.
 *
 * \heading{Description}
 *
 * A lattice of templates for BCV models should include,
 * in addition to the values of \f$(\psi_0, \psi_3)\f$
 * a range of \f$f_{\mathrm{cut}}\f$ -- the cutoff frequency.
 * The right approach would be
 * to compute the metric in the three-dimensional space of
 * \f$(\psi_0, \psi_3, f_{\mathrm{cut}})\f$ and to choose templates as
 * dictated by the metric. However, analytic computation of the
 * metric has not been easy. Therefore, it has become necessary
 * (at least for the time being) to make alternate choice of
 * the cutoff frequencies.
 *
 * In this routine we implement a simple
 * choice based on physical grounds: The post-Newtonian models
 * predict an ending frequency that is larger than, but close to,
 * the Schwarzschild last-stable orbit frequency
 * \f$f_{\mathrm{lso}} = (6^{3/2} \pi M )^{-1}\f$ where \f$M\f$ is the total mass,
 * while the effective one-body model has an ending frequency close
 * to the light-ring, whose Schwarzschild value is
 * \f$f_{\mathrm{lr}} = (3^{3/2} \pi M )^{-1}\f$. It is necessary to know
 * the total mass of the system in both cases.  However, not all
 * pairs of \f$(\psi_0, \psi_3)\f$ can be inverted to get a positive
 * \f$M\f$ but only when \f$\psi_0 > 0\f$ and \f$\psi_3 < 0\f$. Even then
 * it is not guaranteed that the symmetric mass ratio will be
 * less than \f$1/4,\f$ a necessary condition so that the component
 * masses are found to be real. However, we do not demand that the
 * symmetric mass ratio is less than a quarter. If the total mass
 * is non-negative then we compute the \f$(f_{\mathrm{lso}}, f_{\mathrm{lr}})\f$
 * and choose a user specified \c numFcutTemplates number of
 * templates with their cutoff frequency <tt>list-\>fFinal</tt> defined
 * uniformly spaced in the range \f$[f_{\mathrm{lso}},\ f_{\mathrm{lr}}]\f$.
 *
 * Furthermore, this routine discards all templates for which
 * either the mass is not defined or, when defined, <tt>list-\>fFinal</tt> is
 * smaller than the user defined lower frequency cutoff or larger
 * than the Nyquist frequency of templates.
 * Thus, the number of templates returned by this routine could
 * be larger or fewer than the input number of templates.
 *
 * \heading{Algorithm}
 * Given \f$(\psi_0, \psi_3)\f$ one can solve for \f$(M, \eta)\f$ using:
 * \f{equation}{
 * M = \frac{-\psi_3}{16 \pi^2 \psi_0},\ \ \eta = \frac{3}{128 \psi_0 (\pi M)^{5/3}}.
 * \f}
 * Given the total mass compute the last stable orbit and light-ring frequencies using
 * \f{equation}{
 * f_{\mathrm{lso}} = (6^{3/2} \pi M)^{-1},\ \  f_{\mathrm{lr}} = (3^{3/2} \pi M)^{-1}.
 * \f}
 * Divide the range \f$(f_{\mathrm{lso}}, f_{\mathrm{lr}})\f$ so as to have \f$n_{\mathrm{cut}}= \mathtt{numFcutTemplates}\f$
 * templates over this range:
 * \f{equation}{
 * df = f_{\mathrm{lr}} \frac {\left( 1 - 2^{-3/2} \right) }{ (n_{\mathrm{cut}} -1) }.
 * \f}
 * Next, choose templates at \f$f_k = f_\mathrm{lr} - k \times df,\f$ where \f$k=0, \ldots, n_\mathrm{cut}-1\f$.
 * Note that by definition \f$f_0 = f_\mathrm{lr}\f$ and \f$f_{n_\mathrm{cut}-1} = f_\mathrm{lso}\f$;
 * there are exatly \f$n_\mathrm{cut}\f$ templates in the range \f$(f_\mathrm{lso}, f_\mathrm{lr})\f$.
 * We discard a template if either \f$M\f$ is not defined or if \f$f_\mathrm{cut}\f$ is smaller
 * than the lower frequency cutoff specified in  <tt>list[j]-\>fLower</tt>.
 */
void
LALInspiralBCVFcutBank (
    LALStatus            *status,	/**< LAL status pointer */
    InspiralTemplateList **list,	/**< [out,in] an array initially containing the template
                                         * bank with the values of <tt>list[j]-\>psi0, list[j]-\>psi3, list[j]-\>fLower,</tt> specified,
                                         * is replaced on return with a re-sized array specifying also <tt>list-\>fFinal.</tt>
                                         */
    INT4                *NList,		/**< [out,in] the number of templates in the Input bank is replaced by the number of templates in the output bank. */
    InspiralCoarseBankIn coarseIn	/**< UNDOCUMENTED */
    )
{
  UINT4 nf; 	/* number of layers */
  UINT4 nlist; 	/* number of final templates */
  UINT4 j;
  UINT4 ndx;	/* temporary number of templates */
  REAL8 frac;	/* general variable*/
  REAL8 fendBCV;
  REAL4 LowGM;
  REAL4 HighGM;

  INITSTATUS(status);

  nf 		= coarseIn.numFcutTemplates;
  ndx 		= nlist = *NList;

  LowGM 	= coarseIn.LowGM;
  HighGM 	= coarseIn.HighGM;

  /* if we have only one layer, we don't need HighGM.
   * And default value for LowGM is  3GM*/
  if ( nf == 1 )
  {
    frac = 1;
  }
  else
  {
    frac = (1.L - 1.L/pow(HighGM/3., 1.5L)) / (nf-1.L);
  }

  /* for each psi0/psi3 pair, we generate the fcut layers */
  for ( j = 0; j < nlist; ++j )
  {
    UINT4 valid = 0;
    /* let us get the estimated params.fFinal*/
    LALPSItoMasses(status,  &((*list)[j].params), &valid , LowGM);

    if ( valid )
    {
      UINT4 i;
      REAL8 fMax;

      fMax = (*list)[j].params.fFinal;
      /* for each fcut layer */
      for ( i = 0; i < nf; ++i )
      {
		fendBCV = fMax * (1.L - (REAL8) i * frac);

        if ( (*list)[j].params.tSampling <= 0 )
        {
          ABORT( status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
        }
        /* the fFinal must be > flower and less than Nyquist.
         * if not, no template is generated. */
        if ( fendBCV > (*list)[j].params.fLower &&
            fendBCV < (*list)[j].params.tSampling / 2.0 )
        {
          ++ndx;

		  *list = (InspiralTemplateList *)
          LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
          if ( ! *list )
          {
            ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
          }
          memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
          (*list)[ndx-1] 				= (*list)[j];
          (*list)[ndx-1].params.fFinal 	= fendBCV;
          (*list)[ndx-1].metric 		= (*list)[0].metric;
          (*list)[ndx-1].nLayer 		= i;
        }
      }
    }
  }
  /**/
  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }
  /**/
  *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  RETURN( status );
}

/** UNDOCUMENTED */
void
LALPSItoMasses (
    LALStatus			*status,
    InspiralTemplate 	*params,
    UINT4            	*valid,
    REAL4             	HighGM
    )
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  if ( params->psi0 <= 0.L || params->psi3 >= 0.L )
  {
    *valid = 0;
  }
  else
  {
    REAL8 totalMass;
    REAL8 eta;
	/*we estimate the total mass and then fFinal*/
    params->totalMass = -params->psi3/(16.L*LAL_PI * LAL_PI * params->psi0);
    eta = params->eta =
      3.L/(128.L * params->psi0 * pow(LAL_PI*params->totalMass, (5./3.)));
    totalMass = params->totalMass;
    params->fFinal = 1.L/( LAL_PI * pow(HighGM, 1.5L) * params->totalMass );
    params->totalMass /= LAL_MTSUN_SI;
    *valid = 1;

#if 0
    if (params->eta > 0.25L)
    {
      *valid = 0;
    }
    else
    {
      LALInspiralParameterCalc( status->statusPtr, params );
      CHECKSTATUSPTR( status );
      *valid = 1;
    }
#endif

    params->t0 = 5.0 / ( 256.0*eta*pow(totalMass, (5./3.)) *
        pow(LAL_PI * params->fLower, (8./3.)));
    params->t3 = LAL_PI/(8.0*eta*pow(totalMass, (2./3.)) *
        pow(LAL_PI * params->fLower, (5./3.)));
  }
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/** UNDOCUMENTED */
void
LALInspiralBCVBankFcutS3S4 (
    LALStatus            	*status,
    InspiralTemplateList 	**list,
    INT4					*NList,
    InspiralCoarseBankIn 	coarseIn
    )

{
  UINT4 nf; 	/* number of layers */
  UINT4 nlist; 	/* number of final templates */
  UINT4 j;
  UINT4 ndx;	/* temporary number of templates */
  REAL8 frac;	/* general variable*/
  REAL8 fendBCV;
  REAL4 LowGM;
  REAL4 HighGM;

  INITSTATUS(status);

  nf    	= coarseIn.numFcutTemplates;
  ndx   	= nlist = *NList;
  LowGM		=  3.;
  HighGM    = coarseIn.HighGM;

  /* for each template, we get the fcut layers*/
  for ( j = 0; j < nlist; ++j )
  {
    UINT4 valid = 0;
    LALEmpiricalPSItoMassesConversion(status,
    	&((*list)[j].params), &valid , LowGM);

    if (valid)
    {
	  UINT4 i;
	  REAL8 fMax;

	  fMax = (*list)[j].params.fFinal;
	  if ( (*list)[j].params.tSampling <= 0 )
	  {
	    ABORT( status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE );
	  }
        /* the user might request only one layer */
	  if (nf == 1)
	  {
	  	  frac = 1;
	  }
	  else
	  {
	    frac = (1.L - 1.L/pow(HighGM/3., 1.5L)) / (nf-1.L);
	  }

      /* sometimes since fMin is greater than the nyquist frequency, there
       * is no template generated. This is not acceptable. We need at least
       * one frequency at the nyquist frequency otherwise low masses
       * systems are missed. */
      if (((fMax * (1.L - (REAL4) (nf-1) * frac)) >= (*list)[j].params.tSampling/2.0))
      {
        fMax = (*list)[j].params.tSampling/2.0 - 1. ;
        frac = -1;
      }

      /*Similarly, for high masses. */
      /*if (((fMax * (1.L - (REAL4) (nf-1) * frac)) <= (*list)[j].params.fLower * 1.5))
       {
          fMax = (*list)[j].params.fLower * 1.5 ;
        }
        */
      for (i=0; i<nf; i++)
      {
        fendBCV = fMax * (1.L - (REAL4) i * frac);
	    /* we search for valid expression of fendBCV and populate the bank */
	    if ( fendBCV >= (*list)[j].params.fLower * 1.5 &&
		  fendBCV < (*list)[j].params.tSampling / 2.0 )
	    {
		  ++ndx;
		  *list = (InspiralTemplateList *)
		  LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
		  if ( ! *list )
		  {
		    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
		  }
		  memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
		  (*list)[ndx-1] = (*list)[j];
		  (*list)[ndx-1].params.fFinal = fendBCV;
		  (*list)[ndx-1].metric = (*list)[0].metric;
		  (*list)[ndx-1].nLayer = i;
		}
	  }
	}
  }
  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }
  *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }


  RETURN( status );
}

/** UNDOCUMENTED */
void
LALInspiralBCVRegularFcutBank (
    LALStatus            	*status,
    InspiralTemplateList 	**list,
    INT4                	*NList,
    InspiralCoarseBankIn 	coarseIn
    )

{
  /* no restriction of  physical masses.
   * And regular layer of templates in the Frequency dimension */
  UINT4 i,nf, nlist, j, ndx;
  REAL8 fendBCV;

  INITSTATUS(status);

  nf = coarseIn.numFcutTemplates;
  ndx = nlist = *NList;

  for ( j = 0; j < nlist; ++j )
  {
    for ( i = 1; i <=nf; ++i )
    {
	  fendBCV = (*list)[j].params.fLower
		+ i * ((*list)[j].params.tSampling/2.0 - (*list)[j].params.fLower) / nf ;
      ++ndx;

	  *list = (InspiralTemplateList *)
      LALRealloc( *list, ndx * sizeof(InspiralTemplateList) );
      if ( ! *list )
      {
        ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
      memset( *list + ndx - 1, 0, sizeof(InspiralTemplate) );
      (*list)[ndx-1] = (*list)[j];
      (*list)[ndx-1].params.fFinal = fendBCV;
      (*list)[ndx-1].metric = (*list)[0].metric;
      (*list)[ndx-1].nLayer = i;
    }
  }


  for ( j = nlist; j < ndx; ++j )
  {
    (*list)[j-nlist] = (*list)[j];
  }

 *NList = ndx - nlist;
  *list = LALRealloc( *list, *NList * sizeof(InspiralTemplateList) );
  if ( ! *list )
  {
    ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
  }

  RETURN( status );
}

/** UNDOCUMENTED */
void
LALEmpiricalPSItoMassesConversion (
    LALStatus 			*status,
    InspiralTemplate 	*params,
    UINT4            	*valid,
    REAL4             	lightring
    )
{

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  if ( params->psi0 <= 0.L || params->psi3 >= 0.L )
  {
    *valid = 0;
  }
  else
  {
    params->totalMass = -params->psi3/(16.L*LAL_PI * LAL_PI * params->psi0);
    params->totalMass = params->totalMass * 2.  ; /* The factor 2 is purely empiricail and
						      comes from simulaitons. ?It seems indeed
						      tyhat the relation between psi0 and psi3
						      which gives the total mass is not really
						      suitable. Ususally, the total mass is
						      twice as much as the estimated one.
						   */
    params->fFinal = 1.L/( LAL_PI * pow(lightring, 1.5L) * params->totalMass );
    params->totalMass /= LAL_MTSUN_SI; /* it it used later ? */

    *valid = 1;
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/** UNDOCUMENTED */
void
LALInspiralCreateFlatBankS3S4 (
    LALStatus            *status,
    REAL4VectorSequence  *list,
    InspiralBankParams   *bankParams,
    InspiralCoarseBankIn coarseIn
    )

{
  InspiralMetric *metric;
  REAL8 minimalMatch;
  REAL8 x0, x1, dx1=0, dx0=0, x=0, y=0;
  UINT4 nlist = 0;
  INT4 layer  = 1;
  INT4 valid = -1;
  REAL4 xp[8] = {3000, 40000, 100000, 300000, 550000, 550000, 250000,3000};
  REAL4 yp[8] = {0, -3000, -3000, -2000, -1500, -300, -300, 0};
  INT4 npol = 8;


  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /* From the knowledge of the metric and the minimal match
     find the constant increments bankParams->dx0 and
     bankParmams->dx1        */
  metric = bankParams->metric;
  minimalMatch = bankParams->minimalMatch;

  switch (coarseIn.gridSpacing){
    case HybridHexagonal:
    case S2BCV:
    case Hexagonal:
      dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
      dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
      dx0 *=3./2./sqrt(2.);
      dx1 *=sqrt(3./2.);
      break;
    case Square:
      dx0 = sqrt(2.L * (1.L - minimalMatch)/metric->g00 );
      dx1 = sqrt(2.L * (1.L - minimalMatch)/metric->g11 );
      break;
    case  HexagonalNotOriented:
      LALInspiralUpdateParams( status->statusPtr,
  			     bankParams, *metric, minimalMatch );
      CHECKSTATUSPTR( status );
      dx0 = bankParams->dx0 * 3./2./sqrt(2.);
      dx1 = bankParams->dx1 * sqrt(3./2.);
      break;

    case  SquareNotOriented:
      LALInspiralUpdateParams( status->statusPtr,
  			     bankParams, *metric, minimalMatch );
      CHECKSTATUSPTR( status );
      dx0 = bankParams->dx0;
      dx1 = bankParams->dx1;
      break;
  }


  switch (coarseIn.gridSpacing)
  {
    case HybridHexagonal:
    case S2BCV:
    case Hexagonal:
    case HexagonalNotOriented:

    /* x1==psi3 and x0==psi0 */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
      {
	layer++;
	for (x0 = bankParams->x0Min - 1e6 +dx0/2.*(layer%2); x0 <= bankParams->x0Max+1e6; x0 += dx0 )
	{
	  UINT4 ndx = 2 * nlist;
	  if ( coarseIn.gridSpacing == Hexagonal)
	  {
	    x =  x0 *cos(metric->theta) + sin(metric->theta)* x1;
		y =  x0 *sin(metric->theta) - cos(metric->theta)* x1;
	  }
	  else
	  {
		x = x0;
		y = x1;
	  }

	  if ( (x > bankParams->x0Min -dx0/2.) && (y < bankParams->x1Max + dx1/2.) &&
		 (x < bankParams->x0Max +dx0/2.) && (y > bankParams->x1Min - dx1/2.))
	  {

		if (coarseIn.insidePolygon == True)
		{
          LALInsidePolygon(status->statusPtr, xp, yp, npol, x, y, &valid);
		}
        else
		{
		  LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
		}


        if (valid == 1)
        {
          list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
          if ( !list->data )
          {
            ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
          }
          list->data[ndx] = x;
          list->data[ndx + 1] = y;
          ++nlist;
        }
      }
	}
  }
  break;

  case  Square:
  case  SquareNotOriented:

    /* !! dx1 and dx0 are computed in a different way de[pending on the
       value of BANKGRId */
    for (x1 = bankParams->x1Min -1e6;  x1 <= bankParams->x1Max + 1e6; x1 += dx1)
    {
	  for (x0 = bankParams->x0Min - 1e6 ; x0 <= bankParams->x0Max+1e6; x0 += dx0 )
	  {
	    UINT4 ndx = 2 * nlist;

	    if (coarseIn.gridSpacing == Square)
	    {
		  x =  x0 *cos(metric->theta) + sin(metric->theta)* x1 ;
		  y =  x0 *sin(metric->theta) - cos(metric->theta)* x1;
	    }
	    else if (coarseIn.gridSpacing == SquareNotOriented)
	    {
		  x = x0;
		  y = x1;
	    }
	    if ( (x > bankParams->x0Min - dx0/2.) && (y < bankParams->x1Max + dx1/2.) &&
		 (x < bankParams->x0Max + dx0/2.) && (y > bankParams->x1Min - dx1/2.))

	    {

		  if (coarseIn.insidePolygon == True)
		    {
		      LALInsidePolygon(status->statusPtr, xp, yp, npol, x, y, &valid);
            }
		  else
		  {
		    LALExcludeTemplate(status->statusPtr, &valid, bankParams, x, y);
		  }
          if (valid)
          {
            list->data = (REAL4 *) LALRealloc( list->data, (ndx+2) * sizeof(REAL4) );
            if ( !list->data )
            {
              ABORT(status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM);
		    }
            list->data[ndx] = x;
            list->data[ndx + 1] = y;
            ++nlist;
          }
	    }
	  }
     }
  break;
  }

  list->length = nlist;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}


/**
 * Thomas: 31 Aug 2006. This function is redundant with the polygon fit.
 * It was design for BBH and therefore had tight boundary. For a more general
 * purpose, I extended the range to generous values
 */
void
LALExcludeTemplate(
    LALStatus            *status,
    INT4                 *valid,
    InspiralBankParams   UNUSED *bankParams,
    REAL4                 x,
    REAL4                 y)
{
  REAL4 psi0Int = 2500000.;
  REAL4 psi3Int = -10000.;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  if (x > psi0Int && y < psi3Int )
  {
    *valid = 0 ;
  }
  else
  {
    *valid = 1;
  }

  DETATCHSTATUSPTR(status);
  RETURN (status);
}
/*@}*/
