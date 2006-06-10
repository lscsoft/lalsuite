          {
	    NRCSID (CFSLOOPALTIVTAG, "$Id$");

            COMPLEX8 *Xalpha_k = Xalpha + sftIndex;

            realXP=0.0;
            imagXP=0.0;

	    /* VERSION 6 - Altivec/scalar hybrid unrolled */
	    {

	      vector float tsin_v, tcos_v;
	      vector float realXP_v0 = (vector float)(0.0f);
	      vector float imagXP_v0 = (vector float)(0.0f);
	      vector float realXP_v1 = (vector float)(0.0f);
	      vector float imagXP_v1 = (vector float)(0.0f);
	      
	      if (klim > tempF_size) {
		tempf = realloc(tempf, sizeof(REAL4)*klim);
		tempF_size = klim;
	      }
	      
	      REAL8 tempf1 = tempFreq1;
	      for (k=0; k+3<klim; k+=4) {
		tempf[k+0] = tempf1 - 0.0;
		tempf[k+1] = tempf1 - 1.0;
		tempf[k+2] = tempf1 - 2.0;
		tempf[k+3] = tempf1 - 3.0;
		tempf1-=4;
	      }
	      for (; k<klim; k++) {
		tempf[k] = tempf1;
		tempf1--;
	      }
	      
	      tsin_v = vec_ld( 0, &tsin );
	      tsin_v = vec_perm( tsin_v, tsin_v, vec_lvsl( 0, &tsin ) );
	      tsin_v = vec_splat( tsin_v, 0 );
	      
	      tcos_v = vec_ld( 0, &tcos );
	      tcos_v = vec_perm( tcos_v, tcos_v, vec_lvsl( 0, &tcos ) );
	      tcos_v = vec_splat( tcos_v, 0 );
	      
	      vector unsigned char permute_v = vec_lvsl( 0, (float *) Xalpha_k );
	      
	      /* Loop over terms in dirichlet Kernel */
	      for(k=0; k+7 < klim ; k+=8)
		{
		  vector float realP_v0, imagP_v0;
		  vector float realP_v1, imagP_v1;
		  vector float xinv_v0, xinv_v1;
		  vector float Xa_re_v0, Xa_im_v0;
		  vector float Xa_re_v1, Xa_im_v1;
		  vector float temp1, temp2, temp3, temp4, temp5, temp6;
		  vector float tempFreq1_v0 = vec_ld( 0, &tempf[k] );
		  vector float tempFreq1_v1 = vec_ld( 16, &tempf[k] );
		  
		  //Get the reciprocal estimate
		  vector float estimate0 = vec_re( tempFreq1_v0 );
		  vector float estimate1 = vec_re( tempFreq1_v1 );
		  //One round of Newton-Raphson refinement
		  estimate0 = vec_madd( vec_nmsub( estimate0, tempFreq1_v0, (vector float) (1.0) ), estimate0, estimate0 );
		  estimate1 = vec_madd( vec_nmsub( estimate1, tempFreq1_v1, (vector float) (1.0) ), estimate1, estimate1 );
		  xinv_v0 = vec_madd( (vector float)(OOTWOPI), estimate0, (vector float)(0) );
		  xinv_v1 = vec_madd( (vector float)(OOTWOPI), estimate1, (vector float)(0) );
		  //xinv_v0 = vec_div( (vector float)(OOTWOPI), tempFreq1_v0 );
		  //xinv_v1 = vec_div( (vector float)(OOTWOPI), tempFreq1_v1 );
		  
		  temp1 = vec_ld( 0, (float *) Xalpha_k );
		  temp2 = vec_ld( 16, (float *) Xalpha_k );
		  temp3 = vec_ld( 32, (float *) Xalpha_k );
		  temp4 = vec_ld( 48, (float *) Xalpha_k );
		  temp5 = vec_ld( 63, (float *) Xalpha_k );
		  
		  temp1 = vec_perm( temp1, temp2, permute_v );
		  temp2 = vec_perm( temp2, temp3, permute_v );
		  temp3 = vec_perm( temp3, temp4, permute_v );
		  temp4 = vec_perm( temp4, temp5, permute_v );
		  
		  temp5 = vec_mergeh( temp1, temp2 );
		  temp6 = vec_mergel( temp1, temp2 );
		  Xa_re_v0 = vec_mergeh( temp5, temp6 );
		  Xa_im_v0 = vec_mergel( temp5, temp6 );
		  
		  temp5 = vec_mergeh( temp3, temp4 );
		  temp6 = vec_mergel( temp3, temp4 );
		  Xa_re_v1 = vec_mergeh( temp5, temp6 );
		  Xa_im_v1 = vec_mergel( temp5, temp6 );
		  
		  Xalpha_k += 8;
               
		  realP_v0 = vec_madd( tsin_v, xinv_v0, (vector float)(0.0f) );
		  imagP_v0 = vec_madd( tcos_v, xinv_v0, (vector float)(0.0f) );
		  realP_v1 = vec_madd( tsin_v, xinv_v1, (vector float)(0.0f) );
		  imagP_v1 = vec_madd( tcos_v, xinv_v1, (vector float)(0.0f) );
		  
		  /* realXP_v = real_XP_v + Xa_re_v * realP_v - Xa_im_v * imagP_v; */
		  realXP_v0 = vec_madd( Xa_re_v0, realP_v0, realXP_v0 );
		  realXP_v0 = vec_nmsub( Xa_im_v0, imagP_v0, realXP_v0 );
		  imagXP_v0 = vec_madd( Xa_re_v0, imagP_v0, imagXP_v0 );
		  imagXP_v0 = vec_madd( Xa_im_v0, realP_v0, imagXP_v0 );				
		  
		  realXP_v1 = vec_madd( Xa_re_v1, realP_v1, realXP_v1 );
		  realXP_v1 = vec_nmsub( Xa_im_v1, imagP_v1, realXP_v1 );
		  imagXP_v1 = vec_madd( Xa_re_v1, imagP_v1, imagXP_v1 );
		  imagXP_v1 = vec_madd( Xa_im_v1, realP_v1, imagXP_v1 );				
		} /* for k < klim */
	      {
		float realXP_float, imagXP_float;
		vector float re_sum = vec_add( realXP_v0, realXP_v1 );
		vector float im_sum = vec_add( imagXP_v0, imagXP_v1 );
		re_sum = vec_add( re_sum, vec_sld( re_sum, re_sum, 8 ) );
		im_sum = vec_add( im_sum, vec_sld( im_sum, im_sum, 8 ) );
		re_sum = vec_add( re_sum, vec_sld( re_sum, re_sum, 4 ) );
		im_sum = vec_add( im_sum, vec_sld( im_sum, im_sum, 4 ) );
		vec_ste( re_sum, 0, &realXP_float);
		vec_ste( im_sum, 0, &imagXP_float);
		realXP = realXP_float;
		imagXP = imagXP_float;
	      }
	      tempFreq1 = tempFreq1 - k;

	      for(; k < klim ; k++)
		{
		  REAL4 xinv = (REAL4)OOTWOPI / (REAL4)tempFreq1;
		  REAL4 Xa_re = Xalpha_k->re;
		  REAL4 Xa_im = Xalpha_k->im;
		  Xalpha_k ++;
		  tempFreq1 --;
		  
		  realP = tsin * xinv;
		  imagP = tcos * xinv;
		  /* compute P*xtilde */
		  realXP += Xa_re * realP - Xa_im * imagP;
		  imagXP += Xa_re * imagP + Xa_im * realP;
		  
		} /* for k < klim */
	    }
	    
          } /* if x cannot be close to 0 */
        
