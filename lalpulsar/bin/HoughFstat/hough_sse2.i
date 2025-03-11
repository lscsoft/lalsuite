/*
 *  Copyright (C) 2008 Bernd Machenschalk
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
 *
 */

#ifndef _OPT_HOUGH_I686_SSE2_H
#define _OPT_HOUGH_I686_SSE2_H

#ifdef __APPLE__
#define AD_FLOAT ".single "
#define AD_ASCII ".ascii "
#define AD_ALIGN16 ".align 4"
#define AD_ALIGN32 ".align 5"
#define AD_ALIGN64 ".align 6"
#else /* x86 gas */
#define AD_FLOAT ".float "
#define AD_ASCII ".string "
#define AD_ALIGN16 ".align 16"
#define AD_ALIGN32 ".align 32"
#define AD_ALIGN64 ".align 64"
#endif

#define ADDPHMD2HD_WLR_LOOP(_XPIXEL,_YLOWER,_YUPPER,_XSIDEP1,_MAP,_WEIGHT)\
__asm __volatile ( 				      \
	"push %%ebx				\n\t" \
	"mov %[xPixel], %%eax  			\n\t" \
	"mov %[yLower], %%ebx  			\n\t" \
	"lea (%%eax,%%ebx,0x2), %%esi  		\n\t" \
	"mov %[xSideP1], %%edx   		\n\t" \
	"mov %[yUpper] , %%edi  		\n\t" \
	"lea -0x2(%%eax,%%edi,0x2),%%eax  	\n\t" \
	"mov %[map] , %%edi  			\n\t" \
	"mov %%ebx,%%ecx  			\n\t" \
	"imul %%edx, %%ecx  			\n\t" \
	"lea (%%edi, %%ecx, 0x8), %%edi  	\n\t" \
	"movsd %[w],%%xmm0  			\n\t" \
	"cmp  %%eax,%%esi  			\n\t" \
	"jmp  2f 				\n\t" \
	AD_ALIGN32                             "\n" \
	"1:  					\n\t" \
	"movzwl (%%esi),%%ebx			\n\t" \
	"movzwl 2(%%esi),%%ecx			\n\t" \
						       \
	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t" \
	"movsd (%%ebx),%%xmm1  			\n\t" \
	"lea (%%edi,%%edx,0x8) , %%edi  	\n\t" \
	"lea (%%edi,%%ecx,0x8) , %%ecx   	\n\t" \
	"movsd (%%ecx),%%xmm2  			\n\t" \
						       \
	"addsd %%xmm0,%%xmm1  			\n\t"	\
	"movsd %%xmm1,(%%ebx)  			\n\t"	\
	"addsd %%xmm0,%%xmm2	  		\n\t"	\
	"movsd %%xmm2,(%%ecx)  			\n\t"	\
	"lea (%%edi,%%edx,0x8), %%edi   	\n\t" \
						      \
	"lea 4(%%esi) , %%esi   		\n\t" \
	"cmp  %%eax,%%esi       		\n" \
						     \
	"2:	  				\n\t" \
	"jbe 1b	  				\n\t" \
	"add $0x2,%%eax				\n\t" \
	"cmp %%eax,%%esi			\n\t" \
	"jne 3f  				\n\t" \
							   \
	"movzwl (%%esi) , %%ebx  		\n\t" \
	"lea (%%edi, %%ebx, 0x8) , %%ebx  	\n\t" \
	"movsd (%%ebx),%%xmm1  			\n\t" \
	"addsd %%xmm0,%%xmm1  			\n\t" \
	"movsd %%xmm1,(%%ebx)  			\n\t" \
						     \
	"3:  					\n\t" \
	"pop %%ebx				\n\t" \
	: 					 \
	:					 \
	[xPixel]  "m" (_XPIXEL) ,		 \
	[yLower]  "m" (_YLOWER) ,		 \
	[yUpper]  "m" (_YUPPER),		 \
	[xSideP1] "m" (_XSIDEP1) ,		 \
	[map]     "m" (_MAP) ,			 \
	[w]       "m" (_WEIGHT)			 \
	:					 \
	"memory","eax", "ecx", "edx", "esi", "edi", "cc", \
	"xmm0","xmm1","xmm2"	 \
	)
#endif
