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


#ifndef _OPT_HOUGH_AMD64_LINUX_H
#define _OPT_HOUGH_AMD64_LINUX_H

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
__asm __volatile (					\
	"xor %%r11,%%r11			\n\t"	\
	"mov %[yUpper] , %%eax  		\n\t"	\
	"test %%eax,%%eax			\n\t"	\
	"mov %[xPixel], %%R8  			\n\t"   \
	"mov %[yLower], %%r11d  		\n\t" 	\
	"lea (%%R8,%%r11,0x2), %%rsi  		\n\t"	\
	"mov %[xSideP1], %%edx   		\n\t"	\
	"js 3f					\n\t"	\
	"lea -0x2(%%r8,%%rax,0x2),%%R8  	\n\t"	\
	"mov %[map] , %%rdi  			\n\t"	\
	"mov %%r11d,%%eax  			\n\t"   \
	"imul %%edx, %%eax  			\n\t"	\
	"lea (%%rdi, %%rax, 0x8), %%rdi  	\n\t"	\
	"movsd %[w],%%xmm0  			\n\t"	\
	"cmp  %%r8,%%rsi  			\n\t"	\
	"jmp  2f 				\n\t"	\
AD_ALIGN32 "\n"	\
	"1:  					\n\t"	\
	"movzwl (%%rsi),%%r11d			\n\t"	\
	"movzwl 2(%%rsi),%%ecx			\n\t"	\
	"lea (%%rdi, %%r11, 0x8) , %%r9  	\n\t"	\
	"movsd (%%r9),%%xmm1  			\n\t"	\
	"lea (%%rdi,%%rdx,0x8) , %%rdi  	\n\t"	\
	"lea (%%rdi,%%rcx,0x8) , %%r10   	\n\t"	\
	"movsd (%%r10),%%xmm2  			\n\t"	\
	"addsd %%xmm0,%%xmm1  			\n\t"	\
	"movsd %%xmm1,(%%r9)  			\n\t"	\
	"addsd %%xmm0,%%xmm2	  		\n\t"	\
	"movsd %%xmm2,(%%r10)  			\n\t"	\
	"lea (%%rdi,%%rdx,0x8), %%rdi   	\n\t"	\
	"lea 4(%%rsi) , %%rsi   		\n\t"	\
	"cmp  %%r8,%%rsi       			\n"	\
	"2:	  				\n\t"	\
	"jbe 1b	  				\n\t"	\
	"add $0x2,%%r8				\n\t"	\
	"cmp %%r8,%%rsi			\n\t"	\
	"jne 3f  				\n\t"	\
	"movzwl (%%rsi) , %%r11d  		\n\t"	\
	"lea (%%rdi, %%r11, 0x8) , %%r9  	\n\t"	\
	"movsd (%%r9),%%xmm1  			\n\t"	\
	"addsd %%xmm0,%%xmm1  			\n\t"	\
	"movsd %%xmm1,(%%r9)  			\n"	\
	"3:  					\n\t"	\
	"NOP  					\n"	\
	: 	\
	:	\
	[xPixel]  "m" (_XPIXEL) ,		 \
	[yLower]  "m" (_YLOWER) ,		 \
	[yUpper]  "m" (_YUPPER),		 \
	[xSideP1] "m" (_XSIDEP1) ,		 \
	[map]     "m" (_MAP) ,			 \
	[w]       "m" (_WEIGHT)			 \
	:					 \
	"memory","r8","r9","r10","r11","rax", "rcx", "rdx", "rsi", "rdi", "cc",	 \
	"xmm0","xmm1","xmm2"	 \
	)

#endif
