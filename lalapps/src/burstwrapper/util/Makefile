LAL = $(LALLOCATION)
MATLAB = /apps/matlab/

CFLAGS = -Wall -I. -g -O
CC = gcc $(CFLAGS)

all : reab sblit ShellInjection ShellInjection1
	cp reab sblit ShellInjection ShellInjection1 ../

ShellInjection : ShellInjection.c Makefile
	gcc -o ShellInjection ShellInjection.c -lm

ShellInjection1 : ShellInjection1.c Makefile
	gcc -o ShellInjection1 ShellInjection1.c -lm

sblit : sblit.o metaio.o 
	$(CC) -o sblit sblit.o metaio.o -O

gtoc : gtoc.o metaio.o 
	$(CC) -o gtoc gtoc.o metaio.o -O

reab : reab.o ReadBurstBin.o
	$(CC) -o reab reab.o ReadBurstBin.o -O

metaio.o : metaio.c metaio.h Makefile

ReadBurstBin.o : ReadBurstBin.c Makefile


#### MATLAB ####

# for linux only

matlab : reab.mexglx LALBurst.mexglx LALBurstOutput.mexglx
	cp reab.mexglx LALBurst.mexglx LALBurstOutput.mexglx ../matlab

reab.mexglx : reab_m.c ReadBurstBin.c
	mcc -x reab_m.c ReadBurstBin.c -M -Dlinux
	mv reab_m.mexglx reab.mexglx

LALBurst.mexglx: LALBurst.c
	mcc -x LALBurst.c -M '-I$(LAL)/include $(LAL)/lib/liblal.a'

LALBurstOutput.mexglx : LALBurstOutput.c StdBurstSearch.c
	mcc -x LALBurstOutput.c StdBurstSearch.c -M '-I$(LAL)/include $(LAL)/lib/liblal.a'


clean : 
	rm -f *.o reab.mex* LALBurst.mex* LALBurstOutput.mex* reab sblit ShellInjection ShellInjection1
