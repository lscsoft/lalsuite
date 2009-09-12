#Works on Caltech cluster

#Use GNU C compiler
#CC = gcc
CC = condor_compile gcc
CFLAGS = -O3 
#CFLAGS = -O3 -static 
#CFLAGS = -O0 -Wall
#CFLAGS = -g -O0 


#Use Intel C compiler
#CC = icc

#wd: disable: remark #981: operands are evaluated in unspecified order
#CFLAGS = -O3 -traceback
#CFLAGS = -O3 -static -traceback
#CFLAGS = -O3 -static -traceback -Wcheck -Wall -wd981 
#CFLAGS = -O0 -traceback -Wcheck -Wall -wd981 
#CFLAGS = -g -O0  #-traceback -Wcheck -Wall -wd981 



MCMCOBJ = SPINspiral_main.o SPINspiral_mcmc.o SPINspiral_templates.o SPINspiral_routines.o SPINspiral_data.o SPINspiral_parameters.o SPINspiral_signal.o SPINspiral_lal.o SPINspiral_3rdparty.o 
DEF_CFLAGS = -I/usr/include -I./include -I ~/bin/include   -I /opt/lscsoft/lal/include -I /opt/lscsoft/libmetaio/include  -I /opt/lscsoft/libframe/include
DEF_LFLAGS = -L/opt/lscsoft/libframe/lib64 -lFrame -lfftw3 -lfftw3f -lgsl -lgslcblas -lm  		
LAL_LFLAGS = -L /opt/lscsoft/lal/lib64/ -llal -llalmetaio




%.o: %.c Makefile include/SPINspiral.h
	$(CC) $(CFLAGS) $(DEF_CFLAGS) -c $<

all: SPINspiral

exit: 




SPINspiral_3rdparty.o: include/remez.h



SPINspiral: $(MCMCOBJ)
	$(CC) -o SPINspiral $(MCMCOBJ) $(LAL_LFLAGS) $(DEF_LFLAGS)



clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f SPINspiral

cleanoutput:
	rm -f mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

