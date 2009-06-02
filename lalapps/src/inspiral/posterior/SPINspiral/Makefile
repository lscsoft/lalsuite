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



MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o 
DEF_CFLAGS = -I/usr/include -I./include -I ~/bin/include   -I /opt/lscsoft/lal/include  -I /opt/lscsoft/libframe/include
DEF_LFLAGS = -L/opt/lscsoft/libframe/lib64 -lFrame -lfftw3 -lfftw3f -lgsl -lgslcblas -lm  		
LAL_LFLAGS = -L /opt/lscsoft/lal/lib64/ -llal 




%.o: %.c Makefile include/mcmc.h
	$(CC) $(CFLAGS) $(DEF_CFLAGS) -c $<

all: mcmc_nolal mcmc_lal SPINspiral

exit: 




mcmc_lal.o: include/mcmc_lal.h
mcmc_3rdparty.o: include/remez.h



mcmc_nolal: $(MCMCOBJ) mcmc_nolal.o
	$(CC) -o mcmc_nolal $(MCMCOBJ) mcmc_nolal.o $(DEF_LFLAGS)

mcmc_lal: $(MCMCOBJ) mcmc_lal.o
	$(CC) -o mcmc_lal $(MCMCOBJ) mcmc_lal.o $(LAL_LFLAGS) $(DEF_LFLAGS)

SPINspiral: mcmc_lal
	cp -f mcmc_lal SPINspiral



clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f mcmc_nolal mcmc_lal SPINspiral

cleanoutput:
	rm -f mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

