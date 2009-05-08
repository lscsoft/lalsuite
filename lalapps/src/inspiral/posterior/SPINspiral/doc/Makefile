
#Use GNU C compiler
#CC = gcc
#CFLAGS = -O3 
#CFLAGS = -O3 -static 
#CFLAGS = -O0 -Wall
#CFLAGS = -g -O0 


#Use Intel C compiler
CC = icc

#wd: disable: remark #981: operands are evaluated in unspecified order
CFLAGS = -O3 -traceback
#CFLAGS = -O3 -static -traceback
#CFLAGS = -O3 -static -traceback -Wcheck -Wall -wd981 
#CFLAGS = -O0 -traceback -Wcheck -Wall -wd981 
#CFLAGS = -g -O0  #-traceback -Wcheck -Wall -wd981 



MCMCOBJ = mcmc_main.o mcmc_mcmc.o mcmc_templates.o mcmc_routines.o mcmc_data.o mcmc_parameters.o mcmc_signal.o mcmc_3rdparty.o 
DEFCFLAGS = -I/usr/include -I./include -I ~/bin/include   -I ~/local/programs/lscsoft/lal/include/
DEFLFLAGS = -L ~/bin/lib -lFrame -lfftw3 -lgsl -lgslcblas -lm
LALLFLAGS = -L ~/bin/lib -llal 





%.o: %.c Makefile include/mcmc.h
	$(CC) $(CFLAGS) $(DEFCFLAGS) -c $<

mcmc_lal.o: include/mcmc_lal.h
mcmc_3rdparty.o: include/remez.h



mcmc_nolal: $(MCMCOBJ) mcmc_nolal.o
	$(CC) -o mcmc_nolal $(MCMCOBJ) mcmc_nolal.o $(DEFLFLAGS)

mcmc_lal: $(MCMCOBJ) mcmc_lal.o
	$(CC) -o mcmc_lal $(MCMCOBJ) mcmc_lal.o $(DEFLFLAGS) $(LALLFLAGS)




clean:
	rm -f *.o include/*.gch

cleanbin:
	rm -f mcmc 

cleanoutput:
	rm -f mcmc.input.?????? mcmc.output.??????.?? *-signal*.dat

cleanall: clean cleanbin cleanoutput

