CC = gcc

#LAL = ${HOME}/LAL
#LALAPPS = ${HOME}/LALAPPS
#DATACOND = ${HOME}/datacond

LAL = $(LALLOCATION)
LALAPPS = $(LALAPPSLOCATION)
DATACOND = $(DATACONDLOCATION)

LDFLAGS= \
	-L$(DATACOND)/lib \
	-L/ldcg/lib \
	-L$(LAL)/lib \
	-L$(LALAPPS)/lib \
	--rpath=$(DATACOND)/lib \
	--rpath=/ldcg/lib

LIBS=-ldatacondCaller -lstlport_gcc -llal -llalframe -llalapps -lm -ldl

INCLUDES=-I$(DATACOND)/include -I$(LAL)/include -I../ -I/ldcg/include/stlport -I/ldcg/include -I. 

OBJS=getFrameData.o  getNonFrameData.o  InitSearch.o  RunSearch.o Translation.o ConditionData.o ReConditionData.o

CFLAGS+=-g $(INCLUDES) 
CXXFLAGS+=-g $(INCLUDES) 


all : burstdso buildcache Standalone_Driver_Example.tar Condor_Driver_Example.tar LDAS_Driver_Example.tar
	cd util && $(MAKE)

.PHONY : matlab
matlab : 
	cd util && $(MAKE) matlab

burstdso : $(OBJS) burstdso.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o burstdso burstdso.o $(OBJS) $(LIBS)

buildcache : buildcache.o 
	$(CXX) $(CXXFLAGS) -o buildcache buildcache.o 

Standalone_Driver_Example.tar : Standalone_Driver_Example
	tar cvf Standalone_Driver_Example.tar Standalone_Driver_Example

Condor_Driver_Example.tar : Condor_Driver_Example
	tar cvf Condor_Driver_Example.tar Condor_Driver_Example

LDAS_Driver_Example.tar : LDAS_Driver_Example
	tar cvf LDAS_Driver_Example.tar LDAS_Driver_Example

clean :
	rm -f *.o burstdso buildcache *.tar
	cd util && $(MAKE) clean
	rm -f matlab/*.mex*
