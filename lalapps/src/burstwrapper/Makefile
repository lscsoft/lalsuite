CC = gcc

LDFLAGS= \
	-L${HOME}/datacond/lib/datacondAPI \
	-L${HOME}/datacond/lib/genericAPI \
	-L${HOME}/datacond/lib \
	-L/ldcg/lib \
	-L${HOME}/linux/lib \
	--rpath=$${HOME}/datacond/lib/datacondAPI \
	--rpath=$${HOME}/datacond/lib/genericAPI \
	--rpath=$${HOME}/datacond/lib \
	--rpath=/ldcg/lib

LIBS=-ldatacondAPI -lgenericAPI -lstlport_gcc -llal -llalframe -llalapps -lm 

INCLUDES=-I/home/jsylvest/datacond/include -I/home/jsylvest/linux/include -I../ -I/ldcg/include/stlport -I/ldcg/include -I. 

OBJS=getFrameData.o  getNonFrameData.o  InitSearch.o  RunSearch.o Translation.o ConditionData.o ReConditionData.o

CFLAGS+=-g $(INCLUDES) 
CXXFLAGS+=-g $(INCLUDES) 

burstdso : $(OBJS) burstdso.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o burstdso burstdso.o $(OBJS) $(LIBS)
