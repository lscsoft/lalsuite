CC = gcc
CFLAGS = -g -W -Wall
CPPFLAGS = `pkg-config --cflags lalsupport lalframe`
LDFLAGS =
LIBS = `pkg-config --libs lalsupport lalframe`
COMPILE = $(CC) $(CFLAGS) $(CPPFLAGS)
LINK = $(CC) $(CFLAGS) $(LDFLAGS)

all: lalapps_exttrig_crosscorr

.c.o:
	$(COMPILE) -c $<

lalapps_exttrig_crosscorr: exttrig_crosscorr.o crosscorr.o whiten.o inject.o signal.o
	$(LINK) -o $@ $^ $(LIBS)
