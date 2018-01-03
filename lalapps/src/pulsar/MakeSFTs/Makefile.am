include $(top_srcdir)/gnuscripts/lalapps.am

if FRAME
frameprogs = lalapps_MakeSFTs
frameexprogs = lalapps_MakeCleanSFTs
else
frameprogs =
frameexprogs =
endif

## PSS requires FrameLib anyway
if PSS
pss_progs =
pss_src = XLALPSSInterface.c
else
pss_progs =
pssinterface =
endif

bin_PROGRAMS = $(frameprogs) $(pss_progs)
EXTRA_PROGRAMS = $(frameexprogs)

if FRAME
lalapps_MakeSFTs_SOURCES = MakeSFTs.c $(pss_src)
lalapps_MakeCleanSFTs_SOURCES = MakeCleanSFTs.c
endif

lalapps_MakeSFTs_CPPFLAGS = $(AM_CPPFLAGS)
if PSS
lalapps_MakeSFTs_CPPFLAGS += -DPSS_ENABLED
endif
