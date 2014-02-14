## Process this file with automake to produce Makefile.in

# for using lalapps-defines + functions
AM_CPPFLAGS = -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

cfsHOME = @CFS_DIR@
hsHOME = @HS_DIR@
eaHOME = $(hsHOME)/EinsteinAtHome
optHOME = $(cfsHOME)/OptimizedCFS
lalappsDIR = @LALAPPS_DIR@
CLEANFILES = 

NVCC ?= nvcc
if WIN32_CYGWIN_BUILD
.cu.o:
	$(NVCC) --compiler-options /GS- --compiler-bindir '$(MSCDIR)' -c `cygpath -w $<` -o `cygpath -w $@`
else
if MAC_BUILD
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@  --compiler-bindir `which gcc-4.0`
else
.cu.o:
	$(NVCC) $(CPPFLAGS) -c $< -o $@
endif
endif

$(optHOME)/FStatOpenCLKernel.c: $(optHOME)/FStatOpenCLKernel.cl
	echo 'char* opencl_kernel_text="\\' > $@ &&\
	sed 's/\\/\\\\/g;s/"/\\"/g;s/$$/\\n\\/' $< >> $@ &&\
	echo '";' >> $@ ||\
	rm -f $@


BUILT_SOURCES = $(optHOME)/FStatOpenCLKernel.c

## start with some common defaults
eah_sources = 	$(lalappsDIR)/lalapps.c \
		$(hsHOME)/HierarchicalSearch.c $(hsHOME)/StackSlideFstat.c \
		$(hsHOME)/HoughFStatToplist.c \
		$(cfsHOME)/HeapToplist.c $(optHOME)/ComputeFstatREAL4.c \
		$(optHOME)/LocalComputeFstat.c $(eaHOME)/LocalComputeFstatHoughMap.c \
		$(eaHOME)/hs_boinc_extras.c $(eaHOME)/hs_boinc_options.cpp
if USE_CUDA
    eah_sources += $(optHOME)/ComputeFstatREAL4CUDA.c $(optHOME)/FStatCUDA.cu
endif
if USE_OPENCL
    eah_sources += $(optHOME)/FStatOpenCLKernel.c $(optHOME)/ComputeFstatREAL4OpenCL.c
endif
eah_cppflags = -Wall -DEAH_BOINC -DUSE_BOINC -DBOINC_COMPRESS -DBOINC_APIV6 -DTARGET='"@BUILD_TARGET@"'
eah_ldadd = -lgsl -lgslcblas -lboinc_zip -lboinc_api -lboinc
eah_ldflags = 

## ----------
## some Mac-specific stuff to allow setting of app-icon
if MAC_BUILD
   eah_cppflags += -I/Developer/Headers/FlatCarbon
endif

if MINGW32_BUILD
   eah_sources += $(eaHOME)/win_lib.cpp $(eaHOME)/exchndl.c
   eah_ldadd += -lbfd -liberty -lintl
   eah_cppflags += -DHAVE_EXCHNDL=1
endif

if WIN32_CYGWIN_BUILD
   eah_sources += $(eaHOME)/win_lib.cpp
endif

if LINUX_BUILD
   eah_sources += $(eaHOME)/erp_execinfo_plus.c
   eah_ldadd += -lbfd -liberty
endif

## ----------
## There are three possible ways in which we can build the E@H application
## BOINC_GRAPHICS = 0: completely non-graphical (USE_GRAPHICS=false)

eah_cppflags += -DNO_BOINC_GRAPHICS -DBOINC_GRAPHICS=0 @CUDA_CFLAGS@
eah_ldadd += @CUDA_LIBS@

bin_PROGRAMS = cfsBOINC

if USE_GRAPHICS
   bin_PROGRAMS += starsphere
endif

EXTRA_PROGRAMS = cfsSCI cfsBOINC.so

## -------- Einstein@Home Main Code
## build the main-code using the settings determined above
cfsBOINC_SOURCES = $(eah_sources)
cfsBOINC_CPPFLAGS = $(eah_cppflags)
cfsBOINC_LDADD = $(eah_ldadd) @MAC_LIBS@
cfsBOINC_LDFLAGS = $(eah_ldflags)

## standlone CFS science-code for testing
cfsSCI_SOURCES = $(cfsHOME)/ComputeFStatistic.c $(lalappsDIR)/lalapps.c \
	$(cfsHOME)/clusters.c $(cfsHOME)/DopplerScan.c

## Fstats-comparison code (useful for testing using testStandalone.sh)
#compareFstats_SOURCES = $(cfsHOME)/compareFstats.c $(lalappsDIR)/lalapps.c
#compareFstats_LDADD = -lgsl -lgslcblas

## make sure rebuilding uses the right aclocal-switches
ACLOCAL_AMFLAGS = -I m4

starsphere_CPPFLAGS = -Wall -DUSE_BOINC -DBOINC_APP_GRAPHICS -DBOINC_APIV6 -DBOINC_GRAPHICS=0 @GRAPHICS_CFLAGS@ -Igraphics -I/usr/X11R6/include
starsphere_LDADD    = -lboinc_graphics2 -lboinc_api -lboinc @GRAPHICS_LIBS@ @MAC_LIBS@ # -lglut
starsphere_SOURCES  = graphics/starlist.C graphics/snr_list.C graphics/user_text.C graphics/pulsar_list.C graphics/search_info.C graphics/starsphere.C graphics/main.C



## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
