SUBDIRS = lalapps

SUBDIRS += calibration

if LALFRAME
SUBDIRS += frametools
endif

if LALINSPIRAL
SUBDIRS += findchirp
endif

if LALSTOCHASTIC
SUBDIRS += stochastic
endif

if FRAME
if LALFRAME
if LALMETAIO
if LALINSPIRAL
SUBDIRS += ring
endif
endif
endif
endif

if LALFRAME
if LALMETAIO
if LALBURST
SUBDIRS += string
endif
endif
endif

if LALFRAME
if LALMETAIO
if LALBURST
SUBDIRS += power
endif
endif
endif

if LALPULSAR
SUBDIRS += pulsar
endif

if LALFRAME
if LALINSPIRAL
SUBDIRS += inspiral
endif
endif

if FRAME
if LALFRAME
SUBDIRS += online
endif
endif

SUBDIRS += hwinjection
