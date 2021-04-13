SUBDIRS = lalapps

if LALFRAME
SUBDIRS += frametools
endif

if LALINSPIRAL
SUBDIRS += findchirp
endif

if LALFRAME
if LALMETAIO
if LALINSPIRAL
SUBDIRS += ring
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
if LALMETAIO
if LALINSPIRAL
SUBDIRS += inspiral
endif
endif
endif
