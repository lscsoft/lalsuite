# Tools which should be tested early
SUBDIRS = \
	Tools \
	Injections \
	SFTTools \
	FITSTools \
	$(END_OF_LINE)

# Everything else
SUBDIRS += \
	CreateEphemeris \
	CrossCorr \
	Fscan \
	Fstatistic \
	GCT \
	HeterodyneSearch \
	Hough \
	HoughFstat \
	MakeSFTs \
	SidebandSearch \
	TwoSpect \
	Xray \
	$(END_OF_LINE)

# Because many tests in src/pulsar call executables from other subdirectories in
# src/pulsar, it is safest to make sure the whole of src/pulsar is built first
check-recursive: all-recursive
