SUBDIRS = \
	Tools \
	SFTTools \
	FITSTools \
	Injections \
	MakeSFTs \
	Fscan \
	CreateEphemeris \
	Fstatistic \
	HeterodyneSearch \
	Hough \
	HoughFstat \
	GCT \
	CrossCorr \
	SidebandSearch \
	TwoSpect \
	Weave \
	$(END_OF_LINE)

# Because many tests in lalapps/src/pulsar/ call executables from other
# subdirectories in lalapps/src/pulsar/, it is safest to make sure the
# whole of src/pulsar is built first
check-recursive: all-recursive
