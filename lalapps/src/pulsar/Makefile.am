SUBDIRS = \
	MakeSFTs \
	Injections \
	FITSTools \
	Fstatistic \
	hough \
	GCT \
	crosscorr \
	HeterodyneSearch \
	fscan \
	TwoSpect \
	SFTReferenceLibrary \
	SidebandSearch \
	CreateEphemeris \
	Tools \
	Xray \
	$(END_OF_LINE)

# Because many tests in src/pulsar call executables from other subdirectories in
# src/pulsar, it is safest to make sure the whole of src/pulsar is built first
check-recursive: all-recursive
