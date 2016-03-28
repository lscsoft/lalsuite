# General tools which should be tested first
SUBDIRS = \
	Tools \
	Injections \
	SFTTools \
	FITSTools \
	$(END_OF_LINE)

# Everything else
SUBDIRS += \
	CreateEphemeris \
	Fstatistic \
	GCT \
	HeterodyneSearch \
	HoughFstat \
	MakeSFTs \
	SidebandSearch \
	TwoSpect \
	Xray \
	crosscorr \
	fscan \
	hough \
	$(END_OF_LINE)

# Because many tests in src/pulsar call executables from other subdirectories in
# src/pulsar, it is safest to make sure the whole of src/pulsar is built first
check-recursive: all-recursive
