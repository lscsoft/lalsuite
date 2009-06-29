# Configuration Script for auto_resamp.py script. Times are in seconds.
# For parameters with a _min and _max function, if parameter > 0, then it is 
# set as the value, else it is chosen uniformly from _min to _max

# Strength of signal
h0 = 0.1

# Cosine of iota
cosi = 0
cosi_min = 0
cosi_max = 1

# Polarization Angle.
psi = 0.5
psi_min = 0
psi_max = 1

# Initial Phase
phi0 = 0
phi0_min = 0
phi0_max = 1

# Number of Dirichlet terms used.
Dterms = 64

# Interferometer
IFO = 'H2'

# Start Time
t0 = 820000000 

# Reference Time in SSB
refTime = 820000000

# Output Directory
Out= './SFTs'

# Ephemeris Directory 
Ephem = '/Users/ppatel/home/install/lal/share/lal'

# Ephemeris Year
EphemYear = '05-09'

# Noise Sh
Sh = 2

# Duration of Analysis
TSpan = 25000

# SFT time baseline
TSFT = 1800

# Number of SFTs to add
NumSFTs = 10

# Number of Gaps to add
NumGaps = 5

# Alpha (Right Ascension)
Alpha = 0
Alpha_min = 0
Alpha_max = 6.28

# Delta (Declination)
Delta = 0.8
Delta_min = -1.57
Delta_max = 1.57

# Minimum Frequency
Fmin = 4

# Band of Analysis
Band = 0.5

# Injection Frequency
Finj = 4.9
Finj_min = 55.0
Finj_max = 55.1

# Spindown
FDot = 1e-11*0
FDot_min = 1e-10
FDot_max = 1e-11

# Optional debug
debug = 1

# Resolution
Res = 1.0/TSpan/3

# OutputTimeSeries
TimeSeriesOut = 'TSeries'

# TimeStampsFile
#TimeStampsFile = 'TimeStampsFile'
 
