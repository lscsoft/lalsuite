# Configuration Script for auto_resamp.py script. Times are in seconds.
# For parameters with a _min and _max function, if parameter > 0, then it is 
# set as the value, else it is chosen uniformly from _min to _max

# Strength of signal
h0 = 0.2

# Cosine of iota
cosi = -0.3
cosi_min = 0
cosi_max = 1

# Polarization Angle.
psi = 0.6
psi_min = 0
psi_max = 1

# Initial Phase
phi0 = 1.5
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
Sh = 0.1

# Duration of Analysis
TSpan = 36000

# SFT time baseline
TSFT = 1800

# Number of SFTs to add
NumSFTs = 10

# Number of Gaps to add
NumGaps = 3

# Alpha (Right Ascension)
Alpha = 2.0
Alpha_min = 0
Alpha_max = 6.28

# Delta (Declination)
Delta = 0.5
Delta_min = -1.57
Delta_max = 1.57

# Minimum Frequency
#Fmin = 100.12345 - 2e-3
Fmin = 3.0

# Band of Analysis
Band = 2.0

# Injection Frequency
Finj = 4.1
Finj_min = 55.0
Finj_max = 55.1

# Spindown
FDot = 1e-7*0
FDot_min = 1e-10
FDot_max = 1e-11

# Optional debug
debug = 1

# Resolution
Res = 1.0/TSpan/6

# OutputTimeSeries
TimeSeriesOut = 'TSeries'

# TimeStampsFile
TimeStampsFile = 'TimeStampsFileH1'
 
