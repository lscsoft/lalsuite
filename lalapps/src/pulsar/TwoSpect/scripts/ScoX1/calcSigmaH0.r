# Grant David Meadors
# 02014-04-29
# Simple quadrature error calculator
# and unit coverter for asini
# and log10p
# for the Scorpius X-1 MDC
# 02014-04-30 (JD 2456778)

#testOpen <- 'yes'
testOpen <- 'no'


quadErr <- function(errA, errB) {
    return(sqrt((errA)^2 + (errB^2)))
}

write(quadErr(3, 4),stdout())

# Load the data frame called bestSingleIFO
# containing the best single interferometer detections
# Manually load these into the working directory: this is a safeguard
if (testOpen == 'no'){
    # Closed results 
    bestSeen <- load("bestIFOdetectionClosed.Rda")
}
if (testOpen == 'yes'){
    # For the open test
    bestSeen <- load("bestIFOdetectionOpen.Rda")
}
#write(length(bestSingleIFO$freq), stdout())

# We also need to calculate p-values and adjusted the to natural-log
# for the loudest events in the nonDetectedPulsars. I did this by
# looking at the heat maps and reading off their loudest values.

# Their pulsar numbers are:
nonDetectedPulsars <- c(21, 41, 48, 50, 52, 54, 57, 58, 63, 64, 69, 71, 72, 73, 90, 98)
# If we ever needed to access the list out of the 50 closed, they would be
nonDetectedListClosed <- c(11, 18, 21, 22, 24, 25, 26, 27, 32, 33, 38, 39, 40, 41, 48, 50)

pNonDetected <- c(11.3684, 10.8437, 11.954, 12.328, 14.008, 14.4108, 10.8314, 12.3019, 12.2279, 14.309, 10.8493, 11.3576, 12.9252, 12.7191, 14.4123, 34.3318)
pNaturalNonDetected <- -1*log(10)/log(exp(1)) * pNonDetected
pNaturalDetected <- log(10)/log(exp(1)) * bestSingleIFO$Prob

#write('Non-detected natural log p-values', stdout())
#write(pNaturalNonDetected, stdout())
write('Detected natural log p-values', stdout())
write(pNaturalDetected, stdout())

# For non-detected pulsars, the upper limit is (for all)
# 4.23e-25 = 1.11*3.81e-25

# For parameter estimation, we are categorizing anything with a frequency
# >= 1050.0 Hz AND an abs(log10p) < 300 as "Loose", everything else is "Tight"

# In the closed set, pulsars 83, 85 and 95 are "Loose"
# These correspond to number 31, 33, and 34 in the detected closed pulsars
closedLoosePulsars <- c(31,33,34)
# For the open data set, use instead...
openLoosePulsars <- c(27,29,30,31)
if (testOpen == 'yes'){
    loosePulsars <- openLoosePulsars
}
if (testOpen == 'no'){
    loosePulsars <- closedLoosePulsars
}

# The uncertainties in parameter estimation are, for tight
# sigma delta f = 0.0003749 Hz,
# sigma delta asini = 0.01917 s
# for loose,
# sigma delta f = 0.1733 Hz
# sigma delta asini = 0.2752 s

matrixFSigma <- array(0.0003749, length(bestSingleIFO$h0), 1)
matrixFSigma[loosePulsars] <- 0.1733
matrixAsiniSigma <- array(0.01917, length(bestSingleIFO$h0), 1)
matrixAsiniSigma[loosePulsars] <- 0.2752

# Here we output the best frequencies:
write(bestSingleIFO$freq, stdout())

# Despite the name, this field has been overwitten for the bestIFOdetection*
# Rda files and contains asini
write(bestSingleIFO$df, stdout())

# We consistently tested only a single period,
# P = 68023.8259 s
# The uncertainty for this is different for pulsars using a 360 s SFT
# (pulsar 32 and above) than for those
# using a 840 s SFT
# (everything below pulsar 32):
# Pdiff_allowed from the template spacing discussed in the methods paper
# offers an uncertainty of 2.7 * (Tsft/1800) + 1.8, which gives

# Evaluates roughly to Pdiff_840s = 3.06 s
# Evaluates roughly to Pdiff_360s = 2.34 s

# ... except wait, Evan says use equation 23
# delta P = P_0^2 / Pdiff Tobs
# where  P _0 = 68023.8259 s
# T_obs = 31536000 
# He says this gives results of period uncertainty: 
# 48 s for 840 s SFTs
Pdiff840s <- 2.7 * (840/1800) + 1.8
write((68023.8259)^2 / (Pdiff840s * 31536000), stdout())
# 63 s for 360 s SFTs
Pdiff360s <- 2.7 * (360/1800) + 1.8
write((68023.8259)^2 / (Pdiff360s * 31536000), stdout())



# Two points: we may need asini for the table, and we need the quadrature
# sure of errors for h0
# Those errors are, for the intrinsic residual,
# 1.287e-26 is tight,
# 2.501e-26 if loose 
# We also found that a 1.11 rescaling factor needed to be applied to h0
# (already noted in the tight, but not the loose, residual error)
# due to the R->h0 calibration inaccuracy and then further
# a 1.74 factor for random cos iota.
# In addition, the cos iota distribution introduces a +- 37% uncertainty
# (see hEffDist.r)
# at one sigma, which must be counted in quadrature with the intrinsic residual error
# This means that out results are as follows:

rescaleFactor <- 1/0.9008746
#write(rescaleFactor, stdout())
cosiSigma <- 0.37
cosiFactor <- 1.74
tightResSigma <- 1.287e-26*cosiFactor
looseResSigma <- 2.501e-26*cosiFactor*rescaleFactor

#write(tightResSigma, stdout())
#write(looseResSigma, stdout())

matrixResSigma <- array(tightResSigma, length(bestSingleIFO$h0), 1)
matrixResSigma[loosePulsars] <- looseResSigma
#write(matrixResSigma, stdout())

# First, the actual estimates of h0:
rescaledRandomCosiH0 <- cosiFactor*rescaleFactor*bestSingleIFO$h0  
write(rescaledRandomCosiH0, stdout())

# Then the uncertainties, in quadrature error:
uncertainH0 <- quadErr(cosiSigma*rescaledRandomCosiH0, matrixResSigma)
write(uncertainH0, stdout())

# Finally, let us see about the errors in the open data set:



MDCtableReader <- function(filename, headerFlag, seperator) {
    names <- c()
    new.dat <- c()
    x <- readLines(filename)
    xNames <- unlist(strsplit(x[1],'\t'))
    x <- x[-1]
    for (i in 1:length(x)) {
        # The file is tab seperated
        y <- unlist(strsplit(x[i], '\t'))
        # There is an error in the table with the first row, 
        # Where there is an extra space; this fixes that
        if (i == 1){
            y <- y[-6]
        }
        # Convert all except the pulsar number into a numeric data type
        # Excessively general, cannot get it to work
        freqdatum <- as.numeric(y[5])
        h0datum <- as.numeric(y[6])
        asinidatum <- as.numeric(y[9])
        Pdatum <- as.numeric(y[11])
        cosidatum <- as.numeric(y[3])
        tabledata <- rbind(freqdatum, h0datum, asinidatum, Pdatum, cosidatum)
        new.dat <- cbind(new.dat, tabledata)
        #y[1:6] <- as.numeric(y[1:6])
        #y[8:length(y)] <- as.numeric(y[8:length(y)])
        #addition <- rbind(new.dat, y)
        #new.dat <- addition
        #names <- rbind(names, xNames)
    }
    #names(new.dat) <- names
    return(new.dat)
}

tabledata <- MDCtableReader('MDCv6_open_table.dat', FALSE, "\t")
tableMatrix <- matrix(tabledata, length(tabledata)/5, 5, byrow=TRUE)


# We need to note the non-detected open pulsars and remove them from the
# table: f, h0, asini, P, cos iota

nonDetectedPulsarsOpen <- rbind(14, 16, 20, 25, 26, 28, 29, 30, 33, 34, 36, 37, 38, 41, 42, 44, 45, 46, 49)
tableOpenDetected <- tableMatrix[-nonDetectedPulsarsOpen,]
#write(tableOpenDetected,stdout())

if (testOpen == 'yes') {
    # Difference in h0
    diffH0 <- rescaledRandomCosiH0 - tableOpenDetected[,2]
    diffF <- bestSingleIFO$freq - tableOpenDetected[,1] 
    diffAsini <- bestSingleIFO$df - tableOpenDetected[,3]
    # Divided by uncertainty
    relH0 <- diffH0 / uncertainH0
    relF <- diffF / matrixFSigma
    relAsini <- diffAsini / matrixAsiniSigma
    # Test within one sigma:
    withinH0 <- 1* (abs(relH0) <= 1)
    withinF <- 1* (abs(relF) <= 1)
    withinAsini <- 1* (abs(relAsini) <= 1)
    # Fraction within one sigma
    fractH0 <- mean(withinH0)
    fractF <- mean(withinF)
    fractAsini <- mean(withinAsini)

    write(relH0, stdout())
    write(relF, stdout())
    write(relAsini, stdout())

    write(fractH0, stdout())
    write(fractF, stdout())
    write(fractAsini, stdout())
}
