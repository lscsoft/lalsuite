# Grant David Meadors
# 02014-04-28 (JD 2456776)
# g m e a d o r s @ u m i c h . e d u
# Random cos iota distribution and effect on h effective

#write("Testing effective of cos iota distribution on h effective", stdout())
#write("as recovered for given detection efficiency", stdout())
#write("motivated by Scorpius X-1 MDC", stdout())

library(MASS)

simSize <- 2000000

cosi <- runif(simSize, -1, 1)
#meanCosi <- apply(cosi, 1, mean)

png("PlotHEffDistCosi.png")
hist(cosi, probability=FALSE, breaks=100)
dev.off()


h0min <- 3.0e-26
h0max <- 3.0e-24
h0exp <- runif(simSize, 0, 1)
h0true <- h0min * (h0max/h0min)^(h0exp)


png("PlotHEffDistH0.png")
hist(log10(h0true), probability=FALSE, breaks=100)
dev.off()



heffectiveCalc <-function(hValue, cosiValue) {
    # We care using the formula for h-effective that
    # h-eff = sqrt(h_x^2 + h_+^2)
    # May need a square root of two too
    plusSquared = ((1 + cosiValue^2)/2)^2
    crossSquared = (cosiValue)^2
    return (sqrt(plusSquared+crossSquared)*hValue/sqrt(2))
}

h0eff <- heffectiveCalc(h0true, cosi)

png("PlotHEffDistH0All.png")
hist(log10(h0eff), probability=FALSE, breaks=100)
dev.off()

write("Generated this many h0-effective outliers: ",stdout())
write(length(h0eff),stdout())

#detectEff <-function(h0outlier) {
#    
#     ## Ad hoc detection efficiency simulator
#     highThres <- 3.0e-25
#     lowThres <- 1.0e-25
#     allOutliers <- h0outlier[h0outlier >= lowThres]
#     effAllOutliers <- allOutliers[(allOutliers >= runif(length(allOutliers), lowThres, highThres))]
#     return (effAllOutliers)
#} # function detectEff


# Better, we can recover the indices that survive
# so we can use them later
detectEffIndex <-function(h0outlier) {
    
     ## Ad hoc detection efficiency simulator
     highThres <- 3.0e-25
     lowThres <- 1.0e-25
     effAllOutliers <- (h0outlier >= runif(length(h0outlier), lowThres, highThres))
     return (effAllOutliers)
} # function detectEff

detectedIndices <- detectEffIndex(h0eff)
h0trueDetected <- h0true[detectedIndices]
h0detected <- h0eff[detectedIndices]


#h0detected <- h0eff[detectEffIndex(h0eff)]

write("This many candidates survive detection efficiency: ", stdout())
write(length(h0detected), stdout())

png("PlotHEffDistH0Detected.png")
hist(log10(h0detected), probability=FALSE, breaks=100)
dev.off()

breaks <- c(2e-26,1e-25, 1.3e-25,1.5e-25,2.5e-25,3e-25,9e-25)
alternativeBreaks <- (1:200)*3e-27

graphDetectionEff <- function(h0eff, h0detected, breaks){
    #sortedHeff <- h0eff[order(h0eff)]
    #detectedList
    detectionOrNot <- 1* (h0eff %in% h0detected)


    binGroups <- cut(h0eff, breaks, labels=breaks[-1])


    binnedDetection <- tapply(detectionOrNot, binGroups, mean)
    nDetection <- tapply(detectionOrNot, binGroups, length)
    uncertaintyDetection <- sqrt((binnedDetection)*(1-binnedDetection)/nDetection)
    png(paste("PlotHeffDistH0DetectionEfficiency", length(breaks), "breaks.png", sep=""))
    plot(breaks[-1], binnedDetection-uncertaintyDetection ,type="s",xlab='highest h-effective', ylab='Fraction detected', main=paste('(Detection efficiency - uncertainty) for Scorpius X-1 MDC mini'))
    grid()
dev.off()
} # function graphDetectionEff
graphDetectionEff(h0eff, h0detected, breaks)
graphDetectionEff(h0eff, h0detected, alternativeBreaks)

# Finally, let us graph the relationship between the h-effective that
# are detected and the true h0

#library(RColorBrewer)
#my.cols <- rev(brewer.pal(k,"RdY1Bu"))

densityOfDetections <- kde2d(h0trueDetected, h0detected, n=c(125,125))
#freqOfDetections <- range(h0trueDetected) * range(h0detected)
#write(freqOfDetections, stdout())

# The mean of injections will return a vector that equals the injection h0
# corresponding to each row, in ascending order
meansOfInjections <- densityOfDetections$x 
SumsOfRows <- apply(densityOfDetections$z, 2, sum)
cumsumOfRows <- apply(densityOfDetections$z, 2, cumsum)
# The mean of detections will return a vector that equals the detected h0
# corresponding to each column, in ascending order
meansOfDetections <- densityOfDetections$y
SumsOfColumns <- apply(densityOfDetections$z, 1, sum)
cumsumOfColumns <- apply(densityOfDetections$z, 1, cumsum)

#write('Bin of each h0 injected: ',stdout())
#write(meansOfInjections,stdout())

#write('Bin of each h0 detected: ',stdout())
#write(meansOfDetections,stdout())

#write('averages of rows and columns: ', stdout())
#write(apply(densityOfDetections$z, 1, mean), stdout())
#write(apply(densityOfDetections$z, 2, mean), stdout())

#write('Density parameters: ', stdout())
#write(SumsOfRows[25], stdout())
#write(cumsumOfRows[,25], stdout())
#write(SumsOfRows, stdout())
#write(SumsOfColumns, stdout())
#write(meanOfRows, stdout())

# We need to restrict the range of detected h0 that we examine,
# because going higher than that means we will fail to see
# the necessary and complete range of true h0 values that
# correspond to the detected h0
#write('Limiting range of h0 detected search...', stdout())
limitingH0detected <- heffectiveCalc(3e-24, 0)
#write(limitingH0detected, stdout())
limited0DetectedBin <- meansOfDetections[meansOfDetections < limitingH0detected]
# Empirically, there is also a great deal of noise at the lowest
# H0 detected. We should set a minimum something like our
# upper limit noise floor, to avoid bias
minH0detected <- 0e-25
limitedDetectedBin <- limited0DetectedBin[limited0DetectedBin > minH0detected]
#write(limitedDetectedBin,stdout())


# Let us set our thresholds to plus and minus one sigma
plusSigma <- 0.8414
minusSigma <- 1-0.8414

meanOfRows <- c()
meanXofRows <- c()
PSofRows <-c()
MSofRows <-c()

# Now it gets really confusing. The apply function works
# ALONG a given dimension, meaning when given 2, it yields
# the sum FOR each row (summing over columns) 
# and likewhise when given 1, it yields
# the sum FOR each column (summing over rows)
# Recall again that each injection h0 corresponds to a
# column, and
# each detected h0 corresponds to a 
# row

# Now, let us compare the cumulative 
# to the sum. 
# We need to compare each cumsumOfRows[,ii] to
# each SumOfRows[ii]
# The index needs to range for as many bins as there are
for (ii in 1:length(limitedDetectedBin)) {
    
    #write(cumsumOfRows[,ii] <= 0.5*SumsOfRows[ii], stdout())
    meanIndexRow <- cumsumOfRows[,ii] <= 0.5*SumsOfRows[ii]
    #meanOfRows <- cbind(meanOfRows, max(cumsumOfRows[meanIndexRow,ii]))
    meanXofRows <- cbind(meanXofRows, max(densityOfDetections$x[meanIndexRow]))
    PSIndexRow <- cumsumOfRows[,ii] <= plusSigma*SumsOfRows[ii]
    PSofRows <- cbind(PSofRows, max(densityOfDetections$x[PSIndexRow]))
    MSIndexRow <- cumsumOfRows[,ii] <= minusSigma*SumsOfRows[ii]
    MSofRows <- cbind(MSofRows, max(densityOfDetections$x[MSIndexRow]))
}

sigmaDF <- data.frame(detectionBin = limitedDetectedBin, meanData = as.numeric(meanXofRows), PSdata = as.numeric(PSofRows), MSdata = as.numeric(MSofRows))

write('Minus one sigma, mean, and plus one sigma per row: ', stdout())
write(MSofRows, stdout())
write(meanXofRows, stdout())
write(PSofRows, stdout())

png("PlotHEffVsH0True.png")
plot(h0trueDetected, h0detected, pch=20, cex=0.001, xlab = 'true injected h0', ylab='effective recovered, detected h0', main='Distribution of effective h0 for given true injected h0')
dev.off()


# And the most important plot of all, the plot with 1-sigma lines
#write(sigmaDF$meanData,stdout())
#write(sigmaDF$detectionBin,stdout())

# Here, we also try to fit regressions to the means and sigmas
regressMean <- lm(sigmaDF$meanData ~ sigmaDF$detectionBin)
regressMS <- lm(sigmaDF$MSdata ~ sigmaDF$detectionBin)
regressPS <- lm(sigmaDF$PSdata ~ sigmaDF$detectionBin)
regressionTitle <- paste('Regression slope: plus, mean, minus 1 sigma interval \n', ' ', round(coef(regressPS)[2], digits=3), ', ', round(coef(regressMean)[2], digits=3),', ', round(coef(regressMS)[2], digits=3), sep="")

png("PlotHeffVsH0TrueRegressions.png")
plot(sigmaDF$detectionBin, sigmaDF$meanData, pch=19, col="blue", xlab = 'Detected h0 bin', ylab='Injected h0 required to reach statistic', main=regressionTitle)
points(sigmaDF$detectionBin, sigmaDF$PSdata, col="red")
points(sigmaDF$detectionBin, sigmaDF$MSdata, col="green")
#contour(densityOfDetections, drawlabels=FALSE, nlevels=100, add=TRUE)
#densityOfDetections$z <- log(densityOfDetections$z)
abline(0,1)

abline(regressMean, col="blue")
abline(regressPS, col="red")
abline(regressMS, col="green")

legend(8e-25, 8e-25, legend=c('Means', '+sigma', '-sigma', 'h-eff = h-inj', 'fit to mean', 'fit to +sigma', 'fit to -sigma'), pch = c(1,1,1,NA,NA,NA,NA), lty=c(NA,NA,NA,1,1,1,1), col=c("blue","red","green","black","blue","red","green"))
dev.off()

#write('Writing coefficients of regressions for -sigma, mean, +sigma', stdout())
#write(coef(regressMS),stdout())
#write(coef(regressMean),stdout())
#write(coef(regressPS),stdout())


invertRegCoef <- function(someRegression){
    # Linear model: y = a + b x
    # invert: x  = (-a/b) + (1/b) y
    aa <- as.numeric(coef(someRegression)[1])
    bb <- as.numeric(coef(someRegression)[2])
    return( c( (-aa/bb), ( 1 / bb )) )
    #return('Hello')
} # function invertRegCoef

#write('Writing inverse coefficients of regressions for -sigma, mean, +sigma', stdout())
#write(invertRegCoef(regressMS),stdout())
#write(invertRegCoef(regressMean),stdout())
#write(invertRegCoef(regressPS),stdout())

png("PlotHEffVsH0TrueWithLines.png")
plot(h0trueDetected, h0detected, pch=20, cex=0.001, xlab = 'true injected h0', ylab='effective recovered, detected h0', main='Distribution of effective h0 for given true injected h0')
abline(0,1)
abline(invertRegCoef(regressMean), col="blue")
abline(invertRegCoef(regressPS), col="red")
abline(invertRegCoef(regressMS), col="green")
legend(2.0e-24, 7e-25, legend=c('h-eff = h-inj', 'fit to mean', 'fit to +sigma', 'fit to -sigma'), lty=c(1,1,1,1), col=c("black","blue","red","green"))
dev.off()

#Keith recommended the following sanity check:
sigmaInside <- 1
sigmaEst <- 0.37
DiffInEst <- abs(1.74*h0detected - h0trueDetected)
trueSigmaDef <- sd(DiffInEst[h0detected < heffectiveCalc(3e-24,0)]/h0trueDetected[h0detected < heffectiveCalc(3e-24,0)])
sigmaInEst <- 1* (DiffInEst/(sigmaEst*h0trueDetected) <= sigmaInside)
sigma1InEst <- 1* (DiffInEst/(sigmaEst*h0trueDetected) <= 1.5*sigmaInside)
sigma2InEst <- 1* (DiffInEst/(sigmaEst*h0trueDetected) <= 2*sigmaInside)
sigmaTrueInEst <- 1* (DiffInEst/(h0trueDetected) <= trueSigmaDef)
sigma2TrueInEst <- 1* (DiffInEst/(h0trueDetected) <= 2*trueSigmaDef)
write(mean(sigmaInEst), stdout())

alternative2Breaks <- c(1:100)*2.5e-26
binGroups <- cut(h0detected, alternative2Breaks, labels=alternative2Breaks[-1])
sigmaInEstMean <- tapply(sigmaInEst, binGroups, mean)
sigma1InEstMean <- tapply(sigma1InEst, binGroups, mean)
sigma2InEstMean <- tapply(sigma2InEst, binGroups, mean)
sigmaTrueInEstMean <- tapply(sigmaTrueInEst, binGroups, mean)
sigma2TrueInEstMean <- tapply(sigma2TrueInEst, binGroups, mean)

png("PlotSigmaDiffVsH0Eff.png")
plot(alternative2Breaks[-1], sigmaInEstMean, main=paste("Green: 1 sigma, Blue: 1.5 sigma, Red: 2 sigma", "\n sigma = ", sigmaEst), ylim=c(0,1), type="s", xlab="Effective, detected h0", ylab="Mean fraction within sigma", col="green")
lines(alternative2Breaks[-1], sigma1InEstMean, col="blue")
lines(alternative2Breaks[-1], sigma2InEstMean, col="red")
grid()
dev.off()

png("PlotSigmaTrueDiffVsH0Eff.png")
plot(alternative2Breaks[-1], sigmaTrueInEstMean, main=paste("Green: 1 sigma, Red: 2 sigma, sigma =", trueSigmaDef), ylim=c(0,1), type="s", xlab="Effective, detected h0", ylab="Mean fraction within sigma", col="green")
lines(alternative2Breaks[-1], sigma2TrueInEstMean, col="red")
grid()
dev.off()





#write("Evaluating contour plot: ", stdout())
#write(meansOfDetections, stdout())
png("PlotHeffVsH0FilledContour.png")
filled.contour(densityOfDetections, nlevels=100)

#abline(invertRegCoef(regressMean), col="blue")
#abline(invertRegCoef(regressPS), col="red")
#abline(invertRegCoef(regressMS), col="green")

dev.off()
