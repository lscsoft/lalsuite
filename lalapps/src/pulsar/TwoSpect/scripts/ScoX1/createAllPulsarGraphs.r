# Grant David Meadors
# g m e a d o r s @  u m i c h . e d u
# 02014-04-16 (JD 2456764)
# Run the Python graphing script that produces graphs
# for the MDCv6 Sco X-1 pulsars analyzed using TwoSpect
# Generate statistics to characterize detection efficiency

#useMode <- 'aggregate'

useMode <- 'histogram'

#useModeHistogram <- 'collectVerboseWhole'
useModeHistogram <- 'simple'

#useModeRead <- 'wholeSummary'
useModeRead <- 'savedRda'
#useModeRead <- 'testTables'

#useModeCut <- 'performCut'
#useModeCut <- 'noCut'
useModeCut <- 'bypassCut'

useModeProbability <- 'performSort'
#useModeProbability <- 'noSort'

useModeDivide <- 'performDivide'
#useModeDivide <- 'noDivide'

#useModeBlinded <- 'open'
#useModeBlinded <- 'closed'
useModeBlinded <- 'both'

useModeDetectionSelection <- 'detected'
#useModeDetectionSelection <- 'nondetected'
#useModeDetectionSelection <- 'all'

useModeParamDiff <- 'params'
#useModeParamDiff <- 'ULs'
#useModeParamDiff <- 'noDiff'


#useModePulsar <- 'J1751'
useModePulsar <- 'ScoX1'

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
MDCmasterTableReader <- function(filename, headerFlag, seperator) {
    new.dat <- c()
    x <- readLines(filename)
    #write('Experimenting with master table reader',stdout())
    freqdata <- x[grep("Freq",x)]
    h0data <- x[grep("h0",x)]
    asinidata <- x[grep("orbitasini",x)]
    Pdata <- x[grep("orbitPeriod",x)]
    cosidata <- x[grep("cosi",x)]
    numberOfLines <- length(h0data)
    for (i in 1:numberOfLines) {
        freqdatum <- miniTableExtractor(i,freqdata)
        h0datum <- miniTableExtractor(i, h0data)
        asinidatum <- miniTableExtractor(i, asinidata)
        Pdatum <- miniTableExtractor(i, Pdata)
        cosidatum <- miniTableExtractor(i, cosidata)
        tabledata <- rbind(freqdatum, h0datum, asinidatum, Pdatum, cosidatum)
        new.dat <- cbind(new.dat, tabledata)
    }
    #write(new.dat, stdout())
    return(new.dat)
}

miniTableExtractor <- function(i, element){
    charElement <- as.character(element[i])
    charDatum <- as.numeric(substr(charElement, regexpr("= ", charElement)+1, nchar(charElement)))
    return(charDatum)
}

if (useMode == 'aggregate') {
    location <- '/home/gmeadors/ScoX1/2014/04/15-extract-results/12-open-H1/Meadors-ScoX1/'
    write(location, stdout())
    # For a noise test, which I think was on 02013-09-04
    #arguments <- paste(' --elsewhere ', location, '--noiseTest')
    arguments <- paste(' --elsewhere ', location, '--multiTemplateSearch 50')

    mdcv <- '6'
    
    # Full interferometers
    ifoList <- c('H1', 'L1', 'V1')
    # Partial interferometers
    #ifoList <- c('H1')
    if (useModeBlinded == 'closed') {
        # closed pulsars below
        pulsarList <- c('001', '002', '003', '005', '011', '014', '015', '017', '019', '020', '021', '023', '026', '029', '032', '035', '036', '041', '044', '047', '048', '050', '051', '052', '054', '057', '058', '059', '060', '061', '062', '063', '064', '065', '066', '067', '068', '069', '071', '072', '073', '075', '076', '079', '083', '084', '085', '090', '095', '098')
    }
    if (useModeBlinded == 'open') {
        # open pulsars below
        pulsarList <- c('004', '006', '007', '008', '009', '010', '012', '013', '016', '018', '022', '024', '025', '027', '028', '030', '031', '033', '034', '037', '038', '039', '040', '042', '043', '045', '046', '049', '053', '055', '056', '070', '074', '077', '078', '080', '081', '082', '086', '087', '088', '089', '091', '092', '093', '094', '096', '097', '099', '100')
    }
    if (useModeBlinded == 'both') {
        # This mode has not been tested, because to date I have 
        # processed open and closed results in seperate directories
        pulsarList <- c(sprintf("%03d", 1:100))
    }
    # starting frequencies of the 5 Hz bands
    for(i in pulsarList){    
        cat(paste(i), '\n')
        for(j in ifoList){
            system(paste('./createOutputSummary.py ', mdcv, ' ', j, ' ', i, arguments))
        }
    }
    cat(paste('Number of pulsars graphed:', str(length(pulsarList)),' \n'))
} # useMode aggregate
if (useMode == 'histogram') {
    if (useModeBlinded == 'both'){
        tabledata <- MDCmasterTableReader('MDCv6_master.dat', FALSE, "\t")
    } else {
        tabledata <- MDCtableReader('MDCv6_open_table.dat', FALSE, "\t")
    }
    tableMatrix <- matrix(tabledata, length(tabledata)/5, 5, byrow=TRUE)
    dfGuess <- 2 * 3.14159265359 * tableMatrix[,1] * tableMatrix[,3] / tableMatrix[,4]
    asiniCalc <- function(dfValue, freqValue) {
        # Use a hardcoded period of 68023.8259 s, as was assumed for search
        hardcodedP = 68023.8259
        return (dfValue * hardcodedP / (2 * 3.1415926535 * freqValue))
    }
    heffectiveCalc <-function(hValue, cosiValue) {
        # We care using the formula for h-effective that
        # h-eff = sqrt(h_x^2 + h_+^2)
        # May need a square root of two too
        plusSquared = ((1 + cosiValue^2)/2)^2 
        crossSquared = (cosiValue)^2
        return (sqrt(plusSquared+crossSquared)*hValue/sqrt(2))    
    }
    # Print a matrix where the first column is the lower edge and the
    # second is the upper edge of the limits that must be cut
    # due to frequency modulations and Earth doppler shift
    # Buffer the limits by a factor of two to be safe and by ten bins
    dfLimits <- cbind(tableMatrix[,1] *(1-2e-4) - 2*dfGuess -10/360, tableMatrix[,1] *(1+2e-4) + 2*dfGuess+10/360)
    # Optionally, write out wat the limits are that will be cut because
    # they contain open injections
    #write("dfLimits", stdout())
    #write(dfLimits[, 1], stdout())
    #write(dfLimits[, 2], stdout())
    # The table data, for confirmation:
    #write(tableMatrix[,1], stdout())
    #write("true h", stdout())
    #write(tableMatrix[,2], stdout())
    #write("h-effective", stdout())
    #write(heffectiveCalc(tableMatrix[,2], tableMatrix[,5]), stdout())
    #write("true asini", stdout())
    #write(tableMatrix[,3], stdout())
    #write("calculated asini", stdout())
    #write(asiniCalc(dfGuess, tableMatrix[,1]),stdout())
    #write(tableMatrix[,4], stdout())
    #write(tableMatrix[,5], stdout())
    #write(dfGuess, stdout())


    extractTwoSpectOutput <- function (j) {
      if (useModeHistogram == 'collectVerboseWhole'){
          system(paste('cat verbose_summary*', j, '* > whole_summary_', j, '.txt', sep=""))
      }
      if (useModeRead == 'wholeSummary') {
          dataSet <- read.table(paste('whole_summary_', j, '.txt',sep=""))
          numericize <- function (data, s) {
              return(as.numeric(gsub(',','',data[, c(s)],fixed=TRUE)))
          }
          # Results will have the first colum be name, the second frequency,
          # the third modulation depth, fourth R, fifth h0, sixth log10p
          # Note that period is consistently 68023.8259 s for all tests
          results <- data.frame(name = j, freq = numericize(dataSet, 3), df=numericize(dataSet, 9), R = numericize(dataSet, 18), h0=numericize(dataSet, 21), Prob = numericize(dataSet, 24))
      }
      if (useModeRead == 'savedRda') {
        if (useModeBlinded == 'open') {
            # open results (safe backup in /home/gmeadors/ScoX1/2014/04/16-R-backups/open/)
            locationRda <- "/home/gmeadors/ScoX1/2014/04/21/Meadors-ScoX1/"
        }
        if (useModeBlinded == 'closed') {
            # closed results (safe backup in /home/gmeadors/ScoX1/2014/04/16-R-backups/closed/)
            locationRda <- "/home/gmeadors/ScoX1/2014/04/22-closed/Meadors-ScoX1/"
        }
        if (useModeCut == 'performCut' | useModeCut == 'bypassCut') {
            load(paste(locationRda,"results",j,".Rda",sep=""))
        }
        if (useModeCut == 'noCut') {
            load(paste(locationRda,"resultsCut",j,".Rda",sep=""))
        }
      }
      if (useModeRead == 'testTables') {
          write('Testing MDC table reading functionality', stdout())
          results <- data.frame(name = character(), freq = as.numeric(character()), df=as.numeric(character()), R = as.numeric(character()), h0=as.numeric(character()), Prob = as.numeric(character()))
      }
      write(str(results),stdout())
      # Make cut based on proximity to injections
      #cutIndices <- c()
      #for (i in 1:length(results$freq)) { 
      #       betweenLimits <- FALSE
      #   for (j in 1:length(dfLimits[, 1])) {
      #       if (results$freq[i] >= dfLimits[j,1] && results$freq[i] <= dfLimits[j,2]) {
      #           betweenLimits <- TRUE
      #           cutIndices <- cbind(cutIndices, i)
      #       } 
      #   } 
      #}
      #results <- results[-cutIndices, ]
      # Alternate approach using built-in R functions
      write("Length of frequency vector before cut",stdout())
      write(length(results$freq),stdout())
      #write("Head of frequency vector: ", stdout())
      #write(results$freq[1:10],stdout())
      FreqCutResults <- results$freq
      cIndicesList <- c()
      if (useModeCut == 'performCut') {
          for (j in 1:length(dfLimits[, 1])) {
              newCutIndices <- c()
              newCutIndices <- which((FreqCutResults >= dfLimits[j,1]) & (FreqCutResults <= dfLimits[j,2]))
              cIndicesList <- c(cIndicesList, newCutIndices)
          }
      }
      if (useModeRead == 'wholeSummary') {
          save(results,file=paste("results",results$name,".Rda",sep=""))
      }
      # Cut occurs here, and will wipe out all results if
      # cIndicesList is empty
      if (useModeCut == 'performCut') { 
          results <- results[-cIndicesList, ]
          write("Length of cut list", stdout())
          write(length(cIndicesList),stdout())
          write("Length of frequency vector after cut", stdout())
          write(length(results$freq),stdout())
      }
      if (useModeDetectionSelection == 'all') {
          write("Probability minimum", stdout())
          write(min(results$Prob), stdout())
      }
      if (useModeProbability == 'performSort') {
          write("Sorting extreme p-values",stdout())
          # Old way of sorting, only yielded probability
          #write.table(sort(results$Prob)[1:2000], paste('ProbabilityMinima',results$name,'.txt',sep=""), sep="\n")
          # New way of sorting
          #write.table(results[order(results$Prob)[1:2000], ], paste('ProbabilityMinima',results$name,'.txt',sep=""), sep="\n")
          if (useModeDivide == 'noDivide') {
          # How to sort the top 2000 from the whole open data set
              exportResults <- results[order(results$Prob)[1:2000], ]
          }
          if (useModeDivide == 'performDivide') {
              startBandListClosed <- c(50.0, 60.0, 70.0, 90.0, 150.0, 180.0, 190.0, 210.0, 230.0, 240.0, 250.0, 270.0, 300.0, 330.0, 360.0, 390.0, 400.0, 450.0, 480.0, 510.0, 520.0, 540.0, 550.0, 560.0, 590.0, 620.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, 710.0, 730.0, 740.0, 750.0, 760.0, 800.0, 810.0, 820.0, 860.0, 880.0, 930.0, 1080.0, 1100.0, 1110.0, 1190.0, 1320.0, 1370.0)
              startBandListOpen <- c(80.0, 100.0, 110.0, 120.0, 130.0, 140.0, 160.0, 170.0, 200.0, 220.0, 260.0, 280.0, 290.0, 310.0, 320.0, 340.0, 350.0, 370.0, 380.0, 410.0, 420.0, 430.0, 440.0, 460.0, 470.0, 490.0, 500.0, 530.0, 580.0, 600.0, 610.0, 770.0, 830.0, 890.0, 920.0, 960.0, 1040.0, 1050.0, 1120.0, 1160.0, 1170.0, 1180.0, 1250.0, 1280.0, 1290.0, 1300.0, 1330.0, 1360.0, 1390.0, 1450.0)
              # ONLY AFTER HAVING DONE detection calls,
              # we can define this list to select which pulsars were seen
              # and which were not
              # This list corresponds to the position in the lists, not
              # to the "pulsar number"
              wholePulsarIndexList <- c(1:50)
              nonDetectedListOpen <- c(14, 16, 20, 25, 26, 28, 29, 30, 33, 34, 36, 37, 38, 41, 42, 44, 45, 46, 49)
              nonDetectedListClosed <- c(11, 18, 21, 22, 24, 25, 26, 27, 32, 33, 38, 39, 40, 41, 48, 50) 
              if (useModeBlinded == 'closed') {
                  startBandListChosen <- startBandListClosed
                  if (useModeDetectionSelection == 'detected') {
                      startBandListChosen <- startBandListChosen[-nonDetectedListClosed]
                  }
                  if (useModeDetectionSelection == 'nondetected') {
                      startBandListChosen <- startBandListChosen[nonDetectedListClosed]
                  }
              }
              if (useModeBlinded == 'open') {
                  startBandListChosen <- startBandListOpen
                  if (useModeDetectionSelection == 'detected') {
                      startBandListChosen <- startBandListChosen[-nonDetectedListOpen]
                  }
                  if (useModeDetectionSelection == 'nondetected') {
                      startBandListChosen <- startBandListChosen[nonDetectedListOpen]
                  }
              }
              write('How many pulsars will be analyzed:', stdout())
              write(length(startBandListChosen), stdout())
              exportResults <- data.frame(name = character(), freq = as.numeric(character()), df=as.numeric(character()), R = as.numeric(character()), h0=as.numeric(character()), Prob = as.numeric(character()))
              for (k in startBandListChosen) {
                  if (useModeDetectionSelection == 'all') {
                      singleIFOprobThreshold <- -7.75
                      if (k >= 360.0) {
                          singleIFOprobThreshold <- -12.0
                      }
                      subBand <- results[which((results$freq >= k) & (results$freq <= k+5) & (results$Prob <= singleIFOprobThreshold )), ]
                      #subBand <- results[which((results$freq >= 1300.0) & (results$freq <= 1305.0) & (results$Prob <= -14.0 )), ]
                      ## How to sort the top few from given bands
                      #exportResults <- subBand[order(subBand$Prob), ]
                      #subBand <- results[which((results$freq >= 1360.0) & (results$freq <= 1365.0) & (results$Prob <= -14.0 )), ]
                      #exportResults <- rbind(exportResults, subBand[order(subBand$Prob), ])
                      if (length(subBand$Prob) >= 1) {
                          exportResults <- rbind(exportResults, subBand[order(subBand$Prob)[1:min(200,length(subBand$Prob))],] )
                      }
                  }
                  # For detected pulsars, we are interested in the parameters at the peak,
                  # and for non-detected pulsars, in their upper-limit values
                  if ((useModeDetectionSelection == 'detected') | (useModeDetectionSelection == 'nondetected')) {
                      write('Value of starting frequency in band', stdout())    
                      write(k, stdout())
                      if (useModeBlinded == 'open') {
                          openNo <- which(5*floor(as.integer(tableMatrix[,1])/5) == as.integer(k))
                          #write(openNo, stdout())
                          #write('Parameters for the band', stdout())
                          #write(tableMatrix[])
                          # The table data, for confirmation:
                          # write('frequency', stdout())
                          #write(tableMatrix[openNo,1], stdout())
                          #write("true h", stdout())
                          #write(tableMatrix[openNo,2], stdout())
                          #write("h-effective", stdout())
                          #write(heffectiveCalc(tableMatrix[openNo,2], tableMatrix[openNo,5]), stdout())
                          #write("true asini", stdout())
                          #write(tableMatrix[openNo,3], stdout())
                          #write("calculated dfGuess asini", stdout())
                          #write(asiniCalc(dfGuess[openNo], tableMatrix[openNo,1]),stdout())
                          #write('period', stdout())
                          #write(tableMatrix[openNo,4], stdout())
                          #write('cos i', stdout())
                          #write(tableMatrix[openNo,5], stdout())
                          #write('modulation depth', stdout())
                          #write(dfGuess[openNo], stdout())
                      }

                      # Now, take a look at the actual open data
                      #write('Probability minimum in actual data', stdout())
                      subBand <- results[which((results$freq >= k) & (results$freq <= k+5)  ), ]
                      #write(min(subBand$Prob), stdout())
                      peakP <- which(subBand$Prob == min(subBand$Prob))
                      #write('Index, p-value, f, df, h0', stdout())
                      #write(peakP, stdout())
                      #write(subBand$Prob[peakP], stdout())
                      #write(subBand$freq[peakP], stdout())
                      #write(subBand$df[peakP], stdout())
                      #write(subBand$h0[peakP], stdout())
                          
                      # Differences in recovered-at-peak-p and true
                      ProbPeakP <- subBand$Prob[peakP]
                      RPeakP <- subBand$R[peakP] 
                      # The subtraction of true signal assumes open data
                      if (useModeParamDiff == 'params') {
                          if (useModeBlinded == 'closed') {
                              write('Attempted to train parameter estimation on closed data, will not work', stdout())
                          }
                          if (useModeBlinded == 'open') {
                              # Units of Hz difference
                              peakFreq <- subBand$freq[peakP] - tableMatrix[openNo, 1]
                              # Units of asini
                              peakAsini <- asiniCalc(subBand$df[peakP], subBand$freq[peakP]) - tableMatrix[openNo, 3]
                              # Units of h_effective
                              peakh0 <- subBand$h0[peakP] - heffectiveCalc(tableMatrix[openNo,2], tableMatrix[openNo,5])
                          }
                      }
                      if ( (useModeParamDiff == 'ULs') | (useModeParamDiff == 'noDiff')) {
                          # Units of Hz
                          peakFreq <- subBand$freq[peakP]
                          # Units of asini
                          peakAsini <- asiniCalc(subBand$df[peakP], subBand$freq[peakP])
                          # Units of h_effective
                          peakh0 <- subBand$h0[peakP]
                      }
                      write ('Values at peak in actual data: f, asini, h0', stdout())
                      write(peakFreq, stdout())
                      write(peakAsini, stdout())
                      write(peakh0, stdout())
                      exportResults <- rbind(exportResults, data.frame(name = as.character(k), freq = as.numeric(peakFreq), df = as.numeric(peakAsini), R = as.numeric(RPeakP), h0 = as.numeric(peakh0), Prob = as.numeric(ProbPeakP)))

                  }
              }
          }
          exportResults$period <- rep(68023.8259, nrow(exportResults)) 
          exportResults$tfnorm <- rep(0.1, nrow(exportResults)) 
          exportResults$jobnumber <- rep(0, nrow(exportResults)) 
          exportResults$ra <- rep(4.2757, nrow(exportResults)) 
          exportResults$dec <- rep(-0.2730, nrow(exportResults)) 
          exportResults <- exportResults[ , c("freq", "period", "df", "ra", "dec", "R", "h0", "Prob", "tfnorm", "jobnumber")]
          write("Exporting extreme p-values",stdout())
          if (useModeDetectionSelection == 'all') {
              exportHeadName <- 'ProbabilityMinima'
          }
          if (useModeDetectionSelection == 'detected') {
              exportHeadName <- 'DetectedParameters'
          }
          if (useModeDetectionSelection == 'nondetected') {
              exportHeadName <- 'NonDetectedULs'
          }
          if (useModeParamDiff == 'params' ) {
              exportHeadName <- paste(exportHeadName, 'Diff', sep = "")
          }
          if (useModeParamDiff == 'ULs') {
              exportHeadName <- paste(exportHeadName, 'CalUL', sep = "")
          }
          if ( (useModeBlinded == 'open') & (useModeDetectionSelection != 'all')) {
              exportHeadName <- paste(exportHeadName, 'Open', sep = "")
          }
          if ( (useModeBlinded == 'closed') & (useModeDetectionSelection != 'all')) {
              exportHeadName <- paste(exportHeadName, 'Closed', sep = "")
          }
          write.table(exportResults, paste(exportHeadName,results$name,'.txt',sep=""), sep=" ", row.names = FALSE, col.names = FALSE)
          save(exportResults,file=paste(exportHeadName,results$name,".Rda",sep=""))
      }
      if (useModeRead == 'wholeSummary') {
        save(results,file=paste("resultsCut",results$name,".Rda",sep=""))
      }
      if (useModeDetectionSelection == 'all') {
          return(results)
      }
      if (useModeDetectionSelection != 'all') {
          return(exportResults)
      }
    } # extractTwoSpectOutput
    # Use to read output from TwoSpect in the Scorpius X-1
    # MDC, per interferometer, and generate outliers from the
    # searched parameter space based on thresholds
    #Hanford = extractTwoSpectOutput('H1')
    #Livingston = extractTwoSpectOutput('L1')
    #Virgo = extractTwoSpectOutput('V1')
 
    # Make plots
    subhistogrammist <- function(site, quantity, breakNo, mainString, funName) {
        if (quantity == 'Prob') {
            dataHist <- site$Prob
        }
        if (quantity == 'h0') {
            dataHist <- site$h0
        }
        if (quantity == 'R') {
            dataHist <- site$R
        }
        # Use the funName to call a function, such as png or png, to save the file
        get(funName)(paste('StatHist',quantity,site$name[1],'.',funName,sep=""))
        # make the histogram
        if (quantity == 'Prob') {
            dataHist <- exp(log(10) * dataHist)
            histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        }
        if (quantity == 'h0') {
            histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        }
        if (quantity == 'R') {
            histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        }
        # send the histogram to file
        dev.off()
    } # function subhistogrammist
    histogrammist <- function(site) {
        breakNo <- 50
        ProbString <- 'Histogram frequency of exp(log(10) * (log10 probability))'
        subhistogrammist(site, 'Prob', breakNo, ProbString, 'pdf')
        subhistogrammist(site, 'Prob', breakNo, ProbString, 'png')
        h0String <- 'Histogram frequency of recovered h0'
        subhistogrammist(site, 'h0', breakNo, h0String, 'pdf')
        subhistogrammist(site, 'h0', breakNo, h0String, 'png')
        RString <- 'Histogram frequency of R statistic'
        subhistogrammist(site, 'R', breakNo, RString, 'pdf')
        subhistogrammist(site, 'R', breakNo, RString, 'png')
    } # function histogrammist
    # For just determining the cut in p, we do not necessarily need to 
    # see these histograms all the time.
    #histogrammist(Hanford)
    #histogrammist(Livingston)
    #histogrammist(Virgo)

    paramGrapher <- function(openness,selectSeen,paramEst) {
        # This function presupposes that all three observatories have been run
        # through extractTwoSpectOutput
        # with the useMode options
        # bypassCut, performSort, performDivide, open, detected, params
        # This function is used for extracting parameter uncertainties
        # For actual detection assessement, one can use either open or closed,
        # but params must be replaced with noDiff 
        write("Calculating uncertainty in parameter estimation", stdout())
        # For the LIGO Hanford cluster:
        #locationHead <- "/home/gmeadors/ScoX1/2014/04/23/Meadors-ScoX1/"
        # For Gallatin at University of Michigan:
        locationHead <- "/home/gmeadors/Documents/ScoX1/2014/04/28-POB-backup/backup_peaks-of-bands/"
        if (openness == 'open') {
            opennessString <- 'Open'
        }
        if (openness == 'closed') {
            opennessString <- 'Closed'
        }
        if (paramEst == 'params') {
            paramString <- 'Diff'
            if (openness == 'both') {
                paramString <- ''
            }
        }
        if (paramEst == 'noDiff') {
            paramString <- ''
        }
        if (openness == 'both') {
            # Loading in both open and closed data
            # First, load open
            opennessString <- 'Open'
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'H1','.Rda',sep=""))
            H1dataO <- get(loadedData)
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'L1','.Rda',sep=""))
            L1dataO <- get(loadedData)
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'V1','.Rda',sep=""))
            V1dataO <- get(loadedData)
            rm(loadedData)

            # Then, load closed
            opennessString <- 'Closed'
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'H1','.Rda',sep=""))
            H1dataC <- get(loadedData)
            H1dataB <- rbind(H1dataO, H1dataC)
            H1data <- H1dataB[order(H1dataB$f),]
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'L1','.Rda',sep=""))
            L1dataC <- get(loadedData)
            L1dataB <- rbind(L1dataO, L1dataC)
            L1data <- L1dataB[order(L1dataB$f),]
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'V1','.Rda',sep=""))
            V1dataC <- get(loadedData)
            V1dataB <- rbind(V1dataO, V1dataC)
            V1data <- V1dataB[order(V1dataB$f),]
            rm(loadedData)
            rm(H1dataO)
            rm(H1dataC)
            rm(H1dataB)
            rm(L1dataO)
            rm(L1dataC)
            rm(L1dataB)
            rm(V1dataO)
            rm(V1dataC)
            rm(V1dataB)
            opennessString <- 'Both'
            #write('ALL frequencies', stdout())
            #write(H1data$f,stdout())
            #write(L1data$f,stdout())
            #write(V1data$f,stdout())
            #write(H1data$Prob,stdout())
            #write(L1data$Prob,stdout())
            #write(V1data$Prob,stdout())
            constructMaxProb <- c(1:length(H1data$Prob)) 
        }
        else {
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'H1','.Rda',sep=""))
            H1data <- get(loadedData)
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'L1','.Rda',sep=""))
            L1data <- get(loadedData)
            rm(loadedData)
            loadedData <- load(paste(locationHead,'DetectedParameters',paramString,opennessString,'V1','.Rda',sep=""))
            V1data <- get(loadedData)
            rm(loadedData)
            #write(H1data$Prob,stdout())
            #write(L1data$Prob,stdout())
            #write(V1data$Prob,stdout())
            constructMaxProb <- c(1:length(H1data$Prob))
       }
        bestSingleIFO <- data.frame(name = character(), freq = as.numeric(character()), df=as.numeric(character()), R = as.numeric(character()), h0=as.numeric(character()), Prob = as.numeric(character()))
        for (indexProb in constructMaxProb){
            H1Prob <- H1data$Prob[indexProb]
            L1Prob <- L1data$Prob[indexProb]
            V1Prob <- V1data$Prob[indexProb]
            bestIFO <- ''
            if ((H1Prob <= L1Prob) & (H1Prob <= V1Prob)){
                bestIFO <- 'H1'
                bestSingleIFO <- rbind(bestSingleIFO, data.frame(name = as.character(indexProb), freq = as.numeric(H1data$freq[indexProb]), df = as.numeric(H1data$df[indexProb]), R = as.numeric(H1data$R[indexProb]), h0 = as.numeric(H1data$h0[indexProb]), Prob = as.numeric(H1data$Prob[indexProb])))
            }
            if ((L1Prob <= V1Prob) & (L1Prob <= H1Prob)){
                bestIFO <- 'L1'
                bestSingleIFO <- rbind(bestSingleIFO, data.frame(name = as.character(indexProb), freq = as.numeric(L1data$freq[indexProb]), df = as.numeric(L1data$df[indexProb]), R = as.numeric(L1data$R[indexProb]), h0 = as.numeric(L1data$h0[indexProb]), Prob = as.numeric(L1data$Prob[indexProb])))
            }
            if ((V1Prob <= H1Prob) & (V1Prob <= L1Prob)){
                bestIFO <- 'V1'
                bestSingleIFO <- rbind(bestSingleIFO, data.frame(name = as.character(indexProb), freq = as.numeric(V1data$freq[indexProb]), df = as.numeric(V1data$df[indexProb]), R = as.numeric(V1data$R[indexProb]), h0 = as.numeric(V1data$h0[indexProb]), Prob = as.numeric(V1data$Prob[indexProb])))
            }    
        } 
        bestSingleIFO$period <- rep(68023.8259, nrow(bestSingleIFO)) 
        bestSingleIFO$tfnorm <- rep(0.1, nrow(bestSingleIFO)) 
        bestSingleIFO$jobnumber <- rep(0, nrow(bestSingleIFO)) 
        bestSingleIFO$ra <- rep(4.2757, nrow(bestSingleIFO)) 
        bestSingleIFO$dec <- rep(-0.2730, nrow(bestSingleIFO)) 
        bestSingleIFO <- bestSingleIFO[ , c("freq", "period", "df", "ra", "dec", "R", "h0", "Prob", "tfnorm", "jobnumber")]

        # Now we need a variety of modifications for the case where we are
        # comparing the output of both open and closed interferometers
        if (openness == 'both' & paramEst == 'params'){
            write('ESTIMATING PARAMETERS', stdout())
            fListForPlot <- bestSingleIFO$freq
            nonDetectedListBoth <- c(21, 27, 30, 37, 41, 43, 45, 48, 49, 50, 52, 53, 54, 55, 57, 58, 63, 64, 69, 71, 72, 73, 74, 77, 80, 81, 82, 88, 89, 90, 92, 93, 94, 99, 98)
            singleDetectionListBoth <- c(12, 17, 25, 26, 33, 67, 79, 85, 87, 96)
            # Take the differences between recovered and true values
            # Units of Hz difference
            bestSingleIFO$freq <- bestSingleIFO$freq - tableMatrix[-nonDetectedListBoth, 1]
            # Units of asini
            bestSingleIFO$df <- bestSingleIFO$df - tableMatrix[-nonDetectedListBoth, 3]
            #Units of h_effective
            bestSingleIFO$h0 <- bestSingleIFO$h0 - heffectiveCalc(tableMatrix[-nonDetectedListBoth,2], tableMatrix[-nonDetectedListBoth,5])

        }

                write.table(bestSingleIFO, paste('bestIFOdetection', paramString,opennessString,'.txt',sep=""), sep=" ", row.names = FALSE, col.names = FALSE)
        save(bestSingleIFO, file=paste('bestIFOdetection', paramString,opennessString,'.Rda',sep=""))
        #write('Value of best single IFO std of difference in f, asini, h0',stdout())
        #write(sd(bestSingleIFO$freq[-27]),stdout())
        #write(sd(bestSingleIFO$df[-27]),stdout())
        #write(sd(bestSingleIFO$h0[-27]),stdout())

        sdFreq <- sd(bestSingleIFO$freq)
        sdAsini <- sd(bestSingleIFO$df)
        sdh0 <- sd(bestSingleIFO$h0)
        # Noticed that pulsars 87, 96, 97 and 100 were poorly detected
        # Each was above 1150 Hz (indeed, the next lowest open was 1050 Hz)
        # and the only thing in between that was OK was pulsar 91, which
        # had a best p-value of log10p = -388
        # My suspicion is that we should designate log10p = -300 and
        # 1050 Hz as dividing lines; if it is below the former and above
        # the latter, then it gets a loose uncertainty, otherwise tight

        # Checking by hand, we can see the pulsars that fail this criterion
        # are line number 27, 29, 30 and 31 of the detected list

        #startBandListOpen <- c(80.0, 100.0, 110.0, 120.0, 130.0, 140.0, 160.0, 170.0, 200.0, 220.0, 260.0, 280.0, 290.0, 310.0, 320.0, 340.0, 350.0, 370.0, 380.0, 410.0, 420.0, 430.0, 440.0, 460.0, 470.0, 490.0, 500.0, 530.0, 580.0, 600.0, 610.0, 770.0, 830.0, 890.0, 920.0, 960.0, 1040.0, 1050.0, 1120.0, 1160.0, 1170.0, 1180.0, 1250.0, 1280.0, 1290.0, 1300.0, 1330.0, 1360.0, 1390.0, 1450.0)

        if (openness == 'open') {
            LooseTightDivide <- c(27,29,30,31)
            #write(LooseTightDivide, stdout())
        }
        if (openness == 'both') {
            # The pulsars that were poorly detected were numbers 83, 85, 87, 95, 96, 97 and 100
            #write(bestSingleIFO$freq,stdout())
            # In the combined list, sorted by frequency, these are
            LooseTightDivide <- c(56, 58, 60, 62, 63, 64, 65)
        }

        sdFreqLoose <- sd(bestSingleIFO$freq[LooseTightDivide])
        sdAsiniLoose <- sd(bestSingleIFO$df[LooseTightDivide])
        sdh0Loose <- sd(bestSingleIFO$h0[LooseTightDivide])

        sdFreqTight <- sd(bestSingleIFO$freq[-LooseTightDivide])
        sdAsiniTight <- sd(bestSingleIFO$df[-LooseTightDivide])
        sdh0Tight <- sd(bestSingleIFO$h0[-LooseTightDivide])


        
        if (paramEst == 'params') {
            png('ErrorF.png')
            plot(log10(abs(bestSingleIFO$Prob)), log10(abs(bestSingleIFO$freq)),xlab='log10(abs(log10p)) of best single IFO candidate', ylab='log10(abs(error in estimating frequency [Hz]))', main=paste('Error \n (sigma (loose) = ', as.character(sdFreqLoose),' Hz) \n in f vs p-value for Scorpius X-1 MDC'))
            regression <- lm(log10(abs(bestSingleIFO$freq))~log10(abs(bestSingleIFO$Prob)))
            abline(regression, col="red")
            abline(h=log10(sdFreq), col="blue")
            dev.off()
            png('ErrorAsini.png')
            plot(log10(abs(bestSingleIFO$Prob)), log10(abs(bestSingleIFO$df)),xlab='log10(abs(log10p)) of best single IFO candidate', ylab='log10(abs(error in estimating asini [s]))', main=paste('Error \n (sigma (loose) = ', as.character(sdAsiniLoose),' s) \n in asini vs p-value for Scorpius X-1 MDC'))
            regression <- lm(log10(abs(bestSingleIFO$df))~log10(abs(bestSingleIFO$Prob)))
            abline(regression, col="red")
            abline(h=log10(sdAsini), col="blue")
            dev.off()
            png('Errorh0.png')
            plot(log10(abs(bestSingleIFO$Prob)), log10(abs(bestSingleIFO$h0)),xlab='log10p(abs(log10p)) of best single IFO candidate', ylab='log10(abs(error in estimating h0))', main=paste('Error \n (sigma (loose) = ', as.character(sdh0Loose),' )) \n in h0 vs p-value for Scorpius X-1 MDC'))
            regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
            abline(regression, col="red")
            abline(h=log10(sdh0), col="blue")
            dev.off()

            # On a hunch, let us try plotting vs frequency proper
            nonDetectedListOpen <- c(14, 16, 20, 25, 26, 28, 29, 30, 33, 34, 36, 37, 38, 41, 42, 44, 45, 46, 49)
            singleDetectionListOpen <- c(7, 13, 18, 40, 47)
            if (openness == 'open' | openness == 'closed') {
                fListForPlot <- tableMatrix[-nonDetectedListOpen,1]
            }
    
            png('ErrorFvsF.png')
            plot(fListForPlot, log10(abs(bestSingleIFO$freq)),xlab='Injection frequency [Hz]', ylab='log10(abs(error in estimating frequency [Hz]))', main=paste('Error \n (sigma (tight) = ', as.character(sdFreqTight),' Hz) \n in f vs injection f for Scorpius X-1 MDC'))
            #regression <- lm(log10(abs(bestSingleIFO$freq))~log10(abs(bestSingleIFO$Prob)))
            #abline(regression, col="red")
            abline(h=log10(sdFreq), col="blue")
            dev.off()
            png('ErrorAsinivsF.png')
            plot(fListForPlot, log10(abs(bestSingleIFO$df)),xlab='Injection frequency [Hz]', ylab='log10(abs(error in estimating asini [s]))', main=paste('Error \n (sigma (tight) = ', as.character(sdAsiniTight),' s) \n in asini vs injection f for Scorpius X-1 MDC'))
           #regression <- lm(log10(abs(bestSingleIFO$df))~log10(abs(bestSingleIFO$Prob)))
            #abline(regression, col="red")
            abline(h=log10(sdAsini), col="blue")
            dev.off()
            png('Errorh0vsF.png')
            plot(fListForPlot, log10(abs(bestSingleIFO$h0)),xlab='Injection frequency [Hz]', ylab='log10(abs(error in estimating h0 ]))', main=paste('Error \n (sigma (tight) = ', as.character(sdh0Tight),' )) \n in h0 vs injection f for Scorpius X-1 MDC'))
            #regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
            #abline(regression, col="red")
            abline(h=log10(sdh0), col="blue")
            dev.off()


            # Finally, we can plot h-recovered vs h-effective
            if (openness == 'open') {
                heffForPlot <- heffectiveCalc(tableMatrix[,2], tableMatrix[,5])[-nonDetectedListOpen]
            }
            if (openness == 'both') {
                heffForPlot <- heffectiveCalc(tableMatrix[,2], tableMatrix[,5])[-nonDetectedListBoth]
            }
            png('detectedHerrVsHeffective.png')
            plot(heffForPlot, log10(abs(bestSingleIFO$h0)),xlab='heffForPlot', ylab='log10(abs(error in estimating h0))', main=paste('Error \n (sigma = ', as.character(sdh0),' )) \n in h0 vs h-effective for Scorpius X-1 MDC'))
            #regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
            #abline(regression, col="red")
            #abline(h=log10(sdh0), col="blue")
            dev.off()


       }    
   } # function paramGrapher 
   paramGrapher(useModeBlinded,useModeDetectionSelection,useModeParamDiff)
    

    upperLimitGrapher <- function(openness,selectSeen,paramEst) {
        # This function presupposes that all three observatories have been run
        # through extractTwoSpectOutput
        # with the useMode options
        # bypassCut, performSort, performDivide, open, nondetected, ULs
        # This function is used for extracting upper limit uncertainties
        # For actual detection assessement, one can use either open or closed,
        # but ULs must be replaced with noDiff 
        write("Calculating uncertainty in upper limits", stdout())
        # For the LIGO Hanford Cluster:
        #locationHead <- "/home/gmeadors/ScoX1/2014/04/23/Meadors-ScoX1/"
        # For Gallatin at University of Michigan:
        locationHead <- "/home/gmeadors/Documents/ScoX1/2014/04/28-POB-backup/backup_peaks-of-bands/"
        if (openness == 'open') {
            opennessString <- 'Open'
        }
        if (openness == 'closed') {
            opennessString <- 'Closed'
        }
        if (paramEst == 'ULs') {
            paramString <- 'CalUL'
        }
        if (paramEst == 'noDiff') {
            paramString <- ''
        }
        # Load in data for the pulsars that are not seen:
        loadedData <- load(paste(locationHead,'NonDetectedULs',paramString,opennessString,'H1','.Rda',sep=""))
        H1data <- get(loadedData)
        rm(loadedData)
        loadedData <- load(paste(locationHead,'NonDetectedULs',paramString,opennessString,'L1','.Rda',sep=""))
        L1data <- get(loadedData)
        rm(loadedData)
        loadedData <- load(paste(locationHead,'NonDetectedULs',paramString,opennessString,'V1','.Rda',sep=""))
        V1data <- get(loadedData)
        rm(loadedData)
        # Now load in data for the pulsars that are seen:
        loadedData <- load(paste(locationHead,'DetectedParameters',opennessString,'H1','.Rda',sep=""))
        seenH1data <- get(loadedData)
        rm(loadedData)
        loadedData <- load(paste(locationHead,'DetectedParameters',opennessString,'L1','.Rda',sep=""))
        seenL1data <- get(loadedData)
        rm(loadedData)
        loadedData <- load(paste(locationHead,'DetectedParameters',opennessString,'V1','.Rda',sep=""))
        seenV1data <- get(loadedData)
        rm(loadedData)

        constructDEunseen <- c(1:length(H1data$Prob))
        constructDEseen <- c(1:length(seenH1data$Prob))
        constructDEUL <- c(1:(length(H1data$Prob)+length(seenH1data$Prob)))
        DEframe <- data.frame(detection = as.numeric(character()), name = character(), freq = as.numeric(character()), df=as.numeric(character()), R = as.numeric(character()), h0=as.numeric(character()), Prob = as.numeric(character()))
        ULframe <- data.frame(H1h0 = as.numeric(character()), L1h0 = as.numeric(character()), V1h0 = as.numeric(character()), freq = as.numeric(character()), pairs = as.numeric(character()))
        for (indexItem in constructDEunseen){
                DEframe <- rbind(DEframe, data.frame(detection = 0, name = as.character(indexItem), freq = as.numeric(H1data$freq[indexItem]), df = as.numeric(H1data$df[indexItem]), R = as.numeric(H1data$R[indexItem]), h0 = as.numeric(H1data$h0[indexItem]), Prob = as.numeric(H1data$Prob[indexItem])))
        } 
        for (indexItem in constructDEseen){
                DEframe <- rbind(DEframe, data.frame(detection = 1, name = as.character(indexItem), freq = as.numeric(seenH1data$freq[indexItem]), df = as.numeric(seenH1data$df[indexItem]), R = as.numeric(seenH1data$R[indexItem]), h0 = as.numeric(seenH1data$h0[indexItem]), Prob = as.numeric(seenH1data$Prob[indexItem])))
        } 

        for (indexItem in constructDEunseen){
            ULframe <- rbind(ULframe, data.frame(H1h0 = as.numeric(H1data$h0[indexItem]), L1h0 = as.numeric(L1data$h0[indexItem]), V1h0 = as.numeric(V1data$h0[indexItem]), freq = as.numeric(H1data$freq[indexItem]), pairs = as.numeric(0)))        
        }
        for (indexItem in constructDEseen){
            ULframe <- rbind(ULframe, data.frame(H1h0 = as.numeric(seenH1data$h0[indexItem]), L1h0 = as.numeric(seenL1data$h0[indexItem]), V1h0 = as.numeric(seenV1data$h0[indexItem]), freq = as.numeric(seenH1data$freq[indexItem]), pairs = as.numeric(3)))        
        }
        # Arrange by order of frequency, monotonic to pulsar number
        # so that we can specify the correct single-pair detections
        # and attached the true effective injection h0
        ULframe <- ULframe[order(ULframe$freq), ]
        singleDetectionListOpen <- c(7, 13, 18, 40, 47)
        for (indexItem in 1:length(ULframe$freq)) {
            #write(any(indexItem == singleDetectionListOpen),stdout())
            if (any(indexItem == singleDetectionListOpen)){
                ULframe$pairs[indexItem] <- as.numeric(1)
            }
        }
        #write(ULframe$pairs, stdout())
        
  

        DEframe$period <- rep(68023.8259, nrow(DEframe)) 
        DEframe$tfnorm <- rep(0.1, nrow(DEframe)) 
        DEframe$jobnumber <- rep(0, nrow(DEframe)) 
        DEframe$ra <- rep(4.2757, nrow(DEframe)) 
        DEframe$dec <- rep(-0.2730, nrow(DEframe)) 
        DEframe <- DEframe[ , c("detection","freq", "period", "df", "ra", "dec", "R", "h0", "Prob", "tfnorm", "jobnumber")]
        # Arrange by order of frequency, monotonic to pulsar number
        # so that we can attach the true effective injection h0
        DEframe <- DEframe[order(DEframe$freq),] 
        #write(DEframe[order(DEframe$freq),]$freq,stdout())
        #write(DEframe$freq,stdout())
        #write(DEframe$detection,stdout())
        trueHeffInj <- heffectiveCalc(tableMatrix[,2], tableMatrix[,5])
        DEframe <- cbind(DEframe, trueHeffInj)
        ULframe <- cbind(ULframe, trueHeffInj)
        #write(DEframe$trueHeffInj,stdout())
        # Now re-arrange by the true effective injection h0
        DEframe <- DEframe[order(DEframe$trueHeffInj),]
        ULframe <- ULframe[order(ULframe$trueHeffInj),]
        #write(DEframe$detection,stdout())

        # Now bin the data using cut and breaks to calculate detection efficiency
        #breaks <-c(1e-27,1e-25,2e-25,3e-25,4e-25,5e-25)
        breaks <- c(2e-26,1e-25, 1.3e-25,1.5e-25,2.5e-25,3e-25,9e-25)
        binGroups <- cut(DEframe$trueHeffInj, breaks, labels=breaks[-1])
        #binnedDetection <- DEframe$detection
        binnedDetection <- tapply(DEframe$detection, binGroups, mean)
        nDetection <- tapply(DEframe$detection, binGroups, length) 
        uncertaintyDetection <- sqrt((binnedDetection)*(1-binnedDetection)/nDetection)
        write(binnedDetection, stdout())
        write(uncertaintyDetection, stdout())
        png('detectionVsHeffective.png')
        #    histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        mp <- plot(breaks[-1], binnedDetection-uncertaintyDetection ,type="s",xlab='highest h-effective', ylab='Fraction detected', main=paste('(Detection efficiency - uncertainty) for Scorpius X-1 MDC'))
        #segments(mp, binnedDetection - uncertaintyDetection, mp, binnedDetection + uncertaintyDetection, lwd=2)
        #segments(mp-0.1, binnedDetection - uncertaintyDetection, mp+0.1, binnedDetection - uncertaintyDetection, lwd=2)
        #segments(mp-0.1, binnedDetection + uncertaintyDetection, mp+0.1, binnedDetection + uncertaintyDetection, lwd=2)
        #points(break[-1], binnedDetection+uncertaintyDetection,col="green")
        #points(break[-1], binnedDetection-uncertaintyDetection,col="green")
        #regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
        #abline(regression, col="red")
        #abline(h=log10(sdh0), col="blue")
        dev.off()

        png('HrecoveredVsHeffective.png')
        #    histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        plot(DEframe$trueHeffInj, DEframe$h0, xlab='injected h-effective', ylab='recovered h-effective', main=paste('h-recovered vs h-injected for Scorpius X-1 MDC'))
        #regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
        #abline(regression, col="red")
        abline(0, 1)
        dev.off()


        pairNone <- (ULframe$pairs == 0)
        pairOne <- (ULframe$pairs == 1)
        pairThree <- (ULframe$pairs == 3)
        # Find out the 95 percent upper limit for which pulsars were
        # NOT detected
        AllFailedDetections <- rbind(ULframe$H1h0[pairNone], ULframe$L1h0[pairNone], ULframe$V1h0[pairNone])
        SortedAllFailedDetections <- AllFailedDetections[order(AllFailedDetections)]
        Index95PercentConf <- ceiling(0.95*length(SortedAllFailedDetections))
        Limit95PercentConf <- SortedAllFailedDetections[Index95PercentConf]
        write(Limit95PercentConf, stdout())


        ULallTripleInj <-(ULframe$trueHeffInj[pairThree] + ULframe$trueHeffInj[pairThree] + ULframe$trueHeffInj[pairThree])/3  
        ULallTripleRec <-(ULframe$H1h0[pairThree]+ ULframe$L1h0[pairThree] + ULframe$V1h0[pairThree])/3  
        regressHscale <- lm(ULallTripleRec ~ ULallTripleInj)  
        #write(ULallTripleRec,stdout())
        write(coef(regressHscale)[2],stdout())

        rescaleHslope <- coef(regressHscale)[2]
        rescaleHsd <- sd(ULallTripleRec/rescaleHslope - ULallTripleInj)

        png('HrecoveredVsHeffectiveFullUL.png')
        #    histCase <- hist(dataHist, freq=TRUE, col='blue', breakNo, xlab=quantity, main = mainString)
        plot(ULframe$trueHeffInj[pairThree], ULframe$H1h0[pairThree], cex.main=0.8, col='green',xlab='injected h-effective', ylab='recovered h-effective', main=paste('All IFOs, h-recovered vs h-injected for Scorpius X-1 MDC \n', 'Green: 3 pairs, Blue: 1 pair, Red: 0 pairs \n', '95 percent confident detection in effective h0:',Limit95PercentConf, '\n slope, sd of recovered v injected:', rescaleHslope, ',', rescaleHsd))
        points(ULframe$trueHeffInj[pairThree], ULframe$L1h0[pairThree], col='green')
        points(ULframe$trueHeffInj[pairThree], ULframe$V1h0[pairThree], col='green')
        points(ULframe$trueHeffInj[pairOne], ULframe$H1h0[pairOne], col='blue')
        points(ULframe$trueHeffInj[pairOne], ULframe$L1h0[pairOne], col='blue')
        points(ULframe$trueHeffInj[pairOne], ULframe$V1h0[pairOne], col='blue')
        points(ULframe$trueHeffInj[pairNone], ULframe$H1h0[pairNone], col='red')
        points(ULframe$trueHeffInj[pairNone], ULframe$L1h0[pairNone], col='red')
        points(ULframe$trueHeffInj[pairNone], ULframe$V1h0[pairNone], col='red')
        #regression <- lm(log10(abs(bestSingleIFO$h0))~log10(abs(bestSingleIFO$Prob)))
        #abline(regression, col="red")
        abline(0, 1)
        abline(Limit95PercentConf,0, col = "red")
        abline(coef(regressHscale)[1], coef(regressHscale)[2], col = "green")
        dev.off()






       

    } # function upperLimitGrapher
    #upperLimitGrapher(useModeBlinded,useModeDetectionSelection,useModeParamDiff)


} # useMode histogram


tableReaderE <- function(filename, headerFlag) {
    # Following a general design for reading tables proposed by
    # Ethan Obie Romero-Severson at Los Alamos National Labs
    d = read.csv(filename,header=headerFlag)
    names = c()
    new.dat = c()
    for(i in 1:ncol(d)) {
        tok = strsplit(as.character(d[,i]), '=')
        val = sapply(tok, '[[', 2)
        names = c(names, tok[[1]][1])
        new.dat=cbind(new.dat, val)
        new.dat = data.frame(new.dat)
        names(new.dat) = names
        return(new.dat)
    }
} # function tableReaderE
