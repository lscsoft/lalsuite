#########################
#                       #
#  mcmcsummary()        #
#  by Christian Roever  #
#                       #
#########################

# "mcmcsummary" creates some summarising plots & tables for MCMC output,
# organises them in a little directory tree and sets up an HTML
# document from where to access graphs and tables.
#
# The function you need to worry about is "mcmcsummary()".
# You need to have the "xtable" package installed
# (see the help on "install.packages()" for how to do that).
# If the "coda" package is installed (not required), you'll also
# get "highest posterior density intervals" instead of just
# "central posterior intervals".
#
# The arguments of "mcmcsummary()" are:
#
#  - data: a data frame or matrix whose columns contain the MCMC samples of 
#          different variables; rows correspond to MCMC iterations.
#          May also be an array, whose "slices" are then taken to be
#          samples from parallel MCMC chains.
#  - targetdirectory: the directiory the files are supposed to go to.
#          Will be created if necessary and will be cleared
#          if "overwrite=TRUE" (see below).
#  - iteration: vector of iterations corresponding to each row of "data"
#  - burnin: the last iteration of the burn-in phase to be discarded.
#          Set to something below "min(iterations)" for not discarding anything.
#          If "data" is an array, "burnin" may be a vector with each
#          element corresponding to one of the chains. In that case, the burnin
#          may be larger than "max(iterations)" (or just "Inf") for some (!) of
#          the chains, but there need to be some non-discarded samples left
#          for at least one of the chains.
#  - autoburnin: flag to tell the code to figure out the burn-in by itself.
#          You need to also specify the "posterior" parameter below.
#  - posterior: vector (or matrix) of posterior density values corresponding
#          to each sample. Used in order to determine the burn-in for each chain.
#          You can also use likelihood values instead.
#  - threshold: criterion to decide when burn-in is finished; default is
#          square root of the number of parameters (columns of data matrix).
#          Burn-in is assumed finished for each chain once its posterior density
#          comes within <threshold> of the maximum value over all chains.
#  - varnames: character vector giving the names of the variables
#          (corresponding to columns of "data")
#  - truevalue: a vector containing the true values of the variables
#          (again, corresponding to columns of "data").
#          If supplied, this will appear in plots and tables.
#          This vector may also contain missing values ("NA"s).
#  - graphicsformats: vector of graphics formats to be produced.
#          Possible values (by now) are "jpeg", "png", "eps", "pdf".
#          "jpeg" graphics are always produced.
#  - overwrite: (logical) flag indicating whether to clear / overwrite files
#          if the "targetdirectory" (see above) already exists and is not empty.
#  - signifdigits: number of significant digits to be supplied in tables. 
#          For the "summary stats" table, each variable's standard deviation 
#          is given with this accuracy, and all other figures are then aligned
#          with this. Defaults to 3.
#  - plotN: the (maximum) number of samples to use in plots where otherwise ALL
#          samples would appear in (e.g. trace plots). Defaults to 500.
#          Note that too large values result in huge file sizes in the
#          vector graphics formats.
#
# The output is eventually accessible from a file "index.html" in the directory
# given by "targetdirectory".
#
#  /!\  CAUTION!
#       Always be sure you know what you're doing when
#       specifying the "targetdirectory", especially when using
#       the "overwrite=TRUE" option !!
#       Otherwise you might accidentally wipe out
#       the wrong directory...
#
# The code should work under Linux, but I haven't tried it on other platforms.


# The probably simplest possible example is:
#
#   source("mcmcsummary.R")
#   mcmcsummary(iris[,1:4], target="mcmcsummary-irisexample")
#
# The "iris" data has nothing to do with MCMC, but this way at least
# some plots & tables are created.
# Now open the file "mcmcsummary-irisexample/index.html" in a web browser.
# You may also want to have a look at the HTML source code itself.


# Read and display results from a single file:
#
#   post1 <- read.table("/home/user/data/file012.txt", header=TRUE)
#
#   mcmcsummary(post1[,6:11], target="example01", overwrite=TRUE, burnin=c(10000),
#               truevalue=c(1,10,1,-46,-4.6,100), iteration=post1[,1])
#
#
#
# In order to diplay results from multiple parallel runs, the data need to be provided
# in an array (instead of a data frame or matrix). First read data...:
#
#   post1 <- read.table("/home/user/data/file012.txt", header=TRUE)
#   post2 <- read.table("/home/user/data/file013.txt", header=TRUE)
#   post3 <- read.table("/home/user/data/file014.txt", header=TRUE)
#   post4 <- read.table("/home/user/data/file015.txt", header=TRUE)
#
# Conversion from data frame to array goes via matrix and vector...:
#
#   post <- array(c(as.vector(as.matrix(post1[,6:11])),
#                   as.vector(as.matrix(post2[,6:11])),
#                   as.vector(as.matrix(post3[,6:11])),
#                   as.vector(as.matrix(post4[,6:11]))), dim=c(nrow(post1),6,4))
#   colnames(post) <- colnames(post1)[6:11]
#
# Note the "dim" argument: all data frames are of same size, i.e. have the same
# number of rows as the 1st one. There are 4 data frames with 6 columns
# (numbers 6-11) each.
# Now run "mcmcsummary()" with that array as argument:
#
#   mcmcsummary(post, target="example02", overwrite=TRUE,
#               burnin=c(10000,20000,50000,15000),
#               truevalue=c(1,10,1,-46,-4.6,100), iteration=post1[,1])
#

mcmcsummary <- function(data, targetdirectory="mcmcsummary",
                        iteration=seq(from=1, to=nrow(data), by=1),
                        burnin=round(mean(range(iteration))),
                        autoburnin=FALSE, posterior=NULL, threshold=sqrt(ncol(data)),
                        varnames=colnames(data),
                        truevalue=NULL,
                        graphicsformats=c("jpeg", "eps"),
                        overwrite=FALSE,
                        signifdigits=3,
                        plotN=500)
# this is the main function.
{
  # some auxiliary function definitions:
  
  computeKDE <- function(x)
  # a little wrapper function
  # for kernel density estimate computation
  {
    bw <- bw.nrd(x)
    if (!(bw>0)) {
      bw <- 1.06 * sd(x) * length(x)^(-1/5)
      # (this is a copy-paste from "bw.nrd()" to avoid zero bandwidth)
    }
    if (!(bw>0)) {
      bw <- 1
      warning("zero bandwidth argument 'manually' fixed for kernel density estimation.")
    }
    result <- density(x, bw=bw, kern="epanechnikov")
  }

  plotKDE <- function(x, main="", xlab="", ylab="density", trueval=NULL)
  # plotting of kernel density estimates
  # (with grey area under curve etc.)
  {
    xrange <- range(x$x)
    if ((!is.null(trueval)) && (!is.na(trueval))) # expand plotting range if necessary
      xrange <- range(c(xrange,trueval))
    plot(xrange, range(c(0,x$y)), type="n",
         main=main, xlab=xlab, ylab=ylab, axes=FALSE)
    polygon(c(min(x$x),x$x, max(x$x)),
          c(0, x$y, 0), col="lightgrey", border=NA)
    lines(x$x, x$y)
    abline(h=0, col="black")
    if ((!is.null(trueval)) && (!is.na(trueval)))
      abline(v=trueval, lty="dashed", col="darkgrey")
    axis(1); axis(2); box()
    return(invisible())
  }

  traceplot <- function(iter, x, xlab="iteration", ylab="", col=1)
  # function for trace plots
  {
    #plot(iter, x, type="l", xlab=xlab, ylab=ylab)
    matplot(iter, x, type="l", lty="solid", col=col, xlab=xlab, ylab=ylab)
    return(invisible())
  }

  scatterplot <- function(x, y, xlab="x", ylab="y", col=1, trueval=NULL)
  # function for scatter plots
  {
    xrange <- range(x)
    yrange <- range(y)
    if (!is.null(trueval)){ # expand plotting ranges if necessary
      if (!is.na(trueval[1])) xrange <- range(c(xrange,trueval[1]))
      if (!is.na(trueval[2])) yrange <- range(c(yrange,trueval[2]))
    }    
    plot(x, y, pch=".", col=col, xlab=xlab, ylab=ylab, axes=FALSE,
         xlim=xrange, ylim=yrange)
    if (!is.null(trueval)){
      if (!is.na(trueval[1])) abline(v=trueval[1], lty="dashed", col="darkgrey")
      if (!is.na(trueval[2])) abline(h=trueval[2], lty="dashed", col="darkgrey")
    }
    axis(1); axis(2); box()
    return(invisible())
  }

  htmlheader <- function(funcall,
                         title = "MCMC output summary")
  {
    htmlcode <- c('<html>',
                  '<head>',
                  paste('  <title>',title,'</title>', sep=""),
                  '  <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">',
                  "</head>", "",
                  #'<body style="font-family: helvetica,arial,sans-serif">',
                  '<body>',
                  "")
    return(htmlcode)
  }
  
  htmlfooter <- function(bottomline = paste('generated', date(), 'using <a href="http://www.r-project.org" title="R project homepage"><b>R</b></a> and the <a href="http://cran.r-project.org/package=xtable"><code>xtable</code></a> ',ifelse(codainstalled,'and <a href="http://cran.r-project.org/package=coda"><code>coda</code></a> packages.','package.')))
  {
    htmlcode <- c('  <hr>',
                  paste('  <div align="right"><small><i>', bottomline, '</i></small></div>', sep=''),
                  '</body>', '</html>')
    return(htmlcode)
  }

  htmlplot <- function(filename, versions = "jpeg")
                       #width=200, height=200)
  {
    htmlcode <- c(paste('  <img src="jpeg/',filename,'.jpeg','"', #, sep=''), paste('
                        ' alt="',filename,'.jpeg"','>', sep=''))
                  #paste('       width="',width,'" height="',height,'"', '>', sep=''))
    versionsline <- '  <small>(image versions: '
    i <- 1
    while (i <= length(versions)){
      imgfile <- paste(versions[i], '/', filename, '.', versions[i], sep='')
      versionsline <- paste(versionsline, '<a href="', imgfile, '" ',
                            switch(versions[i],
                                   jpeg = 'type="image/jpeg"',
                                   png  = 'type="image/png"',
                                   eps  = 'type="application/postscript"',
                                   pdf  = 'type="application/pdf"')
                            ,'>', versions[i], '</a>', sep='')
      if (i < length(versions))
        versionsline <- paste(versionsline, ', ',sep='')
      i <- i+1
    }
    versionsline <- paste(versionsline, ')</small>',sep='')
    htmlcode <- c(htmlcode, '  <br>', versionsline)#, '  <br>', '')
    return(htmlcode)
  }

  htmlcomment <- function(commentstring="...")
  {
    htmlcode <- c(paste('  <!--', commentstring, '-->'), '')
    return(htmlcode)
  }
  
  htmlhorizline <- function()
  {
    htmlcode <- c('<hr>', '')
    return(htmlcode)
  }

  htmlsection <- function(level=1, heading="section")
  {
    htmlcode <- c(paste('<h', level, '>', heading, '</h', level, '>', sep=''), '')
    return(htmlcode)
  }
  
  htmlfragment <- function(fragmentname)
  {
    htmlcode <- c(paste('<a name="', fragmentname, '"></a>',sep=''))
    return(htmlcode)
  }

  sumstats <- function(x)
  {
    if (codainstalled){
      int <- HPDinterval(mcmc(x), prob=0.95)[1,]
      names(int) <- NULL
    }
    else
      int <- quantile(x, c(0.025,0.975), names=F)
    result <- c("mean"=mean(x),
                "median"=median(x),
                "std. dev."=sd(x),
                "95% lower"=int[1],
                "95% upper"=int[2])
    return(result)
  }

  decplaces <- function(x, signifdigits=3)
  # decimal places (after decimal point)
  {
    return(max(c(0, -(floor(log10(x))-(signifdigits-1)))))
  }
  
  #
 ###  here's where the ACTUAL CODE starts:
  #
  
  ok <- require("xtable")
  if (!ok) {
    warning("missing add-on package 'xtable'")
    return(invisible())
  }
  codainstalled <- require("coda")
  
  # some initial definitions:
  nvar   <- ncol(data)                                       # number of variables
  stopifnot(length(varnames)==nvar)
  nobs   <- nrow(data)                                       # number of samples
  nchain <- ifelse(length(dim(data))==3, dim(data)[3], 1)    # number of parallel chains
  if (nchain==1) {
    data <- as.matrix(data)
    dim(data)[3] <- 1
  }
  if (!is.null(truevalue))
    stopifnot(length(truevalue)==nvar)
  if (autoburnin) { # (try to) figure out burn-ins automatically:
    if (is.null(posterior)) {
      warning("Need to provide 'posterior' argument when using 'autoburnin' option!")
      return(invisible())
    }
    stopifnot((((nchain==1) && (is.vector(posterior) || (ncol(posterior)==1))) ||
               ((nchain>1) && (!is.vector(posterior) && (ncol(posterior)==nchain)))))
    maxpost <- max(posterior)
    abfun <- function(postvec)
    {
      result <- Inf
      if (max(postvec)>= maxpost-threshold)
        result <- iteration[max(c(1, which(postvec >= (maxpost-threshold))[1] - 1))]
      return(result)
      # (first sample to cross threshold is non-burnin)
    }
    burnin <- as.vector(apply(as.matrix(posterior),2,abfun))
    #print(burnin)
  }
  stopifnot(length(iteration)==nobs)
  stopifnot(any(burnin < max(iteration)))
  stopifnot((length(burnin)==1) | (length(burnin)==nchain))
  if ((nchain>1) && (length(burnin)==1))
    burnin <- rep(burnin, nchain)
  #traceobs <- seq(1, nobs, by=ceiling(nobs/plotN))          # observations shown in trace plots
  traceobs <- unique(round(seq(1, nobs, by=nobs/plotN)))     # observations shown in trace plots

  posteriorsample <- NULL  # eventual posterior sample (with chains united and burnins discarded)
  chain <- NULL            # corresponding chain number indicator
  scatterobs <- NULL       # indicator for subset to be used in scatter plots (~plotN from each chain)

  # derive "posteriorsample" matrix:
  for (i in 1:nchain) {
    newN <- sum(iteration>burnin[i])
    if (newN > 0) {
      posteriorsample <- rbind(posteriorsample, data[iteration>burnin[i],,i])
      chain <- c(chain, rep(i, newN))
      newObs <- rep(FALSE, newN)
      #newObs[seq(1, newN, by=ceiling(newN/(plotN/nchain)))] <- TRUE
      newObs[unique(round(seq(1, newN, by=(newN/(plotN/nchain)))))] <- TRUE
      scatterobs <- c(scatterobs, newObs)
    }
  }
  remove(newN, newObs)
  scatterobs <- which(scatterobs)
  scatterobs <- sample(scatterobs, length(scatterobs)) # permutate for nicer plotting
  # 'colvec' defines the plotting colours corresponding to chains (trace plots & scatter plots)
  if (nchain==1)
    colvec <- rgb(0,0,0)
  else
    colvec <- rainbow(nchain)

  # check the specified graphics options:
  knowndevices <- is.element(graphicsformats, c("jpeg", "eps", "png", "pdf"))
  if (any(!knowndevices)) warning("unknown graphics device specified")
  graphicsformats <- graphicsformats[knowndevices]
  if (!is.element("jpeg",graphicsformats)) graphicsformats <- c("jpeg", graphicsformats)
  ngraph <- length(graphicsformats)

  # create directory structure:
  ok <- dir.create(targetdirectory, showWarnings=FALSE)
  if (!ok) { # (problem creating subdirectory)
    if (overwrite){ # try to clear the directory:
      # get files in "targetdirectory"
      rootfiles <- paste(targetdirectory,list.files(targetdirectory),sep="/")
      rootinfo <- file.info(rootfiles)
      # now iterate over the subdirectories:
      for (subdir in rootfiles[rootinfo$isdir]){
        # list & remove files:
        subfiles <- list.files(subdir)
        file.remove(paste(subdir,subfiles,sep="/"))
        file.remove(subdir)
      }
      # finally remove files in "targetdirectory":
      file.remove(rootfiles[! rootinfo$isdir])
    }
    else{ # cancel attempt:
      warning("failed to create target directory '",targetdirectory,"'")
      return(invisible())
    }
  }
  # create subdirectories of "targetdirectory":
  ok <- dir.create(paste(targetdirectory, "txt", sep="/"))
  ok <- (ok && dir.create(paste(targetdirectory, "latex", sep="/")))
  for (subdir in graphicsformats)
    ok <- (ok && dir.create(paste(targetdirectory, subdir, sep="/")))
  if (!ok){
    warning("troubles creating subdirectories")
    return(invisible())
  }
  # directory structure now complete.

  # start the HTML file::
  htmlfilename <- paste(targetdirectory, "index.html", sep="/")  
  write(htmlheader(),
        file=htmlfilename, append=FALSE)
  funcall <- match.call()
  write(c('<!--',
          paste('     ', c('', paste('this summary was generated', date(), 'by calling:'), '', deparse(funcall), ''), sep=''),
          '-->', ''),
        file=htmlfilename, append=TRUE)

  
  write(c("", "<!--", ""),
        file=htmlfilename, append=TRUE)
  write(htmlsection(1, paste("ENTER A MAIN TITLE HERE")),
        file=htmlfilename, append=TRUE)
  write(c('  ENTER MORE DETAILS HERE',
          '  and remove the out-commenting marks above and below.',
          '  <br>', ''),
        file=htmlfilename, append=TRUE)
  write(htmlhorizline(),
        file=htmlfilename, append=TRUE)
  write(c("-->", "", ""),
        file=htmlfilename, append=TRUE)

  
  # compute & tabulate summary stats:
  write(htmlcomment("SUMMARY STATS:"),
        file=htmlfilename, append=TRUE)
  write(htmlsection(1, paste("Summary statistics")),
        file=htmlfilename, append=TRUE)
  write(c('',
          paste('  ', formatC(nobs,digits=0,format="f"), ' MCMC samples of ',nvar,' variables (iterations ',formatC(min(iteration),digits=0,format="f"),'&#150;',formatC(max(iteration),digits=0,format="f"),')<br>',sep='')),
        file=htmlfilename, append=TRUE)
  if (nchain>1)
    write(c(paste('  ','in each of ', nchain, ' parallel chains (',nobs*nchain,' samples total).<br>',sep='')),
          file=htmlfilename, append=TRUE)
#  if (any(discard))
#    write(c(paste('  ', 'Iterations ',min(iteration),'&#150;',burnin,' (',sum(discard),' samples) are discarded as <i>burn-in</i>,<br>',sep=''),
#            paste('  ', sum(!discard), ' samples (iterations ',min(iteration[!discard]),'&#150;',max(iteration),') left for analysis.<br>',sep='')),
#        file=htmlfilename, append=TRUE)
  write(c(paste('  ', formatC(nrow(posteriorsample),digits=0,format="f"), ' samples left after discarding burn-in.<br>',sep='')),
        file=htmlfilename, append=TRUE)

  if (autoburnin)    
    write(ifelse(nchain==1, '  (Burn-in was determined automatically.)<br>',
                            '  (Burn-ins were determined automatically.)<br>'),
          file=htmlfilename, append=TRUE)
  
  if (codainstalled && nchain>1){                # R^p computation only for multiple chains
    if (max(burnin) < iteration[nobs-1]){        # minimum requirement: at least two samples left
      common <- which((iteration > max(burnin))) # the (final) common part of all chains
      codaList <- mcmc.list()
      for (i in 1:nchain)
        codaList[[i]] <- mcmc(data=data[common,,i], start=iteration[common[1]], thin=diff(iteration[1:2]))
      # Brooks and Gelman's multivariate "potential scale reduction factor":
      RHatP <- Re(gelman.diag(codaList)$mpsrf)   
      write(paste('  Convergence diagnostic: R<sup>p</sup> = ',
                  formatC(round(RHatP,max(5,signifdigits)), digits=max(5,signifdigits), format="f"),
                  ' (based on iterations ',formatC(iteration[common[1]],digits=0,format="f"),' and following).<br>', sep=''),
            file=htmlfilename, append=TRUE)
      write(htmlcomment(paste("R^p is Brooks and Gelman's 'potential scale reduction factor',",
                              " and should be close to one.\n",
                              "  See: Brooks & Gelman (1998): General methods for monitoring",
                              " convergence of iterative simulations.", sep="")),
            file=htmlfilename, append=TRUE)
      remove(codaList)
      remove(common)
    }
  }

  write(c('  <br>',''), file=htmlfilename, append=TRUE)
  
  interval <- ifelse(codainstalled, "highest posterior density interval", "confidence bounds")
  if (is.null(truevalue))
    write(c(paste("  mean, median, standard deviation, and 95%",interval,":<br>",sep=""), ""),
          file=htmlfilename, append=TRUE)
  else
    write(c(paste("  mean, median, standard deviation, 95% ",interval,", and true value:<br>",sep=""), ""),
          file=htmlfilename, append=TRUE)

  # compute summary stats:
  # INSERT: number of samples and IACT (?)
  sumtab <- t(apply(posteriorsample, 2, sumstats))
  if (!is.null(truevalue))
    sumtab <- cbind(sumtab, "true"=truevalue)
  row.names(sumtab) <- varnames
  # round figures appropriately:
  vardigitsR <- rep(0, nvar)
  for (i in 1:nvar) {
    vardigitsR[i] <- decplaces(sumtab[i,"std. dev."], signifdigits)
    if (vardigitsR[i]>0)
      sumtab[i,] <- round(sumtab[i,], vardigitsR[i])
    else 
      sumtab[i,] <- signif(sumtab[i,], signifdigits)
  }
  # "vardigitsR" now gives the number of digits to be displayed to the RIGHT of the decimal point.
  vardigitsL <- apply(sumtab, 1, function(x){max(c(1, trunc(log10(max(abs(x))))+1), na.rm=T)})
  # "vardigitsL" gives the number of digits to be displayed to the decimal point's LEFT.

  
  xsumtab <- xtable(sumtab, digits=cbind(rep(0,nvar),matrix(vardigitsR, nrow=nvar, ncol=ncol(sumtab))))

  # create text format table:
  textsumtab <- format(c("", varnames), justify="left")
  colwidths <- sapply(strsplit(c(textsumtab[1], colnames(sumtab)), split=""), length)
  for (j in 2:(ncol(sumtab)+1))
    colwidths[j] <- max(c(colwidths[j], vardigitsL + vardigitsR + 2))
  # first row:
  for (j in 2:(ncol(sumtab)+1))
    textsumtab[1] <- paste(textsumtab[1], format(colnames(sumtab)[j-1], width=colwidths[j], justify="right"), sep="  ")
  # further rows:
  for (i in 2:(nrow(sumtab)+1))
    for (j in 2:(ncol(sumtab)+1))
      textsumtab[i] <- paste(textsumtab[i], formatC(sumtab[i-1,j-1], format="f", digits=vardigitsR[i-1], width=colwidths[j]), sep="  ")
  #print(textsumtab)
  
  # write HTML table to file:
  print(xsumtab, type="html", file=htmlfilename, append=TRUE)
  
  # write LaTeX table to file:
  print(xsumtab, type="latex", file=paste(targetdirectory,'/latex/summarystats.TeX',sep=''), append=FALSE)

  # write TXT table to file:
  write(textsumtab, file=paste(targetdirectory,'/txt/summarystats.txt',sep=''), append=FALSE)

  write(c('  <small>(same table in <a href="txt/summarystats.txt" type="text/plain">text format</a>',
          '          and in <a href="latex/summarystats.TeX" type="text/plain">LaTeX format</a>)</small><br>',
          '  <br> <br>', ''),
        file=htmlfilename, append=TRUE)

  
  # compute correlation matrix:
  write(c("correlation matrix:", ""),
        file=htmlfilename, append=TRUE)
  cormat <- round(cor(posteriorsample), signifdigits)
  dimnames(cormat) <- list(varnames,varnames)
  
  xcormat <- xtable(cormat,digits=signifdigits)
  
  # create text format table:
  textcormat <- format(c("", varnames), justify="left")
  colwidths <- sapply(strsplit(c(textcormat[1], varnames), split=""), length)
  for (j in 2:(nvar+1))
    colwidths[j] <- max(c(colwidths[j], signifdigits + 3))
  # first row:
  for (j in 2:(nvar+1))
    textcormat[1] <- paste(textcormat[1], format(varnames[j-1], width=colwidths[j], justify="right"), sep="  ")
  # further rows:
  for (i in 2:(nvar+1))
    for (j in 2:(nvar+1))
      textcormat[i] <- paste(textcormat[i], formatC(cormat[i-1,j-1], format="f", digits=signifdigits, width=colwidths[j]), sep="  ")
  
  # write HTML table to file:
  print(xcormat, type="html", file=htmlfilename, append=TRUE)
  # write LaTeX table to file:
  print(xcormat, type="latex", file=paste(targetdirectory,'/latex/correlations.TeX',sep=''), append=FALSE)
  # write TXT table to file:
  write(textcormat, file=paste(targetdirectory,'/txt/correlations.txt',sep=''), append=FALSE)
  
  write(c('  <small>(same table in <a href="txt/correlations.txt" type="text/plain">text format</a>',
          '          and in <a href="latex/correlations.TeX" type="text/plain">LaTeX format</a>)</small><br>'),
        file=htmlfilename, append=TRUE)
  write(htmlhorizline(),
        file=htmlfilename, append=TRUE)

  
  # create & insert marginal density plots:
  write(htmlcomment("TRACE PLOTS, KERNEL DENSITY ESTIMATES and HISTOGRAMS:"),
        file=htmlfilename, append=TRUE)
  write(htmlsection(1, paste("Trace plots and marginal density estimates")),
        file=htmlfilename, append=TRUE)
  if (length(traceobs)<nobs)
  write(c(paste("  (only",formatC(length(traceobs),digits=0,format="f"),"out of all",formatC(nobs,digits=0,format="f"),"samples are shown in trace plots.)"), ""),
        file=htmlfilename, append=TRUE)
  # first insert posterior density trace plot -- if posterior was supplied!
  if (!is.null(posterior)) {
    write(htmlfragment("posteriordensitytrace"),
          file=htmlfilename, append=TRUE)
    write(htmlsection(2, "posterior density trace:"),
          file=htmlfilename, append=TRUE)
    filename <- "posteriordensitytrace"
    for (j in 1:ngraph){
      graphheight <- 4
      graphwidth <- graphheight * 1.5
      switch(graphicsformats[j],
             jpeg = bitmap(paste(targetdirectory, "/jpeg/", filename, ".jpeg", sep=""), type="jpeg", width=graphwidth/2, height=graphheight/2, res=144),
             png  = bitmap(paste(targetdirectory, "/png/", filename, ".png", sep=""), type="png16m", width=graphwidth/2, height=graphheight/2, res=144),
             eps  = postscript(paste(targetdirectory, "/eps/", filename, ".eps", sep=""), paper="special", width=graphwidth, height=graphheight),
             pdf  = pdf(paste(targetdirectory, "/pdf/", filename, ".pdf", sep=""), paper="special", width=graphwidth, height=graphheight))
      par(mar=c(4,4,1,1)+0.1)
      traceplot(iteration[traceobs], as.matrix(posterior)[traceobs,], ylab="log(posterior density)", col=colvec)
      abline(v=burnin, lty="dotted", col=colvec)
      box()
      dev.off()
    }
    write(c('  <table><tr>', '  <td align="center">'),
          file=htmlfilename, append=TRUE)
    write(htmlplot(filename, graphicsformats),  #, wi=300, he=300),
          file=htmlfilename, append=TRUE)
    write(c('  </td>', '  </tr></table>', ''),
          file=htmlfilename, append=TRUE)    
  }
  
  # parameter trace plots &c.:
  for (i in 1:nvar){
    write(htmlfragment(paste("parameter1D_", varnames[i], sep="")),
          file=htmlfilename, append=TRUE)
    write(htmlsection(2, paste("parameter: <i>", varnames[i], "</i>", sep="")),
          file=htmlfilename, append=TRUE)
    
    filename <- paste("traceplot-", varnames[i], sep="")
    for (j in 1:ngraph){
      graphheight <- 4
      graphwidth <- graphheight * 1.5
      switch(graphicsformats[j],
             jpeg = bitmap(paste(targetdirectory, "/jpeg/", filename, ".jpeg", sep=""), type="jpeg", width=graphwidth/2, height=graphheight/2, res=144),
             png  = bitmap(paste(targetdirectory, "/png/", filename, ".png", sep=""), type="png16m", width=graphwidth/2, height=graphheight/2, res=144),
             eps  = postscript(paste(targetdirectory, "/eps/", filename, ".eps", sep=""), paper="special", width=graphwidth, height=graphheight),
             pdf  = pdf(paste(targetdirectory, "/pdf/", filename, ".pdf", sep=""), paper="special", width=graphwidth, height=graphheight))
      par(mar=c(4,4,1,1)+0.1)
      traceplot(iteration[traceobs], data[traceobs,i,], ylab=varnames[i], col=colvec)
      if ((!is.null(truevalue)) && (!is.na(truevalue[i]))) abline(h=truevalue[i], lty="dashed", col="darkgrey")
      abline(v=burnin, lty="dotted", col=colvec)
      box()
      dev.off()
    }
    write(c('  <table><tr>', '  <td align="center">'),
          file=htmlfilename, append=TRUE)
    write(htmlplot(filename, graphicsformats),  #, wi=300, he=300),
          file=htmlfilename, append=TRUE)
    write(c('  </td>'),
          file=htmlfilename, append=TRUE)
    
    filename <- paste("kde-", varnames[i], sep="")
    kde <- computeKDE(posteriorsample[,i])
    for (j in 1:ngraph){
      graphheight <- 4
      graphwidth <- graphheight
      switch(graphicsformats[j],
             jpeg = bitmap(paste(targetdirectory, "/jpeg/", filename, ".jpeg", sep=""), type="jpeg", width=graphwidth/2, height=graphheight/2, res=144),
             png  = bitmap(paste(targetdirectory, "/png/", filename, ".png", sep=""), type="png16m", width=graphwidth/2, height=graphheight/2, res=144),
             eps  = postscript(paste(targetdirectory, "/eps/", filename, ".eps", sep=""), paper="special", width=graphwidth, height=graphheight),
             pdf  = pdf(paste(targetdirectory, "/pdf/", filename, ".pdf", sep=""), paper="special", width=graphwidth, height=graphheight))
      par(mar=c(4,4,1,1)+0.1)
      plotKDE(kde, main="", xlab=varnames[i], ylab="density", trueval=truevalue[i])
      #if ((!is.null(truevalue)) && (!is.na(truevalue[i]))) abline(v=truevalue[i], lty="dashed", col="darkgrey")
      #box()
      dev.off()
    }
    write(c('  <td align="center">'),
          file=htmlfilename, append=TRUE)
    write(htmlplot(filename, graphicsformats),  #, wi=300, he=300),
          file=htmlfilename, append=TRUE)
    write(c('  </td>'),
          file=htmlfilename, append=TRUE)

    filename <- paste("histogram-", varnames[i], sep="")
    for (j in 1:ngraph){
      graphheight <- 4
      graphwidth <- graphheight
      switch(graphicsformats[j],
             jpeg = bitmap(paste(targetdirectory, "/jpeg/", filename, ".jpeg", sep=""), type="jpeg", width=graphwidth/2, height=graphheight/2, res=144),
             png  = bitmap(paste(targetdirectory, "/png/", filename, ".png", sep=""), type="png16m", width=graphwidth/2, height=graphheight/2, res=144),
             eps  = postscript(paste(targetdirectory, "/eps/", filename, ".eps", sep=""), paper="special", width=graphwidth, height=graphheight),
             pdf  = pdf(paste(targetdirectory, "/pdf/", filename, ".pdf", sep=""), paper="special", width=graphwidth, height=graphheight))
      par(mar=c(4,4,1,1)+0.1)
      hist(posteriorsample[,i], prob=TRUE, col="lightgrey",
           main="", xlab=varnames[i], ylab="density")
      if ((!is.null(truevalue)) && (!is.na(truevalue[i]))) abline(v=truevalue[i], lty="dashed", col="darkgrey")
      abline(h=0, col="black")
      box()
      dev.off()
    }
    
    write(c('  <td align="center">'),
          file=htmlfilename, append=TRUE)
    write(htmlplot(filename, graphicsformats),  #, wi=300, he=300),
          file=htmlfilename, append=TRUE)
    write(c('  </td>', '  </tr></table>', ''),
          file=htmlfilename, append=TRUE)
  }

  
  write(htmlhorizline(),
        file=htmlfilename, append=TRUE)

  
  # create bivariate distributions:
  write(htmlcomment("BIVARIATE DISTRIBUTION PLOTS:"),
        file=htmlfilename, append=TRUE)
  write(htmlsection(1, paste("Bivariate marginal distributions")),
        file=htmlfilename, append=TRUE)
  if (length(scatterobs)<nrow(posteriorsample))
  write(c(paste("  (only",formatC(length(scatterobs),digits=0,format="f"),"out of all",formatC(nrow(posteriorsample),digits=0,format="f"),"samples are shown in scatter plots.)"), ""),
        file=htmlfilename, append=TRUE)
  graphheight <- 3
  graphwidth <- graphheight
  for (i in 1:nvar){ # loop over all variables
    write(htmlfragment(paste("parameter2D_", varnames[i], sep="")),
          file=htmlfilename, append=TRUE)
    write(htmlsection(2, paste("parameter: <i>", varnames[i], "</i>", sep="")),
          file=htmlfilename, append=TRUE)
    write(c('  <table><tr>'),
          file=htmlfilename, append=TRUE)
    filename1 <- paste("scatterplot-", varnames[i], sep="")

    for (k in ((1:nvar)[-i])) { # loop over variables other than i
      write(c('  <td align="center">'),
            file=htmlfilename, append=TRUE)
      filename2 <- paste(filename1, "-", varnames[k], sep="")
      for (j in 1:ngraph){
        switch(graphicsformats[j],
               jpeg = bitmap(paste(targetdirectory, "/jpeg/", filename2, ".jpeg", sep=""), type="jpeg", width=graphwidth/2, height=graphheight/2, res=144),
               png  = bitmap(paste(targetdirectory, "/png/", filename2, ".png", sep=""), type="png16m", width=graphwidth/2, height=graphheight/2, res=144),
               eps  = postscript(paste(targetdirectory, "/eps/", filename2, ".eps", sep=""), paper="special", width=graphwidth, height=graphheight),
               pdf  = pdf(paste(targetdirectory, "/pdf/", filename2, ".pdf", sep=""), paper="special", width=graphwidth, height=graphheight))
        par(mar=c(4,4,1,1)+0.1)
        #traceplot(iteration[traceobs], data[traceobs,i,], ylab=varnames[i], col=colvec)
        scatterplot(posteriorsample[scatterobs,i], posteriorsample[scatterobs,k],
                    xlab=varnames[i], ylab=varnames[k], col=colvec[chain[scatterobs]],
                    trueval=c(truevalue[i], truevalue[k]))
        #if (!is.null(truevalue)){
        #  if (!is.na(truevalue[i])) abline(v=truevalue[i], lty="dashed", col="darkgrey")
        #  if (!is.na(truevalue[i])) abline(h=truevalue[k], lty="dashed", col="darkgrey")
        #  box()
        #}
        dev.off()
      }
      write(htmlplot(filename2, graphicsformats),
            file=htmlfilename, append=TRUE)
      write(c('  </td>'),
            file=htmlfilename, append=TRUE)  
    }
    
    write(c('  </tr></table>', ''),
          file=htmlfilename, append=TRUE)
  }
  
  # finish the HTML file:
  write(htmlcomment("PAGE FOOTER:"),
        file=htmlfilename, append=TRUE)
  write(htmlfooter(),
        file=htmlfilename, append=TRUE)
  
  return(invisible(posteriorsample))
}
