
# Compute PSD average from SFT list, call as below giving lists of sfts as needed

#R --no-save  -f ./fast_averager.R --slave --args flist.txt [flist2.txt ...]


FilePSD<-function(filename, output) {
	Header<-list(key="real8", gps_sec="int", gps_nsec="int", tbase="real8", first_frequency_index="int", nsamples="int", crc="int64", detector="c2", padding1="c2", comment_length="int")
	
	con<-file(filename, open="rb")
	
	Values<-list()
	for(field in names(Header)) {
		Values[[field]]<-switch(Header[[field]], 
			real8=readBin(con, "double", 1, size=8),
			real4=readBin(con, "double", 1, size=4),
			int=readBin(con, "int", 1, size=4),
			c2=readBin(con, "character", 1, size=2),
			int64=readBin(con, "integer", 1, size=8),
			NA)
		cat(field, "=", Values[[field]], "\n")
		}
	
	#comment<-readBin(con, "character", 1, size=Values$comment_length)
	comment<-readChar(con, Values$comment_length,useBytes=TRUE)
	cat('comment="', comment, '"\n', sep="")
	data<- matrix(readBin(con, "double", Values$nsamples*2, size=4))
	dim(data)<-c(2, Values$nsamples)
	data<-t(data)

	close(con)

	print(data[1:10,])
	
	power<- data[,1]^2+data[,2]^2
	
	step<-round(Values$tbase)
	
	f0<-Values$first_frequency_index/Values$tbase
	i<-0
	con<-file(output, open="w")
	cat(file=con, "Frequency\tMedian\tMean\n")
	while(i<Values$nsamples) {
		k<-power[i:(i+step-1)]
	
		cat(file=con, f0+i*1.0/Values$tbase, "\t",  median(k), "\t", mean(k), "\n", sep="")
		i<-i+step
		}
	close(con)
	return(invisible(NULL))
	}

for(lfilename in commandArgs(trailingOnly=TRUE)) {

	con<-file(lfilename, open="r")

	cat("Scanning", lfilename, "\n")
	while ( 1) {
		filename<-readLines(con, 1)
		if(length(filename)<1) { 
			#cat("End of input\n")
			break
			}
		if(filename=="") { 
			next
			}
	
		cat("\t", filename, "\n")
		FilePSD(filename, paste("psd_", basename(filename), ".txt", sep=""))
		}
	close(con)
	}

