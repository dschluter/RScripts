#!/usr/bin/Rscript

# Run as:   /global/software/R-3.1.2/bin/Rscript bamCoverage.R Solitary*.chrXXI.realigned.bam
# or paste above command into bamCoverage.R.pbs and run: qsub bamCoverage.R.pbs
# Take a set of chrxxx.realigned.bam files and determine mean and median read coverage

arguments <- commandArgs(TRUE) # only arguments will be stored
# arguments <- c("Solitary*.chrXXI.realigned.bam")

# test for regular expression as first argument
fnames <- list.files(pattern = glob2rx(arguments[1])) 
bamfiles <- c(fnames, arguments[-1]) # in case there are additional arguments
nbam <- length(bamfiles)

cat("\nbam files inputted:\n")
for(i in 1:nbam){
	cat(bamfiles[i],"\n")
	}

library(Rsamtools, quietly = TRUE)
mean.coverage <- vector()
median.coverage <- vector()
for( i in 1:nbam ){
	# i <- 1
	# root <- gsub("(.)[.]bam$", "\\1", bamfiles[i]) # when bai files are like "paxb04.recal.chrXXI.bai"
	# baifile <- paste(root, "bai", sep = ".")
	# print(baifile)
	baifile <- paste(bamfiles[i], "bai", sep = ".") # when bai files are like "paxb04.recal.chrXXI.bam.bai"
	# print(baifile)
	x <- scanBam(bamfiles[i], index = baifile)
	x1 <- x[[1]]
  	ind <- !is.na(x1$pos)
	bam <- lapply(x1, function(x) x[ind]) # strips the bam object of unaligned reads
  	ranges1 <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
	cover <- as.vector(coverage(ranges1))
	mean.coverage[i] <- mean(cover)
	median.coverage[i] <- median(cover)
	}
	
print( cbind.data.frame(bamfiles,mean=round(mean.coverage,1),median= median.coverage) )

