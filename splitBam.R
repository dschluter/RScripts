#!/usr/bin/Rscript

# Splits the bamfile by chromosomes and saves the results as indexed bam files.
# Use as "Rscript splitBam.R BAMFILENAME"

# Uses filterBam from Rsamtools:
# filterBam parses records in file. Records satisfying the bamWhich bamFlag and bamSimpleCigar
# criteria of param are accumulated to a default of yieldSize = 1000000 records (change this by
# specifying yieldSize when creating a BamFile instance; see BamFile-class). These records are# then parsed to a DataFrame and made available for further filtering by user-supplied FilterRules.# Functions in the FilterRules instance should expect a single DataFrame argument representing all# information specified by param. Each function must return a logical vector equal to the number of# rows of the DataFrame. Return values are used to include (when TRUE) corresponding records in the# filtered BAM file. The BAM file is created at destination. An index file is created on the destination# when indexDestination=TRUE. It is more space- and time-efficient to filter using bamWhich,# bamFlag, and bamSimpleCigar, if appropriate, than to supply FilterRules.

# BamFile() creates a reference to a BAM file. The reference remains
# 	open across calls to methods, avoiding costly index re-loading.
# However, the connection seems to need repositioning after each run through a loop while file remains open
# So it might not be faster to use after all.

args <- commandArgs(TRUE) # first argument will be bamfile name
# args <- "Marine-Pac-Salmon-01-Sara.recal.bam"

bamFileName <- args[1]
print(bamFileName)

library(Rsamtools, quietly = TRUE)
root <- gsub("[.]bam$", "", bamFileName) # fishname part of bam file name

# Extract chromosome information from the bamfile header
z <- scanBamHeader(bamFileName)[[1]]
maxLen <- max(z$targets) + 1 # Must be a number larger than no. bases in largest chromosome
chrname <- names(z$targets)

# This didn't work, no destination file was created
# for(i in chrname){
	# param <- ScanBamParam(what=scanBamWhat(), which = GRanges(i, IRanges(1, maxLen)))
	# x <- scanBam(bamFileName, param = param, destination = paste(root, i, sep = "."), 
		# indexDestination = TRUE)[[1]]
	# }

# Need to close and reopen file otherwise the connection doesn't rewind to the start.
# Can try using "seek" to fix this.
for(i in chrname){
	# i <- "chrM"
	# i <- "chrII"
	outfilename <- paste(root, i, "bam", sep = ".")
	print(outfilename)
	bamSource <- open(BamFile(bamFileName, index = paste(root, "bai", sep = "."))) #, yieldSize = maxLen))
	param <- ScanBamParam(what=scanBamWhat(), which = GRanges(i, IRanges(1, maxLen)))
	filterBam(bamSource, destination = outfilename, # filter=FilterRules(), 
		indexDestination=TRUE, param = param)
	close(bamSource)
	}

# check that it actually happened
# x <- scanBam("Marine-Pac-Salmon-01-Sara.recal.chrM.bam")[[1]]
