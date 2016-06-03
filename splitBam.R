#!/usr/bin/Rscript

# Takes a single bamfile, reads the CHROM column, and splits the bamfile accordingly, 
#		saving the results as sam files.

args <- commandArgs(TRUE) # first argument will be bamfile name
# args <- "Marine-Pac-Salmon-01-Sara.recal.bam"

bamFileName <- args[1]


library(Rsamtools, quietly = TRUE)
root <- gsub("[.]bam$", "", bamFileName) # part of file name

# Extract the chromosome names
param <- ScanBamParam(what = c("rname"))
x <- scanBam(bamFileName, param = param)[[1]]
chrname <- unique(as.character(x$rname))
chrname <- chrnames[!is.na(chrnames)]

g$splitbam(bamFileName, chrname)

maxLen <- 10^8 # a number larger than the largest chromosome

for(i in chrname){
	param <- ScanBamParam(what=scanBamWhat(), which = GRanges(i, IRanges(1, maxLen)))
	x <- scanBam(bamFileName, param = param, destination = root, indexDestination = TRUE)
	}
