#!/usr/bin/Rscript

gvcffile 	<- NULL
ngroups 	<- NULL
genome 		<- NULL
chunksize 	<- 10000
nFirstLines <- 10000
nLinesAtaTime <- 10000
includeAllContigs <- TRUE

args <- commandArgs(TRUE) 
# args <- c("gvcffile=Pup7.QTL_151.realigned.g.vcf.gz","ngroups=5","genome=pupfish9000scaffolds.fa")
# args <- c("gvcffile=Pup1.QTL_10.g.vcf.gz","ngroups=5","genome=pupfish9000scaffolds.fa")

# Parses the args into a data frame with two columns and then assigns variables 
x <- read.table(text = args, sep = "=", colClasses = "character")
for(i in 1:nrow(x)){assign(x[i,1], x[i,2])}
print(x)
        # V1                              V2
# 1 gvcffile Pup7.QTL_151.realigned.g.vcf.gz
# 2  ngroups                              10
# 3   genome         pupfish9000scaffolds.fa

if(is.null(gvcffile)) stop("Provide gvcffile= in arguments")
if(is.null(ngroups)) stop("Provide ngroups= in arguments")
if(is.null(genome)) stop("Provide genome= in arguments")

ngroups 	<- as.integer(ngroups)
chunksize 	<- as.integer(chunksize)
nFirstLines <- as.integer(nFirstLines)
nLinesAtaTime <- as.integer(nLinesAtaTime)


# g$splitGVCF(gvcffile = gvcffile, ngroups = ngroups, genome = genome, chunksize = chunksize,
	# nFirstLines = nFirstLines, nLinesAtaTime = nLinesAtaTime, includeAllContigs = includeAllContigs)

# Moved the code back to here from genome.r

# g$splitGVCF <- function(
	# gvcffile,
	# ngroups =  5,
	# genome = "pupfish9000scaffolds.fa",
	# chunksize = 10000,
	# nFirstLines = 10000,
	# nLinesAtaTime = 10000,
	# includeAllContigs = TRUE){

	# Function to split a g.vcf.gz file into about equal-sized parts
	# nFirstLines must be large enough to contain the entire header (~10000 lines in pupfish)
	# ngroups is the number of contig groups to split into
	# If includeAllContigs=FALSE, only relevant contig lines are included in each vcf file portion

	gvcfname <- gsub(".g.vcf.gz", "", gvcffile)
	
	# Open file and read nFirstLines lines from vcf file
	INFILE <- file(gvcffile, open="r")
	# close(INFILE)
	
	# read the first lines
	x <- readLines(con = INFILE, n = nFirstLines) 
	
	# get number of header lines ('#')
	nhead <- length(grep('^#', x)) # 9353 in pupfish
	if(nhead == nFirstLines) stop("nFirstLines is too small -- increase to greater than no. header lines")
	
	# Grab the '##contig=' lines, label line, and ref genome line, and remaining header lines
	whichcontig <- grep('^##contig', x)
	contiglines <- x[whichcontig]
	headerlines <- x[1:(whichcontig[1] - 1)]
	labelline <- x[nhead] # line with '#CHROM	POS	ID...'
	refline <- x[nhead - 1] # line with '##reference=file:...'
	
	# Hive off the data portion of the lines
	x <- x[ (nhead+1):length(x) ]
	nlines <- length(x)
	
	# head(contiglines)
	# [1] "##contig=<ID=NW_015150453.1,length=4516015>"
	# [2] "##contig=<ID=NW_015150454.1,length=3735821>"
	# [3] "##contig=<ID=NW_015150455.1,length=3397602>"
	# [4] "##contig=<ID=NW_015150456.1,length=3333070>"
	# [5] "##contig=<ID=NW_015150457.1,length=3291551>"
	# [6] "##contig=<ID=NW_015150458.1,length=3169797>"
	
	# grab the contig lengths from the contig lines
	contignumbers <- as.integer( gsub("##contig=<ID=.*length=([0-9]*)>$", "\\1", contiglines) )
	cumlengths <- cumsum(contignumbers)
	
	# break the contigs into groups of about the same total contig lengths
	maxCL <- max(cumlengths) 
	groupbreaks <- seq(0, maxCL+1, by = round(maxCL/ngroups))
	groupbreaks[length(groupbreaks)] <- maxCL+1 # in case the breaks don't include the highest avlue
	
	# Identify contig bins of about equal cumulative lengths 
	contigBins <- cut(cumlengths, breaks = groupbreaks)

	# table(contigBins) # tabulates number of contigs in each bin
	
	# Split the list of contigs accordingly
	contigsplit <- split(contiglines, contigBins)
	
	# grab the contig names included in each bin
	contiggroups <- lapply(contigsplit, function(z){
		# z <- contigsplit[[1]]
		contiggroups <- gsub("##contig=<ID=(.*),length=[0-9]*>$", "\\1", z)
		})
	
	# Name the new pseudo "chr" groups as part1, part2, etc
	chrnames <- paste("part", 1:ngroups, sep = "")
	
	# Open multiple write connections
	outfiles <- vector("list")
	for(i in chrnames) outfiles[[i]] <- file(paste(gvcfname, i, "g.vcf", sep="."), "w")
	# for(i in chrnames) close(outfiles[[i]])
	
	filenames <- paste(gvcfname, chrnames, "g.vcf", sep=".")

	# Write header lines to files
	for( i in 1:length(chrnames) ){
		writeLines(headerlines, outfiles[[chrnames[i]]])
		if(includeAllContigs)
			writeLines(contiglines, outfiles[[chrnames[i]]]) else
			writeLines(contigsplit[[i]], outfiles[[chrnames[i]]])
		writeLines(refline, outfiles[[chrnames[i]]])
		writeLines(labelline, outfiles[[chrnames[i]]])
		}
	
	# Loop
	while(nlines > 0){
	
		# grab the CHROM column numbers from the vcf lines
		CHROM <- gsub("^([A-z0-9_.-]*)\t.*", "\\1", x)
	
		# print contents to corresponding files
		for(i in 1:length(contiggroups)){
			# i <- 1
			k <- which(is.element(CHROM, contiggroups[[i]]))
			if(length(k) > 0) writeLines(x[k], outfiles[[i]])
			}
	
		x <- readLines(con = INFILE, n = nLinesAtaTime)
		nlines <- length(x)
	
		} # end while loop
	
	close(INFILE)
	
	for(i in chrnames) close(outfiles[[i]])
	
	# Index resulting vcf files (need to bgzip first)
	# Function to be used (make standalone rather than use g$vcfBgzip):
	vcfBgzip <- function(vcfname = "", renameAsGz = FALSE){
		library(VariantAnnotation, quietly = TRUE)
	
		root <- sub("[.]gz$", "", vcfname) # remove the "gz" if present (need not be)
		compressedVcfname <- paste(root, "bgz", sep = ".") 
			
		# compress and save file name
		bgzname <- bgzip(vcfname, compressedVcfname)
		
		# rename file
		if(renameAsGz){
			file.rename(bgzname, vcfname)
			idx <- indexTabix(vcfname, "vcf")
			} else {
			idx <- indexTabix(bgzname, "vcf")
			}
		}

	for(i in 1:length(chrnames)) vcfBgzip(filenames[i])
	
	# Can delete the uncompressed .vcf files
	for(i in 1:length(chrnames)) file.remove(filenames[i])
	
	# }

# ---- old -----
# This version is designed for Jess' RNAseq gvcf files
# TRanscripts numbered TR1 to TR125591
# setwd("/Users/schluter/Documents/Research/genomics/jess salmon RNAseq")

# vcfname <- gsub(".vcf", "", vcffile)

# nFirstLines <- 400000 # Make this >=400K because header is 395175 lines -- all transcript are given as ##contig=
# nLinesAtaTime <- 10000
# maxTR <- 125591   # the largest number corresponding to a transcript
# ngroups <-  30    # the number of transcript groups to split into
# groupbreaks <- seq(0, maxTR+1, by = round(maxTR/ngroups))
# groupbreaks[length(groupbreaks)] <- maxTR+1 # in case the breaks don't include the highest avlue

# # Open file and read nFirstLines lines from vcf file
# INFILE <- file(vcffile, open="r")
# #close(INFILE)

# # according to http://www.r-bloggers.com/faster-files-in-r/ the following reads much faster
# # I'd need to fiddle with the last line if it is not all the way to the "\n" at the end
# # buf = readChar( fname, s, useBytes=T)
# # strsplit( buf,"\n",fixed=T,useBytes=T)[[1]]

# x <- readLines(con = INFILE, n = nFirstLines) 

# # get number of header lines ('#')
# nhead <- length(grep('^#', x))

# # Grab the '##contig=' lines, label line, and ref genome line, and remaining header lines
# whichcontig <- grep('^##contig', x)
# contiglines <- x[whichcontig]
# headerlines <- x[1:(whichcontig[1] - 1)]
# labelline <- x[nhead] # line with '#CHROM	POS	ID...'
# refline <- x[nhead - 1] # line with '##reference=file:...'

# # Hive off the data portion of the lines
# x <- x[ (nhead+1):length(x) ]
# nlines <- length(x)

# # grab the transcript numbers from the contig lines
# contignumbers <- as.integer( gsub("##contig=<ID=TR([0-9]+)[|].*", "\\1", contiglines) )
# contigBins <- cut(contignumbers, breaks = groupbreaks)
# contigsplit <- split(contiglines, contigBins)
# # table(contigBins)
	
# # chrnames <- gsub("^@[SQPG]+\tSN:([A-Za-z0-9.-_]+)\t[LN:0-9]+", "\\1", headerlines[-nhead]) # awkward but faster
# # chrnames <- sapply( strsplit(headerlines[-nhead], split = "[\t:]+"), function(z){z[3]}) # takes twice as long
# chrnames <- paste("part", 1:ngroups, sep = "")

# # Open multiple write connections
# outfiles <- vector("list")
# for(i in chrnames) outfiles[[i]] <- file(paste(vcfname, i, "vcf", sep="."), "w")
# # for(i in chrnames) close(outfiles[[i]])

# # Write header lines to files
# for( i in 1:length(chrnames) ){
	# writeLines(headerlines, outfiles[[chrnames[i]]])
	# writeLines(contigsplit[[i]], outfiles[[chrnames[i]]])
	# writeLines(refline, outfiles[[chrnames[i]]])
	# writeLines(labelline, outfiles[[chrnames[i]]])
	# }

# # Loop
# while(nlines > 0){

	# # grab the transcript numbers from the vcf lines
	# transcriptnumbers <- as.integer( gsub("TR([0-9]+)[|].*", "\\1", x) )
	# transcriptBins <- cut(transcriptnumbers, breaks = groupbreaks)

	# xsplit <- split(x, transcriptBins)

	# # print contents to corresponding files
	# for(i in 1:length(xsplit)) writeLines(xsplit[[i]], outfiles[[i]])

	# x <- readLines(con = INFILE, n = nLinesAtaTime)
	# nlines <- length(x)

	# } # end while loop

# close(INFILE)
# for(i in chrnames) close(outfiles[[i]])
