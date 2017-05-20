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


g$splitGVCF(gvcffile = gvcffile, ngroups = ngroups, genome = genome, chunksize = chunksize,
	nFirstLines = nFirstLines, nLinesAtaTime = nLinesAtaTime, includeAllContigs = includeAllContigs)

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
